proc version_flair {} {
	if {[catch {catch_exec flair --version} version]} {
		set flairdir [findflair]
		set temp [split [file tail $flairdir] -]
		set version [lindex $temp 1]
		return $version
	}
	lindex $version end
}

proc flair_bin {bin} {
	if {[catch {exec which $bin} msg]} {
		return $bin.py
	} else {
		return $bin
	}
}

proc flair_getref {refseq} {
	set refseq [refseq $refseq]
	set refdir [file dir $refseq]
	if {[catch {
		set ref [lindex [bsort [glob $refdir/extra/gencode*.gtf]] end]
	}]} {
		set keeppwd [pwd]
		cd $refdir/hg38/extra
	#	exec wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
	#	exec gunzip gencode.v37.annotation.gtf.gz
		exec wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
		exec gunzip gencode.v39.annotation.gtf.gz
		cd $keeppwd
	}
	set ref [lindex [bsort [glob $refdir/extra/gencode*.gtf]] end]
}

proc findflair {} {
	global flair
	if {![info exists flair]} {
		set flair [searchpath flair flair flair*]
		if {$flair eq ""} {
			set flair [searchpath FLAIR flair flair*]
		}
		set flair [file_resolve $flair]
		if {![file isdir $flair]} {set flair [file dir $flair]}
		set ::env(PATH) $flair:$::env(PATH)
	}
	return $flair
}

proc cg_flair_genecounts {isoformcounts genecounts {genefield associated_gene}} {
	set tempfile [tempfile]
	cg select -overwrite 1 \
		-g $genefield \
		-gc {distinct(structural_category),distinct(subcategory),distinct(chromosome),min(begin),max(end),distinct(strand),ucount(transcript),sum(counts-*),sum(tpm-*)} \
		$isoformcounts $tempfile
	catch {close $f} ; catch {close $o}
	set f [open $tempfile]
	set header [tsv_open $f]
	set oheader {type gene gene_type chromosome begin end strand nrtranscripts}
	foreach field [lrange $header 8 end] {
		regsub sum_ $field {} field
		lappend oheader $field
	}
	set o [open $genecounts w]
	puts $o [deindent {
		#filetype	tsv/genecountsfile
		#fileversion	0.99
		#fields	table
		#fields	field	number	type	description
		#fields	gene	1	String	gene id (usually gene_id from GTF)
		#fields	genetype	1	String	type of gene
		#fields	chromosome	1	String	Chromosome name
		#fields	strand	1	String	+ or - for strand
		#fields	begin	1	Integer	Transcription start position
		#fields	end	1	Integer	Transcription end position
		#fields	nrtranscripts	1	Integer	number of transcripts found
		#fields	counts	1	Integer	Number of reads mapping to gene
		#fields	tpm	1	Float	Number of gene transcripts per million total transcripts
		#fields	type	1	String	type of element
	}]
	puts $o [join $oheader \t]
	while {[gets $f line] != -1} {
		set line [split $line \t]
		foreach {gene category sub} $line break
		set novel [regexp {^novel} $gene]
		if {!$novel && ![regexp _ $gene]} {
			set genetype known
		} elseif {[regexp fusion $category]} {
			set genetype fusion
		} elseif {!$novel} {
			set genetype overlap
		} elseif {[regexp antisense $category]} {
			set genetype novel_antisense
		} elseif {[regexp intergenic $category]} {
			set genetype novel_intergenic
		} elseif {[regexp genic_intron $category]} {
			set genetype novel_genic_intron
		} else {
			set genetype novel
		}
		puts $o gene\t$gene\t$genetype\t[join [lrange $line 3 end] \t]
	}
	close $o
	close $f
}

# set totalcountsfile {}
# foreach {target transcript_classification_file transcripts_genepred_file counts_matrix_file totalcountsfile} $args break
proc cg_flair_mergeresults {target transcript_classification_file transcripts_genepred_file counts_matrix_file {totalcountsfile {}}} {
	#
	# read classifications in classa
	unset -nocomplain classa
	set f [gzopen $transcript_classification_file]
	set classheader [tsv_open $f]
	set classidpos [lsearch $classheader isoform]
	while {[gets $f line] != -1} {
		set line [split $line \t]
		set classid [lindex $line $classidpos]
		set classa($classid) $line
	}
	close $f
	# read counts in counta
	unset -nocomplain counta
	set f [gzopen $counts_matrix_file]
	set cntheader [tsv_open $f]
	set cntheader [list_regsub {_conditionA_batch1$} $cntheader {}]
	set numsamples [expr {[llength $cntheader] -1}]
	set totalcounts [list_fill $numsamples 0]
	while {[gets $f line] != -1} {
		set line [split $line \t]
		set cntid [lindex $line 0]
		set cntid [lindex [split $cntid _] 0]
		set counta($cntid) [lrange $line 1 end]
		set temp {}
		foreach tcnt $totalcounts cnt [lrange $line 1 end] {
			set tcnt [expr {$tcnt + $cnt}]
			lappend temp $tcnt
		}
		set totalcounts $temp
	}
	close $f
	set exproot [file_rootname $target]
	regsub \\.genepred $exproot {} exproot
	if {$totalcountsfile eq ""} {
		set totalcountsfile [file dir $target]/totalcounts-$exproot.tsv
	}
	file_write $totalcountsfile [join [lrange $cntheader 1 end] \t]\n[join $totalcounts \t]\n
	#
	# write combi file
	catch {close $fg} ; catch {close $o} ; 
	set fg [gzopen $transcripts_genepred_file]
	# set gheader [lrange [tsv_open $fg] 1 end]
	set gheader [tsv_open $fg]
	set gposs [tsv_basicfields $gheader 12 0]
	set gposs [list_sub $gposs {0 6 7 8}]
	set namepos [lsearch $gheader name]
	set genepos [lsearch $gheader gene]
	set scatpos [lsearch $classheader structural_category]
	set o [open $target.temp w]
	set newheader $gheader
	set newheader [list_replace $newheader {chrom chromosome txStart begin txEnd end}]
	lappend newheader category
	set pclassheader [list_remove $classheader id isoform chrom strand]
	set poss [list_cor $classheader $pclassheader]
	lappend newheader {*}$pclassheader
	lappend newheader type
	foreach sample [lrange $cntheader 1 end] {
		lappend newheader counts-flair-$sample
	}
	foreach sample [lrange $cntheader 1 end] {
		lappend newheader tpm-flair-$sample
	}
	puts $o [join $newheader \t]
	set emptycntresult [list_fill [expr {2*[llength $totalcounts]}] 0]
	set nr 1
	set gr 0
	array set cattrans {
		full-splice_match known
		intergenic intergenic
		incomplete-splice_match novel_incomplete 
		novel_in_catalog novel_in_catalog novel_not_in_catalog novel_not_in_catalog
		genic novel_genic antisense novel_antisense fusion novel_fusion genic_intron novel_genic_intron
	}
	while 1 {
		if {[gets $fg gline] == -1} break
		set gline [split $gline \t]
		set gid [lindex $gline 0]
		if {![info exists classa($gid)]} {
			error "class data not found for $gid"
		}
		set classline $classa($gid)
		set classid [lindex $classline $classidpos]
		set scat [lindex $classline $scatpos]
		if {[info exists cattrans($scat)]} {set scat $cattrans($scat)}
		if {$scat ne "known"} {
			foreach {chromosome strand exonStarts exonEnds} [list_sub $gline $gposs] break
			set name [iso_name $chromosome $strand $exonStarts $exonEnds]
			lset gline $namepos $name
		}
		set result [join $gline \t]\t$scat
		append result \t[join [list_sub $classline $poss] \t]
		if {[info exists counta($gid)]} {
			append result \ttranscript\t[join $counta($gid) \t]
			foreach cnt $counta($gid) tcnt $totalcounts {
				append result \t[expr {1000000.0*$cnt/$tcnt}]
			}
			puts $o $result
		} else {
			puts stderr "count not find count for $gid at line $nr ($transcripts_genepred_file and $counts_matrix_file)"
			append result \ttranscript\t[join $emptycntresult \t]
			puts $o $result
		}
		incr nr
	}
	close $o
	close $fg
	# write final target
	file_write $target.temp2 [iso_flair_comments]\n
	cg select -s - $target.temp >> $target.temp2
	file rename -force $target.temp2 $target
	file delete -force $target.temp
}

proc iso_flair_comments {} {
	deindent {
		#filetype	tsv/transcriptsfile
		#fileversion	0.99
		#fields	table
		#fields	field	number	type	description
		#fields	name	1	String	Name of gene (usually transcript_id from GTF)
		#fields	chromosome	1	String	Chromosome name
		#fields	strand	1	String	+ or - for strand
		#fields	begin	1	Integer	Transcription start position
		#fields	end	1	Integer	Transcription end position
		#fields	cdsStart	1	Integer	Coding region start
		#fields	cdsEnd	1	Integer	Coding region end
		#fields	exonCount	1	Integer	Number of exons
		#fields	exonStarts	E	Integer	Exon start positions
		#fields	exonEnds	E	Integer	Exon end positions
		#fields	score	1	Float	Score
		#fields	name2	1	String	Alternate name (e.g. gene_id from GTF)
		#fields	cdsStartStat	1	String	Status of CDS start annotation (none, unknown, incomplete, or complete)
		#fields	cdsEndStat	1	String	Status of CDS end annotation (none, unknown, incomplete, or complete)
		#fields	exonFrames	E	Integer	Exon frame offsets {0,1,2}
		#fields	length	1	Integer	isoform length
		#fields	exons	1	Integer	Number of exons
		#fields	structural_category	1	String	one of the SQANTI3 isoform categories based on the best matching reference transcript (https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-output-explanation)
		#fields	associated_gene	1	String	the reference gene name
		#fields	associated_transcript	1	String	the reference transcript name
		#fields	ref_length	1	Integer	reference transcript length
		#fields	ref_exons	1	Integer	reference transcript number of exons
		#fields	diff_to_TSS	1	Integer	distance of query isoform 5' start to reference transcript start end. Negative value means query starts downstream of reference
		#fields	diff_to_TTS	1	Integer	distance of query isoform 3' end to reference annotated end site. Negative value means query ends upstream of reference
		#fields	diff_to_gene_TSS	1	Integer	distance of query isoform 5' start to the closest start end of any transcripts of the matching gene
		#fields	diff_to_gene_TTS	1	Integer	distance of query isoform 3' end to the closest end of any transcripts of the matching gene
		#fields	subcategory	1	String	additional splicing categorization, separated by semi-colons
		#fields	RTS_stage	1	String	TRUE if one of the junctions could be a RT switching artifact
		#fields	all_canonical	1	String	TRUE if all junctions have canonical splice sites
		#fields	min_sample_cov	1	String	sample with minimum coverage
		#fields	min_cov	1	Integer	minimum junction coverage based on short read STAR junction output file. NA if no short read given
		#fields	min_cov_pos	1	Integer	the junction that had the fewest coverage. NA if no short read data given
		#fields	sd_cov	1	Integer	standard deviation of junction coverage counts from short read data. NA if no short read data given
		#fields	FL	1	Integer	FL count associated with this isoform per sample if --fl_count is provided, otherwise NA
		#fields	n_indels	1	Integer	total number of indels based on alignment
		#fields	n_indels_junc	1	Integer	number of junctions in this isoform that have alignment indels near the junction site (indicating potentially unreliable junctions)
		#fields	bite	1	String	TRUE if contains at least one "bite" positive SJ
		#fields	iso_exp	1	Integer	short read expression for this isoform if --expression is provided, otherwise NA
		#fields	gene_exp	1	Integer	short read expression for the gene associated with this isoform (summing over all isoforms) if --expression is provided, otherwise NA
		#fields	ratio_exp	1	Integer	ratio of iso_exp to gene_exp if --expression is provided, otherwise NA
		#fields	FSM_class	1	String	classifies the transcript according to the expression of other isoforms in the gene to which the transcript belongs
		#fields	coding	1	String	coding or non_coding transcript
		#fields	ORF_length	1	Integer	predicted ORF length
		#fields	CDS_length	1	Integer	predicted CDS length
		#fields	CDS_start	1	Integer	CDS start
		#fields	CDS_end	1	Integer	CDS end
		#fields	CDS_genomic_start	1	Integer	genomic coordinate of the CDS start. If on - strand, this coord will be greater than the end
		#fields	CDS_genomic_end	1	Integer	genomic coordinate of the CDS end. If on - strand, this coord will be smaller than the start
		#fields	predicted_NMD	1	String	TRUE if there's a predicted ORF and CDS ends at least 50bp before the last junction; FALSE if otherwise. NA if non-coding
		#fields	perc_A_downstream_TTS	1	Float	percent of genomic "A"s in the downstream 20 bp window. If this number if high (say > 0.8), the 3' end site of this isoform is probably not reliable
		#fields	seq_A_downstream_TTS	1	String	sequence of the downstream 20 bp window
		#fields	dist_to_cage_peak	1	Integer	distance to closest TSS based on CAGE Peak data
		#fields	within_cage_peak	1	String	TRUE if the PacBio transcript start site is within a CAGE Peak
		#fields	dist_to_polya_site	1	Integer	if --polyA_motif_list is given, shows the location of the last base of the hexamer. Position 0 is the putative poly(A) site. This distance is hence always negative because it is upstream
		#fields	within_polya_site	1	String	
		#fields	polyA_motif	1	String	if --polyA_motif_list is given, shows the top ranking polyA motif found within 50 bp upstream of end
		#fields	polyA_dist	1	Integer	: if --polyA_motif_list is given, shows the location of the last base of the hexamer. Position 0 is the putative poly(A) site. This distance is hence always negative because it is upstream
		#fields	ORF_seq	1	String	ORF sequence
		#fields	ratio_TSS	1	Float	Using Short-Read data, we measure the mean coverage of the 100bp upstream and downstream a reported TSS.
		#fields	type	1	String	type of element
		#fields	counts	1	Integer	Number of reads mapping to isoform
		#fields	tpm	1	Float	Transcripts per million (number of reads mapping nomralized to 1m reads total)
	}
}

proc iso_flair_job {args} {
	upvar job_logdir job_logdir
	flair_job {*}$args
}

proc flair_job {args} {
	upvar job_logdir job_logdir
	global appdir
	# putslog [list flair_job {*}$args]
	set cmdline [clean_cmdline cg flair {*}$args]
	set refseq {}
	set skips {}
	set plotgenes {}
	set sqanti 1
	set compar multitranscript
	set reftranscripts {}
	set threads 8
	cg_options flair args {
		-refseq {
			set refseq $value
		}
		-plotgenes {
			set plotgenes $value
		}
		-reftranscripts {
			set reftranscripts $value
		}
		-compar {
			if {$value ni "joint multitranscript"} {
				error "unknown value $value for flair -compar, must be one of: 0 joint multitranscript"
			}
			set compar $value
		}
		-distrreg {
			# this option is not actually supported by flair, 
			# but present for compatibilty with generic call from process_*
		}
		-threads {
			set threads $value
		}
		-skip {
			lappend skips -skip $value
		}
		-cleanup {
			# currently unused
		}
	} {projectdir}
	set projectdir [file_absolute $projectdir]
	if {[file isdir $projectdir]} {
		set sampledirs [glob -nocomplain $projectdir/samples/*]
		if {[llength $sampledirs] == 0} {
			set sampledirs [list $projectdir]
			set compar 0
		}
		# select bam later
		set bam {}
	} else {
		# not a dir, so should be a bamfile
		set bam $projectdir
		set sampledirs [list [file dir $bam]]
		set projectdir [file dir $bam]
		set compar 0
	}
	# cd $projectdir
	set refseq [refseq $refseq]
	if {$reftranscripts eq ""} {
		set reftranscripts [ref_gtftranscripts $refseq]
	} else {
		set reftranscripts [file_absolute $reftranscripts]
	}
	set flairdir [findflair]
	job_logfile $projectdir/flair_[file tail $projectdir] $projectdir $cmdline \
		{*}[versions flair dbdir zstd os]
	# analysis per sample
	set allseq_fasqfiles {}
	foreach sampledir $sampledirs {
		putsvars sampledir
		cd $sampledir
		set sample [file tail $sampledir]
		if {$bam eq ""} {
			set bams [jobglob $sampledir/map-*.bam]
		} else {
			set bams [list $bam]
		}
		if {![llength $bams]} continue
		foreach bam $bams {
			set rootname [file_rootname $bam]
			mkdir flair-$rootname
			job flair_correct-$rootname {*}$skips -skip flair-$rootname/counts_matrix-flair-$rootname.tsv \
			-cores $threads \
			-deps {
				$bam $bam.bai $refseq $reftranscripts
			} -targets {
				flair-$rootname/all_corrected-flair-$rootname.bed
			} -vars {
				bam rootname flairdir refseq reftranscripts threads sample
			} -code {
				analysisinfo_write $bam $target \
					analysis flair-$rootname sample $sample \
					isocaller_reftranscripts [file tail $reftranscripts] \
					isocaller_distrreg 0 \
					isocaller flair isocaller_version [version flair]
				set bed12 [file root $bam].bed12
				exec [flair_bin bam2Bed12] -i $bam > $bed12.temp
				file rename -force $bed12.temp $bed12
				catch_exec [flair_bin flair] correct -t $threads \
					-g $refseq \
					--gtf $reftranscripts \
					-q $bed12 \
					-o flair-$rootname/transcripts-flair-$rootname.temp >@ stdout 2>@ stderr
				file rename flair-$rootname/transcripts-flair-$rootname.temp_all_corrected.bed flair-$rootname/all_corrected-flair-$rootname.bed
				catch {file rename flair-$rootname/transcripts-flair-$rootname.temp_all_inconsistent.bed flair-$rootname/all_inconsistent-flair-$rootname.bed}
				file delete $bed12
			}
			set fastqfiles [jobglob fastq/*.fastq.gz]
			set allseq_fasqfile flair-$rootname/allseq-$rootname.fastq.gz
			lappend allseq_fasqfiles [file_absolute $allseq_fasqfile]
			job flair_allseq-$rootname {*}$skips -skip flair-$rootname/counts_matrix-flair-$rootname.tsv \
			-deps $fastqfiles -targets {
				flair-$rootname/allseq-$rootname.fastq.gz
			} -vars {
				rootname fastqfiles
			} -code {
				# next best on combined data from samples
				set tempfastq flair-$rootname/allseq-$rootname.fastq.gz
				if {![file exists $tempfastq]} {
					puts "Making $tempfastq"
					exec cat {*}$fastqfiles > $tempfastq.temp
					file rename $tempfastq.temp $tempfastq
				}
			}
			job flair_collapse-$rootname {*}$skips -skip flair-$rootname/counts_matrix-flair-$rootname.tsv \
			-cores $threads \
			-deps {
				flair-$rootname/allseq-$rootname.fastq.gz
				flair-$rootname/all_corrected-flair-$rootname.bed
				$refseq $reftranscripts
			} -targets {
				flair-$rootname/transcripts-flair-$rootname.isoforms.gtf
				flair-$rootname/transcripts-flair-$rootname.isoforms.bed
				flair-$rootname/transcripts-flair-$rootname.isoforms.fa
			} -vars {
				rootname refseq reftranscripts threads
			} -code {
				analysisinfo_write flair-$rootname/all_corrected-flair-$rootname.bed $target
				analysisinfo_write flair-$rootname/all_corrected-flair-$rootname.bed flair-$rootname/transcripts-flair-$rootname.isoforms.fa
				puts "collapse -> flair-$rootname-collapse"
				catch_exec [flair_bin flair] collapse \
					-t $threads \
					-g $refseq \
					--gtf $reftranscripts \
					-r flair-$rootname/allseq-$rootname.fastq.gz \
					-q flair-$rootname/all_corrected-flair-$rootname.bed \
					-o flair-$rootname/temptranscripts-flair-$rootname >@ stdout 2>@ stderr
				foreach file [glob flair-$rootname/temptranscripts-flair-$rootname*] {
					file rename -force $file flair-$rootname/[string range [file tail $file] 4 end]
				}
			}
			job flair_quantify-$rootname {*}$skips -cores $threads -deps {
				flair-$rootname/transcripts-flair-$rootname.isoforms.fa
				flair-$rootname/allseq-$rootname.fastq.gz
			} -targets {
				flair-$rootname/counts_matrix-flair-$rootname.tsv
			} -vars {
				rootname sample threads
			} -code {
				analysisinfo_write flair-$rootname/transcripts-flair-$rootname.isoforms.fa $target \
					analysis flair-$rootname sample $sample \
					isocaller flair isocaller_version [version flair]
				set manifestdata {}
				lappend manifestdata [join [list $sample conditionA batch1 flair-$rootname/allseq-$rootname.fastq.gz] \t]
				file_write reads_manifest.tsv [join $manifestdata \n]\n
				puts "quantify -> flair-$rootname/counts_matrix-flair-$rootname.tsv"
				catch_exec [flair_bin flair] quantify \
					--threads $threads \
					-r reads_manifest.tsv \
					-i flair-$rootname/transcripts-flair-$rootname.isoforms.fa \
					-o flair-$rootname/counts_matrix-flair-$rootname.tsv.temp >@ stdout 2>@ stderr
				if {[file exists flair-$rootname/counts_matrix-flair-$rootname.tsv.temp]} {
					# older versions
					file rename -force flair-$rootname/counts_matrix-flair-$rootname.tsv.temp flair-$rootname/counts_matrix-flair-$rootname.tsv
				} else {
					# newer versions
					file rename -force flair-$rootname/counts_matrix-flair-$rootname.tsv.temp.counts.tsv flair-$rootname/counts_matrix-flair-$rootname.tsv
				}
			}
			job flair_sqanti-$rootname {*}$skips -deps {
				flair-$rootname/transcripts-flair-$rootname.isoforms.gtf
				flair-$rootname/counts_matrix-flair-$rootname.tsv
				$reftranscripts
				$refseq
			} -targets {
				isoform_counts-flair-$rootname.tsv
				gene_counts-flair-$rootname.tsv
				totalcounts-flair-$rootname.tsv
				sqanti3-flair-$rootname/sqanti3-flair-${rootname}_classification.txt
				sqanti3-flair-$rootname/sqanti3-flair-${rootname}_junctions.txt
				sqanti3-flair-$rootname/sqanti3-flair-${rootname}_corrected.gtf
				sqanti3-flair-$rootname/sqanti3-flair-${rootname}_corrected.genePred
				sqanti3-flair-$rootname/sqanti3-flair-${rootname}_corrected.genepred.tsv
				sqanti3-flair-$rootname/sqanti3-flair-${rootname}_corrected.fasta
			} -vars {
				rootname sample refseq reftranscripts
			} -code {
				mkdir sqanti3-flair-$rootname
				set extrainfo [list \
					analysis flair-$rootname sample $sample \
					isocaller flair isocaller_version [version flair] \
					isocaller_flair_sqanti3 [version sqanti3_qc.py]
				]
				analysisinfo_write flair-$rootname/transcripts-flair-$rootname.isoforms.gtf \
					isoform_counts-flair-$rootname.tsv \
					{*}$extrainfo
				analysisinfo_write flair-$rootname/transcripts-flair-$rootname.isoforms.gtf \
					gene_counts-flair-$rootname.tsv \
					{*}$extrainfo
				set target isoform_counts-flair-$rootname.tsv
				if {[file size flair-$rootname/transcripts-flair-$rootname.isoforms.gtf] == 0} {
					foreach target $targets {
						file_write $target ""
					}
					file_write isoform_counts-flair-$rootname.tsv [join {chromosome begin end strand exonStarts exonEnds cdsStart cdsEnd transcript gene geneid counts-flair-$sample tpm-flair-$sample} \t]\n
					file_write gene_counts-flair-$rootname.tsv [join [list type gene gene_type chromosome begin end strand nrtranscripts counts-flair-$sample tpm-flair-$sample] \t]\n
					file_write total_counts-flair-$rootname.tsv $sample\n
				} else {
					catch_exec sqanti3_qc.py \
						flair-$rootname/transcripts-flair-$rootname.isoforms.gtf \
						$reftranscripts \
						$refseq \
						-d sqanti3-flair-$rootname \
						-o sqanti3-flair-$rootname \
						--report skip
					cg genepred2tsv \
						sqanti3-flair-$rootname/sqanti3-flair-${rootname}_corrected.genePred \
						sqanti3-flair-$rootname/sqanti3-flair-${rootname}_corrected.genepred.tsv
					# merge results
					cg_flair_mergeresults $target \
						sqanti3-flair-$rootname/sqanti3-flair-${rootname}_classification.txt \
						sqanti3-flair-$rootname/sqanti3-flair-${rootname}_corrected.genepred.tsv \
						flair-$rootname/counts_matrix-flair-$rootname.tsv \
						totalcounts-flair-$rootname.tsv
					cg_flair_genecounts isoform_counts-flair-$rootname.tsv gene_counts-flair-$rootname.tsv
				}
			}
			foreach genename $plotgenes {
				job flair_plotisoforms-$rootname-$genename {*}$skips -deps {
					flair-$rootname/transcripts-flair-$rootname.isoforms.bed 
					flair-$rootname/counts_matrix-flair-$rootname.tsv
				} -targets {
					flair_results/${genename}_isoforms.png
				} -vars {
					rootname
				} -code {
					mkdir flair_results
					puts "plot $rootname"
					catch {
						exec [flair_bin plot_isoform_usage] flair-$rootname/transcripts-flair-$rootname.isoforms.bed flair-$rootname/counts_matrix-flair-$rootname.tsv $genename
					}
					file rename ${genename}_isoforms.png flair_results/${genename}_isoforms.png
				}
			}
		}
	}
	if {$compar eq "0"} return
	# combined analysis
	if {$compar eq "joint"} {
		cd $projectdir
		mkdir compar
		set exproot [file tail $projectdir]
		set bedfiles [jobglob samples/*/flair-*/all_corrected-flair-*.bed]
		job flair_compar-$exproot {*}$skips \
		-cores $threads \
		-deps [list_concat $bedfiles $allseq_fasqfiles] \
		-targets {
			compar/flair-$exproot/counts_matrix-flair-$exproot.tsv
			compar/flair-$exproot/transcripts-flair-$exproot.isoforms.tsv
			compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf
			compar/flair-$exproot/transcripts-flair-$exproot.isoforms.bed
		} -vars {
			bedfiles allseq_fasqfiles exproot sample refseq reftranscripts threads
		} -code {
			analysisinfo_write [lindex $bedfiles 0] $target flair [version flair]
			analysisinfo_write [lindex $bedfiles 0] compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf flair [version flair]
			set flairdir [findflair]
			mkdir compar/flair-$exproot
			exec cat {*}$bedfiles > compar/flair-$exproot/all_corrected-flair-$exproot.bed
			# 
			putslog "collapse transcript info -> trancripts-flair-$exproot.*"
			catch_exec [flair_bin flair] collapse \
				-t $threads \
				-g $refseq \
				--gtf $reftranscripts \
				-r [join $allseq_fasqfiles ,] \
				-q compar/flair-$exproot/all_corrected-flair-$exproot.bed \
				-o compar/flair-$exproot/transcripts-flair-$exproot >@ stdout 2>@ stderr
			file delete compar/flair-$exproot/all_corrected-flair-$exproot.bed
			cg gtf2tsv -separate 1 compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf compar/flair-$exproot/transcripts-flair-$exproot.isoforms.tsv
			cg gtf2tsv -separate 0 compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf compar/flair-$exproot/transcripts-flair-$exproot.isoforms.tsv
			#
			# make manifest (in flair.temp)
			mkdir compar/flair-$exproot/flair.temp
			unset -nocomplain manifestdata
			set condition A
			foreach sample [dirglob samples *] {
				set bam [lindex [glob samples/$sample/map-sminimap*.bam map-*.bam] 0]
				set rootname [file_rootname $bam]
				set fastq samples/$sample/allseq-$rootname.fastq.gz
				if {![file exists $fastq]} {
					puts "Making $fastq"
					exec cat {*}[glob samples/$sample/fastq/*.f*q.gz] > $fastq.temp
					file rename $fastq.temp $fastq
				}
				lappend manifestdata($condition) [join [list $sample ${condition} batch1 $fastq] \t]
			}
			set c {}
			foreach condition [array names manifestdata] {
				append c [join $manifestdata($condition) \n]\n
			}
			file_write compar/flair-$exproot/flair.temp/reads_manifest-flair-$exproot.tsv $c
			#
			puts "quantify -> compar/flair-$exproot/counts_matrix-flair-$exproot.tsv"
			catch_exec [flair_bin flair] quantify \
				--tpm \
				-t $threads \
				-r compar/flair-$exproot/flair.temp/reads_manifest-flair-$exproot.tsv \
				-i compar/flair-$exproot/transcripts-flair-$exproot.isoforms.fa \
				-o compar/flair-$exproot/flair.temp/counts_matrix-flair-$exproot.tsv >@ stdout 2>@ stderr
			# remove condition/batch (A_batch1) from fields in header -> just samples
			file delete compar/flair-$exproot/flair.temp/reads_manifest-flair-$exproot.tsv
			file rename -force compar/flair-$exproot/flair.temp/counts_matrix-flair-$exproot.tsv compar/flair-$exproot/flair.temp/counts_matrix-flair-$exproot.tsv.ori
			catch {close $f} ; set f [open compar/flair-$exproot/flair.temp/counts_matrix-flair-$exproot.tsv.ori]
			catch {close $o} ; set o [open compar/flair-$exproot/flair.temp/counts_matrix-flair-$exproot.tsv w]
			set line [split [gets $f] \t]
			puts $o [join [list_regsub -all {_A_batch1$} $line {}] \t]
			fcopy $f $o
			close $o
			close $f
			foreach file [glob compar/flair-$exproot/flair.temp/*] {
				file rename -force $file compar/flair-$exproot/[file tail $file]
			}
			file delete compar/flair-$exproot/flair.temp
			# puts "diffExp on compar/counts_matrix-flair-$exproot.tsv"
			# file delete -force flair-diffexp-$exproot
			# exec flair diffExp -q compar/counts_matrix-flair-$exproot.tsv -o compar/diffexp-flair-$exproot
			# exec python $flairdir/bin/bin/diff_iso_usage.py compar/diffiso-flair-$exproot.tsv
		}
		if {$sqanti} {
			set transcript_classification_file compar/sqanti3-flair-$exproot/sqanti3-flair-${exproot}_classification.txt
			set transcripts_genepred_file compar/sqanti3-flair-$exproot/sqanti3-flair-${exproot}_corrected.genepred.tsv
			set counts_matrix_file compar/flair-$exproot/counts_matrix-flair-$exproot.tsv
			job flair_sqanti_compar-$exproot {*}$skips -deps {
				compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf
				$reftranscripts
				$refseq
				$counts_matrix_file
			} -targets {
				compar/isoform_counts-flair-$exproot.genepred.tsv
				compar/gene_counts-flair-$exproot.genepred.tsv
				compar/totalcounts-flair-$exproot.tsv
				compar/sqanti3-flair-$exproot/sqanti3-flair-${exproot}_classification.txt
				$transcript_classification_file
				$transcripts_genepred_file
				compar/sqanti3-flair-$exproot/sqanti3-flair-${exproot}_classification.txt
				compar/sqanti3-flair-$exproot/sqanti3-flair-${exproot}_junctions.txt
				compar/sqanti3-flair-$exproot/sqanti3-flair-${exproot}_corrected.gtf
				compar/sqanti3-flair-$exproot/sqanti3-flair-${exproot}_corrected.genePred
				compar/sqanti3-flair-$exproot/sqanti3-flair-${exproot}_corrected.genepred.tsv
				compar/sqanti3-flair-$exproot/sqanti3-flair-${exproot}_corrected.fasta
			} -vars {
				transcript_classification_file transcripts_genepred_file counts_matrix_file
				exproot sample refseq reftranscripts
			} -code {
				analysisinfo_write compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf compar/isoform_counts-flair-$exproot.genepred.tsv sqanti3 [version sqanti3_qc.py]
				analysisinfo_write compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf compar/gene_counts-flair-$exproot.genepred.tsv sqanti3 [version sqanti3_qc.py]
				analysisinfo_write compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf compar/sqanti3-flair-$exproot/sqanti3-${exproot}_classification.txt sqanti3 [version sqanti3_qc.py]
				mkdir compar/sqanti3-flair-$exproot
				if {[catch {
					catch_exec sqanti3_qc.py \
						compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf \
						$reftranscripts \
						$refseq \
						-d compar/sqanti3-flair-$exproot \
						-o sqanti3-flair-$exproot \
						--saturation \
						--report pdf
				}]} {
					catch_exec sqanti3_qc.py \
						compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf \
						$reftranscripts \
						$refseq \
						-d compar/sqanti3-flair-$exproot \
						-o sqanti3-flair-$exproot \
						--saturation \
						--report skip
				}
				cg genepred2tsv compar/sqanti3-flair-$exproot/sqanti3-flair-${exproot}_corrected.genePred $transcripts_genepred_file 
				# merge results
				set isoformcounts compar/isoform_counts-flair-$exproot.genepred.tsv
				cg_flair_mergeresults $isoformcounts \
					$transcript_classification_file \
					$transcripts_genepred_file \
					$counts_matrix_file \
					compar/totalcounts-flair-$exproot.tsv
				#
				# make compar/gene_counts-flair-$exproot.genepred.tsv
				set genecounts compar/gene_counts-flair-$exproot.genepred.tsv
				cg_flair_genecounts $isoformcounts $genecounts
			}
		}
	} else {
		cd $projectdir
		mkdir compar
		set exproot [file tail $projectdir]
		set isoformfiles [bsort [jobglob samples/*/isoform_counts-flair-*.tsv]]
		set totalcountsfiles [bsort [jobglob samples/*/totalcounts-flair-*.tsv]]
		job flair_compar-$exproot {*}$skips \
		-cores $threads \
		-deps [list_concat $isoformfiles $totalcountsfiles] \
		-targets {
			compar/isoform_counts-flair-$exproot.genepred.tsv
			compar/gene_counts-flair-$exproot.genepred.tsv
			compar/totalcounts-flair-$exproot.tsv
		} -vars {
			isoformfiles totalcountsfiles exproot
		} -code {
			analysisinfo_write [lindex $isoformfiles 0] $target
			analysisinfo_write [lindex $isoformfiles 0] compar/gene_counts-flair-$exproot.genepred.tsv
			analysisinfo_write [lindex $totalcountsfiles 0] compar/totalcounts-flair-$exproot.tsv
			cg paste {*}[bsort $totalcountsfiles] > compar/totalcounts-flair-$exproot.tsv
			set isoformcounts compar/isoform_counts-flair-$exproot.genepred.tsv
			cg multitranscript -match . $isoformcounts {*}$isoformfiles
			set genecounts compar/gene_counts-flair-$exproot.genepred.tsv
			cg_flair_genecounts $isoformcounts $genecounts
		}
	}
}

proc cg_flair_plot_isoform_usage {args} {
	set cutnames {}
	cg_options cg_flair_plot_isoform_usage args {
		-cutnames {
			if {$cutnames eq "0"} {
				set cutnames {}
			} else {
				set cutnames $value
			}
		}
	} {flairtranscriptsbed flaircounts_matrix genename resultprefix}
	if {$cutnames ne ""} {
		set tempfile [tempfile]
		catch {close $f} ; set f [open $flaircounts_matrix]
		catch {close $o} ; set o [open $tempfile w]
		set line [split [gets $f] \t]
		set newheader [list [lindex $line 0]]
		foreach el [lrange $line 1 end] {
			set temp [split $el _]
			if {[llength $temp] > 3} {
				set temp [lreplace $temp end-2 end-2]
				set el [join [lrange $temp 0 end-1] _]
			}
			lappend newheader $el
		}
		puts $o [join $newheader \t]
		fcopy $f $o
		close $o
		close $f
		set flaircounts_matrix $tempfile
	}
	set flairdir [findflair]
	if {[info exists ::env(FONTCONFIG_PATH)]} {
		append ::env(FONTCONFIG_PATH) :/etc/fonts
	} else {
		set ::env(FONTCONFIG_PATH) /etc/fonts
	}
	exec [flair_bin plot_isoform_usage] $flairtranscriptsbed $flaircounts_matrix $genename $resultprefix
}

proc cg_flair {args} {
	set args [job_init {*}$args]
	flair_job {*}$args
	job_wait
}

proc cg_iso_flair {args} {
	set args [job_init {*}$args]
	flair_job {*}$args
	job_wait
}

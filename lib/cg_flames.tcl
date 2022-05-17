proc version_flames {} {
	set flamesdir [findflames]
	lindex [split [file tail $flamesdir] -] end
}

proc flames_getref {refseq} {
	set refseq [refseq $refseq]
	set refdir [file dir $refseq]
	if {[catch {
		set ref [lindex [bsort [glob $refdir/extra/gencode*.gtf]] end]
	}]} {
		set keeppwd [pwd]
		cd $refdir/hg38/extra
		exec wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
		exec gunzip gencode.v39.annotation.gtf.gz
		cd $keeppwd
	}
	set ref [lindex [bsort [glob $refdir/extra/gencode*.gtf]] end]
}

proc findflames {} {
	global flames
	if {![info exists flames]} {
		set flames [searchpath flames flames flames-*]
		if {$flames eq ""} {
			set flames [searchpath flames flames flames*]
		}
		set ::env(PATH) $flames:$::env(PATH)
	}
	return $flames
}

proc cg_flames_genecounts {isoformcounts genecounts} {
	set tempfile [tempfile]
	cg select -overwrite 1 \
		-g {associated_gene} \
		-gc {distinct(structural_category),distinct(subcategory),sum(counts-*),sum(tpm-*)} \
		$isoformcounts $tempfile
	catch {close $f} ; catch {close $o}
	set f [open $tempfile]
	set header [tsv_open $f]
	set oheader {gene gene_type}
	foreach field [lrange $header 3 end] {
		regsub sum_ $field {} field
		lappend oheader $field
	}
	set o [open $genecounts w]
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
		puts $o $gene\t$genetype\t[join [lrange $line 3 end] \t]
	}
	close $o
	close $f
}

proc cg_flames_mergeresults {target transcript_classification_file transcripts_genepred_file counts_matrix_file {totalcountsfile {}}} {
	#
	# combined
	# adapted tempfiles
	set tempclassification [tempfile]
	set tempgenepred [tempfile]
	set tempcount [tempfile]
	# class
	set tempfile [tempfile]
	cg select -overwrite 1 -f {id="a$isoform" *} $transcript_classification_file $tempfile
	cg select -s id $tempfile $tempclassification
	# genepred
	cg select -overwrite 1 -f {id="a$name" *} $transcripts_genepred_file $tempfile
	cg select -s id $tempfile $tempgenepred
	# remake ids in countfile
	set tempfile [tempfile]
	catch {close $f} ; catch {close $o}
	set f [gzopen $counts_matrix_file]
	set o [open $tempfile w]
	set cntheader [tsv_open $f]
	set cntheader [list_regsub {_conditionA_batch1$} $cntheader {}]
	puts $o id\t[join $cntheader \t]
	set numsamples [expr {[llength $cntheader] -1}]
	set totalcounts [list_fill $numsamples 0]
	while {[gets $f line] != -1} {
		set line [split $line \t]
		set cntid [lindex $line 0]
		set cntid [lindex [split $cntid _] 0]
		puts $o a$cntid\t[join $line \t]
		set temp {}
		foreach tcnt $totalcounts cnt [lrange $line 1 end] {
			set tcnt [expr {$tcnt + $cnt}]
			lappend temp $tcnt
		}
		set totalcounts $temp
	}
	close $o
	close $f
	cg select -s id $tempfile $tempcount
	set exproot [file_rootname $target]
	regsub \\.genepred $exproot {} exproot
	if {$totalcountsfile eq ""} {
		set totalcountsfile [file dir $target]/totalcounts-$exproot.tsv
	}
	file_write $totalcountsfile [join [lrange $cntheader 1 end] \t]\n[join $totalcounts \t]\n
	#
	# write combi file
 	catch {close $fg} ; catch {close $fclass} ; catch {close $fcnt} ; catch {close $o} ; 
	set fg [gzopen $tempgenepred]
	set gheader [lrange [tsv_open $fg] 1 end]
 	set fclass [gzopen $tempclassification]
	set classheader [lrange [tsv_open $fclass] 1 end]
 	set fcnt [gzopen $tempcount]
	set cntheader [lrange [tsv_open $fcnt] 1 end]
	set o [open $target.temp w]
	set newheader $gheader
	set pclassheader [list_remove $classheader id isoform chrom strand]
	set poss [list_cor $classheader $pclassheader]
	lappend newheader {*}$pclassheader
	foreach sample [lrange $cntheader 1 end] {
		lappend newheader counts-$sample
	}
	foreach sample [lrange $cntheader 1 end] {
		lappend newheader tpm-$sample
	}
	puts $o [join $newheader \t]
	set emptycntresult [list_fill [expr {2*[llength $totalcounts]}] 0]
	set nr 1
	if {[gets $fcnt cntline] == -1} break
	set cntline [split $cntline \t]
	set cntid [list_shift cntline]
	set gr 0
	while 1 {
		if {[gets $fg gline] == -1} break
		set gline [split $gline \t]
		set gid [list_shift gline]
		if {[gets $fclass classline] == -1} break
		set classline [split $classline \t]
		set classid [list_shift classline]
		if {$classid ne $gid} {
			error "difference in ids at line $nr between $transcripts_genepred_file and $transcript_classification_file: $gid vs $classid"
		}
		set result [join $gline \t]
		append result \t[join [list_sub $classline $poss] \t]
		if {$cntid eq $gid} {
			append result \t[join [lrange $cntline 1 end] \t]
			foreach cnt [lrange $cntline 1 end] tcnt $totalcounts {
				append result \t[expr {1000000.0*$cnt/$tcnt}]
			}
			puts $o $result
			set gr [gets $fcnt cntline]
			set cntline [split $cntline \t]
			set cntid [list_shift cntline]
		} else {
			puts stderr "count not found for $gid at line $nr ($transcripts_genepred_file and $counts_matrix_file)"
			append result \t[join $emptycntresult \t]
			puts $o $result
		}
		incr nr
	}
	close $o
	close $fg ; close $fclass ; close $fcnt
	# write final target
	file_write $target.temp2 [deindent {
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
		#fields	structural_category	1	String	one of the isoform categories absed on the best matching reference transcript (https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-output-explanation)
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
		#fields	counts	1	Integer	Number of reads mapping to isoform
		#fields	tpm	1	Float	Transcripts per million (number of reads mapping nomralized to 1m reads total)
	}]\n
	cg select -s {chromosome begin end exonStarts exonEnds strand} $target.temp >> $target.temp2
	file rename -force $target.temp2 $target
	file delete -force $target.temp
}

proc flames_job {args} {
	# putslog [list flames_job {*}$args]
	upvar job_logdir job_logdir
	set cmdline "[list cd [pwd]] \; [list cg flames {*}$args]"
	global appdir
	set refseq {}
	set skips {}
	set genes {}
	set extraopts {}
	set sqanti 1
	set compar joint
	set threads 8
	set resultfile {}
	set reftranscripts {}
	cg_options flames args {
		-refseq {
			set refseq $value
		}
		-compar {
			set compar $value
		}
		-threads {
			set threads $value
		}
		-reftranscripts {
			set reftranscripts $value
		}
		-skip {
			lappend skips -skip $value
		}
	} {fastqdir resultfile} 1 2
	set fastqdir [file_absolute $fastqdir]
	if {$resultfile eq ""} {
		set root [file_rootname [file dir $fastqdir]]
		set resultfile [file dir $fastqdir]/isoform_counts-sqanti3-flames-$root.genepred.tsv
	} else {
		set resultfile [file_absolute $resultfile]
		set root [file_rootname $resultfile]
	}
	set resultfile [file_absolute $resultfile]
	set resulttail [file tail $resultfile]
	set destdir [file dir $resultfile]
	#
	set refseq [refseq $refseq]
	# 
	job_logfile $destdir/flames_[file tail $destdir] $destdir $cmdline \
		{*}[versions flames dbdir zstd os]
	# analysis per sample
	set flamesdir $destdir/flames-$root
	mkdir $flamesdir
	job flames-[file tail $resultfile] {*}$skips -skip flames-$rootname/counts_matrix-flames-$rootname.tsv \
	-cores $threads \
	-deps {
		$fastqdir $refseq $gtfannotation
	} -targets {
		$resultfile
	} -vars {
		fastqdir resultfile reftranscripts root destdir flamesdir refseq gtfannotation threads
	} -code {
		analysisinfo_write [gzfile $fastqdir/*.fastq $fastqdir/*.fq] $resultfile flames [version flames]
		set keeppwd [pwd]
		cd $flamesdir
		file_write config.json [deindent {
			{
			    "comment":"this is the default config for SIRV spike-in data. use splice annotation on alignment.",
			    "pipeline_parameters":{
			        "do_genome_alignment":true,
			        "do_isoform_identification":true,
			        "do_read_realignment":true,
			        "do_transcript_quantification":true
			    },
			    "global_parameters":{
			        "generate_raw_isoform":true,
			        "has_UMI":false
			    },
			    "isoform_parameters":{
			        "MAX_DIST":10,
			        "MAX_TS_DIST":100,
			        "MAX_SPLICE_MATCH_DIST":10,
			        "min_fl_exon_len":40,
			        "Max_site_per_splice":3,
			        "Min_sup_cnt":10,
			        "Min_cnt_pct":0.01,
			        "Min_sup_pct":0.2,
			        "strand_specific":0,
			        "remove_incomp_reads":5
			    },
			    "alignment_parameters":{
			        "use_junctions":true,
			        "no_flank":true
			    },
			    "realign_parameters":{
			        "use_annotation":true
			    },
			    "transcript_counting":{
			        "min_tr_coverage":0.75,
			        "min_read_coverage":0.75
			    }
			}
		}]
		mklink $refseq [file tail $refseq]
		mklink $bamfile [file tail $bamfile]
		mklink $bamfile.bai [file tail $bamfile].bai
		mklink $reftranscripts [file tail $reftranscripts]
		mkdir fastq
		set fastqfiles [bsort [gzfiles $fastqdir/*.fastq $fastqdir/*.fq]]
		set mergedfastqfile fastq/$root.fastq
		if {[llength $fastqfiles] > 1} {
			exec cg zcat {*}$fastqfiles | gzip --fast > $mergedfastqfile.gz.temp
			file rename $mergedfastqfile.gz.temp $mergedfastqfile.gz
		} else {
			set fastqfile [lindex $fastqfiles 0]
			mklink -absolute 0 $fastqfile $mergedfastqfile[gzext $fastqfile]
		}
		exec flames bulk_long_pipeline.py {*}$extraopts \
			--genomefa [file tail $refseq] \
			--gff3 [file tail $reftranscripts] \
			--outdir tempflames \
			--config_file config.json \
			--fq_dir fastq \
			>@ stdout 2>@ stderr
		foreach file [glob tempflames/*] {
			file rename -force $file .
		}
		file delete tempflames
		#
		cd $destdir
		mklink -absolute 0 flames-$root/isoform_annotated.filtered.gff3 transcripts-flames-$root.isoforms.gff3
		mklink -absolute 0 flames-$root/transcript_assembly.fa transcripts-flames-$root.isoforms.fa
		set f [gzopen flames-$root/transcript_count.csv.gz]
		set header [split [gets $f] ,]
		set newheader {transcript_id gene_id}
		foreach field 
		mklink -absolute 0 flames-$root/
		mklink -absolute 0 flames-$root/
	}
	job flames_sqanti-$rootname {*}$skips -deps {
		flames-$rootname/transcripts-flames-$rootname.isoforms.gtf
		flames-$rootname/counts_matrix-flames-$rootname.tsv
		$gtfannotation
		$refseq
	} -targets {
		isoform_counts-sqanti3-flames-$rootname.genepred.tsv
		totalcounts-sqanti3-flames-$rootname.tsv
		sqanti3-flames-$rootname/sqanti3-flames-${rootname}_classification.txt
		sqanti3-flames-$rootname/sqanti3-flames-${rootname}_junctions.txt
		sqanti3-flames-$rootname/sqanti3-flames-${rootname}_corrected.gtf
		sqanti3-flames-$rootname/sqanti3-flames-${rootname}_corrected.genePred
		sqanti3-flames-$rootname/sqanti3-flames-${rootname}_corrected.genepred.tsv
		sqanti3-flames-$rootname/sqanti3-flames-${rootname}_corrected.fasta
	} -vars {
		rootname sample refseq gtfannotation gtfannotation
	} -code {
		mkdir sqanti3-flames-$rootname
		analysisinfo_write flames-$rootname/transcripts-flames-$rootname.isoforms.gtf isoform_counts-sqanti3-flames-$rootname.genepred.tsv sqanti3 [version sqanti3_qc.py]
		analysisinfo_write flames-$rootname/transcripts-flames-$rootname.isoforms.gtf sqanti3-flames-$rootname/sqanti3-${rootname}_classification.txt sqanti3 [version sqanti3_qc.py]
		catch_exec sqanti3_qc.py \
			flames-$rootname/transcripts-flames-$rootname.isoforms.gtf \
			$gtfannotation \
			$refseq \
			-d sqanti3-flames-$rootname \
			-o sqanti3-flames-$rootname \
			--report skip
		cg genepred2tsv \
			sqanti3-flames-$rootname/sqanti3-flames-${rootname}_corrected.genePred \
			sqanti3-flames-$rootname/sqanti3-flames-${rootname}_corrected.genepred.tsv
		# merge results
		set target isoform_counts-sqanti3-flames-$rootname.genepred.tsv
		cg_flames_mergeresults $target \
			sqanti3-flames-$rootname/sqanti3-flames-${rootname}_classification.txt \
			sqanti3-flames-$rootname/sqanti3-flames-${rootname}_corrected.genepred.tsv \
			flames-$rootname/counts_matrix-flames-$rootname.tsv \
			totalcounts-sqanti3-flames-$rootname.tsv
	}
}

if 0 {
	if {$compar eq "0"} return
	# combined analysis
	if {$compar eq "joint"} {
		cd $projectdir
		mkdir compar
		set exproot [file tail $projectdir]
		set bedfiles [jobglob samples/*/all_corrected-flames-*.bed]
		set fastqfiles [jobglob samples/*/fastq/allseq-*.fastq.gz]
		job flames_compar-$exproot {*}$skips \
		-cores $threads \
		-deps [list_concat $bedfiles $allseq_fasqfiles] \
		-targets {
			compar/flames-$exproot/counts_matrix-flames-$exproot.tsv
			compar/flames-$exproot/transcripts-flames-$exproot.isoforms.genepred.tsv
			compar/flames-$exproot/transcripts-flames-$exproot.isoforms.gtf
			compar/flames-$exproot/transcripts-flames-$exproot.isoforms.bed
		} -vars {
			bedfiles allseq_fasqfiles exproot sample refseq gtfannotation threads
		} -code {
			analysisinfo_write [lindex $bedfiles 0] $target flames [version flames]
			analysisinfo_write [lindex $bedfiles 0] compar/flames-$exproot/transcripts-flames-$exproot.isoforms.gtf flames [version flames]
			set flamesdir [findflames]
			mkdir compar/flames-$exproot
			exec cat {*}$bedfiles > compar/flames-$exproot/all_corrected-flames-$exproot.bed
			# 
			putslog "collapse transcript info -> trancripts-flames-$exproot.*"
			catch_exec flames.py collapse \
				-t $threads \
				-g $refseq \
				--gtf $gtfannotation \
				-r [join $allseq_fasqfiles ,] \
				-q compar/flames-$exproot/all_corrected-flames-$exproot.bed \
				-o compar/flames-$exproot/transcripts-flames-$exproot >@ stdout 2>@ stderr
			file delete compar/flames-$exproot/all_corrected-flames-$exproot.bed
			cg gtf2tsv -separate 1 compar/flames-$exproot/transcripts-flames-$exproot.isoforms.gtf compar/flames-$exproot/transcripts-flames-$exproot.isoforms.tsv
			cg gtf2tsv -separate 0 compar/flames-$exproot/transcripts-flames-$exproot.isoforms.gtf compar/flames-$exproot/transcripts-flames-$exproot.isoforms.genepred.tsv
			#
			# make manifest (in flames.temp)
			mkdir compar/flames-$exproot/flames.temp
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
			file_write compar/flames-$exproot/flames.temp/reads_manifest-flames-$exproot.tsv $c
			#
			puts "quantify -> compar/flames-$exproot/counts_matrix-flames-$exproot.tsv"
			catch_exec flames.py quantify \
				--tpm \
				-t $threads \
				-r compar/flames-$exproot/flames.temp/reads_manifest-flames-$exproot.tsv \
				-i compar/flames-$exproot/transcripts-flames-$exproot.isoforms.fa \
				-o compar/flames-$exproot/flames.temp/counts_matrix-flames-$exproot.tsv >@ stdout 2>@ stderr
			# remove condition/batch (A_batch1) from fields in header -> just samples
			file delete compar/flames-$exproot/flames.temp/reads_manifest-flames-$exproot.tsv
			file rename -force compar/flames-$exproot/flames.temp/counts_matrix-flames-$exproot.tsv compar/flames-$exproot/flames.temp/counts_matrix-flames-$exproot.tsv.ori
			catch {close $f} ; set f [open compar/flames-$exproot/flames.temp/counts_matrix-flames-$exproot.tsv.ori]
			catch {close $o} ; set o [open compar/flames-$exproot/flames.temp/counts_matrix-flames-$exproot.tsv w]
			set line [split [gets $f] \t]
			puts $o [join [list_regsub -all {_A_batch1$} $line {}] \t]
			fcopy $f $o
			close $o
			close $f
			foreach file [glob compar/flames-$exproot/flames.temp/*] {
				file rename -force $file compar/flames-$exproot/[file tail $file]
			}
			file delete compar/flames-$exproot/flames.temp
			# puts "diffExp on compar/counts_matrix-flames-$exproot.tsv"
			# file delete -force flames-diffexp-$exproot
			# exec flames.py diffExp -q compar/counts_matrix-flames-$exproot.tsv -o compar/diffexp-flames-$exproot
			# exec python $flamesdir/bin/bin/diff_iso_usage.py compar/diffiso-flames-$exproot.tsv
		}
		if {$sqanti} {
			set transcript_classification_file compar/sqanti3-flames-$exproot/sqanti3-flames-${exproot}_classification.txt
			set transcripts_genepred_file compar/sqanti3-flames-$exproot/sqanti3-flames-${exproot}_corrected.genepred.tsv
			set counts_matrix_file compar/flames-$exproot/counts_matrix-flames-$exproot.tsv
			job flames_sqanti_compar-$exproot {*}$skips -deps {
				compar/flames-$exproot/transcripts-flames-$exproot.isoforms.gtf
				$gtfannotation
				$refseq
				$counts_matrix_file
			} -targets {
				compar/isoform_counts-sqanti3-flames-$exproot.genepred.tsv
				compar/gene_counts-sqanti3-flames-$exproot.genepred.tsv
				compar/totalcounts-sqanti3-flames-$exproot.tsv
				compar/sqanti3-flames-$exproot/sqanti3-flames-${exproot}_classification.txt
				$transcript_classification_file
				$transcripts_genepred_file
				compar/sqanti3-flames-$exproot/sqanti3-flames-${exproot}_classification.txt
				compar/sqanti3-flames-$exproot/sqanti3-flames-${exproot}_junctions.txt
				compar/sqanti3-flames-$exproot/sqanti3-flames-${exproot}_corrected.gtf
				compar/sqanti3-flames-$exproot/sqanti3-flames-${exproot}_corrected.genePred
				compar/sqanti3-flames-$exproot/sqanti3-flames-${exproot}_corrected.genepred.tsv
				compar/sqanti3-flames-$exproot/sqanti3-flames-${exproot}_corrected.fasta
				compar/sqanti3-flames-$exproot/sqanti3-flames-${exproot}_SQANTI3_report.pdf
			} -vars {
				transcript_classification_file transcripts_genepred_file counts_matrix_file
				exproot sample refseq gtfannotation gtfannotation
			} -code {
				analysisinfo_write compar/flames-$exproot/transcripts-flames-$exproot.isoforms.gtf compar/isoform_counts-sqanti3-flames-$exproot.genepred.tsv sqanti3 [version sqanti3_qc.py]
				analysisinfo_write compar/flames-$exproot/transcripts-flames-$exproot.isoforms.gtf compar/gene_counts-sqanti3-flames-$exproot.genepred.tsv sqanti3 [version sqanti3_qc.py]
				analysisinfo_write compar/flames-$exproot/transcripts-flames-$exproot.isoforms.gtf compar/sqanti3-flames-$exproot/sqanti3-${exproot}_classification.txt sqanti3 [version sqanti3_qc.py]
				mkdir compar/sqanti3-flames-$exproot
				catch_exec sqanti3_qc.py \
					compar/flames-$exproot/transcripts-flames-$exproot.isoforms.gtf \
					$gtfannotation \
					$refseq \
					-d compar/sqanti3-flames-$exproot \
					-o sqanti3-flames-$exproot \
					--saturation \
					--report pdf
				cg genepred2tsv compar/sqanti3-flames-$exproot/sqanti3-flames-${exproot}_corrected.genePred $transcripts_genepred_file 
				# merge results
				set isoformcounts compar/isoform_counts-sqanti3-flames-$exproot.genepred.tsv
				cg_flames_mergeresults $isoformcounts \
					$transcript_classification_file \
					$transcripts_genepred_file \
					$counts_matrix_file
				#
				# make compar/gene_counts-sqanti3-flames-$exproot.genepred.tsv
				set genecounts compar/gene_counts-sqanti3-flames-$exproot.genepred.tsv
				cg_flames_genecounts $isoformcounts $genecounts
			}
		}
	} else {
		cd $projectdir
		mkdir compar
		set exproot [file tail $projectdir]
		set isoformfiles [bsort [jobglob samples/*/isoform_counts-sqanti3-flames-*.tsv]]
		set totalcountsfiles [bsort [jobglob samples/*/totalcounts-sqanti3-flames-*.tsv]]
		job flames_compar-$exproot {*}$skips \
		-cores $threads \
		-deps [list_concat $isoformfiles $totalcountsfiles] \
		-targets {
			compar/isoform_counts-sqanti3-flames-$exproot.genepred.tsv
			compar/gene_counts-sqanti3-flames-$exproot.genepred.tsv
			compar/totalcounts-sqanti3-flames-$exproot.tsv
		} -vars {
			isoformfiles totalcountsfiles exproot
		} -code {
			analysisinfo_write [lindex $isoformfiles 0] $target
			analysisinfo_write [lindex $isoformfiles 0] compar/gene_counts-sqanti3-flames-$exproot.genepred.tsv
			analysisinfo_write [lindex $totalcountsfiles 0] compar/totalcounts-sqanti3-flames-$exproot.tsv
			cg paste {*}[bsort $totalcountsfiles] > compar/totalcounts-sqanti3-flames-$exproot.tsv
			set isoformcounts compar/isoform_counts-sqanti3-flames-$exproot.genepred.tsv
			cg multitranscript $isoformcounts {*}$isoformfiles
			set genecounts compar/gene_counts-sqanti3-flames-$exproot.genepred.tsv
			cg_flames_genecounts $isoformcounts $genecounts
		}
	}
}

proc cg_flames {args} {
	set args [job_init {*}$args]
	flames_job {*}$args
	job_wait
}

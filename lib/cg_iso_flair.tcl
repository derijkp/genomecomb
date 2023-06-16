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
	if {![catch {exec which $bin} msg]} {
		return $bin
	}
	if {$bin in "flair"} {
		return $bin.py
	}
	set dir [file dir [file_resolve [exec which flair.py]]]
	regsub {\.py$} $bin {} bin2
	if {[file exists $dir/$bin]} {
		return $dir/$bin
	} elseif {[file exists $dir/$bin.py]} {
		return $dir/$bin.py
	} elseif {[info exists bin2] && [file exists $dir/$bin2]} {
		return $dir/$bin2
	} elseif {[file exists $dir/bin/$bin]} {
		return $dir/bin/$bin
	} elseif {[file exists $dir/bin/bin/$bin]} {
		return $dir/bin/bin/$bin
	} elseif {[info exists bin2] && [file exists $dir/bin/$bin2]} {
		return $dir/bin/$bin2
	} elseif {[info exists bin2] && [file exists $dir/bin/bin/$bin2]} {
		return $dir/bin/bin/$bin2
	} elseif {[info exists bin2] && [file exists $dir/bin/bin/$bin2.py]} {
		return $dir/bin/bin/$bin2.py
	}
	error "flair binary $bin no found (looking in $dir)"
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

proc cg_flair_genecounts {isoformcounts genecounts {genefield gene}} {
	set tempfile [tempfile]
	cg select -overwrite 1 \
		-g $genefield \
		-gc {distinct(chromosome),min(begin),max(end),distinct(strand),ucount(transcript),sum(counts-*)} \
		$isoformcounts $tempfile
	catch {close $f} ; catch {close $o}
	set f [open $tempfile]
	set header [tsv_open $f]
	set oheader {type gene gene_type chromosome begin end strand nrtranscripts}
	foreach field [lrange $header 6 end] {
		regsub sum_ $field {} field
		lappend oheader $field
	}
	set o [open $genecounts.temp w]
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
	while 1 {
		if {[gets $f line] == -1} break
		set line [split $line \t]
		foreach {gene} $line break
		set novel [regexp {^novel} $gene]
		if {[regexp {^novel} $gene]} {
			set genetype novel
		} else {
			set genetype known
		}
		puts $o gene\t$gene\t$genetype\t[join [lrange $line 1 end] \t]
	}
	close $o
	close $f
	cg select -s - $genecounts.temp $genecounts.temp2
	file rename -force $genecounts.temp2 $genecounts
	file delete $genecounts.temp
}

proc cg_flair_mergeresults {rootname tsvreftranscripts flairtranscripts flaircountmatrix out_isoform_counts_file out_gene_counts_file out_total_counts_file} {
	cg gtf2tsv $flairtranscripts [file root $flairtranscripts].tsv
	#
	set basicfields {transcript gene geneid chromosome strand begin end exonStarts exonEnds cdsStart cdsEnd exonCount}
	catch {gzclose $fr} ; set fr [gzopen $tsvreftranscripts]
	set frheader [tsv_open $fr]
	set rposs [list_cor $frheader $basicfields]
	unset -nocomplain refa
	unset -nocomplain refgenea
	while 1 {
		if {[gets $fr line] == -1} break
		set line [list_sub [split $line \t] $rposs]
		set refa([lindex $line 0]) $line
		set refgenea([lindex $line 1]) 1
	}
	gzclose $fr
	catch {gzclose $fc} ; set fc [gzopen $flaircountmatrix]
	set fcheader [tsv_open $fc]
	unset -nocomplain countsa
	set totalcount 0
	while 1 {
		if {[gets $fc line] == -1} break
		set line [split $line \t]
		foreach {transcript count} $line break
		set countsa($transcript) $count
		set totalcount [expr {$totalcount + $count}]
	}
	gzclose $fc
	#
	file_write $out_total_counts_file totalcount-flair-$rootname\n$totalcount\n
	#
	set tempfile [tempfile]
	set o [open $tempfile w]
	puts $o [deindent {
		#filetype	tsv/transcriptsfile
		#fileversion	0.99
		#fields	table
		#fields	field	number	type	description
		#fields	chromosome	1	String	Chromosome name
		#fields	begin	1	Integer	Transcription start position
		#fields	end	1	Integer	Transcription end position
		#fields	type	1	Integer	Type of element (transcript typically)
		#fields	transcript	1	String	Name of transcript (usually transcript_id from GTF)
		#fields	gene	1	String	Alternate name / name of gene (e.g. gene_name or gene_id from GTF)
		#fields	strand	1	String	+ or - for strand
		#fields	cdsStart	1	Integer	Coding region start
		#fields	cdsEnd	1	Integer	Coding region end
		#fields	exonCount	1	Integer	Number of exons
		#fields	exonStarts	E	Integer	Exon start positions
		#fields	exonEnds	E	Integer	Exon end positions
		#fields	source	1	String	Source of data
		#fields	count	1	Integer	number of supporting reads for isoform
	}]
	puts $o [join [list {*}$basicfields type transcripttype counts-flair-$rootname] \t]
	catch {gzclose $ff} ; set ff [gzopen [file root $flairtranscripts].tsv]
	set ffheader [tsv_open $ff]
	set fposs [list_cor $ffheader $basicfields]
	while 1 {
		if {[gets $ff line] == -1} break
		set line [list_sub [split $line \t] $fposs]
		foreach {transcript gene geineid} $line break
		set len [string length $gene]
		# if {[string range $transcript end-$len end] eq "_$gene"} {
		# 	set transcript [string range $transcript 0 end-[expr {$len+1}]]
		# }
		if {[info exists refa($transcript)]} {
			set line $refa($transcript)
			lappend line transcript known
		} else {
			foreach {temp gene geneid chromosome strand begin end exonStarts exonEnds} $line break
			lset line 0 [iso_name $chromosome $strand $exonStarts $exonEnds size]
			if {![info exists refgenea($gene)]} {
				lset line 1 novelg_$gene
			}
			lappend line transcript novel
		}
		lappend line [get countsa(${transcript}_$gene) 0]
		puts $o [join $line \t]
	}
	gzclose $o
	cg select -overwrite 1 -s - $tempfile {*}[compresspipe $out_isoform_counts_file] > $out_isoform_counts_file.temp
	file rename -force $out_isoform_counts_file.temp $out_isoform_counts_file
	cg_flair_genecounts $out_isoform_counts_file $out_gene_counts_file
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
	ref_transcripts_convert $reftranscripts tsvreftranscripts gtfreftranscripts
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
				$bam $bam.bai $refseq $gtfreftranscripts
			} -targets {
				flair-$rootname/all_corrected-flair-$rootname.bed
			} -vars {
				bam rootname flairdir refseq gtfreftranscripts threads sample
			} -code {
				analysisinfo_write $bam $target \
					analysis flair-$rootname sample $sample \
					isocaller_reftranscripts [file tail $gtfreftranscripts] \
					isocaller_distrreg 0 \
					isocaller flair isocaller_version [version flair]
				set bed12 [file root $bam].bed12
				exec [flair_bin bam2Bed12] -i $bam > $bed12.temp
				file rename -force $bed12.temp $bed12
				catch_exec [flair_bin flair] correct -t $threads \
					-g $refseq \
					--gtf $gtfreftranscripts \
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
				$refseq $gtfreftranscripts
			} -targets {
				flair-$rootname/transcripts-flair-$rootname.isoforms.gtf
				flair-$rootname/transcripts-flair-$rootname.isoforms.bed
				flair-$rootname/transcripts-flair-$rootname.isoforms.fa
			} -vars {
				rootname refseq gtfreftranscripts threads
			} -code {
				analysisinfo_write flair-$rootname/all_corrected-flair-$rootname.bed $target
				analysisinfo_write flair-$rootname/all_corrected-flair-$rootname.bed flair-$rootname/transcripts-flair-$rootname.isoforms.fa
				puts "collapse -> flair-$rootname-collapse"
				exec [flair_bin identify_gene_isoform.py] \
					flair-$rootname/all_corrected-flair-$rootname.bed \
					$gtfreftranscripts \
					flair-$rootname/all_corrected_renamed-flair-$rootname.bed
				catch_exec [flair_bin flair] collapse \
					-t $threads \
					-g $refseq \
					--gtf $gtfreftranscripts \
					-r flair-$rootname/allseq-$rootname.fastq.gz \
					-q flair-$rootname/all_corrected_renamed-flair-$rootname.bed \
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
			job flair_convert-$rootname {*}$skips -deps {
				flair-$rootname/transcripts-flair-$rootname.isoforms.gtf
				
				$tsvreftranscripts
				$refseq
			} -targets {
				isoform_counts-flair-$rootname.tsv
				gene_counts-flair-$rootname.tsv
				totalcounts-flair-$rootname.tsv
			} -vars {
				rootname sample refseq gtfreftranscripts reftranscripts tsvreftranscripts
			} -code {
				set flairtranscripts flair-$rootname/transcripts-flair-$rootname.isoforms.gtf
				set flaircountmatrix flair-$rootname/counts_matrix-flair-$rootname.tsv
				set out_isoform_counts_file isoform_counts-flair-$rootname.tsv
				set out_gene_counts_file gene_counts-flair-$rootname.tsv
				set out_total_counts_file totalcounts-flair-$rootname.tsv
				set extrainfo [list \
					analysis flair-$rootname sample $sample \
					isocaller flair isocaller_version [version flair] \
					reftranscripts $reftranscripts \
				]
				analysisinfo_write $flairtranscripts \
					$out_isoform_counts_file \
					{*}$extrainfo
				analysisinfo_write $flairtranscripts \
					$out_gene_counts_file \
					{*}$extrainfo
				set target $out_isoform_counts_file
				if {[file size $flairtranscripts] == 0} {
					foreach target $targets {
						file_write $target ""
					}
					file_write $out_isoform_counts_file [join {chromosome begin end strand exonStarts exonEnds cdsStart cdsEnd transcript gene geneid counts-flair-$sample tpm-flair-$sample} \t]\n
					file_write $out_gene_counts_file [join [list type gene gene_type chromosome begin end strand nrtranscripts counts-flair-$sample tpm-flair-$sample] \t]\n
					file_write total_counts-flair-$rootname.tsv $sample\n
				} else {
					cg_flair_mergeresults $rootname \
						$tsvreftranscripts $flairtranscripts $flaircountmatrix \
						$out_isoform_counts_file $out_gene_counts_file $out_total_counts_file
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
			bedfiles allseq_fasqfiles exproot sample refseq gtfreftranscripts threads
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
				--gtf $gtfreftranscripts \
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
				$gtfreftranscripts
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
				exproot sample refseq gtfreftranscripts
			} -code {
				analysisinfo_write compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf compar/isoform_counts-flair-$exproot.genepred.tsv sqanti3 [version sqanti3_qc.py]
				analysisinfo_write compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf compar/gene_counts-flair-$exproot.genepred.tsv sqanti3 [version sqanti3_qc.py]
				analysisinfo_write compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf compar/sqanti3-flair-$exproot/sqanti3-${exproot}_classification.txt sqanti3 [version sqanti3_qc.py]
				mkdir compar/sqanti3-flair-$exproot
				if {[catch {
					catch_exec sqanti3_qc.py \
						compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf \
						$gtfreftranscripts \
						$refseq \
						-d compar/sqanti3-flair-$exproot \
						-o sqanti3-flair-$exproot \
						--saturation \
						--report pdf
				}]} {
					catch_exec sqanti3_qc.py \
						compar/flair-$exproot/transcripts-flair-$exproot.isoforms.gtf \
						$gtfreftranscripts \
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

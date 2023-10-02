proc version_flair {} {
	if {[catch {catch_exec flair --version} version]} {
		set flairdir [findflair]
		set temp [split [file tail $flairdir] -]
		set version [lindex $temp 1]
	}
	set result [lindex $version end]
	regsub ^v $result {} result
	return $result
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
		if {![catch {exec which flair} temp]} {
			set temp [file_resolve $temp]
			set flair [file dir $temp]
		} else {
			set flair [searchpath flair flair flair*]
			if {$flair eq ""} {
				set flair [searchpath FLAIR flair flair*]
			}
			set flair [file_resolve $flair]
			if {![file isdir $flair]} {set flair [file dir $flair]}
		}
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
	set samples [lrange $fcheader 1 end]
	if {[llength $samples] > 1 && $rootname ne ""} {
		error "more than one sample when rootname given"
	}
	foreach sample $samples {
		set totalcounta($sample) 0
	}
	unset -nocomplain countsa
	while 1 {
		if {[gets $fc line] == -1} break
		set line [split $line \t]
		set transcript [lindex $line 0]
		foreach {count} [lrange $line 1 end] sample $samples {
			set countsa($transcript,$sample) $count
			set totalcounta($sample) [expr {$totalcounta($sample) + $count}]
		}
	}
	gzclose $fc
	#
	set o [open $out_total_counts_file w]
	set newheader {}
	set result {}
	if {$rootname ne ""} {
		lappend newheader totalcount-$rootname
		lappend result $totalcounta($sample)
	} else {
		foreach sample $samples {
			lappend newheader totalcount-$sample
			lappend result $totalcounta($sample)
		}
	}
	puts $o [join $newheader \t]
	puts $o [join $result \t]
	gzclose $o
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
	set newheader [list {*}$basicfields type transcripttype]
	if {$rootname ne ""} {
		lappend newheader counts-$rootname
	} else {
		foreach sample $samples {
			lappend newheader counts-$sample
		}
	}
	puts $o [join $newheader \t]
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
		foreach sample $samples {
			lappend line [get countsa(${transcript}_$gene,$sample) 0]
		}
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
	set analysisname flair
	set cmdline [clean_cmdline cg flair {*}$args]
	set refseq {}
	set skips {}
	set plotgenes {}
	set compar multitranscript
	set reftranscripts {}
	set threads 8
	set resultfile {}
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
		-organelles {
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
	} {projectdir resultfile} 1 2
	set projectdir [file_absolute $projectdir]
	if {[file isdir $projectdir]} {
		if {$resultfile ne ""} {
			error "cannot specify resultfile when giving a directory as input"
		}
		set useresultfile 0
		set bams [glob -nocomplain $projectdir/samples/*/*.bam $projectdir/samples/*/*.cram]
		if {[llength $bams] == 0} {
			set bams [glob -nocomplain $projectdir/*.bam $projectdir/*.cram]
			if {[llength $bams] == 0} {
				error "no bams/crams found in $projectdir/samples/*/*.*am or in $projectdir/*.*am"
			}
			set compar 0
		}
	} else {
		# not a dir, so should be a bamfile
		set useresultfile 1
		set bam $projectdir
		if {$resultfile eq ""} {
			set rootname ${analysisname}-[file_rootname $bam]
			set resultfile [file dir $bam]/isoform_counts-$rootname.tsv
		}			
		set resultfile [file_absolute $resultfile]
		set bams [list $projectdir]
		set compar 0
		if {$resultfile eq ""} {
			set projectdir [file dir $projectdir]
		} else {
			set projectdir [file dir $resultfile]
		}
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
	set allseq_fastqfiles {}
	foreach bam $bams {
		set sampledir [file dir $bam]
		if {$useresultfile} {
			set rootname [file_rootname $resultfile]
		} else {
			set rootname ${analysisname}-[file_rootname $bam]
			set resultfile [file dir $bam]/isoform_counts-$rootname.tsv
		}
		set resultdir [file dir $resultfile]
		set sample [file tail [file dir $resultfile]]
		cd $sampledir
		set workdir $resultdir/$rootname
		set sample [file tail $sampledir]
		mkdir $workdir
		job flair_correct-$rootname {*}$skips -skip $workdir/counts_matrix-$rootname.tsv \
		-cores $threads \
		-deps {
			$bam $bam.bai $refseq $gtfreftranscripts
		} -targets {
			$workdir/all_corrected-$rootname.bed
		} -vars {
			bam rootname workdir refseq gtfreftranscripts threads sample
		} -code {
			analysisinfo_write $bam $target \
				analysis $workdir sample $sample \
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
				-o $workdir/transcripts-$rootname.temp >@ stdout 2>@ stderr
			file rename $workdir/transcripts-$rootname.temp_all_corrected.bed $workdir/all_corrected-$rootname.bed
			catch {file rename $workdir/transcripts-$rootname.temp_all_inconsistent.bed $workdir/all_inconsistent-$rootname.bed}
			file delete $bed12
		}
		set fastqfiles [jobglob fastq/*.fastq.gz]
		lappend allseq_fastqfiles $workdir/allseq-$rootname.fastq.gz
		job flair_allseq-$rootname {*}$skips -skip $rootname/counts_matrix-$rootname.tsv \
		-deps $fastqfiles -targets {
			$workdir/allseq-$rootname.fastq.gz
		} -vars {
			rootname fastqfiles workdir
		} -code {
			# next best on combined data from samples
			if {![llength $fastqfiles]} {error "No fastq files found"}
			set tempfastq $workdir/allseq-$rootname.fastq.gz
			if {![file exists $tempfastq]} {
				puts "Making $tempfastq"
				exec cat {*}$fastqfiles > $tempfastq.temp
				file rename $tempfastq.temp $tempfastq
			}
		}
		job flair_collapse-$rootname {*}$skips -skip $workdir/counts_matrix-$rootname.tsv \
		-cores $threads \
		-deps {
			$workdir/allseq-$rootname.fastq.gz
			$workdir/all_corrected-$rootname.bed
			$refseq $gtfreftranscripts
		} -targets {
			$workdir/transcripts-$rootname.isoforms.gtf
			$workdir/transcripts-$rootname.isoforms.bed
			$workdir/transcripts-$rootname.isoforms.fa
		} -vars {
			rootname refseq gtfreftranscripts threads workdir
		} -code {
			analysisinfo_write $workdir/all_corrected-$rootname.bed $workdir/transcripts-$rootname.isoforms.gtf
			analysisinfo_write $workdir/all_corrected-$rootname.bed $workdir/transcripts-$rootname.isoforms.fa
			puts "collapse -> $rootname-collapse"
			exec [flair_bin identify_gene_isoform.py] \
				$workdir/all_corrected-$rootname.bed \
				$gtfreftranscripts \
				$workdir/all_corrected_renamed-$rootname.bed \
				>@ stdout 2>@ stderr
			catch_exec [flair_bin flair] collapse \
				-t $threads \
				-g $refseq \
				--gtf $gtfreftranscripts \
				-r $workdir/allseq-$rootname.fastq.gz \
				-q $workdir/all_corrected_renamed-$rootname.bed \
				-o $workdir/temptranscripts-$rootname \
				>@ stdout 2>@ stderr
			foreach file [glob $workdir/temptranscripts-$rootname*] {
				file rename -force $file $workdir/[string range [file tail $file] 4 end]
			}
		}
		job flair_quantify-$rootname {*}$skips -cores $threads -deps {
			$workdir/transcripts-$rootname.isoforms.fa
			$workdir/allseq-$rootname.fastq.gz
		} -targets {
			$workdir/counts_matrix-$rootname.tsv
		} -vars {
			rootname sample threads workdir
		} -code {
			analysisinfo_write $workdir/transcripts-$rootname.isoforms.fa $target \
				analysis $rootname sample $sample \
				isocaller flair isocaller_version [version flair]
			set manifestdata {}
			lappend manifestdata [join [list $sample conditionA batch1 $workdir/allseq-$rootname.fastq.gz] \t]
			file_write reads_manifest.tsv [join $manifestdata \n]\n
			puts "quantify -> $workdir/counts_matrix-$rootname.tsv"
			catch_exec [flair_bin flair] quantify \
				--threads $threads \
				-r reads_manifest.tsv \
				-i $workdir/transcripts-$rootname.isoforms.fa \
				-o $workdir/counts_matrix-$rootname.tsv.temp >@ stdout 2>@ stderr
			if {[file exists $workdir/counts_matrix-$rootname.tsv.temp]} {
				# older versions
				file rename -force $workdir/counts_matrix-$rootname.tsv.temp $workdir/counts_matrix-$rootname.tsv
			} else {
				# newer versions
				file rename -force $workdir/counts_matrix-$rootname.tsv.temp.counts.tsv $workdir/counts_matrix-$rootname.tsv
			}
		}
		job flair_convert-$rootname {*}$skips -deps {
			$workdir/transcripts-$rootname.isoforms.gtf
			$workdir/counts_matrix-$rootname.tsv
			$tsvreftranscripts
			$refseq
		} -targets {
			$resultfile
			$resultdir/gene_counts-$rootname.tsv
			$resultdir/totalcounts-$rootname.tsv
		} -vars {
			rootname sample refseq gtfreftranscripts reftranscripts tsvreftranscripts resultfile resultdir workdir
		} -code {
			set flairtranscripts $workdir/transcripts-$rootname.isoforms.gtf
			set flaircountmatrix $workdir/counts_matrix-$rootname.tsv
			set out_isoform_counts_file $resultfile
			set out_gene_counts_file $resultdir/gene_counts-$rootname.tsv
			set out_total_counts_file $resultdir/totalcounts-$rootname.tsv
			set extrainfo [list \
				analysis $rootname sample $sample \
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
				file_write total_counts-$rootname.tsv $sample\n
			} else {
				cg_flair_mergeresults $rootname \
					$tsvreftranscripts $flairtranscripts $flaircountmatrix \
					$out_isoform_counts_file $out_gene_counts_file $out_total_counts_file
			}
		}
		foreach genename $plotgenes {
			job flair_plotisoforms-$rootname-$genename {*}$skips -deps {
				$workdir/transcripts-$rootname.isoforms.bed 
				$workdir/counts_matrix-$rootname.tsv
			} -targets {
				flair_results/${genename}_isoforms.png
			} -vars {
				rootname workdir
			} -code {
				mkdir flair_results
				puts "plot $rootname"
				catch {
					exec [flair_bin plot_isoform_usage] $workdir/transcripts-$rootname.isoforms.bed $workdir/counts_matrix-$rootname.tsv $genename
				}
				file rename ${genename}_isoforms.png flair_results/${genename}_isoforms.png
			}
		}
	}
	if {$compar eq "0"} return
	# combined analysis
	if {$compar eq "joint"} {
		cd $projectdir
		mkdir compar
		set exproot [file tail $projectdir]
		set rootname flair-$exproot
		set workdir compar/$rootname
		set resultdir compar
		set bedfiles [jobglob samples/*/flair-*/all_corrected-flair-*.bed]
		job flair_compar-$exproot {*}$skips \
		-cores $threads \
		-deps [list_concat $bedfiles $allseq_fastqfiles] \
		-targets {
			$workdir/counts_matrix-$rootname.tsv
			$workdir/transcripts-$rootname.isoforms.tsv
			$workdir/transcripts-$rootname.isoforms.gtf
			$workdir/transcripts-$rootname.isoforms.bed
		} -vars {
			workdir rootname bedfiles allseq_fastqfiles exproot sample refseq gtfreftranscripts threads
		} -code {
			analysisinfo_write [lindex $bedfiles 0] $target flair [version flair]
			analysisinfo_write [lindex $bedfiles 0] $workdir/transcripts-$rootname.isoforms.gtf flair [version flair]
			set flairdir [findflair]
			mkdir $workdir
			exec cat {*}$bedfiles > $workdir/all_corrected-$rootname.bed
			# 
			putslog "collapse transcript info -> trancripts-$rootname.*"
			catch_exec [flair_bin flair] collapse \
				-t $threads \
				-g $refseq \
				--gtf $gtfreftranscripts \
				-r [join $allseq_fastqfiles ,] \
				-q $workdir/all_corrected-$rootname.bed \
				-o $workdir/transcripts-$rootname >@ stdout 2>@ stderr
			file delete $workdir/all_corrected-$rootname.bed
			cg gtf2tsv -separate 1 $workdir/transcripts-$rootname.isoforms.gtf $workdir/transcripts-$rootname.isoforms.tsv
			cg gtf2tsv -separate 0 $workdir/transcripts-$rootname.isoforms.gtf $workdir/transcripts-$rootname.isoforms.tsv
			#
			# make manifest (in flair.temp)
			mkdir $workdir/flair.temp
			unset -nocomplain manifestdata
			set condition A
			foreach sample [dirglob samples *] {
				set bam [lindex [glob samples/$sample/map-sminimap*.bam samples/$sample/map-*.bam] 0]
				set fastq samples/$sample/allseq-[file_rootname $bam].fastq.gz
				if {![file exists $fastq]} {
					puts "Making $fastq"
					exec cat {*}[glob samples/$sample/fastq/*.f*q.gz] > $fastq.temp
					file rename $fastq.temp $fastq
				}
				lappend manifestdata($condition) [join [list flair-$sample ${condition} batch1 $fastq] \t]
			}
			set c {}
			foreach condition [array names manifestdata] {
				append c [join $manifestdata($condition) \n]\n
			}
			file_write $workdir/flair.temp/reads_manifest-$rootname.tsv $c
			#
			puts "quantify -> $workdir/counts_matrix-$rootname.tsv"
			catch_exec [flair_bin flair] quantify \
				--tpm \
				-t $threads \
				-r $workdir/flair.temp/reads_manifest-$rootname.tsv \
				-i $workdir/transcripts-$rootname.isoforms.fa \
				-o $workdir/flair.temp/counts_matrix-$rootname.tsv >@ stdout 2>@ stderr
			# remove condition/batch (A_batch1) from fields in header -> just samples
			file delete $workdir/flair.temp/reads_manifest-$rootname.tsv
			if {[file exists $workdir/flair.temp/counts_matrix-$rootname.tsv]} {
				file rename -force $workdir/flair.temp/counts_matrix-$rootname.tsv $workdir/flair.temp/counts_matrix-$rootname.tsv.ori
			} else {
				file rename -force $workdir/flair.temp/counts_matrix-$rootname.tsv.counts.tsv $workdir/flair.temp/counts_matrix-$rootname.tsv.ori
			}
			catch {close $f} ; set f [open $workdir/flair.temp/counts_matrix-$rootname.tsv.ori]
			catch {close $o} ; set o [open $workdir/flair.temp/counts_matrix-$rootname.tsv w]
			set line [split [gets $f] \t]
			puts $o [join [list_regsub -all {_A_batch1$} $line {}] \t]
			fcopy $f $o
			close $o
			close $f
			foreach file [glob $workdir/flair.temp/*] {
				file rename -force $file $workdir/[file tail $file]
			}
			file delete $workdir/flair.temp
			# puts "diffExp on compar/counts_matrix-$rootname.tsv"
			# file delete -force flair-diffexp-$exproot
			# exec flair diffExp -q compar/counts_matrix-$rootname.tsv -o compar/diffexp-$rootname
			# exec python $flairdir/bin/bin/diff_iso_usage.py compar/diffiso-$rootname.tsv
		}
		job flair_convert-compar-$exproot {*}$skips -deps {
			$workdir/transcripts-$rootname.isoforms.gtf
			$workdir/counts_matrix-$rootname.tsv
			$tsvreftranscripts
			$refseq
		} -targets {
			$resultdir/isoform_count-$rootname.tsv
			$resultdir/gene_counts-$rootname.tsv
			$resultdir/totalcounts-$rootname.tsv
		} -vars {
			workdir rootname sample refseq gtfreftranscripts reftranscripts tsvreftranscripts resultdir workdir
		} -code {
			set flairtranscripts $workdir/transcripts-$rootname.isoforms.gtf
			set flaircountmatrix $workdir/counts_matrix-$rootname.tsv
			set out_isoform_counts_file $resultdir/isoform_count-$rootname.tsv
			set out_gene_counts_file $resultdir/gene_counts-$rootname.tsv
			set out_total_counts_file $resultdir/totalcounts-$rootname.tsv
			set extrainfo [list \
				analysis $rootname sample $sample \
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
				file_write total_counts-$rootname.tsv $sample\n
			} else {
				cg_flair_mergeresults {} \
					$tsvreftranscripts $flairtranscripts $flaircountmatrix \
					$out_isoform_counts_file $out_gene_counts_file $out_total_counts_file
			}
			file delete -force $workdir
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

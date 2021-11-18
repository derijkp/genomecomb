proc version_flair {} {
	set flairdir [findflair]
	lindex [split [file tail $flairdir] -] end
}

proc flair_getref {} {
	cd /complgen/z/hg38/extra
	exec wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
	exec gunzip gencode.v37.annotation.gtf.gz
}

proc findflair {} {
	global flair
	if {![info exists flair]} {
		set flair [searchpath flair flair flair*]
		if {$flair eq ""} {
			set flair [searchpath FLAIR flair flair*]
		}
		set ::env(PATH) $flair/bin:$::env(PATH)
	}
	return $flair
}

proc flair_job {args} {
	# putslog [list flair_job {*}$args]
	set cmdline "[list cd [pwd]] \; [list cg flair {*}$args]"
	global appdir
	set refseq {}
	set skips {}
	set genes {}
	set sqanti 1
	set compar 0
	upvar job_logdir job_logdir
	cg_options flair args {
		-refseq {
			set refseq $value
		}
		-genes {
			set genes $value
		}
		-compar {
			set compar $value
		}
		-sqanti {
			set sqanti $value
		}
		-skip {
			lappend skips -skip $value
		}
	} {projectdir}
	set projectdir [file_absolute $projectdir]
	set refseq [refseq $refseq]
	set gtfannotation [glob [file dir $refseq]/extra/gencode*.gtf]
	set flairdir [findflair]
	set samples [glob samples/*]
	cd $projectdir
	job_logfile $projectdir/flair_[file tail $projectdir] $projectdir $cmdline \
		{*}[versions flair dbdir zstd os]
	# analysis per sample
	foreach sample $samples {
		putsvars sample
		cd $projectdir/$sample
		set bam [lindex [glob map-sminimap*.bam map-*.bam] 0]
		set rootname [file_rootname $bam]
		job flair_correct-[file tail $sample] {*}$skips -skip counts_matrix-flair-$rootname.tsv -deps {
			$bam $refseq $gtfannotation
		} -targets {
			all_corrected-flair-$rootname.bed
		} -vars {
			bam rootname flairdir refseq gtfannotation
		} -code {
			set bed12 [file root $bam].bed12
			exec $flairdir/bin/bin/bam2Bed12.py -i $bam > $bed12.temp
			file rename -force $bed12.temp $bed12
			exec flair.py correct -t 8 \
				-g $refseq \
				--gtf $gtfannotation \
				-q $bed12 \
				-o transcripts-flair-$rootname.temp >@ stdout 2>@ stderr
			file rename transcripts-flair-$rootname.temp_all_corrected.bed all_corrected-flair-$rootname.bed
			catch {file rename transcripts-flair-$rootname.temp_all_inconsistent.bed all_inconsistent-flair-$rootname.bed}
			file delete $bed12
		}
		set fastqfiles [glob fastq/*.fastq.gz]
		job flair_allseq-[file tail $sample] {*}$skips -skip counts_matrix-flair-$rootname.tsv \
		-deps $fastqfiles -targets {
			allseq-$rootname.fastq.gz
		} -vars {
			rootname fastqfiles
		} -code {
			# next best on combined data from samples
			set tempfastq allseq-$rootname.fastq.gz
			if {![file exists $tempfastq]} {
				puts "Making $tempfastq"
				set o [open [list | cg gzip -compressionlevel 1 > $tempfastq.temp] w]
				foreach tfile $fastqfiles {
					catchchildkilled_exec {*}[gzcat $tfile] $tfile >@ $o
				}
				close $o
				file rename $tempfastq.temp $tempfastq
			}
		}
		job flair_collapse-[file tail $sample] {*}$skips -skip counts_matrix-flair-$rootname.tsv -deps {
			allseq-$rootname.fastq.gz all_corrected-flair-$rootname.bed $refseq $gtfannotation
		} -targets {
			transcripts-flair-$rootname.isoforms.fa transcripts-flair-$rootname.isoforms.bed
		} -vars {
			rootname refseq gtfannotation
		} -code {
			puts "collapse -> flair-$rootname-collapse"
			exec flair.py collapse \
				-g $refseq \
				--gtf $gtfannotation \
				-r allseq-$rootname.fastq.gz \
				-q all_corrected-flair-$rootname.bed \
				-o temptranscripts-flair-$rootname-collapse >@ stdout 2>@ stderr
			foreach file [glob temptranscripts-flair-$rootname-collapse*] {
				file rename -force $file [string range $file 4 end]
			}
		}
		job flair_quantify-[file tail $sample] {*}$skips -deps {
			transcripts-flair-$rootname.isoforms.fa allseq-$rootname.fastq.gz
		} -targets {
			counts_matrix-flair-$rootname.tsv
		} -vars {
			rootname sample
		} -code {
			set manifestdata {}
			lappend manifestdata [join [list [file tail $sample] conditionA batch1 allseq-$rootname.fastq.gz] \t]
			file_write reads_manifest.tsv [join $manifestdata \n]\n
			puts "quantify -> counts_matrix-flair-$rootname.tsv"
			exec flair.py quantify \
				-r reads_manifest.tsv \
				-i transcripts-flair-$rootname.isoforms.fa \
				-o counts_matrix-flair-$rootname.tsv.temp >@ stdout 2>@ stderr
			file rename -force counts_matrix-flair-$rootname.tsv.temp counts_matrix-flair-$rootname.tsv
		}
		foreach genename $genes {
			job flair_plotisoforms-[file tail $sample]-$genename {*}$skips -deps {
				transcripts-flair-$rootname.isoforms.bed counts_matrix-flair-$rootname.tsv
			} -targets {
				flair_results/${genename}_isoforms.png
			} -vars {
				rootname sample
			} -code {
				mkdir flair_results
				puts "plot $rootname"
				set plot_isoform_usage $flairdir/bin/plot_isoform_usage.py
				catch {
					exec python $plot_isoform_usage transcripts-flair-$rootname.isoforms.bed counts_matrix-flair-$rootname.tsv $genename
				}
				file rename ${genename}_isoforms.png flair_results/${genename}_isoforms.png
			}
		}
	}
	if {!$compar} return
	# combined analysis
	mkdir compar
	set exproot [file tail $projectdir]
	set bedfiles [jobglob samples/*/all_corrected-flair-*-.bed]
	job flair_compar-[file tail $sample] {*}$skips -deps {
		$bedfiles
	} -targets {
		compar/counts_matrix-flair-$exproot.tsv
		compar/transcripts-flair-$exproot.isoforms.genepred.tsv
		compar/transcripts-flair-$exproot.isoforms.gtf
		compar/transcripts-flair-$exproot.isoforms.bed
	} -vars {
		exproot sample refseq gtfannotation
	} -code {
		set tempfastq [tempfile].fastq.gz
		set fastqfiles [glob samples/*/fastq/*.fastq.gz]
		set o [wgzopen $tempfastq]
		foreach file $fastqfiles {
			exec zcat $file >@ $o
		}
		gzclose $o
		if {![file exists compar/transcripts-flair.isoforms.bed]} {
			exec cat {*}[glob samples/*/all_corrected-flair-*.bed] > compar/all_corrected-flair-$exproot.bed
			putslog "collapse -> trancripts-flair-$exproot"
			exec flair.py collapse \
				-g $refseq \
				--gtf $gtfannotation \
				-r $tempfastq \
				-q compar/all_corrected-flair-$exproot.bed \
				-o compar/transcripts-flair-$exproot >@ stdout 2>@ stderr
			file delete compar/all_corrected-flair-$exproot.bed
		}
		cg gtf2tsv -separate 1 compar/transcripts-flair-$exproot.isoforms.gtf compar/transcripts-flair-$exproot.isoforms.tsv
		cg gtf2tsv -separate 0 compar/transcripts-flair-$exproot.isoforms.gtf compar/transcripts-flair-$exproot.isoforms.genepred.tsv
		unset -nocomplain manifestdata
		set condition A
		foreach sample [dirglob samples *] {
			set bam [lindex [glob samples/$sample/map-sminimap*.bam map-*.bam] 0]
			set rootname [file_rootname $bam]
			set fastq samples/$sample/allseq-$rootname.fastq.gz
			if {![file exists $fastq]} {
				puts "Making $fastq"
				set o [open [list | cg gzip -compressionlevel 1 > $tempfastq.temp] w]
				foreach tfile [glob samples/$sample/fastq/*.f*q.gz] {
					catchchildkilled_exec {*}[gzcat $tfile] $tfile >@ $o
				}
				close $o
				file rename $tempfastq.temp $tempfastq
			}
			lappend manifestdata($condition) [join [list $sample ${condition} batch1 $fastq] \t]
		}
		set c {}
		foreach condition [array names manifestdata] {
			append c [join $manifestdata($condition) \n]\n
		}
		file_write compar/reads_manifest-flair-$exproot.tsv $c
		puts "quantify -> compar/counts_matrix-flair-$exproot.tsv"
		mkdir compar/flair.temp
		exec flair.py quantify \
			-r compar/flair.temp/reads_manifest-flair-$exproot.tsv \
			-i compar/flair.temp/transcripts-flair-$exproot.isoforms.fa \
			-o compar/flair.temp/counts_matrix-flair-$exproot.tsv >@ stdout 2>@ stderr
		# remove condition/batch (A_batch1) from fields in header -> just samples
		file rename compar/flair.temp/counts_matrix-flair-$exproot.tsv compar/flair.temp/counts_matrix-flair-$exproot.tsv.ori
		catch {close $f} ; set f [open compar/flair.temp/counts_matrix-flair-$exproot.tsv.ori]
		catch {close $o} ; set o [open compar/flair.temp/counts_matrix-flair-$exproot.tsv w]
		set line [split [gets $f] \t]
		puts $o [join [list_regsub -all {_A_batch1$} $line {}] \t]
		fcopy $f $o
		close $o
		close $f
		foreach file [glob compar/flair.temp/*] {
			file rename $file compar/[file tail $file]
		}
		# puts "diffExp on compar/counts_matrix-flair-$exproot.tsv"
		# file delete -force flair-diffexp-$exproot
		# exec flair.py diffExp -q compar/counts_matrix-flair-$exproot.tsv -o compar/diffexp-flair-$exproot
		# exec python $flairdir/bin/bin/diff_iso_usage.py compar/diffiso-flair-$exproot.tsv
	}
	if {$sqanti} {
		job flair_sqanti_compar-[file tail $sample] {*}$skips -deps {
			compar/transcripts-flair-$exproot.isoforms.gtf
			$gtfannotation
			$refseq$exproot-hg38s-GRN_cDNA_lblast_brain.tsv
			compar/counts_matrix-flair-$exproot.tsv
		} -targets {
			compar/sqanti3-$exproot/transcripts-sqanti3-flair-{$exproot}.tsv
			compar/transcripts-sqanti-flair-$exproot.isoforms.gtf
		} -vars {
			exproot sample refseq gtfannotation gtfannotation
		} -code {
			mkdir compar/sqanti3-$exproot
#			# commented out because data does not get into the results (different ids needed?)
#			set flfile [tempfile].fl
#			catch {close $f} ; set f [open compar/counts_matrix-flair-$exproot.tsv]
#			catch {close $o} ; set o [open $flfile w]
#			set header [split [gets $f] \t]
#			puts $o superPBID\t[join [lrange $header 1 end] \t]
#			while {[gets $f line] != -1} {
#				set line [split $line \t]
#				set new [list [lindex $line 0]]
#				puts $o [lindex $line 0]\t[join [list_regsub -all \\.0 [lrange $line 1 end] {}] \t]
#			}
#			close $o
#			close $f
			exec sqanti3_qc.py \
				compar/transcripts-flair-$exproot.isoforms.gtf \
				$gtfannotation \
				$refseq \
				-d compar/sqanti3-$exproot \
				-o sqanti3-$exproot \
				--saturation \
				--report html
			mklink compar/sqanti3-$exproot/sqanti3-${exproot}_classification.txt compar/transcripts-sqanti3-flair-$exproot.tsv
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
	set plot_isoform_usage $flairdir/bin/bin/plot_isoform_usage.py
	if {[info exists ::env(FONTCONFIG_PATH)]} {
		append ::env(FONTCONFIG_PATH) :/etc/fonts
	} else {
		set ::env(FONTCONFIG_PATH) /etc/fonts
	}
	exec python $plot_isoform_usage $flairtranscriptsbed $flaircounts_matrix $genename $resultprefix
}

proc cg_flair {args} {
	set args [job_init {*}$args]
	flair_job {*}$args
	job_wait
}

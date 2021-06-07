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
	upvar job_logdir job_logdir
	cg_options flair args {
		-refseq {
			set refseq $value
		}
		-genes {
			set genes $value
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
	foreach sample $samples {
		putsvars sample
		cd $projectdir/$sample
		set bam [lindex [glob map-sminimap*.bam] 0]
		set rootname [file_rootname $bam]
		job flair_correct-[file tail $sample] {*}$skips -skip flair-$rootname-counts_matrix.tsv -deps {
			$bam $refseq $gtfannotation
		} -targets {
			flair-$rootname-all_corrected.bed
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
				-o flair-$rootname-corrected.temp >@ stdout 2>@ stderr
			file rename flair-$rootname-corrected.temp_all_corrected.bed flair-$rootname-all_corrected.bed
			catch {file rename flair-$rootname-corrected.temp_all_inconsistent.bed flair-$rootname-all_inconsistent.bed}
			file delete $bed12
		}
		set fastqfiles [glob fastq/*.fastq.gz]
		job flair_allseq-[file tail $sample] {*}$skips -skip flair-$rootname-counts_matrix.tsv \
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
		job flair_collapse-[file tail $sample] {*}$skips -skip flair-$rootname-counts_matrix.tsv -deps {
			allseq-$rootname.fastq.gz flair-$rootname-all_corrected.bed $refseq $gtfannotation
		} -targets {
			flair-$rootname-collapse.isoforms.fa flair-$rootname-collapse.isoforms.bed
		} -vars {
			rootname refseq gtfannotation
		} -code {
			puts "collapse -> flair-$rootname-collapse"
			exec flair.py collapse \
				-g $refseq \
				--gtf $gtfannotation \
				-r allseq-$rootname.fastq.gz \
				-q flair-$rootname-all_corrected.bed \
				-o tempflair-$rootname-collapse >@ stdout 2>@ stderr
			foreach file [glob tempflair-$rootname-collapse*] {
				file rename -force $file [string range $file 4 end]
			}
		}
		job flair_quantify-[file tail $sample] {*}$skips -deps {
			flair-$rootname-collapse.isoforms.fa allseq-$rootname.fastq.gz
		} -targets {
			flair-$rootname-counts_matrix.tsv
		} -vars {
			rootname sample
		} -code {
			set manifestdata {}
			lappend manifestdata [join [list [file tail $sample] conditionA batch1 allseq-$rootname.fastq.gz] \t]
			file_write reads_manifest.tsv [join $manifestdata \n]\n
			puts "quantify -> flair-$rootname-counts_matrix.tsv"
			exec flair.py quantify \
				-r reads_manifest.tsv \
				-i flair-$rootname-collapse.isoforms.fa \
				-o flair-$rootname-counts_matrix.tsv.temp >@ stdout 2>@ stderr
			file rename -force flair-$rootname-counts_matrix.tsv.temp flair-$rootname-counts_matrix.tsv
		}
		foreach genename $genes {
			job flair_plotisoforms-[file tail $sample]-$genename {*}$skips -deps {
				flair-$rootname-collapse.isoforms.bed flair-$rootname-counts_matrix.tsv
			} -targets {
				flair_results/${genename}_isoforms.png
			} -vars {
				rootname sample
			} -code {
				mkdir flair_results
				puts "plot $rootname"
				set plot_isoform_usage $flairdir/bin/plot_isoform_usage.py
				catch {
					exec python $plot_isoform_usage flair-$rootname-collapse.isoforms.bed flair-$rootname-counts_matrix.tsv flair_results/$genename
				}
				file rename flair_results/${genename}_isoforms.png flair_results/${genename}_isoforms.png
			}
		}
	}
	
}

proc cg_flair {args} {
	set args [job_init {*}$args]
	flair_job {*}$args
	job_wait
}

proc cg_bams2crams {args} {
	set args [job_init {*}$args]
	job_logfile [pwd]/bams2crams [pwd] [list cg bams2crams $args]
	set refseq {}
	set threads 1
	set handlebam old
	set links rename
	cg_options cg_bam2cram args {
		-refseq {
			set refseq [refseq $value]
		} 
		-handlebam {
			if {$value ni "rm old keep"} {
				error "wrong value for -handlebam, must be one of: rm, old, keep"
			}
			set handlebam $value
		}
		-links {
			if {$value ni "ignore rename convert error"} {
				error "wrong value for -links, must be one of: ignore rename convert error"
			}
			set links $value
		}
		-threads {set threads $value}
	} {bamfile} 1 ...
	set bamfiles [list $bamfile {*}$args]
	foreach bam $bamfiles {
		set target [file root $bam].cram
		job bams2crams-[file tail $bam] -cores $threads -deps {
			$bam
		} -targets {
			$target
		} -vars {
			bam refseq threads handlebam links
		} -code {
			cg_bam2cram -index 1 -links $links -refseq $refseq -handlebam $handlebam -threads $threads $bam $target
		} 
	}
	job_wait
}

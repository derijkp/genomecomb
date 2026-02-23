proc ubam_split_job {args} {
	upvar job_logdir job_logdir
	set infile -
	set numseq 1000000
	set maxparts .
	set threads 1
	cg_options ubam_split args {
		-numseq {
			set numseq $value
		}
		-maxparts {
			set maxparts $value
		}
		-parts {
			set parts $value
		}
		-threads {
			set threads $value
		}
	} {infile outtemplate} 2 2
	if {![info exists outtemplate]} {
		set outtemplate $infile
		set infile -
	}
	job_logfile [file dir $outtemplate]/ubam_split_[file tail $outtemplate] [file dir $outtemplate]
	set outdir [file dir $outtemplate]
	set outfile [file tail $outtemplate]
	file mkdir $outdir
	if {[info exists parts]} {
		set files {}
		set skips {}
		for {set part 1} {$part <= $parts} {incr part} {
			lappend files $outdir/p${part}_$outfile.temp
		}
		set totalnumseq [exec samtools view --no-PG $infile | countlines]
		set numseq [expr {$totalnumseq / $parts}]
		set maxparts $parts
	}
	set header [exec samtools view --no-PG -H $infile | grep -v ^@PG]
	exec samtools view --no-PG $infile | splitubam $outdir $outfile.temp $numseq $header\n $threads $maxparts
	if {![info exists files]} {
		set files [glob $outdir/p*_$outfile.temp]
	}
	set result {}
	foreach file $files {
		set target [file root $file]
		file rename -- $file $target
		lappend result $target
	}
	return $result
}

proc cg_ubam_split {args} {
	set result [ubam_split_job {*}$args]
	return $result
}

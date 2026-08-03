proc ubam_split_job {args} {
	upvar job_logdir job_logdir
	set infile -
	set numseq 1000000
	set maxparts .
	set threads 1
	set aligned 0
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
		-aligned {
			set aligned $value
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
			lappend files $outdir/p${part}_$outfile
		}
		job [job_relfile2name ubam_split- $infile] -deps {
			$infile
		} -targets $files -vars {
			infile outdir outfile parts threads aligned files
		} -code {
			if {$aligned} {
				set workdir [shadow_workdir $outfile]
				catch_exec samtools collate -l 1 -@ $threads --no-PG $infile $workdir/collate.bam
				catch_exec samtools reset --output-fmt BAM,level=1 $workdir/collate.bam > $workdir/ubam.bam
				set infile $workdir/ubam.bam
			}
			set totalnumseq [exec samtools view -c $infile]
			set numseq [expr {$totalnumseq / $parts}]
			set maxparts $parts
			set header [exec samtools view --no-PG -H $infile | grep -v ^@PG]
			exec samtools view --no-PG $infile | splitubam $outdir $outfile.temp $numseq $header\n $threads $maxparts
			set result {}
			foreach file $files {
				file rename -- $file.temp $file
			}
		}
		set result $files
	} else {
		if {$aligned} {
			set tempfile [tempfile]
			samtools collate -l 1 -@ 4 --no-PG $infile $tempfile.collate.bam
			samtools reset --output-fmt BAM,level=1 $tempfile.collate.bam > $tempfile.ubam.bam
			set infile $tempfile.ubam.bam
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
	}
	return $result
}

proc cg_ubam_split {args} {
	set result [ubam_split_job {*}$args]
	return $result
}

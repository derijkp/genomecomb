proc fastq_split_job {args} {
	set infile -
	set numseq 1000000
	set maxparts .
	set threads 1
	cg_options fastq_split args {
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
	} {infile outtemplate} 1 2
	if {![info exists outtemplate]} {
		set outtemplate $infile
		set infile -
	}
	job_logfile_set [file dir $outtemplate]/fastq_split_[file tail $outtemplate] [file dir $outtemplate]
	set compressed [gziscompressed $outtemplate]
	set ext [gzext $outtemplate]
	set outdir [file dir $outtemplate]
	set outfile [file tail $outtemplate]
	file mkdir $outdir
	if {$infile eq "-"} {
		if {[info exists parts]} {
			error "Cannot user -parts option on stdin"
		}
		set o [open "| splitfastq $outdir $outfile.temp $numseq $maxparts" w+]
		fcopy stdin $o
		set files [read $o]
		close $o
		set files [split [string trim $files] \n]
	} elseif {[info exists parts]} {
		set files {}
		for {set part 1} {$part <= $parts} {incr part} {
			lappend files $outdir/p${part}_$outfile.temp
		}
		job [job_relfile2name fastq_split- $infile] -deps {
			$infile
		} -targets $files -vars {
			infile outdir outfile parts
		} -code {
			set totalnumseq [exec {*}[gzcat $infile] $infile | countlines]
			set numseq [expr {$totalnumseq / (4 * $parts)}]
			set maxparts $parts
			exec {*}[gzcat $infile] $infile | splitfastq $outdir $outfile.temp $numseq $maxparts
		}
	} else {
		set files [exec {*}[gzcat $infile] $infile | splitfastq $outdir $outfile.temp $numseq $maxparts]
		set files [split [string trim $files] \n]
	}
	if {$compressed} {
		set result {}
		foreach file $files {
			set target [file root $file]
			job [job_relfile2name fastq_split_compres- $file] -cores $threads -deps {
				$file
			} -targets {
				$target
			} -vars {
				file ext threads
			} -code {
				compress $file $target 1 0 $threads
			}
		}
	} else {
		set result {}
		foreach file $files {
			set target [file root $file]
			job [job_relfile2name fastq_split_rename- $file] -deps {$file} -targets {[file root $file]} -code {
				file rename -- $dep $target
			}
			lappend result $target
		}
	}
	return $result
}

proc cg_fastq_split {args} {
	set args [job_init {*}$args]
	set result [fastq_split_job {*}$args]
	job_wait
	return $result
}

proc cg_fastq_split.tcl {args} {
	set infile -
	set numseq 1000000
	set partslimit 4294967296
	set threads {}
	cg_options fastq_split args {
		-numseq {
			set numseq $value
		}
		-maxparts {
			set partslimit [expr {$value+1}]
		}
		-threads {
			set threads $value
		}
	} {infile outtemplate} 1 2
	if {![info exists outtemplate]} {
		set outtemplate $infile
		set infile -
	}
	set compressed [gziscompressed $outtemplate]
	set compresspipe [compresspipe $outtemplate -1 $threads]
	set outdir [file dir $outtemplate]
	file mkdir $outdir
	if {$infile eq "-"} {
		set f stdin
	} else {
		set f [gzopen $infile]
	}
	set count 0
	set part 1
	set outlist {}
	while {[gets $f l] != -1} {
		if {$count == 0 && $part < $partslimit} {
			catch {close $o}
			set out $outdir/p${part}_[file tail $outtemplate]
			lappend outlist $out
			if {$compressed} {
				set o [open [list {*}$compresspipe > $out.temp] w]
			} else {
				set o [open $out.temp w]
			}
			file rename -- $out.temp $out
			incr part
			set count $numseq
		}
		puts $o $l
		for {set i 0} {$i < 3} {incr i} {
			set l [gets $f]
			puts $o $l
		}
		incr count -1
	}
	catch {close $o}
}


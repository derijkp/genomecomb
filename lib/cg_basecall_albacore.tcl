# basecaller_albacore_mvresults $name $tempdir/tempfastq $resultdir $barcoding

proc basecaller_albacore_mvresults {barcoding name srcdir resultdir} {
	if {![llength $barcoding]} {
		exec cat {*}[glob $srcdir/workspace/pass/*.fastq] | gzip > $resultdir/pass_$name.fastq.gz
		exec cat {*}[glob $srcdir/workspace/fail/*.fastq] | gzip > $resultdir/fail_$name.fastq.gz
	} else {
		foreach dir [glob -nocomplain $srcdir/workspace/pass/*] {
			set barcode [file tail $dir]
			exec cat {*}[glob $dir/*.fastq] | gzip > $resultdir/${barcode}_pass_$name.fastq.gz
		}
		foreach dir [glob -nocomplain $srcdir/workspace/fail/*] {
			set barcode [file tail $dir]
			exec cat {*}[glob $dir/*.fastq] | gzip > $resultdir/${barcode}_fail_$name.fastq.gz
		}
		foreach barcode $barcoding {
			if {![file exists $resultdir/${barcode}_pass_$name.fastq.gz]} {
				file_write $resultdir/${barcode}_pass_$name.fastq.gz ""
			}
			if {![file exists $resultdir/${barcode}_fail_$name.fastq.gz]} {
				file_write $resultdir/${barcode}_fail_$name.fastq.gz ""
			}
		}
	}
	file copy $srcdir/sequencing_summary.txt $resultdir/sequencing_summary_$name.txt
}

proc basecaller_albacore_job {args} {
	upvar job_logdir job_logdir
	set keepargs $args
	set numreads 4000
	set subdirs 1
	set config {}
	set flowcell {}
	set kit {}
	set opts {}
	set threads 1
	set barcoding {}
	cg_options var_sam args {
		-n - -numreads {
			set numreads $value
		}
		-s - -subdirs {
			set subdirs $value
		}
		-c - -config {
			set config $value
		}
		-f - -flowcell {
			set flowcell $value
		}
		-k - -kit {
			set kit $value
		}
		-filter {
			if {![true $value]} {lappend opts --disable_filtering}
		}
		-b - -barcoding {
			if {$value eq "1"} {
				set bsarcoding {}
				for {set i 1} {$i <= 48} {incr i} {
					lappend barcoding barcode[format %02d $i]
				}
			} else {
				set barcoding $value
			}
			if {$barcoding ne ""} {
				lappend opts --barcoding
			}
		}
		-t - -threads {
			set threads $value
		}
		default {
			lappend opts $key $value
		}
	} {resultdir sourcedir} 2 ... {
		basecall nanopore reads using albacore
	}
	set resultdir [file_absolute $resultdir]
	set sourcedirs [list [file_absolute $sourcedir]]
	foreach dir $args {
		lappend sourcedirs [file_absolute $dir]
	}
	if {$config eq ""} {
		if {$kit eq "" || $flowcell eq ""} {
			error "Either -config or -flowcell and -kit must be specified"
		}
		lappend opts -kit $kit -flowcell $flowcell
	} else {
		lappend opts --config $config
	}
	if {![info exists job_logdir]} {
		job_logdir $resultdir/log_jobs
	}
	set name [file tail $resultdir]
	if {$name eq "fastq"} {
		set name [file tail [file dir $resultdir]]
	}
	job_logfile $resultdir/basecall_albacore_$name $resultdir \
		[list cg basecaller_albacore {*}$keepargs] \
		{*}[versions albacore]

	if {$subdirs} {
		foreach sourcedir $sourcedirs {
			set sourcename [file tail $sourcedir]
			# find subdirs
			set dirs [lsort -dict [dirglob $sourcedir fast5/*]]
			file mkdir $resultdir
			file_write $resultdir/info.txt "source\t$sourcedir\nsubdirs\t[llength $dirs]\n"
			foreach dir $dirs {
				set files [dirglob $sourcedir/$dir *.fast5]
				set len [llength $files]
				if {$subdirs == 1 && $len < $numreads} {
					puts "skipping $dir: not full yet"
					continue
				}
				set name ${sourcename}_[file tail $dir]
				if {![llength $barcoding]} {
					set target1 $resultdir/pass_$name.fastq.gz
					set target2 $resultdir/fail_$name.fastq.gz
					set target3 $resultdir/sequencing_summary_$name.txt
					set targets [list $target1 $target2 $target3]
				} else {
					set targets {}
					foreach barcode $barcoding {
						lappend targets $resultdir/${barcode}_pass_$name.fastq.gz
						lappend targets $resultdir/${barcode}_fail_$name.fastq.gz
					}
					lappend targets $resultdir/${barcode}_sequencing_summary_$name.txt
				}
				
				job albacore-[file tail $sourcedir]-$name -cores $threads \
				-deps {$sourcedir/$dir} \
				-targets $targets \
				-vars {sourcedir from to todo opts threads barcoding resultdir name} \
				-code {
					albacoreinit
					set tempdir [scratchdir]
					exec read_fast5_basecaller.py {*}$opts \
						--files_per_batch_folder 0 \
						--reads_per_fastq_batch 0 \
						-t $threads -o fastq \
						-i $dep \
						-s $tempdir/tempfastq >@ stdout 2>@ stderr
					basecaller_albacore_mvresults $name $tempdir/tempfastq $resultdir $barcoding
#					if {![llength $barcoding]} {
#						exec cat {*}[glob $tempdir/tempfastq/workspace/pass/*.fastq] | gzip > $target1
#						exec cat {*}[glob $tempdir/tempfastq/workspace/fail/*.fastq] | gzip > $target2
#					} else {
#						set target3 [list_pop targets]
#						foreach dir [glob -nocomplain $tempdir/tempfastq/workspace/pass/*] {
#							set barcode [file tail $dir]
#							exec cat {*}[glob $dir/*.fastq] | gzip > $resultdir/${barcode}_pass_$name.fastq.gz
#						}
#						foreach dir [glob -nocomplain $tempdir/tempfastq/workspace/fail/*] {
#							set barcode [file tail $dir]
#							exec cat {*}[glob $dir/*.fastq] | gzip > $resultdir/${barcode}_fail_$name.fastq.gz
#						}
#						foreach barcode $barcoding {
#							if {![file exists $resultdir/${barcode}_pass_$name.fastq.gz]} {
#								file_write $resultdir/${barcode}_pass_$name.fastq.gz ""
#							}
#							if {![file exists $resultdir/${barcode}_fail_$name.fastq.gz]} {
#								file_write $resultdir/${barcode}_fail_$name.fastq.gz ""
#							}
#						}
#					}
					file copy $tempdir/tempfastq/sequencing_summary.txt $target3
					file delete -force $tempdir/fast5 $tempdir/tempfastq
				}
			}
		}
	# old method of selecting 1 directory and finding all the *.fast5
	} else {
		# find fast5 files
		set files [lsort -dict [dirglob $sourcedir *.fast5]]
		set len [llength $files]
		if {$len == 0} {error "no fast5 files found (use -subdirs?)"}
		file mkdir $resultdir
		file_write $resultdir/info.txt "source\t$sourcedir\nsize\t$len\n"
		for {set from 0} {$from < $len} {incr from $numreads} {
			set to [expr {$from+$numreads-1}]
			if {$to > $len} {set to $len}
			set todo [lrange $files $from $to]
			set name [file tail $resultdir]-$from-$to
			if {![llength $barcoding]} {
				set targets [list \
					$resultdir/pass_$name.fastq.gz \
					$resultdir/fail_$name.fastq.gz \
					$resultdir/sequencing_summary_$name.txt]
			} else {
				set targets {}
				foreach barcode $barcoding {
					lappend targets $resultdir/${barcode}_pass_$name.fastq.gz
					lappend targets $resultdir/${barcode}_fail_$name.fastq.gz
				}
				lappend targets $resultdir/${barcode}_sequencing_summary_$name.txt
			}
			job albacore-[file tail $sourcedir]-$from-$to -cores $threads \
			-deps $sourcedir \
			-targets $targets \
			-vars {name sourcedir resultdir barcoding from to todo opts threads} \
			-code {
				albacoreinit
				set tempdir [scratchdir]
				file mkdir $tempdir/fast5
				foreach file $todo {
					mklink $sourcedir/$file $tempdir/fast5/$file 1
				}
				exec read_fast5_basecaller.py {*}$opts \
					--files_per_batch_folder 0 \
					--reads_per_fastq_batch 0 \
					-t $threads -o fastq \
					-i $tempdir/fast5 \
					-s $tempdir/tempfastq >@ stdout 2>@ stderr
				basecaller_albacore_mvresults $name $tempdir/tempfastq $resultdir $barcoding
				file delete -force $tempdir/fast5 $tempdir/tempfastq
			}
		}
	}

}

proc cg_basecaller_albacore {args} {
	set args [job_init {*}$args]
	set result [basecaller_albacore_job {*}$args]
	job_wait
	return $result
}

proc fastq_clipadapters_job {args} {
	upvar job_logdir job_logdir
	set skips {}
	set resultdir {}
	set adapterfile {}
	set paired 1
	set removeskew 1
	set compress gz
	cg_options fastq_clipadapters args {
		-adapterfile {
			set adapterfile $value
		}
		-resultdir {
			set resultdir [file_absolute $value]
		}
		-paired {
			set paired $value
		}
		-skips {
			set skips $value
		}
		-skip {
			lappend skips -skip $value
		}
		-removeskew {
			set removeskew $value
		}
		-compress {
			set compress $value
		}
	} {fastqfile} 1 ... {
		Use fastq-mcf to clip adaptors from fastqs, results are in a dir
		given by -resultdir or named fastq.clipped in the same dir as the (first) fastq
	}
	set fastqfiles [list $fastqfile {*}$args]
	set fastqfiles [bsort $fastqfiles]
	set adapterfile [adapterfile $adapterfile]
	if {$resultdir eq ""} {
		set resultdir [file dir [lindex $fastqfiles 0]]/fastq.clipped
	}
	file mkdir $resultdir
	set resultfastqs {}
	set resultanalysisinfo {}
	foreach file $fastqfiles {
		set file [file_absolute [gzroot $file]]
		set root [file root $file]
		if {$compress eq "gz"} {
			lappend resultfastqs $resultdir/[file tail $root].clipped.fastq.gz
		} elseif {$compress eq "zst"} {
			lappend resultfastqs $resultdir/[file tail $root].clipped.fastq.zst
		} else {
			lappend resultfastqs $resultdir/[file tail $root].clipped.fastq
		}
		lappend resultanalysisinfo [analysisinfo_file $resultdir/[file tail [gzroot $root]].clipped.fastq]
	}
	# paired files need to be clipped together!
	set temptargets {}
	if {[llength $fastqfiles] == 1 || !$paired} {
		foreach {dep} $fastqfiles {target} $resultfastqs a1 $resultanalysisinfo {
			set name [file tail [file dir [file dir $target]]]__[file tail $target]
			job clip-$name {*}$skips -deps {
				$dep
			} -targets {
				$target $a1
			} -vars {
				resultfastqs adapterfile paired removeskew ::analysisinfo compress
			} -code {
				file mkdir [file dir $target]
				analysisinfo_write $dep $target clipping fastq-mcf clipping_version [version fastq-mcf] clipping_cg_version [version genomecomb]
				set tempout1 [filetemp $target]
				if {$compress eq "gz"} {
					set compresspipe "| gzip --fast"
				} elseif {$compress eq "zst"} {
					set compresspipe "| cg zst -compressionlevel 1"
				} else {
					set compresspipe ""
				}
				exec fastq-mcf -k $removeskew -a $adapterfile $dep {*}$compresspipe > $tempout1 2>@ stderr
				file rename -force -- $tempout1 $target
			}
		}
	} else {
		foreach {dep dep2} $fastqfiles {target target2} $resultfastqs {a1 a2} $resultanalysisinfo {
			set name [file tail [file dir $target]]__[file tail $target]
			job clip-$name {*}$skips -deps {
				$dep $dep2
			} -targets {
				$target $target2 $a1 $a2
			} -vars {
				resultfastqs adapterfile paired removeskew ::analysisinfo compress
			} -code {
				file mkdir [file dir $target]
				file mkdir [file dir $target2]
				analysisinfo_write $dep $target clipping fastq-mcf clipping_version [version fastq-mcf] clipping_cg_version [version genomecomb]
				analysisinfo_write $dep2 $target2 clipping fastq-mcf clipping_version [version fastq-mcf] clipping_cg_version [version genomecomb]
				set tempfile1 [tempfile]
				set tempfile2 [tempfile]
				set tempout1 [filetemp $target]
				set tempout2 [filetemp $target2]
				exec fastq-mcf -k $removeskew -a -o $tempfile1 -o $tempfile2 $adapterfile $dep $dep2 2>@ stderr
				if {$compress eq "gz"} {
					exec gzip --fast -c $tempfile1 > $tempout1
					exec gzip --fast -c $tempfile2 > $tempout2
					file delete $tempfile1
					file delete $tempfile2
					file rename -force -- $tempout1 $target
					file rename -force -- $tempout2 $target2
				} elseif {$compress eq "zst"} {
					exec cg zst -k 0 -compressionlevel 1 -o $tempout1 $tempfile1
					exec cg zst -k 0 -compressionlevel 1 -o $tempout2 $tempfile2
					file rename -force -- $tempout1 $target
					file rename -force -- $tempout2 $target2
				} else {
					file rename -force -- $tempfile1 $target
					file rename -force -- $tempfile2 $target2
				}
			}
		}
	}
	return $resultfastqs
}

proc cg_fastq_clipadapters {args} {
	set args [job_init {*}$args]
	fastq_clipadapters_job {*}$args
	job_wait
}

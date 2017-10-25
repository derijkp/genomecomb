proc fastq_clipadapters {files targets args} {
	set adapterfile {}
	set removeskew 1
	set paired 1
	foreach {key value} $args {
		if {$key eq "-adapterfile"} {
			set adapterfile $value
		} elseif {$key eq "-paired"} {
			set paired $value
		} elseif {$key eq "-removeskew"} {
			set removeskew $value
		} else {
			lappend opts $key $value
		}
	}
	set adapterfile [adapterfile $adapterfile]
	# clip primers, quality
	set temptargets {}
	if {[llength $files] == 1 || !$paired} {
		foreach {f1} $files {t1} $targets {
			file mkdir [file dir $t1]
			set tempout1 [filetemp $t1]
			exec fastq-mcf -k $removeskew -a -o $tempout1 $adapterfile $f1 2>@ stderr
			lappend temptargets $tempout1
		}
	} else {
		foreach {f1 f2} $files {t1 t2} $targets {
			file mkdir [file dir $t1]
			file mkdir [file dir $t2]
			set tempout1 [filetemp $t1]
			set tempout2 [filetemp $t2]
			exec fastq-mcf -k $removeskew -a -o $tempout1 -o $tempout2 $adapterfile $f1 $f2 2>@ stderr
			lappend temptargets $tempout1 $tempout2
		}
	}
	foreach target $targets temptarget $temptargets {
		file mkdir [file dir $target]
		file rename -force $temptarget $target
	}
}

proc fastq_clipadapters_job {files args} {
	upvar job_logdir job_logdir
	set targets {}
	set skips {}
	set adapterfile {}
	set paired 1
	set removeskew 1
	set files [ssort -natural $files]
	foreach {key value} $args {
		if {$key eq "-adapterfile"} {
			set adapterfile $value
		} elseif {$key eq "-paired"} {
			set paired $value
		} elseif {$key eq "-skips"} {
			set skips $value
		} elseif {$key eq "-removeskew"} {
			set removeskew $value
		} else {
			lappend opts $key $value
		}
	}
	set adapterfile [adapterfile $adapterfile]
	foreach file $files {
		set file [file_absolute [gzroot $file]]
		set root [file root $file]
		file mkdir [file dir $root].clipped
		lappend targets [file dir $root].clipped/[file tail $root].clipped.fastq
	}
	# paired files need to be clipped together!
	job clip-[file dir [file dir $root]] -deps $files -targets $targets \
	-vars {adapterfile paired removeskew ::analysisinfo} {*}$skips -code {
		fastq_clipadapters $deps $targets -removeskew $removeskew -adapterfile $adapterfile -paired $paired
		foreach dep $deps $target $targets {
			analysisinfo_write $dep $target clipping fastq-mcf clipping_version [version fastq-mcf] clipping_cg_version [version genomecomb]
		}
	}
	return $targets
}


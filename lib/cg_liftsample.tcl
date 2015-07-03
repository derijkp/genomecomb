proc liftsample_job {args} {
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-split - -s {
				set split $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {([llength $args] < 3)} {
		errorformat liftover
		exit 1
	}
	foreach {srcdir destdir liftoverfile} $args break
	if {[file exists $destdir] && ![file isdir $destdir]} {
		error "$destdir already exists and is not a directory"
	}
	if {![file isdir $srcdir]} {
		error "$srcdir is not a (sample) directory"
	}
	job_logdir $destdir/log_jobs
	unset -nocomplain infoa
	array set infoa [fileinfo $srcdir]
	if {[file exists $srcdir/sampleinfo.tsv]} {
		array set infoa [infofile_read $srcdir/sampleinfo.tsv]
	}
	if {[file exists $srcdir/info.txt]} {
		set c [split [file_read $srcdir/info.txt] \n]
		foreach line $c {
			foreach {key value} [split [string range $line 1 end] \t] break
			set infoa($key) $value
		}
	}
	if {[info exists infoa(split)]} {
		if {$infoa(split) ne $split} {
			error "split option $split given, but sample is $infoa(split) according to sampleinfo.tsv"
		}
	} elseif {![info exists split]} {
		set split 1
	}
	set infoa(split) $split
	lappend infoa(liftover) $liftoverfile
	file mkdir $destdir
	infofile_write $destdir/sampleinfo.tsv [array get infoa]
	foreach file [jobglob $srcdir/var-*.tsv $srcdir/fannotvar-*.tsv] {
		set destfile $destdir/[file tail $file]
		set regionfile [findregionfile $file]
		job liftvar-[file tail $file] -deps {$file $liftoverfile ($regionfile)} \
		-vars {liftoverfile regionfile split} \
		-targets {$destfile} -code {
			if {![catch {file link $dep} link]} {
				putslog "Copying link $dep"
				file copy -force $dep $target
			} else {
				putslog "converting $dep"
				file delete $target.temp
				if {[jobglob $regionfile] ne ""} {
					cg liftover -regionfile $dep2 -split $split $dep $target.temp $liftoverfile 2>@ stderr
				} else {
					cg liftover -split $split $dep $target.temp $liftoverfile 2>@ stderr
				}
				file rename -force $target.temp $target
				catch {file rename -force $target.temp.unmapped $target.unmapped}
			}
		}
	}
	foreach file [jobglob $srcdir/sreg-*.tsv $srcdir/reg_*.tsv] {
		set target $destdir/[file tail $file]
		job liftreg-[file tail $file] -deps {$file $liftoverfile} \
		-vars {file liftoverfile} \
		-targets {$target} -code {
			if {![catch {file link $file} link]} {
				putslog "Copying link $file"
				file copy -force $file $target
			} else {
				putslog "converting region $file"
				cg liftregion $file $target.temp $liftoverfile
				file rename -force $target.temp $target
				catch {file rename -force $target.temp.unmapped $target.unmapped}
			}
		}
	}
	foreach file [jobglob $srcdir/cgcnv-*.tsv $srcdir/cgsv-*.tsv] {
		set target $destdir/[file tail $file]
		job liftreg-[file tail $file] -deps {$file $liftoverfile} \
		-vars {file liftoverfile} \
		-targets {$target} -code {
			if {![catch {file link $file} link]} {
				putslog "Copying link $file"
				file copy -force $file $target
			} else {
				putslog "converting $file"
				cg liftover $file $target.temp $liftoverfile 2>@ stderr
				file rename -force $target.temp $target
				catch {file rename -force $target.temp.unmapped $target.unmapped}
			}
		}
	}
}

proc cg_liftsample {args} {
	set args [job_init {*}$args]
	liftsample_job {*}$args
	job_wait
}

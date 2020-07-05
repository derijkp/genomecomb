proc liftsample_job {args} {
	cg_options liftsample args {
		-split - -s {
			set split $value
		}
	} {srcdir destdir liftoverfile} 3 3
	if {[file exists $destdir] && ![file isdir $destdir]} {
		error "$destdir already exists and is not a directory"
	}
	if {![file isdir $srcdir]} {
		error "$srcdir is not a (sample) directory"
	}
	if {![info exists job_logdir]} {
		set_job_logdir $destdir/log_jobs
	}
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
	set newsample [file tail $destdir]
	infofile_write $destdir/sampleinfo.tsv [array get infoa]
	foreach file [jobglob $srcdir/var-*.tsv $srcdir/fannotvar-*.tsv] {
		regsub -- {-[^-]+.tsv$} [file tail $file] -$newsample.tsv temp
		set target $destdir/$temp
		if {[regexp {(.*)fannotvar-(.*)\.tsv$} $file temp temp1 temp2]} {
			set regionfile [findregionfile $temp1/var-cg-cg-$temp2.tsv]
		} else {
			set regionfile [findregionfile $file]
		}
		job [job_relfile2name liftvar- $file] -deps {$file $liftoverfile ($regionfile)} \
		-vars {file liftoverfile regionfile split newsample} \
		-targets {$target} -code {
			if {![catch {file link $dep} link]} {
				putslog "Copying link $dep"
				regsub -- {-[^-]+.tsv$} $link -$newsample.tsv newlink
				file_link $target $newlink
			} else {
				putslog "converting $file"
				file delete $target.temp
				set regionfile [jobglob $regionfile]
				if {$regionfile ne ""} {
					cg liftover -regionfile $regionfile -split $split $file $target.temp $liftoverfile 2>@ stderr
				} else {
					cg liftover -split $split $file $target.temp $liftoverfile 2>@ stderr
				}
				file rename -force -- $target.temp $target
				catch {file rename -force -- $target.temp.unmapped $target.unmapped}
			}
		}
	}
	foreach file [jobglob $srcdir/sreg-*.tsv $srcdir/reg_*.tsv] {
		regsub -- {-[^-]+.tsv$} [file tail $file] -$newsample.tsv temp
		set target $destdir/$temp
		job [job_relfile2name liftreg- $file] -deps {$file $liftoverfile} \
		-vars {file liftoverfile newsample} \
		-targets {$target} -code {
			if {![catch {file link $file} link]} {
				putslog "Copying link $file"
				regsub -- {-[^-]+.tsv$} $link -$newsample.tsv newlink
				file_link $target $newlink
			} else {
				putslog "converting region $file"
				cg liftregion $file $target.temp $liftoverfile
				file rename -force -- $target.temp $target
				catch {file rename -force -- $target.temp.unmapped $target.unmapped}
			}
		}
	}
	foreach file [jobglob $srcdir/cgcnv-*.tsv $srcdir/cgsv-*.tsv] {
		regsub -- {-[^-]+.tsv$} [file tail $file] -$newsample.tsv temp
		set target $destdir/$temp
		job [job_relfile2name liftreg- $file] -deps {$file $liftoverfile} \
		-vars {file liftoverfile newsample} \
		-targets {$target} -code {
			if {![catch {file link $file} link]} {
				putslog "Copying link $file"
				regsub -- {-[^-]+.tsv$} $link -$newsample.tsv newlink
				file_link $target $newlink
			} else {
				putslog "converting $file"
				cg liftover $file $target.temp $liftoverfile 2>@ stderr
				file rename -force -- $target.temp $target
				catch {file rename -force -- $target.temp.unmapped $target.unmapped}
			}
		}
	}
}

proc cg_liftsample {args} {
	set args [job_init {*}$args]
	liftsample_job {*}$args
	job_wait
}

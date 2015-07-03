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
		if {[file exists $destfile]} continue
		if {![catch {file link $file} link]} {
			putslog "Copying link $file"
			file copy -force $file $destfile
		} else {
			putslog "converting $file"
			set regionfile [findregionfile $file]
			if {[jobglob $regionfile] ne ""} {
				cg liftover -regionfile $regionfile -split $split $file $destfile.temp $liftoverfile 2>@ stderr
			} else {
				cg liftover -split $split $file $destfile.temp $liftoverfile 2>@ stderr
			}
			file rename $destfile.temp $destfile
			catch {file rename $destfile.temp.unmapped $destfile.unmapped}
		}
	}
	foreach file [jobglob $srcdir/sreg-*.tsv $srcdir/reg_*.tsv] {
		set destfile $destdir/[file tail $file]
		if {[file exists $destfile]} continue
		if {![catch {file link $file} link]} {
			putslog "Copying link $file"
			file copy -force $file $destfile
		} else {
			putslog "converting region $file"
			cg liftregion $file $destfile.temp $liftoverfile
			file rename $destfile.temp $destfile
			catch {file rename $destfile.temp.unmapped $destfile.unmapped}
		}
	}
	foreach file [jobglob $srcdir/cgcnv-*.tsv $srcdir/cgsv-*.tsv] {
		set destfile $destdir/[file tail $file]
		if {[file exists $destfile]} continue
		if {![catch {file link $file} link]} {
			putslog "Copying link $file"
			file copy -force $file $destfile
		} else {
			putslog "converting $file"
			cg liftover $file $destfile.temp $liftoverfile 2>@ stderr
			file rename $destfile.temp $destfile
			catch {file rename $destfile.temp.unmapped $destfile.unmapped}
		}
	}
}

proc cg_liftsample {args} {
	set args [job_init {*}$args]
	liftsample_job {*}$args
	job_wait
}

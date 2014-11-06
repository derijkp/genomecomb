proc indexdir {mainfile} {
	global configdata
	set file [file normalize [gzroot $mainfile]]
	if {[info exists configdata(indexdir,$file)]} {return $configdata(indexdir,$file)}
	set configdir [configdir]
	set tail [file tail $file]
	set dirs [glob -nocomplain $configdir/indexdirs/${tail}_*]
	set num 0
	foreach indexdir $dirs {
		set normfilename [file_read $indexdir/info.normfilename]
		if {$normfilename eq $file} {
			set configdata(indexdir,$file) $indexdir
			return $indexdir
		}
		regexp {[0-9]+$} $indexdir dirnum
		if {$dirnum > $num} {set num $dirnum}
	}
	incr num
	set indexdir $configdir/indexdirs/${tail}_$num
	file mkdir $indexdir
	file_write $indexdir/info.normfilename $file
	if {$file eq $mainfile} {
		file_write $indexdir/info.filesize [file size $file]
	}
	set configdata(indexdir,$file) $indexdir
}

proc indexdir_clean {} {
	set configdir [configdir]
	set dirs [glob -nocomplain $configdir/indexdirs/*]
	foreach indexdir $dirs {
		set normfilename [file_read $indexdir/info.normfilename]
		if {![file exists $normfilename]} {
			file delete -force $indexdir
		}
	}
}

proc cg_indexclean {} {
	indexdir_clean
}

proc indexdir_cache_check {mainfile indexdir indexfile} {
	set depfile $indexdir/$indexfile
	if {![file exists $depfile]} {return 0}
	if {![file readable $depfile]} {return 0}
	if {[file mtime $depfile] < [file mtime $mainfile]} {return 0}
	if {[file exists $indexdir/info.filesize]} {
		if {[file mtime $depfile] < [file mtime $indexdir/info.filesize]} {return 0}
		set filesize [file_read $indexdir/info.filesize]
		if {[file size $mainfile] != $filesize} {
			return 0
		}
	}
	return 1
}

# use this if you want to write to the index file, regardless of whether an up to date version exists
proc indexdir_filewrite {mainfile indexfile} {
	set indexdir [file normalize [gzroot $mainfile]].index
	# first try sib indexdir
	if {[file exists $indexdir/$indexfile]} {
		# if we can delete the obsolete file, we can also write the new one
		if {![catch {
			file delete $indexdir/$indexfile
		}]} {
			return $indexdir/$indexfile
		}
	} else {
		# check if we can write in the sib indexdir
		if {![catch {
			file mkdir [file dir $indexdir/$indexfile.temp]
			file delete $indexdir/$indexfile.temp
			file_write $indexdir/$indexfile.temp {}
			file delete $indexdir/$indexfile.temp
		}]} {
			return $indexdir/$indexfile
		}
	}
	# not writable in sib indexdir, return user indexdir
	set userindexdir [indexdir $mainfile]
	return $userindexdir/$indexfile
}

# returns the indexfile in the indexdir.
# The var ok will be 1 if an up to date indexfile is present
# if ok is 0 (indexfile is not up to date), the returned file location is writable
# if okVar is not given, only an existing, up to date indefile will be returned,
# if there is none, the result will be empty
proc indexdir_file {mainfile indexfile {okVar {}}} {
	if {$okVar ne ""} {upvar $okVar ok}
	set ok 1
	set indexdir [file normalize [gzroot $mainfile]].index
	# is a good file existing in either user indexdir or sib indexdir (ok is 1 if we do find one)
	set userindexdir [indexdir $mainfile]
	if {[indexdir_cache_check $mainfile $userindexdir $indexfile]} {
		return $userindexdir/$indexfile
	}
	if {[indexdir_cache_check $mainfile $indexdir $indexfile]} {
		return $indexdir/$indexfile
	}
	if {$okVar eq ""} {return {}}
	# nothing found, see if we can create in sib indexdir
	set ok 0
	if {[file exists $indexdir/$indexfile]} {
		# if we can delete the obsolete file, we can also write the new one
		if {![catch {
			file delete $indexdir/$indexfile
		}]} {
			return $indexdir/$indexfile
		}
	} else {
		# check if we can write in the sib indexdir
		if {![catch {
			file mkdir [file dir $indexdir/$indexfile.temp]
			file delete $indexdir/$indexfile.temp
			file_write $indexdir/$indexfile.temp {}
			file delete $indexdir/$indexfile.temp
		}]} {
			return $indexdir/$indexfile
		}
	}
	# sib indexdir does not work, use user local indexdir
	return $userindexdir/$indexfile
}


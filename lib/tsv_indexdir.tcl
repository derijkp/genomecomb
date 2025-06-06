proc indexdirfind {mainfile} {
	global configdata
	set file [file_absolute [gzroot $mainfile]]
	if {[info exists configdata(indexdir,$file)]} {
		 if {[file exists $configdata(indexdir,$file)]} {
			return $configdata(indexdir,$file)
		} else {
			unset configdata(indexdir,$file)
		}
	}
	set configdir [configdir]
	set tail [file tail $file]
	set dirs [glob -nocomplain $configdir/indexdirs/${tail}/*]
	set num 0
	foreach indexdir $dirs {
		if {![file exists $indexdir/info.normfilename]} {
			continue
		}
		set normfilename [file_read $indexdir/info.normfilename]
		if {$normfilename eq $file} {
			set configdata(indexdir,$file) $indexdir
			return $indexdir
		}
		regexp {[0-9]+$} $indexdir dirnum
		if {$dirnum > $num} {set num $dirnum}
	}
	return {}
}

proc indexdir {mainfile} {
	global configdata
	set file [file_absolute [gzroot $mainfile]]
	# if {[info exists configdata(indexdir,$file)]} {return $configdata(indexdir,$file)}
	set configdir [configdir]
	set tail [file tail $file]
	set dirs [glob -nocomplain $configdir/indexdirs/${tail}/*]
	set num 0
	foreach indexdir $dirs {
		if {![file exists $indexdir/info.normfilename]} {
			continue
		}
		set normfilename [file_read $indexdir/info.normfilename]
		if {$normfilename eq $file} {
			set configdata(indexdir,$file) $indexdir
			return $indexdir
		}
		regexp {[0-9]+$} $indexdir dirnum
		if {$dirnum > $num} {set num $dirnum}
	}
	incr num
	set indexdir $configdir/indexdirs/${tail}/$num
	mkdir $indexdir
	file_write $indexdir/info.normfilename $file
	if {$file eq $mainfile} {
		if {[file exists $file]} {
			file_write $indexdir/info.filesize [file size $file]
		}
	}
	set configdata(indexdir,$file) $indexdir
}

proc indexdir_clean {} {
	set configdir [configdir]
	set dirs [glob -nocomplain $configdir/indexdirs/*/*]
	foreach indexdir $dirs {
		if {[file exists $indexdir/info.normfilename]} {
			set normfilename [file_read $indexdir/info.normfilename]
		} else {
			set normfilename {}
		}
		if {![file exists $normfilename]} {
			file delete -force $indexdir
			set dir [file dir $indexdir]
			if {[catch {glob $dir/*}]} {
				file delete $dir
			}
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
	if {![file exists $mainfile] || [file mtime $depfile] < [file mtime $mainfile]} {return 0}
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
	set indexdir [file_absolute [gzroot $mainfile]].index
	# first try sib indexdir
	if {[file exists $indexdir/$indexfile]} {
		# if it is a dir, check if we can write in it
		if {[file isdir $indexdir/$indexfile]} {
			if {![catch {
				file_write $indexdir/$indexfile/indexdir_filewrite.test {}
			}]} {
				file delete $indexdir/$indexfile/indexdir_filewrite.test
				return $indexdir/$indexfile
			}
			
		} else {
			# if we can delete the obsolete file, we can also write the new one
			if {![catch {
				file delete $indexdir/$indexfile
			}]} {
				return $indexdir/$indexfile
			}
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
	set mainfile [gzfile $mainfile]
	set ok 1
	set indexdir [file_absolute [gzroot $mainfile]].index
	# if a good indexfile exists in indexdir, return it  (ok is 1 if we do find one)
	if {[indexdir_cache_check $mainfile $indexdir $indexfile]} {
		return $indexdir/$indexfile
	}
	if {$okVar eq ""} {
		# try to avoid using userindexdir (slow if it grows too much)
		# -> prefer making a new one in indexdir
		# only check for existing indexfile in userindexdir if creating a new indexdir is not allowed
		set userindexdir [indexdirfind $mainfile]
		if {$userindexdir ne "" && [indexdir_cache_check $mainfile $userindexdir $indexfile]} {
			return $userindexdir/$indexfile
		}
		return {}
	}
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
	# sib indexdir does not work, use user local indexdir, create if not found
	set userindexdir [indexdir $mainfile]
	if {$userindexdir ne "" && [indexdir_cache_check $mainfile $userindexdir $indexfile]} {
		return $userindexdir/$indexfile
	}
	return $userindexdir/$indexfile
}


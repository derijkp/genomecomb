proc compress {file {ext .lz4}} {
	set file [file_absolute $file]
	if {[file exists $file$ext]} {file delete $file$ext}
	if {[inlist {.rz} $ext]} {
		exec razip $file
	} elseif {[inlist {.lz4} $ext]} {
		exec lz4c -q -9 -c $file > $file.lz4.temp
		file rename -force $file.lz4.temp $file.lz4
	} elseif {[inlist {.gz} $ext]} {
		exec gzip $file
	} elseif {[inlist {.bgz} $ext]} {
		exec bgzip $file
	} elseif {[inlist {.bz2} $ext]} {
		exec bzip2 $file
	} else {
		error "Unknown extension $ext"
	}
}

proc gzopen {file {pos -1}} {
	if {![file exists $file]} {
		exiterror "Error: couldn't open \"$file\": no such file or directory"
	}
	set file [file_absolute $file]
	set ext [file extension $file]
	if {[inlist {.rz} $ext]} {
		if {$pos == -1} {
			set f [open "| razip -d -c $file"]
		} else {
			set f [open "| razip -d -c -b $pos $file"]
		}
	} elseif {[inlist {.lz4} $ext]} {
		if {$pos == -1} {
			set f [open "| lz4c -d -c $file"]
		} else {
			set f [open "| lz4ra $file $pos"]
		}
	} elseif {[inlist {.bgz .gz} $ext]} {
		if {$pos == -1} {
			set f [open "| zcat $file"]
		} else {
			error "positioning not supported in (b)gz files"
		}
	} elseif {[inlist {.bz2} $ext]} {
		if {$pos == -1} {
			set f [open "| bzcat $file"]
		} else {
			error "positioning not supported in bz2 files"
		}
	} else {
		set f [open $file]
		if {$pos != -1} {
			seek $f $pos
		}
	}
	return $f
}

proc gzclose {f} {
	if {$f in {stdin stdout}} return
	if {[catch {close $f} error]} {
		if {$error eq "child killed: write on pipe with no readers"} return
		if {[regexp {Successfully decoded [0-9]+ bytes} $error]} return
		error $error
	}
}

proc gzcatch {cmd} {
	if {[catch {uplevel $cmd} error]} {
		if {$error eq "child killed: write on pipe with no readers"} return
		if {[regexp {Successfully decoded [0-9]+ bytes} $error]} return
		error $error
	}
}

proc gzroot filename {
	if {[inlist {.rz .lz4 .gz .bgz .bz2} [file extension $filename]]} {
		return [file root $filename]
	} else {
		return $filename
	}
}

proc gziscompressed filename {
	if {[inlist {.rz .lz4 .gz .bgz .bz2} [file extension $filename]]} {
		return 1
	} else {
		return 0
	}
}

proc gzexists {filename {checkcompressed 1}} {
	if {$checkcompressed} {
		expr {[file exists $filename] || [file exists $filename.rz] || [file exists $filename.lz4] || [file exists $filename.gz] ||[file exists $filename.bgz] || [file exists $filename.bz2]}
	} else {
		file exists $filename
	}
}

proc gzfile {args} {
	foreach filename $args {
		if {![catch {glob $filename $filename.rz $filename.lz4 $filename.bgz $filename.gz $filename.bz2} list]} {
			return [lindex $list 0]
		}
	}
	return [lindex $args 0]
}

proc gzfile_multi {filelist} {
	set result {}
	foreach filename $filelist {
		if {![catch {glob $filename $filename.rz $filename.lz4 $filename.bgz $filename.gz $filename.bz2} list]} {
			lappend result [lindex $list 0]
		} else {
			lappend result $filename
		}
	}
	return $result
}

proc gzfiles {args} {
	set result {}
	foreach filename $args {
		if {![catch {glob $filename $filename.rz $filename.lz4 $filename.bgz $filename.gz $filename.bz2} list]} {
			lappend result {*}$list
		}
	}
	return $result
}

proc gzarraynames {aVar pattern} {
	upvar $aVar a
	set result [lsort [list_remdup [list_concat [array names a $pattern] [array names a $pattern.rz] [array names a $pattern.lz4] [array names a $pattern.gz] [array names a $pattern.bgz] [array names a $pattern.bz2]]]]
	return $result
}

proc gzcat {filename} {
	switch [file extension $filename] {
		.rz {set cat "razip -d -c"}
		.lz4 {set cat "lz4c -d -c"}
		.gz - .bgz {set cat zcat}
		.bz2 {set cat bzcat}
		default {set cat cat}
	}
	return $cat
}

proc gztemp {filename} {
	set ext [file extension $filename]
	switch $ext {
		.gz {
			set tempfile [scratchfile get]
			exec gunzip -d -c $filename > $tempfile
			set ::gztemp_files($tempfile) 1
			return $tempfile
		}
		.rz {
			set tempfile [scratchfile get]
			exec razip -d -c $filename > $tempfile
			set ::gztemp_files($tempfile) 1
			return $tempfile
		}
		.lz4 {
			set tempfile [scratchfile get]
			exec lz4c -q -d -c $filename > $tempfile
			set ::gztemp_files($tempfile) 1
			return $tempfile
		}
		.bz2 {
			set tempfile [scratchfile get]
			exec bzcat $filename > $tempfile
			set ::gztemp_files($tempfile) 1
			return $tempfile
		}
		default {
			set ::gztemp_files($filename) 0
			return $filename
		}
	}
}

proc gzrmtemp {args} {
	foreach filename $args {
		if {$::gztemp_files($filename)} {
			file delete $filename
		}
	}
}

proc decompress {file args} {
	if {[llength $args]} {
		set resultfile [lindex $args 0]
	} else {
		set resultfile [file root $file]
	}
	if {[file extension $file] eq ".lz4"} {
		set error [catch {exec lz4c -q -d $file > $resultfile.temp} result]
	} else {
		set error [catch {exec zcat $file > $resultfile.temp} result]
	}
	if $error {
		if {![regexp "decompression OK, trailing garbage ignored" $result] && ![regexp {Successfully decoded} $result]} {
			error $result
		}
	}
	file rename -force $resultfile.temp $resultfile
	if {![llength $args]} {
		file delete $file
	}
	return $result
}

proc gunzip {file args} {
	if {[llength $args]} {
		set resultfile [lindex $args 0]
	} else {
		set resultfile [file root $file]
	}
	set error [catch {exec zcat $file > $resultfile.temp} result]
	if $error {
		if {![regexp "decompression OK, trailing garbage ignored" $result] && ![regexp {Successfully decoded} $result]} {
			error $result
		}
	}
	file rename -force $resultfile.temp $resultfile
	if {![llength $args]} {
		file delete $file
	}
	return $result
}

proc razip_job {file args} {
	set deps [gzfiles $file {*}$args]
	uplevel [list job razip-$file -checkcompressed 0 -deps $deps -targets $file.rz -rmtargets $file -code {
		if {![file exists $dep]} {error "error compressing: file $dep does not exist"}
		cg_razip $dep
	}]
}

proc lz4_job {file args} {
	upvar job_logdir job_logdir
	set deps [jobglob $file]
	job lz4-$file -checkcompressed 0 -deps $deps -targets $file.lz4 -rmtargets $file -vars args -code {
		if {![file exists $dep]} {error "error compressing: file $dep does not exist"}
		cg_lz4 -keep 0 {*}$args $dep
	}
}

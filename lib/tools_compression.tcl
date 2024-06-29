proc setdefcompressionlevel {level} {
	set ::defcompressionlevel $level
}

proc defcompressionlevel {{default 8}} {
	if {[info exists ::defcompressionlevel]} {
		return $::defcompressionlevel
	} else {
		return $default
	}
}

proc compressionlevel {{compressionlevel {}} {default 8} {min 1} {max 9}} {
	if {$compressionlevel in {{} -1}} {
		if {[info exists ::defcompressionlevel]} {
			set compressionlevel $::defcompressionlevel
		} else {
			set compressionlevel $default
		}
	}
	if {$compressionlevel < $min} {return $min}
	if {$compressionlevel > $max} {return $max}
	return $compressionlevel
}

proc setdefcompressionthreads {level} {
	set ::defcompressionthreads $level
}

proc defcompressionthreads {{default 1}} {
	if {[info exists ::defcompressionthreads]} {
		return $::defcompressionthreads
	} else {
		return $default
	}
}

proc compressionthreads {{compressionthreads {}} {default 1}} {
	if {![isint $compressionthreads] || $compressionthreads < 1} {
		if {[info exists ::defcompressionthreads]} {
			set compressionthreads $::defcompressionthreads
		} else {
			set compressionthreads $default
		}
	}
	return $compressionthreads
}

proc compress {file {destfile {}} {index 1} {keep 1} {threads {}} {compressionlevel {}} {blocksize {}} args} {
	# putsvars file destfile index keep compressionlevel blocksize args
	set ext [file extension $destfile]
	set method [string range $ext 1 end]
	if {$compressionlevel eq ""} {set compressionlevel [defcompressionlevel]}
	if {[auto_load compress_$method]} {
		compress_$method $file $destfile $index $keep $threads $compressionlevel $blocksize
	} elseif {[gziscompressed $file]} {
		set temp [filetemp $destfile]
		exec {*}[gzcat $file] $file > $temp
		file rename -force -- $temp $destfile
		if {!$keep} {file delete $file}
	} elseif {$keep} {
		set temp [filetemp $destfile]
		file copy -force $file $temp
		file rename -force -- $temp $destfile
	} else {
		file rename -force -- $file $destfile
	}
}

proc compress_template {file destfile method cmd {index 1} {keep 1}} {
	# putsvars file destfile index keep
	if {$destfile eq ""} {set destfile $file.$method}
	if {$file eq "-"} {
		if {$destfile eq "-"} {
			exec {*}$cmd <@ stdin >@ stdout 2>@ stderr
		} else {
			set temp [filetemp $destfile]
			exec {*}$cmd <@ stdin > $temp 2>@ stderr
		}
	} elseif {[gziscompressed $file]} {
		if {$destfile eq "-"} {
			exec {*}[gzcat $file] | {*}$cmd >@ stdout
		} else {
			set temp [filetemp $destfile]
			exec {*}[gzcat $file] $file | {*}$cmd > $temp
		}
		if {!$keep} {file delete $file}
	} else {
		set temp [filetemp $destfile]
		exec {*}$cmd $file > $temp
		if {!$keep} {file delete $file}
	}
	if {$destfile ne "-"} {
		file rename -force -- $temp $destfile
		if {$index} {
			catch {index_$method $destfile}
		}
	}
}

proc wgzopen {file {compressionlevel -1} {threads {}} {pipe {}}} {
	set root [file root [gzroot $file]]
	if {$root eq "-"} {
		if {![gziscompressed $file]} {
			if {$pipe eq ""} {
				return stdout
			} else {
				set o [open [list {*}$pipe >@ stdout" w]
				set ::genomecomb_gzopen_info($o) $pipe
				return $o
			}
		} else {
			set o [open [list {*}$pipe {*}[compresspipe $file $compressionlevel $threads] >@ stdout] w]
			set ::genomecomb_gzopen_info($o) $pipe
			return $o
		}
	} elseif {![gziscompressed $file]} {
		if {$pipe eq ""} {
			set o [open $file w]
			set ::genomecomb_gzopen_info($o) $file
			return $o
		} else {
			set o [open [list {*}$pipe > $file] w]
			set ::genomecomb_gzopen_info($o) $pipe
			return $o
		}
	} elseif {$pipe eq ""} {
		set o [open [list {*}[compresspipe $file $compressionlevel $threads] > $file] w]
		set ::genomecomb_gzopen_info($o) $file
		return $o
	} else {
		set o [open [list {*}$pipe {*}[compresspipe $file $compressionlevel $threads] > $file] w]
		set ::genomecomb_gzopen_info($o) $pipe
		return $o
	}
}

proc gzopen {file {pos -1}} {
	set root [file root [gzroot $file]]
	set ext [file extension $file]
	if {$root eq "-"} {
		if {$pos != -1} {
			error "cannot provide random access on stdin"
		}
		set in "<@ stdin"
	} else {
		if {![file exists $file]} {
			error "Error: couldn't open \"$file\": no such file or directory"
		}
		set file [file_absolute $file]
		if {[file size $file] == 0} {
			# we sometimes make empty files with compression extension
			# avoid these giving errors
			set f [open $file]
			if {$pos != -1} {
				seek $f $pos
			}
			set ::genomecomb_gzopen_info($f) $file
			return $f
		}
		set in [list $file]
	}
	if {[inlist {.rz} $ext]} {
		if {$pos == -1} {
			set f [open "| razip -d -c $in"]
		} else {
			set f [open "| razip -d -c -b $pos $in"]
		}
	} elseif {[inlist {.zst} $ext]} {
		if {$pos == -1} {
			set f [open "| zstd-mt -T 1 -k -d -c $in"]
		} else {
			set f [open "| zstdra [list $file] $pos"]
		}
	} elseif {[inlist {.lz4} $ext]} {
		if {$pos == -1} {
			set f [open "| lz4 -d -c $in"]
		} else {
			set f [open "| lz4ra [list $file] $pos"]
		}
	} elseif {[inlist {.gz} $ext]} {
		if {$pos == -1} {
			set f [open "| zcat $in"]
		} else {
			error "positioning not supported in gz files"
		}
	} elseif {[inlist {.bgz} $ext]} {
		if {$pos == -1} {
			set f [open "| zcat $in"]
		} else {
			set f [open "| bgzip -d -c -b $pos [list $file]"]
		}
	} elseif {[inlist {.bz2} $ext]} {
		if {$pos == -1} {
			set f [open "| bzcat $in"]
		} else {
			error "positioning not supported in bz2 files"
		}
	} elseif {$root eq "-"} {
		return stdin
	} else {
		set f [open $file]
		if {$pos != -1} {
			seek $f $pos
		}
	}
	set ::genomecomb_gzopen_info($f) $file
	return $f
}

proc gzclose {f} {
	if {$f in {stdin stdout stderr}} return
	if {[catch {close $f} error]} {
		# if {$error eq "child process exited abnormally"} return
		if {[join [list_remdup [split $error \n]] \n] eq "child killed: write on pipe with no readers"} return
		if {$error eq "error writing \"stdout\": broken pipe"} return
		if {[regexp {Successfully decoded [0-9]+ bytes} $error]} return
		error "error closing file \"[get ::genomecomb_gzopen_info($f) $f]\": $error"
	}
	unset -nocomplain ::genomecomb_gzopen_info($f)
}

proc gzcloseout {f} {
	if {$f in {stdin stdout stderr}} return
	if {[catch {close $f} error]} {
		regsub "\n?child process exited abnormally\$" $error {} error
		error $error
	}
	unset -nocomplain ::genomecomb_gzopen_info($f)
}

proc gzcatch {cmd} {
	if {[catch {uplevel $cmd} error]} {
		if {$error eq "child killed: write on pipe with no readers"} return
		if {[regexp {Successfully decoded [0-9]+ bytes} $error]} return
		error $error
	}
}

array set ::gzexts {
	.zst zst
	.rz razip
	.lz4 lz4
	.gz gzip
	.bgz bgzip
	.bz2 bzip2
}

proc isgzext ext {
	info exists ::gzexts($ext)
}

proc gzext file {
	set ext [file extension $file]
	if {[isgzext $ext]} {
		return $ext
	} else {
		return ""
	}
}

proc gzroot filename {
	if {[isgzext [file extension $filename]]} {
		return [file root $filename]
	} else {
		return $filename
	}
}

proc gziscompressed filename {
	if {[isgzext [file extension $filename]]} {
		return 1
	} else {
		return 0
	}
}

proc gzexists {filename {checkcompressed 1}} {
	if {$checkcompressed} {
		expr {[file exists $filename] || [file exists $filename.zst] || [file exists $filename.rz] || [file exists $filename.lz4] || [file exists $filename.gz] ||[file exists $filename.bgz] || [file exists $filename.bz2]}
	} else {
		file exists $filename
	}
}

proc gzfile {args} {
	foreach filename $args {
		if {$filename eq ""} {return {}}
		if {![catch {glob $filename $filename.zst $filename.rz $filename.lz4 $filename.bgz $filename.gz $filename.bz2} list]} {
			return [lindex $list 0]
		}
	}
	return [lindex $args 0]
}

proc gzfile_multi {filelist} {
	set result {}
	foreach filename $filelist {
		if {![catch {glob $filename $filename.zst $filename.rz $filename.lz4 $filename.bgz $filename.gz $filename.bz2} list]} {
			lappend result [lindex $list 0]
		} else {
			lappend result $filename
		}
	}
	return $result
}

#proc gzfiles {args} {
#	foreach filename $args {
#		set list [glob -nocomplain $filename $filename.zst $filename.lz4 $filename.rz $filename.bgz $filename.gz $filename.bz2]
#		foreach file $list {
#			set root [gzroot $file]
#			if {[info exists a($root)]} continue
#			set a($root) $file
#		}
#	}
#	set result {}
#	foreach file [array names a] {
#		lappend result $a($file)
#	}
#	return $result
#}

# list files matching any of the given patterns
# including if they are compressed (give pattern for uncompressed filename)
# if the same file is present with multiple compression, only one version is given (in gzorder)
# results are sorted:
#   first all files matching first given pattern, then the second pattern, etc.
#   per pattern files are sorted (natural using bsort)
proc gzfiles {args} {
	set result {}
	unset -nocomplain a
	set gzorder {{} .zst .lz4 .rz .bgz .gz .bz2}
	foreach filename $args {
		if {[string first * $filename] != -1} {
			set list [glob -nocomplain $filename $filename.zst $filename.lz4 $filename.rz $filename.bgz $filename.gz $filename.bz2]
			unset -nocomplain thisa
			foreach file $list {
				set root [gzroot $file]
				if {[info exists a($root)]} continue
				if {[info exists thisa($root)]} {
					if {[lsearch $gzorder [gzext $file]] < [lsearch $gzorder [gzext $thisa($root)]]} {
						set thisa($root) $file
					}
				} else {
					set thisa($root) $file
				}
			}
			foreach root [bsort [array names thisa]] {
				set a($root) $thisa($root)
				lappend result $thisa($root)
			}
		} else {
			foreach ext {{} .zst .lz4 .rz .bgz .gz .bz2} {
				set file $filename$ext
				if {[file exists $file]} {
					set root [gzroot $file]
					set a($root) $file
					lappend result $file
					break
				}
			}
		}
	}
	return $result
}

proc gzarraynames {aVar pattern} {
	upvar $aVar a
	set result [lsort [list_remdup [list_concat [array names a $pattern.zst] [array names a $pattern] [array names a $pattern.zst] [array names a $pattern.zst] [array names a $pattern.rz] [array names a $pattern.lz4] [array names a $pattern.gz] [array names a $pattern.bgz] [array names a $pattern.bz2]]]]
	return $result
}

proc gzcat {filename {checkexists 1}} {
	if {$checkexists && ![file exists $filename]} {error "file $filename does not exist"}
	switch [file extension $filename] {
		.zst {return "zstd-mt -T 1 -k -q -d -c"}
		.rz {return "razip -d -c"}
		.lz4 {return "lz4 -q -d -c"}
		.gz - .bgz {return zcat}
		.bz2 {return bzcat}
		default {return cat}
	}
}

proc gzcatra {filename {pos 0} {checkexists 1}} {
	if {$checkexists && ![file exists $filename]} {error "file $filename does not exist"}
	if {$pos == 0} {
		list {*}[gzcat $filename] $filename
	} else {
		switch [file extension $filename] {
			.zst {
				return [list zstdra $filename $pos]
			}
			.rz {
				return [list razip -d -c -b $pos $file]
			}
			.lz4 {
				return [list lz4ra $filename $pos]
			}
			.gz {
				return [list zcat $filename | tail -c +[expr {$pos + 1}]]
			}
			.bgz {
				return [list bgzip -d -c -b $pos $filename]
			}
			.bz2 {
				return [list bzcat $filename | tail -c +[expr {$pos + 1}]]
			}
			default {
				return [list tail -c +[expr {$pos + 1}] $filename]
			}
		}
	}
}


proc compresscmd {target {threads {}} {compressionlevel {}} {blocksize {}}} {
	set method [string range [file extension $target] 1 end]
	if {[auto_load compresscmd_$method]} {
		compresscmd_$method $threads $compressionlevel $blocksize
	} else {
		error "compression $method not supported"
	}
}

proc compresspipe {target {compressionlevel {}} {threads {}} {blocksize {}}} {
	if {![catch {
		compresscmd [file extension $target] $threads $compressionlevel $blocksize
	} compresscmd]} {
		return "| $compresscmd"
	} else {
		return {}
	} 
}

proc gztemp {filename} {
	set ext [file extension $filename]
	switch $ext {
		.zst {
			set tempfile [scratchfile get]
			exec zstd-mt -T 1 -k -q -d -c $filename > $tempfile
			set ::gztemp_files($tempfile) 1
			return $tempfile
		}
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
			exec lz4 -q -d -c $filename > $tempfile
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
	set error [catch {
		exec {*}[gzcat $file] $file > $resultfile.temp
	} result]
	if $error {
		if {![regexp "decompression OK, trailing garbage ignored" $result] && ![regexp {Successfully decoded} $result]} {
			error $result
		}
	}
	file rename -force -- $resultfile.temp $resultfile
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
	file rename -force -- $resultfile.temp $resultfile
	if {![llength $args]} {
		file delete $file
	}
	return $result
}

proc razip_job {file args} {
	set deps [gzfiles $file {*}$args]
	uplevel [list job razip-$file -checkcompressed 0 -deps $deps -targets {$file.rz} -rmtargets {$file} -code {
		if {![file exists $dep]} {error "error compressing: file $dep does not exist"}
		cg razip $dep
	}]
}

proc lz4_job {file args} {
	upvar job_logdir job_logdir
	set deps [jobglob -checkcompressed 1 $file]
	set defcompressionlevel [defcompressionlevel]
	job lz4-$file -checkcompressed 0 -deps $deps -targets {$file.lz4} -rmtargets {$file} -vars {
		args defcompressionlevel
	} -code {
		defcompressionlevel $defcompressionlevel
		if {![file exists $dep]} {error "error compressing: file $dep does not exist"}
		cg_lz4 -keep 0 {*}$args $dep
	}
}

proc lz4index_job {args} {
	upvar job_logdir job_logdir
	set skips {}
	cg_options lz4index args {
		-skip {
			lappend skips -skip $value
		}
	} {file} 1
	job lz4index-$file {*}$skips -checkcompressed 0 -deps {$file} -targets {$file.lz4i} -code {
		if {![file exists $dep]} {error "error indexing: file $dep does not exist"}
		cg_lz4index $dep
	}
}

proc zst_job {file args} {
	upvar job_logdir job_logdir
	set deps [jobglob -checkcompressed 1 $file]
	set defcompressionlevel [defcompressionlevel]
	job zst-$file -checkcompressed 0 -deps $deps -targets {$file.zst} -rmtargets {$file} -vars {
		args defcompressionlevel
	} -code {
		defcompressionlevel $defcompressionlevel
		if {![file exists $dep]} {error "error compressing: file $dep does not exist"}
		zst -keep 0 {*}$args $dep
	}
}

proc zst {args} {
	cg_compress_job -method zst {*}$args
}

proc zstindex_job {args} {
	upvar job_logdir job_logdir
	set skips {}
	cg_options zstindex args {
		-skip {
			lappend skips -skip $value
		}
	} {file} 1
	job zstindex-$file {*}$skips -checkcompressed 0 -deps {$file} -targets {$file.zsti} -code {
		if {![file exists $dep]} {error "error indexing: file $dep does not exist"}
		zstindex $dep
	}
}

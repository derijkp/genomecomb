proc bcol_indexlines {file indexfile} {
	set time [file mtime $file]
	set ext [file extension $file]
	if {[inlist {.rz .bgz .gz} $ext]} {set compressed 1} else {set compressed 0}
	if {![file exists $indexfile] || [file mtime $indexfile] < $time} {
		file mkdir [file dir $indexfile]
		if {$compressed} {
			progress start 1 "uncompressing $file for indexing, please be patient"
			progress message "uncompressing $file for indexing, please be patient (no progress shown)"
			if {$ext eq ".rz"} {
				set indexdir [gzroot $file].index
				set tempfile $indexdir/[file root [file tail $file]]
				if {![file exists $tempfile]} {
					gunzip $file $tempfile
				}
			} else {
				gunzip $file
				set file [file root $file]
				set tempfile $file
			}
			progress stop
		} else {
			set tempfile $file
		}
		set f [open $tempfile]
		set header [tsv_open $f comment]
		close $f
		if {[catch {tsv_basicfields $header 4} poss]} {
			if {![catch {tsv_basicfields $header 3} poss]} {
				lappend poss -1
			} else {
				set poss {-1 -1 -1 -1}
			}
		}
		if {![file exists $indexfile] || [file mtime $indexfile] < $time} {
			progress start [file size $tempfile] "Indexing $file, please be patient"
			progress message "Indexing $file, please be patient"
			exec bcol_indexfile $tempfile $indexfile.temp $indexfile.bin.temp {*}$poss
			file rename -force $indexfile.bin.temp $indexfile.bin
			file rename -force $indexfile.temp $indexfile
			progress stop
		}
		if {$compressed} {
			if {$ext eq ".rz"} {
				file delete $tempfile
			} else {
				cg razip $file
				set file $file.rz
			}
		}
	}
}

proc bcol_open indexfile {
	set result [dict create objtype bcol file $indexfile binfile $indexfile.bin compressedbin 0]
	set f [open $indexfile]
	set header [tsv_open $f comment]
	set comment [split $comment \n]
	if {[lindex $comment 0] ne "# binary column"} {error "file \"$indexfile\" is not a binary column file"}
	if {$header ne "begin type offset"} {error "file \"$indexfile\" has an incorrect table header"}
	foreach line $comment {
		set line [string range $line 1 end]
		dict set result [lindex $line 0] [lindex $line 1]
	}
	set table {}
	while {![eof $f]} {
		set line [gets $f]
		if {![llength $line]} continue
		lappend table [split $line \t]
	}
	close $f
	set tabled [dict create]
	list_foreach {num type offset} $table {
		if {$type eq "end"} break
		dict set tabled $num $offset
	}
	dict set result table $table
	dict set result tabled $tabled
	dict set result max [lindex $table end 0]
	if {[file exists $indexfile.bin]} {
		set fi [open $indexfile.bin]
		fconfigure $fi -encoding binary -translation binary
		dict set result fi $fi
	} elseif {[file exists $indexfile.bin.rz]} {
		dict set result binfile $indexfile.bin.rz
		dict set result compressedbin 1
	} else {
		exiterror "binfile $indexfile.bin not found"
	}
	return $result
}

proc bcol_size bcol {
	expr {[dict get $bcol max]+1}
}

proc bcol_close bcol {
	catch {close [dict get $bcol fi]}
}

proc bcol_get {bcol start end} {
	if {[dict get $bcol objtype] ne "bcol"} {error "This is not a bcol object: $bcol"}
	set type [dict get $bcol type]
	if {$type eq "lineindex"} {set type iu}
	set btype [string index $type 0]
	array set typesizea {c 1 s 2 i 4 w 8 f 4 d 8}
	set typesize $typesizea($btype)
	set max [dict get $bcol max]
 	set table [dict get $bcol table]
	set tabled [dict get $bcol tabled]
	if {[catch {dict get $bcol default} default]} {
		set default 0
	}
	if {$start > $max} {return [list_fill [expr {$end-$start+1}] $default]}
	if {$end > $max} {set uend $max} else {set uend $end}
	set offset 0
	set result {}
	if {[llength $table] > 2} {
		list_foreach {num type noffset} $table {
			if {$start < $num} {
				break
			}
			set offset $noffset
		}
	} else {
		set num [lindex $table 0 0]
		while {$start < $num} {
			if {$start == $end} {return $result}
			lappend result $default
			incr start
			incr curpos
		}
	}
	if {$start < $num} {
		set result [list_fill [expr {$num-$start}] $default
		set start $num
	}
	set len [expr {$uend-$start+1}]
	if {$len <= 0} {return $result}
	set pos [expr {$typesize*($start-$num)}]
	set binfile [dict get $bcol binfile]
	if {![dict get $bcol compressedbin]} {
		set f [open $binfile]
		seek $f $pos
		fconfigure $f -encoding binary -translation binary
		set b [read $f [expr {$typesize*$len}]]
	} else {
		# set f [open "| razip -d -c -s [expr {$typesize*$len}] -b $pos $binfile"]
		set b [exec razip -d -c -s [expr {$typesize*($len+1)}] -b $pos $binfile]
	}
	binary scan $b $type$len sresult
	catch {close $f}
	foreach v $sresult {
		catch {set offset [dict get $tabled $start]}
		incr start
		lappend result [expr {$v + $offset}]
	}
	if {$uend < $end} {lappend result {*}[list_fill [expr {$end-$uend}] $default]}
	return $result
}

proc bcol_table {bcol {start {}} {end {}}} {
	if {[dict get $bcol objtype] ne "bcol"} {error "This is not a bcol object: $bcol"}
	set type [dict get $bcol type]
	if {$type eq "lineindex"} {set type iu}
	set btype [string index $type 0]
	array set typesizea {c 1 s 2 i 4 w 8 f 4 d 8}
	set typesize $typesizea($btype)
	set max [dict get $bcol max]
 	set table [dict get $bcol table]
	set tabled [dict get $bcol tabled]
	if {[catch {dict get $bcol default} default]} {
		set default 0
	}
	if {$start eq ""} {set start 0}
	if {$end eq ""} {set end $max}
	if {$start > $max} {return [list_fill [expr {$end-$start+1}] $default]}
	if {$end > $max} {set uend $max} else {set uend $end}
	set offset 0
	set curpos $start
	set o stdout
	set name [file root [file tail [dict get $bcol file]]]
	puts $o "pos\tvalue"
	if {[llength $table] > 2} {
		list_foreach {num type noffset} $table {
			if {$start < $num} {
				break
			}
			set offset $noffset
		}
	} else {
		set num [lindex $table 0 0]
	}
	while {$start < $num} {
		if {$start == $end} {return}
		puts $o $curpos\t$default
		incr start
		incr curpos
	}
	set len [expr {$uend-$start+1}]
	if {$len <= 0} {return {}}
	set pos [expr {$typesize*($start-$num)}]
	set binfile [dict get $bcol binfile]
	if {![dict get $bcol compressedbin]} {
		set f [open $binfile]
		seek $f $pos
		fconfigure $f -encoding binary -translation binary
	} else {
		set f [open "| razip -d -c -b $pos $binfile"]
		fconfigure $f -encoding binary -translation binary
	}
	while {$curpos <= $end} {
		catch {set offset [dict get $tabled $start]}
		incr start
		set b [read $f $typesize]
		binary scan $b $type value
		set value [expr {$value + $offset}]
		puts $o $curpos\t$value
		incr curpos
	}
	catch {close $f}
	while {$uend < $end} {
		puts $o $curpos\t$default
		incr curpos
		incr uend
	}
}

array set bcol_typea {c,mn -127 s,mn -32768 i,mn -2147483648 w,mn -9223372036854775808 cu,mn 0 su,mn 0 iu,mn 0}
array set bcol_typea {c,mx  127 s,mx  32767 i,mx  2147483647 w,mx  9223372036854775807 cu,mx 255 su,mx 65535 iu,mx 4294967295}

proc bcol_make {args} {
	global bcol_typea
	set type iu
	set offsetcol {}
	set defaultvalue 0
	set pos 0
	set distribute 0
	set start 0
	foreach {key value} $args {
		switch -- $key {
			-t - --type {
				set type $value
			}
			-p - --poscol {
				set offsetcol $value
			}
			-d - --default {
				set defaultvalue $value
			}
			-n - --nan {
				set nanvalue $value
			}
			-c - --chromosomecol {
				set chromosomecol $value
				set distribute 1
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {[llength $args] != 2} {
		exiterror "wrong # args: should be \"cg bcol make ?options? bcolprefix column\""
	}
	foreach {prefix valuecolumn} $args break
	set result $prefix.bcol
	if {[info exists bcol_typea($type,mx)]} {
		set max $bcol_typea($type,mx)
		set min $bcol_typea($type,mn)
	} elseif {![inlist {f d} $type]} {
		exiterror "error: type $type not supported (must be one of: c s i w cu su iu f d)"
	}
	if {[info exists max]} {
		if {$defaultvalue > $max} {
			exiterror "default value $v too large for type $type"
		} elseif {$defaultvalue < $min} {
			exiterror "default value $v too small for type $type"
		}
	}
	set btype [string index $type 0]
	# putslog "Making $result"
	set f stdin
	set header [tsv_open $f comment]
	set colpos [lsearch $header $valuecolumn]
	if {$colpos == -1} {
		exiterror "error: valuecolumn $valuecolumn not found"
	}
	if {$offsetcol eq ""} {
		set offsetpos -1
	} else {
		set offsetpos [lsearch $header $offsetcol]
		if {$offsetpos == -1} {
			exiterror "error: pos column $offsetcol not found"
		}
	}
	if {$distribute} {
		set chrompos [lsearch $header $chromosomecol]
		if {$chrompos == -1} {
			exiterror "error: chromosome column $chromosomecol not found"
		}
	}
	set bo [open $result.bin.temp w]
	fconfigure $bo -encoding binary -translation binary
	set prevchr {}
	set len 0
	set poffset -1
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		if {$distribute} {
			set chr [lindex $line $chrompos]
			if {$chr ne $prevchr} {
				close $bo
				if {$prevchr ne ""} {
					set size [expr {$start+$len-1}]
					set o [open $prefix-$prevchr.bcol.temp w]
					puts $o "# binary column"
					puts $o "# type $type"
					puts $o "# default $defaultvalue"
					puts $o [join {begin type offset} \t]
					puts $o $start\t$type\t0
					puts $o $size\tend\t0
					close $o
					file rename -force $prefix-$prevchr.bcol.bin.temp $prefix-$prevchr.bcol.bin
					file rename -force $prefix-$prevchr.bcol.temp $prefix-$prevchr.bcol
				}
				set prevchr $chr
				set poffset -1
				set len 0
				# putslog "Making $prefix-$chr.bcol"
				set bo [open $prefix-$chr.bcol.bin.temp w]
				fconfigure $bo -encoding binary -translation binary
			}
		}
		if {$offsetpos != -1} {
			set offset [lindex $line $offsetpos]
			if {$poffset == -1} {
				set start $offset
			} elseif {$poffset != $offset} {
				set size [expr {$offset-$poffset}]
				if {$size < 0} {
					exiterror "error: cannot make position based bcol on unsorted data (sort on position first)"
				}
				while {$poffset < $offset} {
					puts -nonewline $bo [binary format $btype $defaultvalue]
					incr poffset
				}
				incr len $size
			}
			set poffset [incr offset]
		}
		set v [lindex $line $colpos]
		if {![isdouble $v] || $v eq "N"} {
			if {[info exists nanvalue]} {
				set v $nanvalue
			} else {
				exiterror "value $v is not a number"
			}
		} elseif {[info exists max]} {
			if {$v > $max} {
				exiterror "value $v too large for type $type"
			} elseif {$v < $min} {
				exiterror "value $v too small for type $type"
			}
		}
		puts -nonewline $bo [binary format $btype $v]
		incr len
	}
	close $bo
	if {$distribute} {
		set size [expr {$start+$len-1}]
		set o [open $prefix-$prevchr.bcol.temp w]
	} else {
		set size [expr {$start+$len-1}]
		set o [open $result.temp w]
	}
	puts $o "# binary column"
	puts $o "# type $type"
	puts $o "# default $defaultvalue"
	puts $o [join {begin type offset} \t]
	puts $o $start\t$type\t0
	puts $o $size\tend\t0
	close $o
	if {$distribute} {
		file rename -force $prefix-$prevchr.bcol.bin.temp $prefix-$prevchr.bcol.bin
		file rename -force $prefix-$prevchr.bcol.temp $prefix-$prevchr.bcol
	} else {
		file rename -force $result.bin.temp $result.bin
		file rename -force $result.temp $result
	}
}

proc cg_bcol {cmd args} {
	switch $cmd {
		get {
			if {[llength $args] != 3} {
				exiterror "wrong # args: should be \"cg bcol get file start end\""
			}
			foreach indexfile $args break
			set bcol [bcol_open $indexfile]
			set result [bcol_get $bcol {*}[lrange $args 1 end]]
			bcol_close $bcol
			puts $result
			return $result
		}
		table {
			if {[llength $args] < 1 || [llength $args] > 3} {
				exiterror "wrong # args: should be \"cg bcol table file ?start? ?end?\""
			}
			foreach {indexfile start end} $args break
			set bcol [bcol_open $indexfile]
			set result [bcol_table $bcol {*}[lrange $args 1 end]]
			bcol_close $bcol
			puts $result
			return $result
		}
		size {
			if {[llength $args] != 1} {
				exiterror "wrong # args: should be \"cg bcol size file\""
			}
			foreach indexfile $args break
			set bcol [bcol_open $indexfile]
			set size [bcol_size $bcol]
			bcol_close $bcol
			puts $size
			return $size
		}
		make {
			bcol_make {*}$args
		}
		default {
			error "unknown subcommand $cmd of bcol, must be one of: get, size"
		}
	}
}

proc cg_index file {
	set time [file mtime $file]
	set indexdir [gzroot $file].index
	set ext [file extension $file]
	if {[inlist {.rz .bgz .gz} $ext]} {set compressed 1} else {set compressed 0}
	file mkdir $indexdir
	set indexfile $indexdir/lines.bcol
	bcol_indexlines $file $indexfile
	if {![file exists $indexdir/info.tsv] || [file mtime $indexdir/info.tsv] < $time} {
		set f [gzopen $file]
		set header [tsv_open $f]
		catch {close $f}
		set bcol [bcol_open $indexfile]
		set size [bcol_size $bcol]
		bcol_close $bcol
		set f [open $indexdir/info.tsv w]
		puts $f key\tvalue
		puts $f file\t$file
		puts $f lineindexfile\t[file tail $indexfile]
		puts $f header\t$header
		puts $f size\t$size
		close $f
	}
	return $indexfile
}

proc cg_size file {
	set indexfile [cg_index $file]
	set bcol [bcol_open $indexfile]
	set size [bcol_size $bcol]
	bcol_close $bcol
	puts $size
	return $size
}

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
			# putslog "bcol_indexfile $tempfile $indexfile.temp $indexfile.bin.temp {*}$poss"
			exec bcol_indexfile $tempfile $indexfile.temp $indexfile.bin.temp {*}$poss 2> /dev/null
			file rename -force $indexfile.bin.temp $indexfile.bin
			file rename -force $indexfile.temp $indexfile
			progress stop
		}
		if {$compressed} {
			if {$ext eq ".rz"} {
				file delete $tempfile
			} else {
				exec razip -c $file > $file.rz
			}
		}
	}
}

proc bcol_open indexfile {
	set result [dict create objtype bcol file $indexfile binfile $indexfile.bin compressedbin 0]
	set f [open $indexfile]
	set header [tsv_open $f comment]
	set comment [split $comment \n]
	if {[lindex $comment 0] ne "# binary column"} {
		close $f
		error "file \"$indexfile\" is not a binary column file"
	}
	if {$header ne "begin type offset"} {
		close $f
		error "file \"$indexfile\" has an incorrect table header"
	}
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
		dict set tabled [expr {$num+1}] $offset
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
		set tempfile [tempfile]
		exec razip -d -c $indexfile.bin.rz > $tempfile
		dict set result tempfile $tempfile
		set fi [open $tempfile]
		fconfigure $fi -encoding binary -translation binary
		dict set result fi $fi
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
	dict set bcol fi {}
	if {[dict exists $bcol tempfile]} {
		file delete [dict get $bcol tempfile]
	}
}

# Warning:
# while cg bcol get follows the half-open convention, end position is not included
# the procedure bcol_get follows the Tcl lrange convention: end position is included!
proc bcol_get {bcol start {end {}}} {
	if {[dict get $bcol objtype] ne "bcol"} {error "This is not a bcol object: $bcol"}
	if {$end eq ""} {
		set end $start
	}
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
	set len [expr {$end-$start+1}]
	if {$len <= 0} {return {}}
	if {$start > $max} {return [list_fill $len $default]}
	if {$end > $max} {set uend $max} else {set uend $end}
	set result {}
	set offset 0
	if {[llength $table] > 2} {
		list_foreach {num temp noffset} $table {
			if {$start <= $num} {
				break
			}
			set offset $noffset
		}
		set num [lindex $table 0 0]
	} else {
		set num [lindex $table 0 0]
	}
	if {$start < $num} {
		if {$end < $num} {
			return [list_fill [expr {$end-$start+1}] $default]
		} else {
			set result [list_fill [expr {$num-$start}] $default]
		}
		set start $num
	}
	set len [expr {$uend-$start+1}]
	if {$len <= 0} {return $result}
	set pos [expr {$typesize*($start-$num)}]
	set binfile [dict get $bcol binfile]
	# get data from already opened file
	set f [dict get $bcol fi]
	seek $f $pos
	fconfigure $f -encoding binary -translation binary
	set b [read $f [expr {$typesize*$len}]]
	binary scan $b $type$len sresult
	foreach v $sresult {
		catch {set offset [dict get $tabled $start]}
		incr start
		lappend result [expr {$v + $offset}]
	}
	if {$uend < $end} {lappend result {*}[list_fill [expr {$end-$uend}] $default]}
	return $result
}

# Warning:
# while cg table get follows the half-open convention, end position is not included
# the procedure bcol_table follows the Tcl lrange convention: end position is included!
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
	set len [expr {$end-$start+1}]
	if {$len <= 0} {return {}}
	if {$start > $max} {return [list_fill $len $default]}
	if {$end > $max} {set uend $max} else {set uend $end}
	set offset 0
	set curpos $start
	set o stdout
	set name [file root [file tail [dict get $bcol file]]]
	puts $o "pos\tvalue"
	if {[llength $table] > 2} {
		list_foreach {num type noffset} $table {
			if {$start <= $num} {
				break
			}
			set offset $noffset
		}
	} else {
		set num [lindex $table 0 0]
	}
	while {$start < $num} {
		if {$start > $end} {return}
		puts $o $curpos\t$default
		incr start
		incr curpos
	}
	set len [expr {$uend-$start+1}]
	if {$len <= 0} {return {}}
	set pos [expr {$typesize*($start-$num)}]
	# get data from already opened file
	set f [dict get $bcol fi]
	seek $f $pos
	fconfigure $f -encoding binary -translation binary
	while {$curpos <= $end} {
		catch {set offset [dict get $tabled $start]}
		incr start
		set b [read $f $typesize]
		binary scan $b $type value
		set value [expr {$value + $offset}]
		puts $o $curpos\t$value
		incr curpos
	}
	while {$uend < $end} {
		puts $o $curpos\t$default
		incr curpos
		incr uend
	}
}

array set bcol_typea {c,mn -127 s,mn -32768 i,mn -2147483648 w,mn -9223372036854775808 cu,mn 0 su,mn 0 iu,mn 0}
array set bcol_typea {c,mx  127 s,mx  32767 i,mx  2147483647 w,mx  9223372036854775807 cu,mx 255 su,mx 65535 iu,mx 4294967295}

proc cg_bcol_make {args} {
	global bcol_typea
	set type iu
	set chromosomecol coverage
	set chrompos -1
	set offsetcol {}
	set defaultvalue 0
	set distribute 0
	set start 0
	set pos 0
	set header 1
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
			-h - --header {
				set header $value
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
	set prefix [file normalize $prefix]
	set tail [file tail $prefix]
	file mkdir $prefix.bcol.temp
	if {[info exists bcol_typea($type,mx)]} {
		set max $bcol_typea($type,mx)
		set min $bcol_typea($type,mn)
	} elseif {![inlist {c s i w cu su iu wu f d} $type]} {
		exiterror "error: type $type not supported (must be one of: c s i w cu su iu wu f d)"
	}
	if {[info exists max]} {
		if {$defaultvalue > $max} {
			exiterror "default value $v too large for type $type"
		} elseif {$defaultvalue < $min} {
			exiterror "default value $v too small for type $type"
		}
	}
	set btype [string index $type 0]
	set f stdin
	if {$header} {
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
	} else {
		set offsetpos $offsetcol
		set colpos $valuecolumn
	}
# puts "bcol_make [list $prefix.bcol.temp/$tail] $type $colpos $chrompos $offsetpos $defaultvalue"
	set pipe [open "| bcol_make [list $prefix.bcol.temp/$tail] $type $colpos $chrompos $offsetpos $defaultvalue >@ stdout 2>@ stderr" w]
	fconfigure $f -encoding binary -translation binary
	fconfigure $pipe -encoding binary -translation binary
	fcopy $f $pipe
	if {[catch {close $pipe} err]} {
		regsub {child process exited abnormally} $err {} err
		puts stderr $err
		exit 1
	}
	foreach file [glob $prefix.bcol.temp/${tail}*.bin] {
		exec razip -c $file > $file.rz
		file delete $file
	}
	file rename {*}[glob $prefix.bcol.temp/*] [file dir $prefix]
	file delete -force $prefix.bcol.temp
}

proc cg_bcol_get {args} {
	if {[llength $args] != 3} {
		exiterror "wrong # args: should be \"cg bcol get file start end\""
	}
	foreach {indexfile begin end} $args break
	set bcol [bcol_open $indexfile]
	if {[isint $end]} {incr end -1}
	set result [bcol_get $bcol $begin $end]
	bcol_close $bcol
	puts $result
	return $result
}

proc cg_bcol_table {args} {
	if {[llength $args] < 1 || [llength $args] > 3} {
		exiterror "wrong # args: should be \"cg bcol table file ?start? ?end?\""
	}
	foreach {indexfile begin end} $args break
	set bcol [bcol_open $indexfile]
	if {[isint $end]} {incr end -1}
	set result [bcol_table $bcol $begin $end]
	bcol_close $bcol
	puts $result
	return $result
}

proc cg_bcol_size {args} {
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

proc cg_bcol_histo {args} {
	set pos 0
	set namecol name
	foreach {key value} $args {
		switch -- $key {
			-n - --namecol {
				set namecol $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {[llength $args] != 3} {
		exiterror "wrong # args: should be \"cg bcol histo regionfile bcolfile/prefix intervals\""
	}
	foreach {regionfile bcolfile intervals} $args break
	set f [gzopen $regionfile]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	set namepos [lsearch $header $namecol]
	lappend poss $namepos
	if {![file exists $bcolfile]} {
		set prefix 1
		set pchr {}
	} else {
		set prefix 0
		set bcol [bcol_open $bcolfile]
	}
	set biv <[lindex $intervals 0]
	set iv $biv
	set header [list name]
	set p {}
	set tota($iv) 0
	foreach limit $intervals {
		set tota($limit) 0
	}
	set a($iv) 0
	foreach limit $intervals {
		set a($limit) 0
		lappend header r$p<$limit
		set p $limit
	}
	lappend header r${limit}<
	puts [join $header \t]
	set line [getline $f]
	set prevname [lindex $line 3]
	while {1} {
		foreach {chr begin end name} [list_sub $line $poss] break
		if {$prefix && $chr ne $pchr} {
			catch {bcol_close $bcol}
			set file [lindex [gzfile $bcolfile-$chr-*.bcol $bcolfile-chr$chr-*.bcol] 0]
			set bcol [bcol_open $file]
		}
		if {[eof $f] || $name ne $prevname} {
			incr tota($biv) $a($biv)
			set result [list $a($biv)]
			set iv $biv
			set a($iv) 0
			foreach limit $intervals {
				incr tota($limit) $a($limit)
				lappend result $a($limit)
				set a($limit) 0
			}
			puts $prevname\t[join $result \t]
			set prevname $name
			if {[eof $f]} break
		}
		incr end -1
		set data [bcol_get $bcol $begin $end]
		foreach v $data {
			set iv $biv
			foreach limit $intervals {
				if {$v < $limit} break
				set iv $limit
			}
			incr a($iv)
		}
		set line [getline $f]
	}
	close $f
	bcol_close $bcol
	puts ----------
	set result [list $tota($biv)]
	foreach limit $intervals {
		lappend result $tota($limit)
	}
	puts Total\t[join $result \t]
	set tot [lmath_sum $result]
	set presult [list [format %.2f [expr {100*$tota($biv)/$tot}]]]
	foreach limit $intervals {
		lappend presult [format %.2f [expr {100*$tota($limit)/$tot}]]
	}
	puts Totalpercent\t[join $presult \t]
}

proc cg_bcol {cmd args} {
	if {![llength [info commands cg_bcol_$cmd]]} {
		set list [info commands cg_bcol_*]
		set list [list_regsub ^cg_bcol_ $list {}]
		exiterror "cg bcol has no subcommand $cmd, must be one of: [join $list ,]"
	}
	cg_bcol_$cmd {*}$args
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
		catch {file delete $indexdir/info.tsv}
		set f [gzopen $file]
		set header [tsv_open $f]
		catch {close $f}
		set bcol [bcol_open $indexfile]
		set size [bcol_size $bcol]
		bcol_close $bcol
		set f [open $indexdir/info.tsv.temp w]
		puts $f key\tvalue
		puts $f file\t$file
		puts $f lineindexfile\t[file tail $indexfile]
		puts $f header\t$header
		puts $f size\t$size
		close $f
		file rename $indexdir/info.tsv.temp $indexdir/info.tsv
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

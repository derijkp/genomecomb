proc bcol_progress {args} {
puts "progress: $args"
	set pos [lindex [lindex $args end] end]
	if {![isint $pos]} {
		progress message $args
		return
	}
	if {[catch {
		progress set $pos
	}]} {
		global bgexechandle
		puts error
		Extral::bgexec_cancel $bgexechandle
	}
}

proc bcol_indexlines {file indexfile {colinfo 0}} {
	set file [file_absolute $file]
	set indexfile [file_absolute $indexfile]
	set time [file mtime $file]
	set ext [file extension $file]
	set compressed [gziscompressed $file]
	if {![file exists $indexfile] || [file mtime $indexfile] < $time} {
		file mkdir [file dir $indexfile]
		if {$compressed} {
			progress start 1 "uncompressing $file for indexing, please be patient"
			progress message "uncompressing $file for indexing, please be patient (no progress shown)"
			set file [gztemp $file]
			progress stop
		}
		set f [open $file]
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
			progress start [file size $file] "Indexing $file, please be patient"
			progress message "Indexing $file, please be patient"
			# putslog "bcol_indexfile $file $indexfile.temp $indexfile.bin.temp {*}$poss"
			if {$colinfo} {
				set indexdir [file dir $indexfile]
				file mkdir $indexdir/colinfo.temp/
				# exec bcol_indexfile_all $file $indexfile.temp $indexfile.bin.temp {*}$poss $indexdir/colinfo.temp/ $header 2>@ stderr
				if {[catch {
					Extral::bgexec -progresscommand bcol_progress -no_error_redir -channelvar bgexechandle \
						bcol_indexfile_all $file $indexfile.temp $indexfile.bin.temp {*}$poss $indexdir/colinfo.temp/ $header 2>@1
					catch {file delete -force $indexdir/colinfo}
					file rename -force $indexdir/colinfo.temp/ $indexdir/colinfo
				}]} {
					Extral::bgexec -progresscommand bcol_progress -no_error_redir -channelvar bgexechandle \
						bcol_indexfile $file $indexfile.temp $indexfile.bin.temp {*}$poss 2>@1
				}
			} else {
				Extral::bgexec -progresscommand bcol_progress -no_error_redir -channelvar bgexechandle \
					bcol_indexfile $file $indexfile.temp $indexfile.bin.temp {*}$poss 2>@1
			}
			file rename -force $indexfile.bin.temp $indexfile.bin
			file rename -force $indexfile.temp $indexfile
			progress stop
		}
		if {$compressed} {
			gzrmtemp $file
		}
	}
}

# if ra == 1 , a compressed bin file will not be unzipped and opened (fi == "")
proc bcol_open {indexfile {ra 0}} {
	set result [dict create objtype bcol file $indexfile binfile $indexfile.bin compressedbin 0]
	set f [open $indexfile]
	set header [tsv_open $f comment]
	set comment [split $comment \n]
	if {[lindex $comment 0] ne "# binary column"} {
		close $f
		error "file \"$indexfile\" is not a binary column file"
	}
	# offset support was removed (complicates things, and was never actually used)
	if {$header eq "begin type offset"} {
		set version 0
	} elseif {$header eq "chromosome begin end"} {
		set version 1
	} else {
		close $f
		error "file \"$indexfile\" has an incorrect table header"
	}
	dict set result version $version
	foreach line $comment {
		set line [string range $line 1 end]
		dict set result [lindex $line 0] [lindex $line 1]
	}
	if {[catch {dict get $bcol multi} multi]} {
		set multi {}
	}
	set multilist [split $multi ,]
	set multilen [llength $multilist]
	dict set bcol multilen $multilen
	set table {}
	while {![eof $f]} {
		set line [gets $f]
		if {![llength $line]} continue
		lappend table [split $line \t]
	}
	close $f
	if {$version == 0} {
		set begin [lindex $table 0 0]
		set end [lindex $table end 0]
		# end is not included for half open
		incr end
		set table [list [list {} $begin $end]]
	}
	set newtable {}
	set totalnum 0
	list_foreach {chromosome begin end} $table {
		dict set tablechr $chromosome [list $begin $end]
		lappend newtable [list $chromosome $begin $end $totalnum]
		set totalnum [expr {$totalnum + $end - $begin}]
	}
	set table $newtable
	dict set result table $table
	dict set result totalnum $totalnum
	set tablechr [dict create]
	foreach line $table {
		set chromosome [lindex $line 0]
		dict set tablechr $chromosome $line
	}
	dict set result tablechr $tablechr
	set binfile [lindex [glob -nocomplain $indexfile.bin [file root $indexfile].bin $indexfile.bin.lz4 [file root $indexfile].bin.lz4 $indexfile.bin.rz [file root $indexfile].bin.rz] 0]
	if {$binfile eq ""} {set binfile [lindex [gzfiles $indexfile.bin [file root $indexfile].bin] 0]}
	if {$binfile eq ""} {exiterror "binfile $indexfile.bin not found"}
	dict set result binfile $binfile
	set ext [file extension $binfile]
	if {![gziscompressed $binfile]} {
		set fi [open $binfile]
		fconfigure $fi -encoding binary -translation binary
		dict set result fi $fi
	} elseif {$ra} {
		dict set result fi {}
	} else {
		set tempfile [gztemp $binfile]
		dict set result tempfile $tempfile
		set fi [open $tempfile]
		fconfigure $fi -encoding binary -translation binary
		dict set result fi $fi
		dict set result compressedbin 1
	}
	return $result
}

proc bcol_size bcol {
	dict get $bcol totalnum
}

proc bcol_first bcol {
	lindex [dict get $bcol table] 0 1
}

proc bcol_last bcol {
	lindex [dict get $bcol table] 0 2
}

proc bcol_close bcol {
	set fi [dict get $bcol fi]
	if {$fi ne ""} {
		catch {close $fi}
	}
	dict set bcol fi {}
	if {[dict exists $bcol tempfile]} {
		file delete [dict get $bcol tempfile]
	}
}

# Warning:
# while cg bcol get follows the half-open convention, end position is not included
# the procedure bcol_get follows the Tcl lrange convention: end position is included!
# chromosome "" will take the empty chromosome, or else the first in the list
proc bcol_get {bcol start {end {}} {chromosome {}}} {
	if {[dict get $bcol objtype] ne "bcol"} {error "This is not a bcol object: $bcol"}
	if {$end eq ""} {
		set end $start
	}
	set type [dict get $bcol type]
	if {$type eq "lineindex"} {set type iu}
	set btype [string index $type 0]
	array set typesizea {c 1 s 2 i 4 w 8 f 4 d 8}
	set typesize $typesizea($btype)
	set chrline [lindex [bcol_chrlines $bcol $chromosome] 0]
	foreach {chr min max chrstart} $chrline break
	# chrline is in half open, we are using (Tcl style) closed interval here
	incr max -1
	if {[catch {dict get $bcol default} default]} {
		set default 0
	}
	set len [expr {$end-$start+1}]
	if {$len <= 0} {return {}}
	if {$start > $max} {return [list_fill $len $default]}
	if {$end > $max} {set uend $max} else {set uend $end}
	set result {}
	if {$start < $min} {
		if {$end < $min} {
			return [list_fill [expr {$end-$start+1}] $default]
		} else {
			set result [list_fill [expr {$min-$start}] $default]
		}
		set start $min
	}
	set len [expr {$uend-$start+1}]
	if {$len <= 0} {return $result}
	set pos [expr {$typesize*($chrstart + $start-$min)}]
	set binfile [dict get $bcol binfile]
	# get data from already opened file
	set f [dict get $bcol fi]
	if {$f ne ""} {
		seek $f $pos
		fconfigure $f -encoding binary -translation binary
		set b [read $f [expr {$typesize*$len}]]
	} else {
		set tf [gzopen [gzfile $binfile] $pos]
		fconfigure $tf -encoding binary -translation binary
		set b [read $tf [expr {$typesize*$len}]]
		catch {close $tf}
	}
	binary scan $b $type$len sresult
	if {![info exists sresult]} {
		error "error getting data from $binfile for position $start-$end"
	}
	foreach v $sresult {
		incr start
		lappend result $v
	}
	if {$uend < $end} {lappend result {*}[list_fill [expr {$end-$uend}] $default]}
	return $result
}

proc bcol_chrlines {bcol chromosome} {
	set tablechr [dict get $bcol tablechr]
	if {$chromosome eq "all"} {
		set chrlines [dict get $bcol table]
	} elseif {[dict exists $tablechr $chromosome]} {
		set chrlines [list [dict get $bcol tablechr $chromosome]]
	} elseif {$chromosome eq ""} {
		set table [dict get $bcol table]
		set chrlines [list [lindex $table 0]]
	} else {
		error "chromosome \"$chromosome\" not found in bcol file [dict get $bcol file]"
	}
	return $chrlines
}

# bcol_table now also follows the half-open convention, end position is not included
# offset support was removed (complicates things, and was never actually used)
proc bcol_table {bcol {start {}} {end {}} {chromosome {}} {showchr 1} {byrownum 0} {precision {}}} {
	if {[dict get $bcol objtype] ne "bcol"} {error "This is not a bcol object: $bcol"}
	set type [dict get $bcol type]
	if {[catch {dict get $bcol multi} multi]} {
		set multi {}
	}
	set multilist [split $multi ,]
	set multilen [llength $multilist]
	if {$type eq "lineindex"} {set type iu}
	set btype [string index $type 0]
	array set typesizea {c 1 s 2 i 4 w 8 f 4 d 8}
	set typesize $typesizea($btype)
	set table [dict get $bcol table]
	set tablechr [dict get $bcol tablechr]
	set chrlines [bcol_chrlines $bcol $chromosome]
	set o stdout
	if {!$multilen} {
		if {$showchr} {
			puts $o "chromosome\tpos\tvalue"
		} else {
			puts $o "pos\tvalue"
		}
	} else {
		if {$showchr} {
			puts $o "chromosome\tpos\talt\tvalue"
		} else {
			puts $o "pos\talt\tvalue"
		}
	}
	set ostart $start
	set oend $end
	foreach chrline $chrlines {
		foreach {chr min max chrstart} $chrline break
		if {[catch {dict get $bcol default} default]} {
			set default 0
		}
		if {$start eq ""} {set start $min}
		if {$end eq ""} {set end $max}
		set len [expr {$end-$start}]
		if {$len <= 0} {return {}}
		if {$end > $max} {set uend $max} else {set uend $end}
		set curpos $start
		while {$start < $min} {
			if {$start >= $end} {return}
			if {$showchr} {
				puts $o $chr\t$curpos\t$default
			} {
				puts $o $curpos\t$default
			}
			incr start
			incr curpos
		}
		set len [expr {$uend-$start}]
		if {$len <= 0} {return {}}
		if {!$multilen} {
			set pos [expr {$typesize*($chrstart + $start-$min)}]
		} else {
			set pos [expr {$typesize*$multilen*($chrstart + $start-$min)}]
		}
		# get data from already opened file
		set f [dict get $bcol fi]
		if {$f ne ""} {
			seek $f $pos
		} else {
			set f [gzopen [gzfile $binfile] $pos]
		}
		fconfigure $f -encoding binary -translation binary
		if {!$multilen} {
			while {$curpos < $uend} {
				incr start
				set b [read $f $typesize]
				binary scan $b $type value
				if {$showchr} {
					puts $o $chr\t$curpos\t$value
				} else {
					puts $o $curpos\t$value
				}
				incr curpos
			}
		} else {
			set readsize [expr {$typesize*$multilen}]
			while {$curpos < $uend} {
					incr start
					set b [read $f $readsize]
					binary scan $b $type$multilen value
					if {$precision ne ""} {
						set tempvalue {}
						foreach v $value {
							lappend tempvalue [format %.${precision}f $v $precision]
						}
						set value $tempvalue
					}
					if {$showchr} {
						puts $o $chr\t$curpos\t$multi\t[join $value ,]
					} else {
						puts $o $curpos\t$multi\t[join $value ,]
					}
					incr curpos
				}
		}
		if {[dict get $bcol fi] eq ""} {
			close $f
		}
		while {$curpos < $end} {
			puts $o $curpos\t$default
			incr curpos
			incr uend
		}
		set start $ostart
		set end $oend
	}
}

array set bcol_typea {c,mn -127 s,mn -32768 i,mn -2147483648 w,mn -9223372036854775808 cu,mn 0 su,mn 0 iu,mn 0}
array set bcol_typea {c,mx  127 s,mx  32767 i,mx  2147483647 w,mx  9223372036854775807 cu,mx 255 su,mx 65535 iu,mx 4294967295}

proc cg_bcol_update {newbcol oldbcol args} {
	set newbinfile $newbcol.bin
	set bcol [bcol_open $oldbcol 1]
	set type [dict get $bcol type]
	set default [dict get $bcol default]
	bcol_close $bcol
	set o [open $newbcol w]
	puts $o "\# binary column"
	puts $o "\# type $type"
	puts $o "\# default $default"
	puts $o "chromosome\tbegin\tend"
	set args [list $oldbcol {*}$args]
	foreach file $args {
		set bcol [bcol_open $file 1]
		set table [dict get $bcol table]
		set version [dict get $bcol version]
		set binfile [dict get $bcol binfile]
		exec {*}[gzcat $binfile] $binfile >> $newbinfile
		if {$version == 0} {
			set line [lindex $table 0]
			if {[regexp -- {-c?h?r?([0-9XYM]+)-[^-]+.bcol$} [gzfile $file] temp chr]
			||[regexp -- {-([^-]+).bcol$} [gzfile $file] temp chr]} {
				lset line 0 [chr_clip $chr]
			}
			puts $o [join $line \t]
		} else {
			foreach line $table {
				puts $o [join $line \t]
			}
		}
	}
	cg lz4 -c 9 $newbinfile
}

proc cg_bcol_make {args} {
	global bcol_typea
	set type iu
	set compress 9
	set chromosomecol coverage
	set chrompos -1
	set offsetcol {}
	set defaultvalue 0
	set multicol {}
	set multilist {}
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
			-c - --chromosomecol {
				set chromosomecol $value
				set distribute 1
			}
			-co - --compress {
				set compress $value
			}
			-h - --header {
				set header $value
			}
			-m - --multicol {
				set multicol $value
			}
			-l - --multilist {
				set multilist $value
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
		errorformat bcol_make
		exit 1
	}
	foreach {bcolfile valuecolumn} $args break
	set bcolfile [file_absolute $bcolfile]
	set tail [file tail $bcolfile]
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
			if {[isint $chromosomecol]} {
				set chrompos $chromosomecol
			} else {
				set chrompos [lsearch $header $chromosomecol]
			}
			if {$chrompos == -1} {
				exiterror "error: chromosome column $chromosomecol not found"
			}
		}
		if {$multicol ne ""} {
			set multipos [lsearch $header $multicol]
			if {$multipos == -1} {
				exiterror "error: multicolumn $multicol not found"
			}
		}
	} else {
		if {$distribute && [isint $chromosomecol]} {
			set chrompos $chromosomecol
		}
		set offsetpos $offsetcol
		set colpos $valuecolumn
		set multipos $multicol
	}
	if {$compress} {
		set compresspipe "| lz4c -$compress -c "
		set tempbinfile $bcolfile.temp.bin.lz4
	} else {
		set compresspipe ""
		set tempbinfile $bcolfile.temp.bin
	}
	if {$multicol eq ""} {
		# puts "bcol_make $bcolfile.temp $type $colpos $chrompos $offsetpos $defaultvalue"
		set pipe [open "| bcol_make [list $bcolfile.temp] $type $colpos $chrompos $offsetpos $defaultvalue $compresspipe > [list $tempbinfile] 2>@ stderr" w]
	} else {
		# putsvars bcolfile type colpos multipos multilist chrompos offsetpos defaultvalue
		# puts "bcol_make_multi $bcolfile.temp $type $multipos $multilist $colpos $chrompos $offsetpos $defaultvalue"
		set pipe [open "| bcol_make_multi [list $bcolfile.temp] $type $multipos $multilist $colpos $chrompos $offsetpos $defaultvalue $compresspipe > [list $tempbinfile] 2>@ stderr" w]
	}
	fconfigure $f -encoding binary -translation binary
	fconfigure $pipe -encoding binary -translation binary
	fcopy $f $pipe
	if {[catch {close $pipe} err]} {
		regsub {child process exited abnormally} $err {} err
		puts stderr $err
		exit 1
	}
	if {$compress} {
		file rename -force $tempbinfile $bcolfile.bin.lz4
	} else {
		file rename -force $tempbinfile $bcolfile.bin
	}
	file rename -force $bcolfile.temp $bcolfile
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
	set chromosome all
	set showchr 1
	set byrownum 0
	set pos 0
	set precision {}
	foreach {key value} $args {
		switch -- $key {
			-c - --chromosome {
				set chromosome $value
			}
			-s - --showchromosome {
				set showchr $value
			}
			-r - --byrownum {
				set byrownum $value
			}
			-p - --precision {
				set precision $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {[llength $args] < 1 || [llength $args] > 3} {
		exiterror "wrong # args: should be \"cg bcol table ?-c chromosome? ?-s showchr? file ?start? ?end?\""
	}
	foreach {indexfile begin end} $args break
	set bcol [bcol_open $indexfile]
	bcol_table $bcol $begin $end $chromosome $showchr $byrownum $precision
	bcol_close $bcol
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
		set bcol [bcol_open $bcolfile 1]
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
	set totsize 0
	set totsum 0
	set totmin {}
	set totmax {}
	set size 0
	set sum 0
	set mins {}
	set maxs {}
	lappend header size avg min max
	puts [join $header \t]
	set line [getline $f]
	set prevname [lindex $line 3]
	while {1} {
		foreach {chr begin end name} [list_sub $line $poss] break
		set chr [chr_clip $chr]
		if {$prefix && $chr ne $pchr} {
			catch {bcol_close $bcol}
			set file [lindex [gzfile $bcolfile-chr$chr-*.bcol $bcolfile-chr$chr.bcol $bcolfile-$chr-*.bcol $bcolfile-$chr.bcol] 0]
			set bcol [bcol_open $file 1]
			set getchr {}
		} else {
			set getchr $chr
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
			set min [lmath_min $mins]
			set max [lmath_max $maxs]
			lappend result $size
			if {$size == 0} {
				lappend result ""
			} else {
				lappend result [format %.2f [expr {$sum/double($size)}]]
			}
			lappend result $min
			lappend result $max
			puts $prevname\t[join $result \t]
			set prevname $name
			incr totsize $size
			incr totsum $sum
			if {$totmin eq "" || $min < $totmin} {set totmin $min}
			if {$totmax eq "" || $max > $totmax} {set totmax $max}
			set size 0
			set sum 0
			set mins {}
			set maxs {}
			if {[eof $f]} break
		}
		incr end -1
		set data [bcol_get $bcol $begin $end $getchr]
		foreach v $data {
			set iv $biv
			foreach limit $intervals {
				if {$v < $limit} break
				set iv $limit
			}
			incr a($iv)
		}
		incr size [llength $data]
		set sum [expr {$sum + round([lmath_sum $data])}]
		lappend mins [lmath_min $data]
		lappend maxs [lmath_max $data]
		set line [getline $f]
	}
	close $f
	bcol_close $bcol
	puts ----------
	set result [list $tota($biv)]
	foreach limit $intervals {
		lappend result $tota($limit)
	}
	set tot [lmath_sum $result]
	set presult [list [format %.2f [expr {100*$tota($biv)/$tot}]]]
	foreach limit $intervals {
		lappend presult [format %.2f [expr {100*$tota($limit)/$tot}]]
	}
	lappend result $totsize [format %.2f [expr {$totsum/double($totsize)}]] $totmin $totmax
	puts Total\t[join $result \t]
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

proc cg_size file {
	set indexfile [cg_index $file]
	set bcol [bcol_open $indexfile 1]
	set size [bcol_size $bcol]
	bcol_close $bcol
	puts $size
	return $size
}
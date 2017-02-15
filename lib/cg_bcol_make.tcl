array set bcol_typea {c,mn -127 s,mn -32768 i,mn -2147483648 w,mn -9223372036854775808 cu,mn 0 su,mn 0 iu,mn 0}
array set bcol_typea {c,mx  127 s,mx  32767 i,mx  2147483647 w,mx  9223372036854775807 cu,mx 255 su,mx 65535 iu,mx 4294967295}

proc cg_bcol_make {args} {
	global bcol_typea
	set type iu
	set compress 9
	set chrompos -1
	set chromosomename {}
	set defaultvalue 0
	set multicol {}
	set multilist {}
	set start 0
	set precision -1
	set header 1
	set srcfile {}
	cg_options bcol_make args {
		-t - --type {
			set type $value
		}
		-p - --poscol {
			set offsetcol $value
		}
		-e - --endcol {
			set endcol $value
		}
		-d - --default {
			set defaultvalue $value
		}
		-c - --chromosomecol {
			set chromosomecol $value
		}
		-n - --chromosomename {
			set chromosomename $value
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
		--precision {
			set precision $value
		}
	} {bcolfile valuecolumn srcfile} 2 3
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
	if {$srcfile eq ""} {
		set f stdin
	} else {
		set f [open $srcfile]
	}
	if {$header} {
		set header [tsv_open $f comment]
		set poss [tsv_basicfields $header 3 0]
		set colpos [lsearch $header $valuecolumn]
		if {$colpos == -1} {
			exiterror "error: valuecolumn $valuecolumn not found"
		}
		if {![info exists offsetcol]} {
			set offsetpos [lindex $poss 1]
		} elseif {$offsetcol eq ""} {
			set offsetpos -1
		} else {
			set offsetpos [lsearch $header $offsetcol]
			if {$offsetpos == -1} {
				exiterror "error: pos column $offsetcol not found"
			}
		}
		if {![info exists endcol]} {
			set endpos [lindex $poss 2]
		} elseif {$endcol eq ""} {
			set endpos -1
		} else {
			if {$offsetpos == -1} {error "Cannit use -e (endcol) without -p (poscol)"}
			if {$multicol ne ""} {error "Cannit use -e (endcol) with -m (multicol)"}
			set endpos [lsearch $header $endcol]
			if {$endpos == -1} {
				exiterror "error: pos column $endcol not found"
			}
		}
		if {![info exists chromosomecol]} {
			set chrompos [lindex $poss 0]
		} elseif {$chromosomecol ne ""} {
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
		if {[isint $chromosomecol]} {
			set chrompos $chromosomecol
		} else {
			set chromosomecol {}
		}
		if {![info exists offsetcol]} {set offsetcol {}}
		if {![info exists endcol]} {set endcol {}}
		set offsetpos $offsetcol
		set colpos $valuecolumn
		set multipos $multicol
		if {$endpos != -1} {
			if {$offsetpos == -1} {error "Cannot use -e (endcol) without -p (poscol)"}
			if {$multicol ne ""} {error "Cannot use -e (endcol) with -m (multicol)"}
		}
		set endpos $endcol
	}
	set tempfile [filetemp $bcolfile]
	if {$compress} {
		set compresspipe "| lz4c -$compress -c "
		set tempbinfile $tempfile.bin.lz4
	} else {
		set compresspipe ""
		set tempbinfile $tempfile.bin
	}
	if {$multicol eq ""} {
		# puts "bcol_make [list $tempfile] $type $colpos $chrompos [list $chromosomename] $offsetpos $endpos [list $defaultvalue] $precision"
		set pipe [open "| bcol_make [list $tempfile] $type $colpos $chrompos [list $chromosomename] $offsetpos $endpos [list $defaultvalue] $precision $compresspipe > [list $tempbinfile] 2>@ stderr" w]
	} else {
		# putsvars bcolfile type colpos multipos multilist chrompos offsetpos defaultvalue
		# puts "bcol_make_multi $tempfile $type $multipos $multilist $colpos $chrompos [list $chromosomename] $offsetpos $endpos [list $defaultvalue] $precision"
		set pipe [open "| bcol_make_multi [list $tempfile] $type $multipos $multilist $colpos $chrompos [list $chromosomename] $offsetpos $endpos [list $defaultvalue] $precision $compresspipe > [list $tempbinfile] 2>@ stderr" w]
	}
	fconfigure $f -encoding binary -translation binary
	fconfigure $pipe -encoding binary -translation binary
	fcopy $f $pipe
	if {[catch {close $pipe} err]} {
		regsub {child process exited abnormally} $err {} err
		error $err
	}
	if {$srcfile ne ""} {
		close $f
	}
	if {$compress} {
		exec lz4index $tempbinfile
		file rename -force $tempbinfile $bcolfile.bin.lz4
		file rename -force $tempbinfile.lz4i $bcolfile.bin.lz4.lz4i
	} else {
		file rename -force $tempbinfile $bcolfile.bin
	}
	file rename -force $tempfile $bcolfile
}

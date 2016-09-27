proc cg_cat {args} {
	set fields {}
	set force 0
	set merge 0
	set addcomment 1
	set namefield {}
	set sort 0
	cg_options cat args {
		-f - --force {
			if {$value in "1 0"} {
				set force $value
			} else {
				set force 1
				incr pos -1
			}
		}
		-m - --merge {
			if {$value in "1 0"} {
				set merge $value
			} else {
				set merge 1
				incr pos -1
			}
		}
		-s - --sort {
			if {$value in "1 0"} {
				set sort $value
			} else {
				set sort 1
				incr pos -1
			}
		}
		-c - --comments {
			set addcomment $value
		}
		-n - --fieldname {
			set namefield $value
		}
		-fields - --fields {
			set fields $value
		}
	} 1
	if {$merge} {set force m} elseif {$force} {set force f} else {set force ""}
	if {[llength $args] == 1} {
		set f [gzopen [lindex $args 0]]
		fcopy $f stdout
		exit 0
	}
	if {$sort} {set args [lsort -dict $args]}
	set headers {}
	set comments {}
	foreach file $args {
		set f [gzopen $file]
		set header [tsv_open $f comment]
		if {-1 in [list_cor $header $fields]} {
			error "$file does not have all fields: $fields"
		}
		lappend headers $header
		gzclose $f
		if {$addcomment eq "1"} {
			lappend comments "# ++++ $file ++++"
			if {[inlist {f m} $force]} {
				lappend comments "# ++ $header"
			}
		}
		if {$addcomment in {f 0 1}} {
			foreach line [split $comment \n] {
				if {[string index $line 0] ne "#"} continue
				lappend comments $line
			}
			if {$addcomment eq "f"} {set addcomment n}
		}
	}
	if {$fields ne ""} {
		set header $fields
	} else {
		set header [lindex $headers 0]
		set hlen [llength $header]
		foreach testheader $headers {
			set cor [list_cor $header $testheader]
			if {[llength $testheader] != $hlen} {
				if {$force eq ""} {
					puts stderr "headers do not match, use -f to force or -m to merge"
					exit 1
				}
			}
			set poss [list_find $cor -1]
			if {[llength $poss]} {
				if {$force eq ""} {
					puts stderr "headers do not match, use -f to force or -m to merge"
					exit 1
				}
				lappend header {*}[list_sub $testheader $poss]
			}
		}
	}
	if {$namefield ne ""} {
		lappend header $namefield
	}
	if {[llength $comments]} {
		puts [join $comments \n]
	}
	puts [join $header \t]
	set defcor [list_cor $header $header]
	set lh [llength $header]
	foreach testheader $headers file $args {
		set fname [file tail $file]
		set f [gzopen $file]
		set testheader [tsv_open $f]
		if {($force eq "f" || $force eq "") && $fields eq ""} {
			if {$namefield ne ""} {
				while {![eof $f]} {
					set line [gets $f]
					if {$line eq ""} continue
					puts $line\t$fname
				}
			} else {
				fcopy $f stdout
			}
		} else {
			# merge
			set cor [list_cor $testheader $header]
			if {$cor eq $defcor && [llength $testheader] == $lh} {
				fcopy $f stdout
			} else {
				while {![eof $f]} {
					set line [split [gets $f] \t]
					if {![llength $line]} continue
					set line [list_sub $line $cor]
					if {$namefield ne ""} {lappend line $fname}
					puts [join $line \t]
				}
			}
		}
		gzclose $f
	}
}

if 0 {
	set args {data/reg1.tsv data/reg2.tsv data/reg3.tsv data/reg4.tsv}
}

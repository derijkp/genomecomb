proc cg_cat {args} {
	if {[lindex $args 0] eq "-f"} {
		set force f
		set args [lrange $args 1 end]
	} elseif {[lindex $args 0] eq "-m"} {
		set force m
		set args [lrange $args 1 end]
	} else {
		set force ""
	}
	if {[llength $args] == 0} {
		errorformat cat
		exit 1
	}
	if {[llength $args] == 1} {
		set f [gzopen [lindex $args 0]]
		fcopy $f stdout
		exit 0
	}
	set headers {}
	set comments {}
	foreach file $args {
		set f [gzopen $file]
		set header [tsv_open $f comment]
		lappend headers $header
		close $f
		lappend comments "# ++++ $file ++++"
		if {$force eq "f"} {
			lappend comments "# ++ $header"
		}
		foreach line [split $comment \n] {
			if {[string index $line 0] ne "#"} continue
			lappend comments $line
		}
	}
	set header [lindex $headers 0]
	foreach testheader $headers {
		set cor [list_cor $header $testheader]
		set poss [list_find $cor -1]
		if {[llength $poss]} {
			if {$force eq ""} {
				puts stderr "headers do not match, use -f to force or -m to merge"
				exit 1
			}
			lappend header {*}[list_sub $testheader $poss]
		}
	}
	if {[llength $comments]} {
		puts [join $comments \n]
	}
	puts [join $header \t]
	foreach testheader $headers file $args {
		set f [gzopen $file]
		tsv_open $f
		if {$force eq "f" || $force eq ""} {
			fcopy $f stdout
		} else {
			# merge
			set cor [list_cor $testheader $header]
			if {[lsearch $cor -1] == -1} {
				fcopy $f stdout
			} else {
				while {![eof $f]} {
					set line [split [gets $f] \t]
					if {![llength $line]} continue
					puts [join [list_sub $line $cor] \t]
				}
			}
		}
	}
}

if 0 {
	set args {data/reg1.tsv data/reg2.tsv data/reg3.tsv data/reg4.tsv}
}

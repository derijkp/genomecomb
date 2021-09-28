proc cg_cat {args} {
	set fields {}
	set force 0
	set merge 0
	set namefield {}
	set sort 0
	unset -nocomplain addcomment
	cg_options cat args {
		-f - -force {
			if {$value in "1 0"} {
				set force $value
			} else {
				set force 1
				incr pos -1
			}
		}
		-m - -merge {
			if {$value in "1 0"} {
				set merge $value
			} else {
				set merge 1
				incr pos -1
			}
		}
		-s - -sort {
			if {$value in "1 0"} {
				set sort $value
			} else {
				set sort 1
				incr pos -1
			}
		}
		-c - -comments {
			set addcomment $value
		}
		-n - -fieldname {
			set namefield $value
		}
		-fields {
			set fields $value
		}
	} {} 1
	if {$merge} {set force m} elseif {$force} {set force f} else {set force ""}
	if {![info exists addcomment]} {
		set addcomment m
	}
	if {[llength $args] == 1} {
		set f [gzopen [lindex $args 0]]
		fcopy $f stdout
		gzclose $f
		exit 0
	}
	if {$sort} {set args [bsort $args]}
	set headers {}
	set comments {}
	set files {}
	set diffcomments 0
	unset -nocomplain commenta
	unset -nocomplain firstcomment
	foreach file $args {
		if {[file size $file] == 0} continue
		set f [gzopen $file]
		set header [tsv_open $f comment]
		if {-1 in [list_cor $header $fields]} {
			error "$file does not have all fields: $fields"
		}
		lappend headers $header
		lappend files $file
		gzclose $f
		if {$addcomment eq "m"} {
			if {![info exists firstcomment]} {
				set firstcomment $comment
			} elseif {$comment ne $firstcomment} {
				if {!$diffcomments} {
					tsv_comment2var $firstcomment commenta
					set diffcomments 1
				}
				tsv_comment2var $comment a
				foreach key [array names a] {
					if {![info exists commenta($key)]} {
						set commenta($key) $a($key)
					} else {
						list_addnew commenta($key) {*}$a($key)
					}
				}
			}
		} else {
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
	}
	if {$fields ne ""} {
		set header $fields
	} else {
		set header [lindex $headers 0]
		set hlen [llength $header]
		foreach testheader $headers file $files {
			set cor [list_cor $header $testheader]
			if {[llength $testheader] != $hlen} {
				if {$force eq ""} {
					error "headers do not match, use -f to force or -m to merge (at file $file)"
				}
			}
			set poss [list_find $cor -1]
			if {[llength $poss]} {
				if {$force eq ""} {
					error "headers do not match, use -f to force or -m to merge (at file $file)"
				}
				lappend header {*}[list_sub $testheader $poss]
			}
		}
	}
	if {$addcomment eq "m"} {
		if {!$diffcomments} {
			puts $firstcomment
		} else {
			set commenta(catfiles) $files
			puts [tsv_var2comment commenta]
		}
	} elseif {[llength $comments]} {
		puts [join $comments \n]
	}
	set defcor [list_cor $header $header]
	set lh [llength $header]
	if {$namefield ne ""} {
		puts [join $header \t]\t$namefield
	} else {
		puts [join $header \t]
	}
	foreach testheader $headers file $args {
		if {[file size $file] == 0} continue
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

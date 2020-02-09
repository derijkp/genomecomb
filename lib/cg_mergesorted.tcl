proc cg_mergesorted {args} {
	# job_logdir
	upvar job_logdir job_logdir
	set maxopenfiles {}
	set sortfields {}
	set addcomment 1
	cg_options sam_catmerge args {
		-sortfields {
			set sortfields $value
		}
		-c - -comments {
			set addcomment $value
		}
	} {} 1 ... {
		merge sorted tsv files into a a correctly sorted resultfile
	}
	if {[llength $args] == 1} {
		set f [gzopen [lindex $args 0]]
		fcopy $f stdout
		exit 0
	}
	set headers {}
	set comments {}
	set files {}
	foreach file $args {
		if {[file size $file] == 0} continue
		set f [gzopen $file]
		set theader [tsv_open $f comment]
		if {![info exists header]} {
			set header $theader
			if {![llength $sortfields]} {
				set poss [tsv_basicfields $header 6 0]
				set poss [list_remove $poss -1]
			} else {
				set poss [list_cor $header $sortfields]
				if {-1 in $poss} {
					error "some sortfields not found: [join [list_sub $sortfields [list_find $poss -1]] ,]"
				}
			}
		} elseif {$theader ne $header} {
			error "mismatched headers in files"
		}
		gzclose $f
	}
	exec mergesorted \# 1 {} $poss {*}$args >@ stdout 2>@stderr
}

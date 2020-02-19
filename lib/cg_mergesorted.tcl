proc cg_mergesorted {args} {
	# job_logdir
	upvar job_logdir job_logdir
	set maxopenfiles {}
	set sortfields {}
	set addcomment 1
	set commentchar \#
	set headerline 1
	set header {}
	set sortposs {}
	cg_options mergesorted args {
		-sortfields {
			set sortfields $value
		}
		-commentchar {
			set commentchar $value
		}
		-headerline {
			set headerline $value
		}
		-header {
			set header $value
		}
		-sortpos {
			set sortposs $value
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
	if {!$headerline && ![llength $sortposs]} {
		error "if -headerline is 0, -sortposs must be given"
	}
	set headers {}
	set comments {}
	set files {}
	if {[llength $sortposs]} {
		set poss $sortposs
	} else {
		set poss {}
	}
	if {$headerline} {
		foreach file $args {
			if {[file size $file] == 0} continue
			set f [gzopen $file]
			set theader [tsv_open $f comment]
			if {![info exists fileheader]} {
				set fileheader $theader
				if {![llength $poss]} {
					if {![llength $sortfields]} {
						set poss [tsv_basicfields $fileheader 6 0]
						set poss [list_remove $poss -1]
					} else {
						set poss [list_cor $fileheader $sortfields]
						if {-1 in $poss} {
							error "some sortfields not found: [join [list_sub $sortfields [list_find $poss -1]] ,]"
						}
					}
				}
			} elseif {$theader ne $fileheader} {
				error "mismatched headers in files"
			}
			gzclose $f
		}
	}
	# puts stderr [list mergesorted $commentchar $headerline $header $poss {*}$args]
	exec mergesorted $commentchar $headerline $header $poss {*}$args >@ stdout 2>@stderr
}

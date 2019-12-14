#proc cg_convsolid {args} {
#	global projectfile
#	cgmake_clear
#	cgmakelib_convsolid
#	unset -nocomplain ::projectsamples
#	set args [cgmake_args {*}$args]
#	if {[llength $args] != 1} {
#		error "wrong # args: should be cg convsolid ?options? projectfile"
#	}
#	foreach {projectfile} $args break
#	set projectfile [file_absolute $projectfile]
#	cd [file dir $projectfile]
#	cgmake solid_conv
#}

proc cg_convsolid {args} {
	global projectfile
	unset -nocomplain projectsamples
	if {[llength $args] != 1} {
		error "wrong # args: should be cg convsolid ?options? projectfile"
	}
	foreach {projectfile} $args break
	set projectfile [file_absolute $projectfile]
	cd [file dir $projectfile]

	# get projectsamples from projectfile
	if {![file exists $projectfile]} {error "projectfile $projectfile does not exist"}
	set c [split [string trim [file_read $projectfile]] \n]
	set header [split [list_shift c] \t]
	set poss [list_cor $header {project experiment sample source barcode}]
	if {[lsearch $poss -1] != -1} {
		error "project-$project.txt must have columns: project experiment sample source barcode"
	}
	set c [list_subindex $c $poss]
	set list {}
	set sources {}
	list_foreach {p experiment sample source barcode} $c {
		lappend list $experiment/$sample
		set projectsamplesource($experiment/$sample) $source
		lappend sources $source
		set projectsamplebarcode($experiment/$sample) $barcode
	}
	set projectsamples $list
	set sources [list_remdup $sources]

	# split on barcodes and distribute
	foreach target1 $projectsamples {
		# split on barcodes
		set dep $projectsamplesource($target1)
		file mkdir xsq_split
		file mkdir $target1/xsq
		if {![file exists xsq_split/finished_[file tail $dep]]} {
			set files [glob $dep/result/lane*/*.xsq]
			foreach file $files {
				if {[regexp {\.xsq_} $file]} continue
				putslog "splitting $file"
				exec convertFromXSQ.sh -s $file -o xsq_split
			}
			file_write xsq_split/finished_[file tail $dep] ""
		}
		exec cp -l {*}[glob xsq_split/*_Library$projectsamplebarcode($target1).xsq] $target1/xsq
		file_write $target1/xsq/finished ""
	}
	
	foreach target1 $projectsamples {
		# convert to csfasta
		file mkdir $target1/csfasta
		set files [ssort -natural [glob $target1/xsq/*.xsq]]
		foreach file $files {
			putslog "convert to $target1/csfasta from $file"
			exec convertFromXSQ.sh -c $file -o $target1/csfasta >> log
		}
		file rename -force -- {*}[glob $target1/csfasta/*/*/*] $target1/csfasta
		file delete -force -- $target1/csfasta/Libraries
		file_write $target1/csfasta/finished ""
	}

}


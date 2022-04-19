proc project_transfer {src dest {type soft} {force 0}} {
	if {[file exists $dest]} {
		if {$force} {
			file delete $dest
		} else {
			error "$dest already exists"
		}
	}
	switch $type {
		soft {mklink $src $dest 1}
		rel {mklink $src $dest 0}
		hard {hardcopy $src $dest}
		copy {file copy $src $dest
		}
	}
}

proc cg_project_addsample {args} {
	set targetfile {}
	set amplicons {}
	set transfertype soft
	set force 0
	cg_options addsample args {
		-targetfile {
			set targetfile [file_absolute $value]
			if {$value ne "" && ![jobfileexists $targetfile]} {error "target file $targetfile does not exists"}
		}
		-amplicons {
			set amplicons [file_absolute $value]
			if {$value ne "" && ![jobfileexists $amplicons]} {error "amplicons file $amplicons does not exists"}
		}
		-transfer {
			if {$value ni {soft rel hard copy}} {error "unknown option $value for -transfer, must be one of: soft, rel hard, copy"}
			set transfertype $value
		}
		-force {
			set force [true $value]
		}
	} {projectdir samplename} 2
	if {[regexp -- - [file tail $projectdir]]} {error "[file tail $projectdir] contains a - character (which is not allowed in a projectdir)"}
	if {[regexp -- - $samplename]} {error "$samplename contains a - character (which is not allowed in a samplename)"}
	file mkdir $projectdir/samples/$samplename
	set dir [lindex $args 0]
	if {[llength $args] == 1 && [file isdir $dir] && $dir ne "ASM"} {
		project_transfer $dir $projectdir/samples/$samplename/ori $transfertype $force
	} else {
		file mkdir $projectdir/samples/$samplename/ori
		foreach file $args {
			if {![file exists $file]} {error "file $file does not exist"}
			project_transfer $file $projectdir/samples/$samplename/ori/[file tail $file] $transfertype $force
		}
	}
	if {$amplicons ne ""} {
		set temp [lindex [split [file root [gzroot [file tail $amplicons]]] -] end]
		set target $projectdir/samples/$samplename/reg_amplicons-$temp.tsv[gzext $amplicons]
		project_transfer $amplicons $target $transfertype $force
	}
	if {$targetfile ne ""} {
		set temp [file root [gzroot [file tail $targetfile]]]
		regsub {^reg_[^_]+_exome_} $temp {} temp
		set temp [lindex [split $temp -] end]
		set target $projectdir/samples/$samplename/reg_targets-$temp.tsv[gzext $targetfile]
		project_transfer $targetfile $target $transfertype $force
	}
}

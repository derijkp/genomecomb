proc project_transfer {src dest {type soft} {force 0}} {
	if {[file exists $dest] && [catch {file link $dest}]} {
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

proc cg_make_project {args} {
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
			if {$value ni {soft rel hard copy}} {error "unknown option $value for -transfer, must be one of: soft, rel, hard, copy"}
			set transfertype $value
		}
		-force {
			set force [true $value]
		}
	} {projectdir samplesheet} 2
	set projectdir [file_absolute $projectdir]
	if {[regexp -- - [file tail $projectdir]]} {error "[file tail $projectdir] contains a - character (which is not allowed in a projectdir)"}
	if {[regexp -- , [file tail $projectdir]]} {error "[file tail $projectdir] contains a , character (which is not allowed in a projectdir)"}
	mkdir $projectdir
	set f [gzopen $samplesheet]
	set header [tsv_open $f]
	set poss [list_cor $header {sample seqfiles}]
	if {[lindex $poss 1] == -1} {
		lset poss 1 [lsearch $header fastq]
		if {[lindex $poss 1] == -1} {
			lset poss 1 [lsearch $header fastqs]
		}
	}
	if {-1 in [lrange $poss 0 1]} {
		error "fields sample and seqfiles must be present in samplesheet $samplesheet"
	}
	set optionfields [list_sub $header -exclude $poss]
	unset -nocomplain optionsa
	set options {}
	mkdir $projectdir/samples
	while 1 {
		if {[gets $f line] == -1} break
		set line [split $line \t]
		foreach {samplename seqfiles} [list_sub $line $poss] break
		# if {[regexp -- - $samplename]} {error "$samplename contains a - character (which is not allowed in a samplename)"}
#		if {$preset ne ""} {
#			lappend options $samplename\tpreset\t$preset
#		}
		file mkdir $projectdir/samples/$samplename
		file mkdir $projectdir/samples/$samplename/ori
		if {[file isdir $seqfiles]} {
			set files [gzfiles $seqfiles/*\\.f*q]
			if {![llength $files]} {
				set files [gzfiles $seqfiles/*\\.bam]
			}
		} else {
			set files [gzfiles $seqfiles]
		}
		foreach file $files {
			project_transfer $file $projectdir/samples/$samplename/ori/[file tail $file] $transfertype $force
			if {[file extension [gzroot $file]] in ".fastq .fq"} {
				file mkdir $projectdir/samples/$samplename/fastq
				project_transfer $file $projectdir/samples/$samplename/fastq/[file tail $file] $transfertype $force
			}
			if {[file extension [gzroot $file]] in ".bam"} {
				file mkdir $projectdir/samples/$samplename/ubam
				project_transfer $file $projectdir/samples/$samplename/ubam/[file tail $file] $transfertype $force
			}
		}
		foreach field $optionfields value [list_sub $line -exclude $poss] {
			if {$value eq ""} continue
			if {$value eq "-"} {set value {}}
			set optionsa($samplename,$field) $value
		}
	}
	set keys [bsort [array names optionsa]]
	if {[llength $keys]} {
		set o [open $projectdir/options.tsv w]
		puts $o [join {sample option value} \t]
		foreach key $keys {
			foreach {sample option} [split $key ,] break
			puts $o $sample\t$option\t$optionsa($key)
		}
	}
	set dir [lindex $args 0]
	if {$amplicons ne ""} {
		set temp [lindex [split [file root [gzroot [file tail $amplicons]]] -] end]
		set target $projectdir/reg_amplicons-$temp.tsv[gzext $amplicons]
		project_transfer $amplicons $target $transfertype $force
	}
	if {$targetfile ne ""} {
		set temp [file root [gzroot [file tail $targetfile]]]
		regsub {^reg_[^_]+_exome_} $temp {} temp
		set temp [lindex [split $temp -] end]
		set target $projectdir/reg_targets-$temp.tsv[gzext $targetfile]
		project_transfer $targetfile $target $transfertype $force
	}
}

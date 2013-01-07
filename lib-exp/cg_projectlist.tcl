proc cg_makeprojects {{projectsfile /complgen/project-work/cgdata.tab}} {
	set dir [file dir $projectsfile]
	cd $dir
	set list [split [cg select -f {CGI_Sample_ID SampleDir Project Family users reference} $projectsfile] \n]
	list_shift list
	unset -nocomplain a
	foreach line $list {
		foreach {cgsample sampledir project fam users reference} [split $line \t] break
		if {$project eq ""} continue
		lappend a($project) $line
	}
	set projects [array names a]
	foreach project $projects {
		cd $dir
		puts =====$project=====
		file mkdir $project
		set line [split [lindex $a($project) 0] \t]
		set users [string trim [lindex $line 4] \"]
		set reference [string trim [lindex $line 5] \"]
		set result [list "users\t[list $users]" "build\t$reference" {}]
		file delete -force $dir/$project/oricg
		file mkdir $dir/$project/oricg
		cd $dir/$project/oricg
		list_foreach {cgsample sampledir project fam users reference} $a($project) {
			set oridir [glob -nocomplain ../../../oricg/*/*/$cgsample]
			if {![llength $oridir]} {
				set oridir [glob -nocomplain /complgen/oricg/*/*/$cgsample]
			}
			if {[llength $oridir]} {
				set oridir [lindex $oridir 0]
				file link -symbolic [file tail $oridir] $oridir
			} else {
				set oridir /complgen/oricg/$cgsample
			}
			lappend result $cgsample\t$sampledir
		}
		cd $dir
		if {[file exists $dir/$project/$project.cgprj] && ![file exists $dir/$project/$project.cgprj.save]} {
			file delete $dir/$project/$project.cgprj.save
			file copy $dir/$project/$project.cgprj $dir/$project/$project.cgprj.save
		}
		file delete $dir/$project/$project.cgprj
		file_write $dir/$project/$project.cgprj [join $result \n]\n
	}
}

if 0 {
# set projectsfile /complgen/projects/cgdata.tab
# set projects {cmt154 cmtW pn211 ep136 ep80}

	# check for differences
	diff -r /complgen/projects-work/ /complgen/projects
	# copy differences, making hardlinks
	cp -alf /complgen/projects-work/* /complgen/projects


	cd /complgen
	cp -al projects projects-work
	cp -alf projects-work/* projects

	rm tmp2/test
	echo 'abcde' > tmp1/test
	echo 'tmp2 12345' > tmp2/test
	cat tmp1/test
	cat tmp2/test
	cp -alf tmp2/* tmp1
	ls -i tmp1/test
	ls -i tmp2/test
}

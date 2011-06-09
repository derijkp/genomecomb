proc process_compare_coverage {dir1 dir2 dbdir resultsdir {force 0}} {
	set dir1 [file normalize $dir1]
	set dir2 [file normalize $dir2]
	set dbdir [file normalize $dbdir]
	file mkdir $resultsdir
	cd $resultsdir
	# region files
	set name1 [file tail $dir1]
	set name2 [file tail $dir2]
	set target ${name1}_${name2}_regcommon.tsv.covered
	if {$force || ![file exists $target]} {
		putslog "Making $target"
		cg regsubtract $dir1/sreg-${name1}.tsv $dir2/sreg-${name2}.tsv > temp.tsv
		cg regsubtract $dir1/sreg-${name1}.tsv temp.tsv > ${name1}_${name2}_regcommon.tsv
		cg covered ${name1}_${name2}_regcommon.tsv > $target
	}
	# region from each sample
	set filters {}
	lappend filters refcons [list $dir1/reg_refcons-${name1}.tsv $dir2/reg_refcons-${name2}.tsv]
	lappend filters ns [list $dir1/reg_ns-${name1}.tsv $dir2/reg_ns-${name2}.tsv]
	lappend filters str [list $dbdir/regdb-microsatelite.tsv {}]
	lappend filters trf [list $dbdir/regdb-simple_repeats.tsv {}]
	lappend filters cluster [list $dir1/reg_cluster-${name1}.tsv $dir2/reg_cluster-${name2}.tsv]
	lappend filters rtg [list $dir1/reg_rtg-${name1}.tsv $dir2/reg_rtg-${name2}.tsv]
	lappend filters segdup [list $dbdir/regdb-segdups.tsv {}]
	lappend filters b10 [list $dir1/reg_below10-${name1}.tsv $dir2/reg_below10-${name2}.tsv]
	lappend filters b20 [list $dir1/reg_below20-${name1}.tsv $dir2/reg_below20-${name2}.tsv]
	lappend filters a100 [list $dir1/reg_above100-${name1}.tsv $dir2/reg_above100-${name2}.tsv]
	lappend filters lowscore [list $dir1/reg_lowscore-${name1}.tsv $dir2/reg_lowscore-${name2}.tsv]
	lappend filters selfchain [list $dbdir/regdb-selfchain.tsv {}]
	lappend filters repeat [list $dbdir/regdb-repeatmasker.tsv {}]
	set todo [list_unmerge $filters]
	foreach dbname $todo {
		foreach {rfile1 rfile2} [dict get $filters $dbname] break
		set target ${name1}_${name2}_regcommon-$dbname.tsv.covered
		if {$force || ![file exists $target]} {
			putslog "Making $target"
			cg regsubtract ${name1}_${name2}_regcommon.tsv $rfile1 > temp.tsv
			if {$rfile2 ne ""} {
				file rename -force temp.tsv temp-old.tsv
				cg regsubtract temp-old.tsv $rfile2 > temp.tsv
			}
			file rename -force temp.tsv ${name1}_${name2}_regcommon-$dbname.tsv
			cg covered ${name1}_${name2}_regcommon-$dbname.tsv > temp.covered
			file rename -force temp.covered $target
		} else {
			putslog "$target already there"
		}
	}
	# subsequent subtractions (filters)
	set orders {
		{ns trf str cluster rtg segdup b20 a100 lowscore selfchain repeat}
		{ns trf str cluster rtg segdup b10 a100 lowscore selfchain repeat}
		{ns trf str segdup cluster lowscore b20 a100 rtg}
		{ns trf str segdup cluster lowscore b10 a100 rtg}
		{ns trf str segdup cluster lowscore rtg}
		{ns b20 a100 cluster lowscore trf str segdup rtg}
		{ns b10 a100 cluster lowscore trf str segdup rtg}
		{ns cluster lowscore trf str segdup rtg}
		{ns b20 a100 lowscore cluster trf str segdup rtg}
		{ns b10 a100 lowscore cluster trf str segdup rtg}
		{ns lowscore cluster trf str segdup rtg}
	}
	foreach order $orders {
		set newname ${name1}_${name2}_regcommon-refcons
		foreach dbname $order {
			foreach {rfile1 rfile2} [dict get $filters $dbname] break
			set oldname $newname
			set newname $oldname-$dbname
			set target $newname.tsv.covered
			if {$force || ![file exists $target]} {
				putslog "Making $target"
				cg regsubtract $oldname.tsv  $rfile1 > temp.tsv
				if {$rfile2 ne ""} {
					file rename -force temp.tsv temp-old.tsv
					cg regsubtract temp-old.tsv $rfile2 > temp.tsv
				}
				file rename -force temp.tsv $newname.tsv
				cg covered $newname.tsv > temp.covered
				file rename -force temp.covered $target
			} else {
				putslog "$target already there"
			}
		}
	}
	exec grep total *.covered > summary_genomecoverage_${name1}_${name2}.tsv
}

proc cg_process_compare_coverage {args} {
	if {[llength $args] == 2} {
		set force 0
		foreach {dbdir resultsdir} $args break
		set dbdir [file normalize $dbdir]
		set resultsdir [file normalize $resultsdir]
		foreach {dir1 dir2} [lrange [split $resultsdir _] end-1 end] break
		set dir1 [file dir $resultsdir]/$dir1
		set dir2 [file dir $resultsdir]/$dir2
		set resultsdir $resultsdir/genomecoverage
		file mkdir $resultsdir
	} else {
		if {([llength $args] < 4) || ([llength $args] > 5)} {
			puts stderr "format is: $scriptname $action sampledir1 sampledir2 dbdir resultsdir ?force?"
			puts stderr " - calculate coverages for 2 compared sample directories."
			puts stderr " - By default, only files that are not present already will be created."
			puts stderr " -When force is given as a parameter, everything will be recalculated and overwritten."
			exit 1
		}
		foreach {dir1 dir2 dbdir resultsdir force} $args break
		switch $force {
			force {set force 1}
			"" {set force 0}
			default {error "unrecognized option $force"}
		}
	}
	process_compare_coverage $dir1 $dir2 $dbdir $resultsdir $force
}

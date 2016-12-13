proc hsmetrics_tsv2interval {regionfile resultfile bamfile} {
	if {[file extension $regionfile] eq ".bed"} {
		set tsvfile [tempfile]
		cg bed2sft $regionfile $tsvfile
	} else {
		set tsvfile $regionfile
	}
	set f [gzopen $tsvfile]
	set header [tsv_open $f]
	gzclose $f
	set poss [tsv_basicfields $header 3]
	set fields {}
	foreach {field} [list_sub $header $poss] nfield {chromosome begin end} {
		lappend fields "$nfield=\$$field"
	}
	lappend fields {strand="+"}
	set pos [lsearch $header name]
	if {$pos == -1} {
		set pos [lsearch $header info]
	}
	if {$pos != -1} {lappend fields "name=\$[lindex $header $pos]"} else {lappend fields {name="all"}}
	#remove comment columns & add strand info - due to lack of correct strand info take + as default
	set bamheader [split [exec samtools view -H $bamfile] \n]
	file_write $resultfile.temp [join [list_sub $bamheader [list_find -regexp $bamheader ^@SQ]] \n]\n
	cg select -sh /dev/null -f $fields $tsvfile >> $resultfile.temp
	file rename -force $resultfile.temp $resultfile
	file delete $resultfile.temp.pre
}

proc cg_hsmetrics {args} {
	cg_options hsmetrics args {
		-s - --sample {
			set sample $value
		}
		-b - --baitfile {
			set baitfile $value
		}
	} {bamfile targetfile resultfile} 3 3
	if {![info exists sample]} {
		set sample [file tail [file root $bamfile]]
		regsub ^map- $sample {} sample
	}
	set target_intervals [tempfile]
	hsmetrics_tsv2interval $targetfile $target_intervals $bamfile
	if {![info exists baitfile]} {
		# We have to give a bait interval file, or CalculateHsMetrics wont run
		# if we do not have actual bait regions, use target regions.
		set bait_intervals $target_intervals
	} else {
		set bait_intervals [tempfile]
		hsmetrics_tsv2interval $baitfile $bait_intervals $bamfile
	}
	picard CalculateHsMetrics BAIT_INTERVALS=$bait_intervals TARGET_INTERVALS=$target_intervals I=$bamfile O=$resultfile.temp 2>@ stderr
	cg select -f [list sample=\"$sample\" *] $resultfile.temp $resultfile.temp2
	file rename -force $resultfile.temp2 $resultfile
	file delete $target_intervals
}

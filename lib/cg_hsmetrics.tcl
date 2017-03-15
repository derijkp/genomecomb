proc calculate_hsmetrics_job {bamfile targetsfile {optional 1}} {
	#calculate hsmetrics with picard tools (=coverage statistics): input bamfile & targetsfile
	upvar job_logdir job_logdir
	set bamfile [file_absolute $bamfile]
	set dir [file dir $bamfile]
	set file [file tail $bamfile]
	set root [join [lrange [split [file root $file] -] 1 end] -]
	set target $dir/$root.hsmetrics
	job calc_hsmetrics-$root -optional $optional -deps {$bamfile $bamfile.bai $targetsfile} -targets [list $target] -vars {bamfile targetsfile} -code {
		cg_hsmetrics $bamfile $targetsfile $target
	}
	return $target
}

proc make_hsmetrics_report_job {destdir files {optional 1}} {
	upvar job_logdir job_logdir
	set experiment [file tail $destdir]
	if {![llength $files]} return
	job calc_hsmetrics-$experiment -optional $optional -deps $files -targets $destdir/${experiment}_hsmetrics_report.tsv -code {
		cg cat -c 0 {*}$deps > $target.temp
		cg select -rc 1 $target.temp $target.temp2
		file rename -force $target.temp2 $target
		file delete $target.temp
	}
}

proc hsmetrics_tsv2interval {regionfile bamfile resultfile} {
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
	set ffields [list_sub $header $poss]
	lappend fields "chromosome=\$[lindex $ffields 0]"
	lappend fields "begin=\$[lindex $ffields 1] + 1"
	lappend fields "end=\$[lindex $ffields 2]"
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
		-s - -sample {
			set sample $value
		}
		-b - -baitfile {
			set baitfile $value
		}
	} {bamfile targetfile resultfile} 3 3
	if {![info exists sample]} {
		set sample [file tail [file root $bamfile]]
		regsub ^map- $sample {} sample
	}
	set target_intervals [tempdir]/[file root [file tail [gzroot $targetfile]]].intervals
	hsmetrics_tsv2interval $targetfile $bamfile $target_intervals
	if {![info exists baitfile]} {
		# We have to give a bait interval file, or CalculateHsMetrics wont run
		# if we do not have actual bait regions, use target regions.
		set bait_intervals $target_intervals
	} else {
		set bait_intervals [tempdir]/[file root [file tail [gzroot $baitfile]]].intervals
		hsmetrics_tsv2interval $baitfile $bamfile $bait_intervals
	}
	set temp [tempfile]
	set temptarget [filetemp $resultfile]
	picard CalculateHsMetrics BAIT_INTERVALS=$bait_intervals TARGET_INTERVALS=$target_intervals I=$bamfile O=$temp
	cg select -f [list sample=\"$sample\" *] $temp $temptarget
	file rename -force $temptarget $resultfile
	file delete $target_intervals
}

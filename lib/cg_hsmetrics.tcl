proc calculate_hsmetrics_job {bamfile targetsfile {optional 1}} {
	#calculate hsmetrics with picard tools (=coverage statistics): input bamfile & targetsfile
	upvar job_logdir job_logdir
	set bamfile [file_absolute $bamfile]
	set dir [file dir $bamfile]
	set root [file_rootname $bamfile]
	set target $dir/$root.hsmetrics
	set indexext [indexext $bamfile]
	job calc_hsmetrics-$root -optional $optional -deps {$bamfile $bamfile.$indexext $targetsfile} -targets [list $target] -vars {bamfile targetsfile} -code {
		cg_hsmetrics $bamfile $targetsfile $target
	}
	return $target
}

proc make_hsmetrics_report_job {destdir files {optional 1}} {
	upvar job_logdir job_logdir
	set experiment [file tail $destdir]
	if {![llength $files]} return
	job calc_hsmetrics-$experiment -optional $optional -deps $files -targets {$destdir/${experiment}_hsmetrics_report.tsv} -code {
		cg cat -c 0 {*}$deps > $target.temp
		cg select -rc 1 $target.temp $target.temp2
		file rename -force $target.temp2 $target
		file delete $target.temp
	}
}

proc hsmetrics_tsv2interval {regionfile bamfile resultfile} {
	if {[file extension $regionfile] eq ".bed"} {
		set tsvfile [tempfile]
		cg bed2tsv $regionfile $tsvfile
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
	set num [lindex [cg select -g all $targetfile] end]
	if {$num == 0} {
		# no regions -> write "dummy" with all 0
		set f [open $resultfile.temp w]
		puts $f [join {sample	BAIT_SET	GENOME_SIZE	BAIT_TERRITORY	TARGET_TERRITORY	BAIT_DESIGN_EFFICIENCY	TOTAL_READS	PF_READS	PF_UNIQUE_READS	PCT_PF_READS	PCT_PF_UQ_READS	PF_UQ_READS_ALIGNED	PCT_PF_UQ_READS_ALIGNED	PF_UQ_BASES_ALIGNED	ON_BAIT_BASES	NEAR_BAIT_BASES	OFF_BAIT_BASES	ON_TARGET_BASES	PCT_SELECTED_BASES	PCT_OFF_BAIT	ON_BAIT_VS_SELECTED	MEAN_BAIT_COVERAGE	MEAN_TARGET_COVERAGE	PCT_USABLE_BASES_ON_BAIT	PCT_USABLE_BASES_ON_TARGET	FOLD_ENRICHMENT	ZERO_CVG_TARGETS_PCT	FOLD_80_BASE_PENALTY	PCT_TARGET_BASES_2X	PCT_TARGET_BASES_10X	PCT_TARGET_BASES_20X	PCT_TARGET_BASES_30X	HS_LIBRARY_SIZE	HS_PENALTY_10X	HS_PENALTY_20X	HS_PENALTY_30X	AT_DROPOUT	GC_DROPOUT	SAMPLE	LIBRARY	READ_GROUP} \t]
		set bait [file tail [file root $targetfile]]
		set temp [list_fill 39 {}]
		foreach i {27 28 29 30} {
			lset temp $i 0
		}
		puts $f $sample\t$bait[join $temp \t]
		close $f
		file rename -force $resultfile.temp $resultfile
		return
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
	cg select -overwrite 1 -f [list sample=\"$sample\" *] $temp $temptarget
	file rename -force $temptarget $resultfile
	file delete $target_intervals
}

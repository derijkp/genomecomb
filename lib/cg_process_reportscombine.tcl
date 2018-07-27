proc proces_reportscombine_job {args} {
	cg_options process_reportscombine args {
		-experimentname {
			set experimentname $value
		}
	} {destdir reportsdir} 2 ... {
		combines reports directories (per sample) into one overview reports directory
	}
	set destdir [file_absolute $destdir]
	set reportstodo [list $reportsdir {*}$args]
	upvar job_logdir job_logdir
	if {![info exists job_logdir]} {
		job_logdir $destdir/log_jobs
	}
	if {![info exists experimentname]} {
		set experimentname [file tail $destdir]
		if {$experimentname eq "reports"} {
			set experimentname [file tail [file dir $destdir]]
		}
	}
	set deps {}
	foreach dir $reportstodo {
		lappend deps {*}[jobglob $dir/hsmetrics-*.hsmetrics $dir/reports/]
	}
	if {[llength $deps]} {
		set target $destdir/report_hsmetrics-${experimentname}.tsv
		set target2 $destdir/${experimentname}_hsmetrics_report.tsv
		job reportscombine_hsmetrics-$experimentname -deps $deps -targets {$target $target2} -code {
			file mkdir [file dir $target]
			cg cat -c 0 {*}$deps > $target.temp
			cg select -rc 1 $target.temp $target.temp2
			file rename -force $target.temp2 $target
			file delete $target.temp
			mklink $target $target2
		}
	}
	set deps {}
	foreach dir $reportstodo {
		lappend deps {*}[jobglob $dir/report_*.tsv]
	}
	if {[llength $deps]} {
		set target $destdir/report_stats-${experimentname}.tsv
		set target2 $destdir/report_summarytable-${experimentname}.tsv
		job reportscombine_stats-$experimentname -deps $deps -targets {$target $target2} -code {
			# make report_stats
			file mkdir [file dir $target]
			cg cat -m -c 0 {*}[lsort -dict $deps] > $target.temp
			cg select -rc 1 $target.temp $target.temp2
			file rename -force $target.temp2 $target
			file delete $target.temp
			# make report_summarytable
			unset -nocomplain data
			unset -nocomplain a
			set f [open $target]
			set header [tsv_open $f]
			while {[gets $f line] != -1} {
				foreach {analysis source parameter value} [split $line \t] break
				set a($analysis) 1
				set parameter [string trim $parameter]
				set data($analysis,$source,$parameter) $value
			}
			close $f
			set o [open $target2.temp w]
			puts $o [join {
				analysis sample numreads 
				pf_reads pct_pf_reads pf_unique_reads pct_pf_unique_reads pf_mapped pct_pf_aligned_reads
				targetbases	pct_target_bases_2X pct_target_bases_10X pct_target_bases_20X pct_target_bases_30X
				covered_total qvars qvars_target qvars_refcoding
			} \t]
			set analyses [ssort -natural [array names a]]
			foreach analysis $analyses {
				# only report for "endpoints" analysis (data of partial will be included in endpoints)
				if {[regexp {^cov5-} $analysis]} continue
				if {[llength [array names a *-$analysis]]} continue
				set sample [lindex [split $analysis -] end]
				set resultline [list $analysis $sample]
				# numreads
				set fw_numreads [get data($sample,fastq-stats,fw_numreads) {}]
				set rev_numreads [get data($sample,fastq-stats,rev_numreads) {}]
				if {[isint $fw_numreads] && [isint $rev_numreads]} {
					set numreads [expr {$fw_numreads+$rev_numreads}]
				} else {
					set numreads {}
				}
				lappend resultline $numreads
				# bam based stats
				set bamname [join [lrange [split $analysis -] end-1 end] -]
				set pf_reads [get "data($bamname,flagstat_reads,in total)" {}]
				set pf_duplicates [get "data($bamname,flagstat_reads,duplicates)" {}]
				set pf_mapped [get "data($bamname,flagstat_reads,mapped)" {}]
				set pf_properlypaired [get "data($bamname,flagstat_reads,properly paired)" {}]
				if {[isint $numreads] && [isint $pf_reads]} {
					set pct_pf_reads [format %.2f [expr {$pf_reads*100.0/$numreads}]]
					set pf_unique_reads [expr {($pf_reads - $pf_duplicates)}]
					set pct_pf_unique_reads [format %.2f [expr {$pf_unique_reads*100.0/$pf_reads}]]
					set pct_pf_aligned_reads [format %.2f [expr {($pf_mapped)*100.0/$pf_reads}]]
				} else {
					set pct_pf_reads {}
					set pf_unique_reads {}
					set pct_pf_unique_reads {}
					set pct_pf_aligned_reads {}
				}
				lappend resultline $pf_reads $pct_pf_reads $pf_unique_reads $pct_pf_unique_reads $pf_mapped $pct_pf_aligned_reads
				# rest
				lappend resultline [get data($bamname,histodepth,targetbases) {}]
				lappend resultline [get data($bamname,histodepth,pct_target_bases_2X) {}]
				lappend resultline [get data($bamname,histodepth,pct_target_bases_10X) {}]
				lappend resultline [get data($bamname,histodepth,pct_target_bases_20X) {}]
				lappend resultline [get data($bamname,histodepth,pct_target_bases_30X) {}]
				lappend resultline [get data($analysis,genomecomb,covered_total) {}]
				lappend resultline [get data($analysis,genomecomb,qvars) {}]
				lappend resultline [get data($analysis,genomecomb,qvars_target) {}]
				lappend resultline [get data($analysis,genomecomb,qvars_refcoding) {}]
				puts $o [join $resultline \t]
			}
			close $o
			file rename -force $target2.temp $target2
		}
	}
}

proc cg_process_reportscombine {args} {
	set args [job_init {*}$args]
	proces_reportscombine_job {*}$args
	job_wait
}

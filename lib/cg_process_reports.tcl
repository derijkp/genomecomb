proc process_reports_job {args} {
	set reports all
	cg_options process_reports args {
		-dbdir {
			set dbdir $value
		}
		-r - -reports {
			set reports $value
		}
	} {sampledir dbdir reports} 1 3
	set dbdir [dbdir $dbdir]
	set sampledir [file_absolute $sampledir]
	set sample [file tail $sampledir]
	job_logdir $sampledir/job_logdir
	set ref [dbdir_ref $dbdir]
	if {$reports eq "all"} {
		set reports {flagstats fastqstats fastqc vars hsmetrics covered}
	}
	set bamfiles [jobglob $sampledir/*.bam]
	file mkdir $sampledir/reports
	foreach bamfile $bamfiles {
		set sample [file root [file tail [gzroot $bamfile]]]
		regsub ^map- $sample {} sample
		if {[inlist $reports flagstats]} {
			set dep $bamfile
			set target $sampledir/reports/flagstat-$sample.flagstat
			set target2 $sampledir/reports/report_bam-$sample.tsv
			job reports_flagstats-[file tail $bamfile] -deps {$dep} -targets {$target $target2} -vars {sample} -code {
				exec samtools flagstat $dep > $target.temp
				file rename -force -force $target.temp $target
				set o [open $target2.temp w]
				puts $o [join {sample source parameter value value_qcfail} \t]
				set f [open $target]
				while {[gets $f line] != -1} {
					if {![regexp {^([0-9]+) \+ ([0-9]+) ([^()]*)} $line temp value value_qcfail parameter]} continue
					puts $o "$sample\tflagstat\t$parameter\t$value\t$value_qcfail"
				}
				close $f
				close $o
				file rename -force $target2.temp $target2
			}
		}
		if {[inlist $reports hsmetrics]} {
			set capturefile [lindex [jobglob $sampledir/info_capture.txt] 0]
			set targetfile $sampledir/reg_${ref}_targets.tsv
			if {$capturefile ne ""} {
				job reports_targetfile -optional 1 -deps {$capturefile} -targets {$targetfile} -vars {sample dbdir ref} -code {
					set capture [string trim [file_read $dep]]
					set oritargetfile $dbdir/extra/reg_${ref}_exome_$capture.tsv
					if {![file exists $oritargetfile]} {
						array set transa {seqcapv3 SeqCap_EZ_v3 sure4 SureSelectV4 sure5 SureSelectV5 sure5utr SureSelectV5UTR}
						set capture [get transa($capture) $capture]
						set oritargetfile $dbdir/extra/reg_${ref}_exome_$capture.tsv
					}
					mklink $oritargetfile $target 1
				}
			}
			set dep1 $bamfile
			set dep2 $targetfile
			set target $sampledir/reports/hsmetrics-$sample.hsmetrics
			set target2 $sampledir/reports/report_hsmetrics-$sample.tsv
			job reports_hsmetrics-[file tail $bamfile] -optional 1 -deps {$dep1 $dep2} -targets {$target $target2} -vars {sample bamfile targetfile} -code {
				cg_hsmetrics --sample $sample $bamfile $targetfile $target
				set f [open $target]
				set header [tsv_open $f]
				set data [split [gets $f] \t]
				close $f
				set target2temp [filetemp target2]
				set o [open $target2temp w]
				puts $o [join {sample source parameter value} \t]
				foreach key [lrange $header 1 end] value [lrange $data 1 end] {
					puts $o "$sample\thsmetrics\t$key\t$value"
				}
				close $o
				file rename -force $target2temp $target2
			}
		}
	}
	set fastqfiles [ssort -natural [jobglob $sampledir/fastq/*]]
	if {[inlist $reports fastqc]} {
		foreach fastqfile $fastqfiles {
			set name [file root [file tail [gzroot $fastqfile]]]-$sample
			set dep $fastqfile
			set outdir $sampledir/reports/fastqc
			file mkdir $outdir
			set target $sampledir/reports/fastqc/${name}.fq_fastqc
			job reports_fastqc-[file tail $fastqfile] -deps {$dep} -targets {$target} -vars {outdir} -code {
				exec fastqc -q -o $outdir $dep 2>@ stderr >@ stdout
				file delete $target.zip
			}
		}
	}
	if {[inlist $reports fastqstats]} {
		set fastqfiles_fw [list_unmerge $fastqfiles 1 fastqfiles_rev]
		foreach deps [list $fastqfiles_fw $fastqfiles_rev] dir {fw rev} {
			set target $sampledir/reports/report_fastq_$dir-$sample.tsv
			set target2 $sampledir/reports/fastq_stats_$dir-$sample.txt
			set target3 $sampledir/reports/fastx_$dir-$sample.tsv
			job reports_fastq-stats-$dir -deps $deps -targets {$target $target2 $target3} -vars {sample dir} -code {
				set gzcat [gzcat [lindex $deps 0]]
				exec {*}$gzcat {*}$deps | fastq-stats -x $target3 > $target2
				set o [open $target.temp w]
				puts $o [join {sample source parameter value} \t]
				set f [open $target2]
				set dups {}
				while {[gets $f line] != -1} {
					set line [split $line \t]
					foreach {parameter value} $line break
					set parameter [string_change [string trim $parameter] {% pct_}]
					regsub -all {[^a-zA-z]} $parameter _ parameter
					switch $parameter {
						reads {
							puts $o [join [list $sample fastq-stats ${dir}_numreads $value] \t]
						}
						dup_seq {
							lappend dups [lindex $line 2] [lindex $line 3]
						}
						default {
							puts $o [join [list $sample fastq-stats ${dir}_$parameter $value] \t]
						}
					}
				}
				close $f
				puts $o [join [list $sample fastq-stats ${dir}_dups $dups] \t]
				close $o
				file rename -force $target.temp $target
			}
		}
	}
	if {[inlist $reports vars]} {
		foreach dep [jobglob $sampledir/var-*.tsv] {
			set sample [file root [file tail [gzroot $dep]]]
			regsub ^var- $sample {} sample
			set target $sampledir/reports/report_vars-$sample.tsv
			job reports_vars-$sample -deps {$dep} -targets {$target} -vars sample -code {
				set qvars [lindex [exec cg select -q {$coverage >= 20 and $quality >= 50} -g all $dep] end]
				set vars [lindex [exec cg select -q {$sequenced == "v"} -g all $dep] end]
				set f [open $target.temp w]
				puts $f [join {sample source parameter value} \t]
				puts $f $sample\tgenomecomb\tvars\t$vars
				puts $f $sample\tgenomecomb\tqvars\t$qvars
				close $f
				file rename -force $target.temp $target
			}
		}
	}
	if {[inlist $reports covered]} {
		foreach dep [jobglob $sampledir/sreg-*.tsv] {
			set sample [file root [file tail [gzroot $dep]]]
			regsub ^sreg- $sample {} sample
			set target $sampledir/reports/report_covered-$sample.tsv
			job reports_vars-$sample -deps {$dep} -targets {$target} -vars sample -code {
				set temp [split [exec cg covered $dep] \n]
				set f [open $target.temp w]
				puts $f [join {sample source parameter value} \t]
				foreach line [lrange $temp 1 end] {
					foreach {chr cov} $line break
					puts $f $sample\tgenomecomb\tcovered_$chr\t$cov
				}
				close $f
				file rename -force $target.temp $target
			}
		}
	}
}

proc proces_reportscombine_job {destdir reportstodo} {
	upvar job_logdir job_logdir
	set experiment [file tail $destdir]
	set deps {}
	foreach dir $reportstodo {
		lappend deps {*}[jobglob $dir/hsmetrics-*.hsmetrics]
	}
	if {[llength $deps]} {
		set target $destdir/reports/report_hsmetrics-${experiment}.tsv
		set target2 $destdir/${experiment}_hsmetrics_report.tsv
		job reportscombine_hsmetrics-$experiment -deps $deps -targets {$target $target2} -code {
			file mkdir [file dir $target]
			cg cat -c 0 {*}$deps > $target.temp
			cg select -rc 1 $target.temp $target.temp2
			file rename -force $target.temp2 $target
			file delete $target.temp
			mklink $target $target2 1
		}
	}
	set deps {}
	foreach dir $reportstodo {
		lappend deps {*}[jobglob $dir/report_*.tsv]
	}
	if {[llength $deps]} {
		set target $destdir/reports/report_stats-${experiment}.tsv
		job reportscombine_stats-$experiment -deps $deps -targets {$target} -code {
			file mkdir [file dir $target]
			cg cat -m -c 0 {*}$deps > $target.temp
			cg select -rc 1 $target.temp $target.temp2
			file rename -force $target.temp2 $target
			file delete $target.temp
		}
	}
}

proc cg_process_reports {args} {
	set args [job_init {*}$args]
	if {[llength $args] < 1} {
		errorformat process_reports
	}
	process_reports_job {*}$args
	job_wait
}
proc targetfile_job {sampledir {dbdir {}}} {
	upvar job_logdir job_logdir
	set dbdir [dbdir $dbdir]
	set ref [dbdir_ref $dbdir]
	set targetfile $sampledir/reg_${ref}_targets.tsv
	if {[jobfileexists $targetfile]} {
		return $targetfile
	}
	# get targetfile (if possible)
	set capturefile $sampledir/info_capture.txt
	if {![jobfileexists $capturefile]} {
		if {[file tail [file dir $sampledir]] eq "samples"} {
			set take2 [file dir [file dir $sampledir]]/reg_${ref}_targets.tsv
			set capture [file dir [file dir $sampledir]]/info_capture.tsv
		} else {
			set take2 [file dir $sampledir]/reg_${ref}_targets.tsv
			set capture [file dir $sampledir]/info_capture.tsv
		}
		if {[jobfileexists $take2]} {
			set take2 [find_link $take2 1]
			mklink_job $take2 $targetfile
			return $targetfile
		}
	}
	if {[jobfileexists $capturefile]} {
		job reports_targetfile -deps {$capturefile} -targets {$targetfile} -vars {sample dbdir ref} -code {
			set capture [string trim [file_read $dep]]
			set oritargetfile $dbdir/extra/reg_${ref}_exome_$capture.tsv
			if {![file exists $oritargetfile]} {
				array set transa {seqcapv3 SeqCap_EZ_v3 sure4 SureSelectV4 sure5 SureSelectV5 sure5utr SureSelectV5UTR}
				set capture [get transa($capture) $capture]
				set oritargetfile $dbdir/extra/reg_${ref}_exome_$capture.tsv
			}
			mklink $oritargetfile $target 1
		}
		return $targetfile
	}
	return {}
}

proc process_reports_job {args} {
	set reports all
	cg_options process_reports args {
		-dbdir {
			set dbdir $value
		}
		-r - -reports {
			set reports $value
		}
	} {sampledir dbdir reports} 1 3 {
		Calculates a number of statistics on a sample in the reports subdir
	}
	set dbdir [dbdir $dbdir]
	set sampledir [file_absolute $sampledir]
	set sample [file tail $sampledir]
	job_logdir $sampledir/log_jobs
	set ref [dbdir_ref $dbdir]
	if {$reports eq "all"} {
		set reports {flagstats fastqstats fastqc vars hsmetrics covered}
	}
	set bamfiles [jobglob $sampledir/*.bam]
	set targetfile [targetfile_job $sampledir $dbdir]
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
		if {[inlist $reports hsmetrics] && $targetfile ne ""} {
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
	if {[inlist $reports fastqc] && [llength $fastqfiles]} {
		foreach fastqfile $fastqfiles {
			set name [file root [file tail [gzroot $fastqfile]]]
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
	if {[inlist $reports fastqstats] && [llength $fastqfiles]} {
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
		set refcodingfile $dbdir/extra/reg_hg19_refcoding.tsv
		foreach varfile [jobglob $sampledir/var-*.tsv] {
			set sample [file root [file tail [gzroot $varfile]]]
			regsub ^var- $sample {} sample
			set target $sampledir/reports/report_vars-$sample.tsv
			set deps [list $varfile]
			lappend deps "($refcodingfile)"
			if {$targetfile ne ""} {lappend deps $targetfile}
			job reports_vars-$sample -deps $deps -targets {$target} -vars {sample varfile targetfile refcodingfile dbdir build} -code {
				set tempfile [tempfile]
				set fields {}
				set header [cg select -h $varfile]
				if {"coverage" in $header && "quality" in $header} {
					set usehq 1
					lappend fields {hq=if($coverage >= 20 and $quality >= 50,1,0)}
				} else {
					set usehq 0
					lappend fields {hq=0}
				}
				set annotfiles {}
				if {[file exists $targetfile]} {
					lappend annotfiles $targetfile
					lappend fields {target=if($targets ne "",1,0)}
				} else {
					lappend fields {target=0}
				}
				set refcodingfile [gzfile $refcodingfile]
				if {[file exists $refcodingfile]} {
					lappend annotfiles $refcodingfile
					lappend fields {refcoding=if($refcoding ne "",1,0)}
				} else {
					lappend fields {refcoding=0}
				}
				if {[llength $annotfiles]} {
					cg annotate $varfile $tempfile {*}$annotfiles
				} else {
					mklink $varfile $tempfile
				}
				if {[inlist $header sequenced]} {
					set query {$sequenced == "v"}
				} elseif {[inlist $header zyg]} {
					set query {$zyg in "m t c"}
				} else {
					set query {}
				}
				set temp [exec cg select -q $query \
					-f $fields -g {hq {} target {} refcoding {} type} $tempfile]
				set vars 0 ; set qvars 0 ; set qvars_target 0; set qvars_refcoding 0
				set vars_snp 0 ; set qvars_snp 0 ; set qvars_target_snp 0; set qvars_refcoding_snp 0
				foreach line [lrange [split [string trim $temp] \n] 1 end] {
					foreach {hq ontarget refcoding type count} [split $line \t] break
					if {$hq} {
						incr qvars $count
						if {$type eq "snp"} {incr qvars_snp $count}
						if {$ontarget} {
							incr qvars_target $count
							if {$type eq "snp"} {incr qvars_target_snp $count}
						}
						if {$refcoding} {
							incr qvars_refcoding $count
							if {$type eq "snp"} {incr qvars_refcoding_snp $count}
						}
					}
					incr vars $count
					if {$type eq "snp"} {incr vars_snp $count}
				}
				set f [open $target.temp w]
				puts $f [join {sample source parameter value} \t]
				puts $f $sample\tgenomecomb\tvars\t$vars
				puts $f $sample\tgenomecomb\tvars_snp\t$vars_snp
				puts $f $sample\tgenomecomb\tqvars\t$qvars
				puts $f $sample\tgenomecomb\tqvars_snp\t$qvars_snp
				if {[file exists $targetfile]} {
					puts $f $sample\tgenomecomb\tqvars_target\t$qvars_target
					puts $f $sample\tgenomecomb\tqvars_target_snp\t$qvars_target_snp
				}
				if {[file exists $refcodingfile]} {
					puts $f $sample\tgenomecomb\tqvars_refcoding\t$qvars_refcoding
					puts $f $sample\tgenomecomb\tqvars_refcoding_snp\t$qvars_refcoding_snp
				}
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

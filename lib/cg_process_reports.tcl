proc process_reports_job {args} {
	set reports basic
	set dbdir {}
	set allreports {fastqstats fastqc flagstat_reads histodepth hsmetrics vars covered histo predictgender}
	set basicreports {fastqstats fastqc flagstat_reads histodepth hsmetrics vars covered histo predictgender}
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
		set reports $allreports
	} elseif {[string index $reports 0] eq "-"} {
		set reports [list_lremove $allreports [string range $reports 1 end]]
	} elseif {"basic" in $reports} {
		set pos [lsearch $reports all]
		set reports [list_remdup [lreplace $reports $pos $pos {*}$basicreports]]
	}
	# find regionfile indicating target of sequencing (used by hsmetrics, histodepth, vars, so needs to be here)
	set targetfile [targetfile_job $sampledir $dbdir]
	# logfile
	set cmdline [list cg process_reports]
	foreach option {
		dbdir reports
	} {
		if {[info exists $option]} {
			lappend cmdline -$option [get $option]
		}
	}
	lappend cmdline $sampledir $dbdir
	job_logfile $sampledir/process_reports_$sample $sampledir $cmdline \
		{*}[versions dbdir fastqc fastq-stats fastq-mcf bwa bowtie2 samtools gatk gatk3 picard java gnusort8 lz4 os]
	# start
	set bamfiles [jobglob $sampledir/*.bam]
	set ampliconsfile [ampliconsfile $sampledir $ref]
	file mkdir $sampledir/reports
	foreach bamfile $bamfiles {
		set bamroot [file root [file tail [gzroot $bamfile]]]
		regsub ^map- $bamroot {} bamroot
		if {[inlist $reports flagstat_alignments]} {
			set dep $bamfile
			set target $sampledir/reports/flagstat_alignments-$bamroot.flagstat
			set target2 $sampledir/reports/report_flagstat_alignments-$bamroot.tsv
			job reports_flagstat_alignments-[file tail $bamfile] -deps {$dep} -targets {$target $target2} -vars {bamroot} -code {
				exec samtools flagstat $dep > $target.temp
				file rename -force -force $target.temp $target
				set o [open $target2.temp w]
				puts $o [join {sample source parameter value value_qcfail} \t]
				set f [open $target]
				while {[gets $f line] != -1} {
					if {![regexp {^([0-9]+) \+ ([0-9]+) ([^()]*)} $line temp value value_qcfail parameter]} continue
					puts $o "$bamroot\tflagstat_alignments\t[string trim $parameter]\t$value\t$value_qcfail"
				}
				close $f
				close $o
				file rename -force $target2.temp $target2
			}
		}
		if {[inlist $reports flagstat_reads]} {
			set dep $bamfile
			set target $sampledir/reports/flagstat_reads-$bamroot.flagstat
			set target2 $sampledir/reports/report_flagstat_reads-$bamroot.tsv
			job reports_flagstat_reads-[file tail $bamfile] -deps {$dep} -targets {$target $target2} -vars {bamroot} -code {
				if {[catch {
					exec samtools view -F 256 -h -b $dep | samtools flagstat - > $target.temp
				} msg]} {
					if {$msg ne "\[bam_header_read\] EOF marker is absent. The input is probably truncated."} {error $msg}
				}
				file rename -force -force $target.temp $target
				set o [open $target2.temp w]
				puts $o [join {sample source parameter value value_qcfail} \t]
				set f [open $target]
				while {[gets $f line] != -1} {
					if {![regexp {^([0-9]+) \+ ([0-9]+) ([^()]*)} $line temp value value_qcfail parameter]} continue
					puts $o "$bamroot\tflagstat_reads\t[string trim $parameter]\t$value\t$value_qcfail"
				}
				close $f
				close $o
				file rename -force $target2.temp $target2
			}
		}
		if {[inlist $reports hsmetrics] && [jobfileexists $targetfile]} {
			set dep1 $bamfile
			set dep2 $targetfile
			set target $sampledir/reports/hsmetrics-$bamroot.hsmetrics
			set target2 $sampledir/reports/report_hsmetrics-$bamroot.tsv
			job reports_hsmetrics-[file tail $bamfile] -optional 1 -deps {$dep1 $dep2} -targets {$target $target2} -vars {bamroot bamfile targetfile} -code {
				cg_hsmetrics --sample $bamroot $dep1 $dep2 $target
				set f [open $target]
				set header [tsv_open $f]
				set data [split [gets $f] \t]
				close $f
				set target2temp [filetemp $target2]
				set o [open $target2temp w]
				puts $o [join {sample source parameter value} \t]
				foreach key [lrange $header 1 end] value [lrange $data 1 end] {
					puts $o "$bamroot\thsmetrics\t$key\t$value"
				}
				close $o
				file rename -force $target2temp $target2
			}
		}
		if {[inlist $reports histodepth]} {
			set dep1 $bamfile
			set dep2 $targetfile
			set target $sampledir/reports/histodepth-$bamroot.tsv
			set target2 $sampledir/reports/report_histodepth-$bamroot.tsv
			job reports_histodepth-[file tail $bamfile] -optional 1 -deps {$dep1 ($dep2) $dep1.bai} -targets {$target $target2} -vars {bamroot bamfile} -code {
				set targettemp [filetemp $target]
				if {![file exists $dep2]} {
					set targetfile {}
					set tottarget 0
				} else {
					set targetfile $dep2
					set tottarget [lindex [cg covered $targetfile] end]
				}
				cg depth_histo -max 1000 -q 0 -Q 0 $dep1 $targetfile > $targettemp
				file rename -force $targettemp $target
				set c [split [string trim [file_read $target]] \n]
				set header [split [list_shift c] \t]
				if {$header ne "depth ontarget offtarget"} {
					error "$target has wrong format (should have fields: depth ontarget offtarget)"
				}
				set c [list_reverse $c]
				set totontarget 0
				set totofftarget 0
				set result {}
				foreach line $c {
					foreach {depth ontarget offtarget} [split $line \t] break
					incr totontarget $ontarget
					incr totofftarget $offtarget
					if {$depth in {30 20 10 2 1}} {
						if {$tottarget} {
							lappend result "$bamroot\thistodepth\tontarget_bases_${depth}X\t$totontarget"
							lappend result "$bamroot\thistodepth\tpct_target_bases_${depth}X\t[format %.4f [expr {(100.0*$totontarget)/$tottarget}]]"
						}
						lappend result "$bamroot\thistodepth\tofftarget_bases_${depth}X\t$totofftarget"
					}
				}
				set out [join {sample source parameter value} \t]\n
				if {$tottarget} {
					append out "$bamroot\thistodepth\ttargetbases\t$tottarget\n"
				}
				append out [join [ssort -natural $result] \n]
				append out \n
				set target2temp [filetemp $target2]
				file_write $target2temp $out
				file rename -force $target2temp $target2
			}
		}
		if {[inlist $reports histo] && $ampliconsfile ne ""} {
			set dep1 $bamfile
			set dep2 $ampliconsfile
			set target $sampledir/reports/$bamroot.histo
			job reports_histo-[file tail $bamfile] -optional 1 -deps {$dep1 $dep2 $dep1.bai} -targets {$target} -vars {bamroot bamfile} -code {
				set tempfile [filetemp $target]
				set regionfile [filetemp $target]
				cg regcollapse $dep2 > $regionfile
				cg bam_histo $regionfile $dep {1 5 10 20 50 100 200 500 1000} > $tempfile
				file rename -force $tempfile $target
			}
		}
	}
	if {[inlist $reports predictgender]} {
		set target $sampledir/reports/report_predictgender-$sample.tsv
		job predictgender-$sample -deps {$sampledir/map-*.bam $sampledir/var-*.tsv} -vars {dbdir sampledir} -targets {$target} -code {
			cg predictgender -dbdir $dbdir $sampledir $target
		}
	}
	set fastqfiles [ssort -natural [jobglob $sampledir/fastq/*.fastq.gz $sampledir/fastq/*.fastq $sampledir/fastq/*.fq.gz $sampledir/fastq/*.fq]]
	set fastqfiles_fw [list_unmerge $fastqfiles 1 fastqfiles_rev]
	if {[inlist $reports fastqc] && [llength $fastqfiles]} {
		foreach deps [list $fastqfiles_fw $fastqfiles_rev] dir {fw rev} {
			set target $sampledir/reports/fastqc_$dir-$sample.fastqc
			job reports_fastqc-$dir-$sample -deps $deps -targets {$target} -code {
				file mkdir $target.temp
				set gzcat [gzcat [lindex $deps 0]]
				exec -ignorestderr {*}$gzcat {*}$deps | fastqc -o $target.temp stdin 2>1
				exec unzip $target.temp/stdin_fastqc.zip -d $target.temp
				file rename -force {*}[glob $target.temp/stdin_fastqc/*] $target.temp
				file delete $target.temp/stdin_fastqc
				file delete $target.temp/stdin_fastqc.zip
				file delete $target.temp/stdin_fastqc.html
				file rename -force $target.temp $target
			}
		}
	}
	if {[inlist $reports fastqstats] && [llength $fastqfiles]} {
		foreach deps [list $fastqfiles_fw $fastqfiles_rev] dir {fw rev} {
			set target $sampledir/reports/report_fastq_$dir-$sample.tsv
			set target2 $sampledir/reports/fastq_stats_$dir-$sample.txt
			set target3 $sampledir/reports/fastx_$dir-$sample.tsv
			job reports_fastq-stats-$dir-$sample -deps $deps -targets {$target $target2 $target3} -vars {sample dir} -code {
				set gzcat [gzcat [lindex $deps 0]]
				exec -ignorestderr {*}$gzcat {*}$deps | fastq-stats -x $target3 > $target2
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
			job reports_vars-$sample -deps $deps -targets {
				$target
			} -vars {
				sample varfile targetfile refcodingfile dbdir build
			} -code {
				cg_report_vars -sample $sample -targetfile $targetfile -refcodingfile $refcodingfile $dep $target
			}
		}
	}
	if {[inlist $reports covered]} {
		foreach dep [jobglob $sampledir/sreg-*.tsv] {
			set sample [file root [file tail [gzroot $dep]]]
			regsub ^sreg- $sample {} sample
			set target $sampledir/reports/report_covered-$sample.tsv
			job reports_covered-$sample -deps {$dep} -targets {$target} -vars sample -code {
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

proc cg_process_reports {args} {
	set args [job_init {*}$args]
	if {[llength $args] < 1} {
		errorformat process_reports
	}
	process_reports_job {*}$args
	job_wait
}

#set bamroot hlongshot-sminimap2_splice-rr_NA05_055_v4.0.11+f1071ce_189859_hg38s
#set dir /work/rr/ontr/hg38s/exp157484-hg38s-cDNA-ONT-pilot-FUS/samples/rr_NA05_055_v4.0.11+f1071ce_189859_hg38s/reports
#set option unaligned
#
#set bamfile /work/rr/ontr/hg38s/exp157484-hg38s-cDNA-ONT-pilot-FUS/samples/rr_NA05_055_v4.0.11+f1071ce_189859_hg38s/map-hlongshot-sminimap2_splice-rr_NA05_055_v4.0.11+f1071ce_189859_hg38s.bam
#set target $sampledir/reports/samstats-$bamroot.stats.zst
#set target2 $sampledir/reports/report_samstats_summary-$bamroot.tsv.zst
#
proc reports_samstats {bamfile {option {}} {resultdir {}} {threads 1}} {
	upvar job_logdir job_logdir
	set bamroot [file root [file tail [gzroot $bamfile]]]
	regsub ^map- $bamroot {} bamroot
	if {$resultdir eq ""} {
		set resultdir [file dir $bamfile]/reports
	}
	if {![file exists $resultdir]} {file mkdir $resultdir}
	set dep $bamfile
	set target $resultdir/${option}samstats-$bamroot.stats.zst
	set target2 $resultdir/report_${option}samstats_summary-$bamroot.tsv.zst
	set targets [list $target $target2]
	set sections {
		{FFQ {} {First fragment qualities}}
		{LFQ {} {Last fragment qualities}}
		{GCF {gc_pct	count} {GC content of first fragments}}
		{GCL {gc_pct	count} {GC content of last fragments}}
		{GCC {cycle A_pct C_pxt G_pct T_pct N_pct O_pct} {ACGT content per cycle}}
		{GCT {cycle A_pct C_pxt G_pct T_pct} {ACGT content per cycle, read oriented}}
		{FBC {cycle A_pct C_pxt G_pct T_pct N_pct O_pct} {ACGT content per cycle for first fragments only}}
		{FTC {numA numC numG numT numN} {ACGT raw counters for first fragments}}
		{LBC {cycle A_pct C_pxt G_pct T_pct N_pct O_pct} {ACGT content per cycle for last fragments only}}
		{LTC {numA numC numG numT numN} {ACGT raw counters for last fragments}}
		{IS {insert_size pairs_total inward_oriented_pairs outward_oriented_pairs other_pairs} {Insert sizes}}
		{RL {read_length count} {Read lengths}}
		{FRL {read_length count} {Read lengths for first fragments only}}
		{LRL {read_length count} {Read lengths for last fragments only}}
		{ID {length number_of_insertions number_of_deletions} {Indel size distribution}}
		{IC {cycle number_of_insertions_fwd number_of_insertions_fwd_rev number_of_deletions_fwd number_of_deletions_rev} {Indels per cycle}}
		{COV {range start count} {Coverage (depth) distribution}}
		{GCD {GC_pct unique_sequence_percentiles depth_percentile_10 depth_percentile_25 depth_percentile_50 depth_percentile_75 depth_percentile_90} {GC-depth}}
	}
	foreach section [list_subindex $sections 0] {
		lappend targets $resultdir/${option}samstats_${section}-$bamroot.tsv.zst
	}
	lappend targets $resultdir/${option}samstats_FFQs-$bamroot.tsv.zst
	lappend targets $resultdir/${option}samstats_LFQs-$bamroot.tsv.zst
	if {$option eq ""} {set cores $threads} else {set cores 1}
	job [job_relfile2name reports_${option}samstats- $bamfile] -time 6:00:00 -cores $cores -deps {
		$dep
	} -targets $targets -vars {
		bamroot sections option threads
	} -code {
		list_foreach {key fields descr} $sections {
			set sectionsa($key) $fields
			set sectionsdescra($key) $descr
		}
		set reportsdir [file dir $target]
		if {![file exists $target]} {
			putslog "Making $target"
			analysisinfo_write $dep $target ${option}samstats_tool samtools ${option}samstats_version [version samtools]
			if {$option eq ""} {
				exec samtools stats -@ $threads $dep > $target.temp
			} elseif {$option eq "unaligned"} {
				exec samtools view -b -u -f 4 $dep | samtools stats > $target.temp
			} elseif {$option eq "aligned"} {
				exec samtools view -b -u -F 4 $dep | samtools stats > $target.temp
			} else {
				error "unknown option $option"
			}
			cg zst $target.temp
			file rename -force -- $target.temp.zst $target
		}
		putslog "Making $target2"
		catch {close $f} ; catch {close $o}
		set f [gzopen $target]
		set curtarget $target2
		set o [wgzopen $target2.temp.zst w]
		puts $o [join {sample source parameter value} \t]
		# do SN (summary)
		while {[gets $f line] != -1} {
			if {[regexp ^SN $line]} break
		}
		if {$option ne ""} {set aoption ${option}_} else {set aoption ""}
		while 1 {
			if {[regexp {^#} $line]} break
			set line [split $line \t]
			foreach {key parameter value} $line break
			if {$key ne "SN"} break
			set parameter [string_change [string trim $parameter] {{ } _ : {} \( {} \) {} % pct}]
			puts $o "$bamroot\t${option}samstats_summary\t${aoption}$parameter\t$value"
			if {[gets $f line] == -1} break
		}
		putslog "Making rest of samstats"
		# do rest
		while 1 {
			if {[regexp {^#} $line]} {
				gzclose $o
				file rename -force -- $curtarget.temp.zst $curtarget
				if {![regexp {grep \^([A-Z]+)} $line temp key]} {
					error "$target not expected format; could not extract key from: $line"
				}
				set curtarget $reportsdir/${option}samstats_${key}-$bamroot.tsv.zst
				set o [wgzopen $curtarget.temp.zst w]
				puts $o $line
				if {[gets $f line] == -1} break
				if {[regexp {^#} $line]} {
					if {[regexp {grep \^([A-Z]+)} $line temp key]} continue
					puts $o $line
				}
				if {[string range $key 0 2] in "BCC CRC OXC RXC"} {
					set sectionsa($key) {cycle A_pct C_pxt G_pct T_pct N_pct O_pct}
				}
				if {[string range $key 0 2] in "FFQ LFQ QTQ CYQ BZQ QXQ"} {
					if {[gets $f line] == -1} break
					if {[regexp {^#} $line]} continue
					set header [list cycle]
					set len [expr {[llength [split $line \t]] - 2}]
					for {set i 0} {$i < $len} {incr i} {
						lappend header q$i
					}
					puts $o [join $header \t]
				} else {
					if {[regexp {^#} $line]} continue
					if {[info exists sectionsa($key)]} {
						set header $sectionsa($key)
					} else {
						set header {{# unknown header}}
					}
					puts $o [join $header \t]
				}
			}
			set line [split $line \t]
			puts $o [join [lrange $line 1 end] \t]
			if {[gets $f line] == -1} break
		}
		gzclose $o
		close $f
		file rename -force -- $curtarget.temp.zst $curtarget
		foreach key {FFQ LFQ} {
			set f [gzopen $reportsdir/${option}samstats_${key}-$bamroot.tsv.zst]
			set header [tsv_open $f]
			if {[eof $f]} {
				file_write $reportsdir/${option}samstats_${key}s-$bamroot.tsv.zst {}
				continue
			}
			unset -nocomplain a
			foreach el $header {
				set a($el) 0.0
			}
			while {[gets $f line] != -1} {
				foreach v [split $line \t] q $header {
					set a($q) [expr {$a($q) + $v}]
				}
			}
			close $f
			set o [wgzopen $reportsdir/${option}samstats_${key}s-$bamroot.tsv.temp.zst]
			puts $o quality\tcount
			foreach el [lrange $header 1 end] {
				regsub q $el {} x
				puts $o $x\t$a($el)
			}
			gzclose $o
			file rename $reportsdir/${option}samstats_${key}s-$bamroot.tsv.temp.zst $reportsdir/${option}samstats_${key}s-$bamroot.tsv.zst
		}
		foreach file [lrange $targets 2 end] {
			if {![file exists $file]} {file_write $file ""}
		}
	}
}

proc reports_singlecell {sampledir} {
	set sample [file tail $sampledir]
	set rootname $sample
	if {[file exists $sampledir/reports/report_fastq_fw-$rootname.tsv]} {
		set total_reads [lindex [exec grep fw_numreads $sampledir/reports/report_fastq_fw-$rootname.tsv] end]
	} else {
		set total_reads NA
	}
	set countfile [gzfile $sampledir/umis_per_cell_raw-*$rootname.tsv $sampledir/umis_per_cell_raw*.tsv]
	if {[file exists $countfile]} {
		set cellbarcoded_reads [lindex [cg select -g all -gc sum(count) $countfile] end]
		set pct_cellbarcoded [format %.2f [expr {100.0*$cellbarcoded_reads/$total_reads}]]
	} else {
		set cellbarcoded_reads NA
		set pct_cellbarcoded NA
	}
	set countfile [gzfile $sampledir/sc_gene_counts_raw-*$rootname.tsv]
	if {[file exists $countfile]} {
		set temp [cg select -g all -gc sum(count) $countfile]
		set rawgenecount_reads [lindex $temp end]
		set pct_rawgenecount_reads [format %.2f [expr {100.0*$rawgenecount_reads/$total_reads}]]
	} else {
		set rawgenecount_reads NA
		set pct_rawgenecount_reads NA
	}
	set countsfile [gzfile $sampledir/sc_gene_counts_filtered-*$rootname.tsv]
	if {[file exists $countsfile]} {
		set temp [cg select -g all -gc sum(count) $countsfile]
		set filteredgenecount_reads [lindex $temp end]
		set pct_filteredgenecount_reads [format %.2f [expr {100.0*$filteredgenecount_reads/$total_reads}]]
		set genes_percell [cg select -g {cell * gene *} $countsfile | cg select -g cell]
		set temp [lrange [list_subindex [split $genes_percell \n] 1] 1 end]
		set mean_genes_percell [lmath_average $temp]
		set mediangenes [median $temp]
	} else {
		set filteredgenecount_reads NA
		set pct_filteredgenecount_reads NA
		set mean_genes_percell NA
		set mediangenes NA
	}
	set countsfile [gzfile $sampledir/sc_isoform_counts_filtered-*$rootname.tsv]
	if {[file exists $countsfile]} {
		set temp [cg select -g all -gc sum(counts_weighed) $countsfile]
		set filteredisoformcount_reads [lindex $temp end]
		set pct_filteredisoformcount_reads [format %.2f [expr {100.0*$filteredisoformcount_reads/$total_reads}]]
		set isoforms_percell [cg select -g {cell * transcript *} $countsfile | cg select -g cell]
		set temp [lrange [list_subindex [split $isoforms_percell \n] 1] 1 end]
		set mean_isoforms_percell [lmath_average $temp]
		set median_isoforms_percell [median $temp]
	} else {
		set filteredisoformcount_reads NA
		set pct_filteredisoformcount_reads NA
		set mean_isoforms_percell NA
		set median_isoforms_percell NA
	}
	set countsfile [gzfile $sampledir/sc_cellinfo_filtered-*$rootname.tsv sc_cellinfo_filtered-*.tsv]
	if {[file exists $countsfile]} {
		set readcounts [cg select -sh /dev/null -q {$is_cell} -f readcount $countsfile]
		set mean_readcounts [lmath_average $readcounts]
		set median_readcounts [median $readcounts]
	} else {
		set mean_readcounts NA
		set median_readcounts NA
	}
	#
	set countsfile [gzfile $sampledir/sc_gene_counts_filtered-*$rootname.tsv sc_gene_counts_filtered-*.tsv]
	set nrcells [lindex [cg select -g cell $countsfile | cg select -g all ] end]
	# write file
	set o [open $sampledir/reports/report_singlecell-$sample.tsv w]
	puts $o [join {sample source parameter value} \t]
	puts $o [join [list $sample genomecomb sc_total_reads $total_reads] \t]
	puts $o [join [list $sample genomecomb sc_cellbarcoded_reads $cellbarcoded_reads] \t]
	puts $o [join [list $sample genomecomb sc_pct_cellbarcoded $pct_cellbarcoded] \t]
	puts $o [join [list $sample genomecomb sc_rawgenecount_reads $rawgenecount_reads] \t]
	puts $o [join [list $sample genomecomb sc_pct_rawgenecount_reads $pct_rawgenecount_reads] \t]
	puts $o [join [list $sample genomecomb sc_filteredgenecount_reads $filteredgenecount_reads] \t]
	puts $o [join [list $sample genomecomb sc_pct_filteredgenecount_reads $pct_filteredgenecount_reads] \t]
	puts $o [join [list $sample genomecomb sc_filteredisoformcount_reads $filteredisoformcount_reads] \t]
	puts $o [join [list $sample genomecomb sc_pct_filteredisoformcount_reads $pct_filteredisoformcount_reads] \t]
	puts $o [join [list $sample genomecomb sc_nrcells $nrcells] \t]
	puts $o [join [list $sample genomecomb sc_mean_readcounts $mean_readcounts] \t]
	puts $o [join [list $sample genomecomb sc_median_readcounts $median_readcounts] \t]
	puts $o [join [list $sample genomecomb sc_mean_genes_percell $mean_genes_percell] \t]
	puts $o [join [list $sample genomecomb sc_median_genes_percell $mediangenes] \t]
	puts $o [join [list $sample genomecomb sc_mean_isoforms_percell $mean_isoforms_percell] \t]
	puts $o [join [list $sample genomecomb sc_median_isoforms_percell $median_isoforms_percell] \t]
	close $o
	set o [open $sampledir/reports/report_singlecell-cellgrouping-$sample.tsv w]
	puts $o [join {sample source parameter group value} \t]
	foreach groupfile [gzfiles $sampledir/sc_group-*.tsv] {
		set data [string trim [cg select -g group $groupfile | cg select -sh /dev/null]]
		set root [join [lrange [split [file_rootname $groupfile] -] 0 end-1] -]
		foreach {type count} [split $data \n\t] {
			puts $o [join [list $sample $root sc_nrcells $type $count] \t]
		}
	}
	close $o
	set o [open $sampledir/reports/singlecell-pseudobulk_isoforms-$sample.tsv w]
	puts $o [join {sample source group nrisoforms count} \t]
	foreach pbfile [gzfiles $sampledir/pb_isoform_counts-*.tsv] {
		set data [cg long $pbfile \
			| cg select -q {$counts_weighed > 0} \
			| cg select -g {sample * gene *} \
			| cg select -g {sample * count *} \
		]
		set root [join [lrange [split [file_rootname $pbfile] -] 0 end-1] -]
		foreach line [lrange [split $data \n] 1 end] {
			foreach {sample nrisoforms count} [split $line \t] break
			set group [lindex [split $sample -] 0]
			set sample [lindex [split $sample -] end]
			puts $o [join [list $sample $root $group $nrisoforms $count] \t]
		}
	}
	close $o
}

proc reports_expand {reports} {
	set allreports {fastqstats fastqc flagstat_reads samstats alignedsamstats unalignedsamstats histodepth hsmetrics vars covered histo predictgender}
	set basicreports {fastqstats fastqc flagstat_reads samstats histodepth hsmetrics vars covered histo predictgender}
	if {$reports eq "all"} {
		set reports $allreports
	} elseif {[string index $reports 0] eq "-"} {
		set reports [list_lremove $allreports [string range $reports 1 end]]
	} elseif {"basic" in $reports} {
		set pos [lsearch $reports all]
		set reports [list_remdup [lreplace $reports $pos $pos {*}$basicreports]]
	} else {
		return $reports
	}
}

proc process_reports_job {args} {
	set cmdline [clean_cmdline cg process_reports {*}$args]
	set reports basic
	set dbdir {}
	set resultbamfile {}
	set paired 1
	set depth_histo_max 1000
	set threads 1
	cg_options process_reports args {
		-dbdir {
			set dbdir $value
		}
		-r - -reports {
			set reports $value
		}
		-paired {
			set paired $value
		}
		-depth_histo_max {
			set depth_histo_max $value
		}
		-resultbamfile {
			set resultbamfile [file_absolute $value]
		}
		-threads {
			set threads $value
		}
	} {sampledir dbdir reports} 1 3 {
		Calculates a number of statistics on a sample in the reports subdir
	}
	set dbdir [dbdir $dbdir]
	set sampledir [file_absolute $sampledir]
	set sample [file tail $sampledir]
	set_job_logdir $sampledir/log_jobs
	set ref [dbdir_ref $dbdir]
	set reports [reports_expand $reports]
	# find regionfile indicating target of sequencing (used by hsmetrics, histodepth, vars, so needs to be here)
	set targetfile [targetfile_job $sampledir $dbdir]
	# logfile
	job_logfile $sampledir/process_reports_$sample $sampledir $cmdline \
		{*}[versions dbdir fastq-stats samtools gnusort8 zst os]
	# start
	# get bamfiles; do not include hlongshot- ones if source exists (are the same anyway except one tag that is not used here)
	unset -nocomplain a
	foreach file [bsort [jobglob $sampledir/*.bam $sampledir/*.cram]] {
		set a([file tail $file]) $file
	}
	set bamfiles {}
	foreach file [array names a] {
		if {[regexp ^map-hlongshot- $file]} {
			regsub ^map-hlongshot- $file map- temp
			if {[info exists a($temp)]} continue
		}
		lappend bamfiles $a($file)
	}
	if {$resultbamfile eq ""} {
		set bamfiles [bsort $bamfiles]
		set resultbamfile [lindex $bamfiles end]
	} else {
		list_addnew bamfiles $resultbamfile
	}
	set ampliconsfile [ampliconsfile $sampledir $ref]
	file mkdir $sampledir/reports
	foreach bamfile $bamfiles {
		set bamroot [file root [file tail [gzroot $bamfile]]]
		regsub ^map- $bamroot {} bamroot
		if {[inlist $reports flagstat_alignments]} {
			set dep $bamfile
			set target $sampledir/reports/flagstat_alignments-$bamroot.flagstat
			set target2 $sampledir/reports/report_flagstat_alignments-$bamroot.tsv
			job [job_relfile2name reports_flagstat_alignments- $bamfile] -time 5:00:00 -cores $threads -deps {
				$dep
			} -targets {
				$target $target2
			} -vars {
				bamroot threads
			} -code {
				analysisinfo_write $dep $target flagstat_tool samtools flagstat_version [version samtools]
				exec samtools flagstat -@ $threads $dep > $target.temp
				file rename -force -- $target.temp $target
				set o [open $target2.temp w]
				puts $o [join {sample source parameter value value_qcfail} \t]
				set f [open $target]
				while {[gets $f line] != -1} {
					if {![regexp {^([0-9]+) \+ ([0-9]+) ([^()]*)} $line temp value value_qcfail parameter]} continue
					puts $o "$bamroot\tflagstat_alignments\t[string trim $parameter]\t$value\t$value_qcfail"
				}
				close $f
				close $o
				file rename -force -- $target2.temp $target2
			}
		}
		if {[inlist $reports alignedsamstats]} {
			reports_samstats $bamfile aligned $sampledir/reports $threads
		}
		if {[inlist $reports unalignedsamstats]} {
			reports_samstats $bamfile unaligned $sampledir/reports
		}
		if {[inlist $reports samstats]} {
			reports_samstats $bamfile {} $sampledir/reports
		}
		if {[inlist $reports flagstat_reads]} {
			set dep $bamfile
			set target $sampledir/reports/flagstat_reads-$bamroot.flagstat
			set target2 $sampledir/reports/report_flagstat_reads-$bamroot.tsv
			job [job_relfile2name reports_flagstat_reads- $bamfile] -time 5:00:00 -deps {
				$dep
			} -targets {
				$target $target2
			} -vars {bamroot} -code {
				analysisinfo_write $dep $target flagstat_tool samtools flagstat_version [version samtools]
				if {[catch {
					exec samtools view --no-PG -F 256 -h -b $dep | samtools flagstat - > $target.temp
				} msg]} {
					if {$msg ne "\[bam_header_read\] EOF marker is absent. The input is probably truncated."} {error $msg}
				}
				file rename -force $target.temp $target
				set o [open $target2.temp w]
				puts $o [join {sample source parameter value value_qcfail} \t]
				set f [open $target]
				while {[gets $f line] != -1} {
					if {![regexp {^([0-9]+) \+ ([0-9]+) ([^()]*)} $line temp value value_qcfail parameter]} continue
					puts $o "$bamroot\tflagstat_reads\t[string trim $parameter]\t$value\t$value_qcfail"
				}
				close $f
				close $o
				file rename -force -- $target2.temp $target2
			}
		}
		if {[inlist $reports hsmetrics] && [jobfileexists $targetfile]} {
			set dep1 $bamfile
			set dep2 $targetfile
			set dep3 $bamfile.[indexext $bamfile]
			set target $sampledir/reports/hsmetrics-$bamroot.hsmetrics
			set target2 $sampledir/reports/report_hsmetrics-$bamroot.tsv
			job [job_relfile2name reports_hsmetrics- $bamfile] -optional 1 -deps {
				$dep1 $dep2 $dep3
			} -targets {
				$target $target2
			} -vars {
				bamroot bamfile targetfile dbdir
			} -code {
				analysisinfo_write $dep $target hsmetrics_tool picard hsmetrics_version [version picard]
				cg_hsmetrics -refseq $dbdir --sample $bamroot $dep1 $dep2 $target
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
				file rename -force -- $target2temp $target2
			}
		}
		if {[inlist $reports histodepth]} {
			set dep1 $bamfile
			set dep2 $targetfile
			set target $sampledir/reports/histodepth-$bamroot.tsv
			set target2 $sampledir/reports/report_histodepth-$bamroot.tsv
			set indexext [indexext $bamfile]
			job [job_relfile2name reports_histodepth- $bamfile] -optional 1 -time 5:00:00 -deps {
				$dep1 ($dep2) $dep1.$indexext
			} -targets {
				$target $target2
			} -vars {
				bamroot bamfile depth_histo_max
			} -code {
				analysisinfo_write $dep $target histodepth_tool genomecomb histodepth_version [version genomecomb]
				set targettemp [filetemp $target]
				if {![file exists $dep2]} {
					set targetfile {}
					set tottarget 0
				} else {
					set targetfile $dep2
					set tottarget [lindex [cg covered $targetfile] end]
				}
				cg depth_histo -max $depth_histo_max -q 0 -Q 0 $dep1 $targetfile > $targettemp
				file rename -force -- $targettemp $target
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
				append out [join [bsort $result] \n]
				append out \n
				set target2temp [filetemp $target2]
				file_write $target2temp $out
				file rename -force -- $target2temp $target2
			}
		}
		if {[inlist $reports histo] && $ampliconsfile ne ""} {
			set dep1 $bamfile
			set dep2 $ampliconsfile
			set target $sampledir/reports/$bamroot.histo
			set indexext [indexext $bamfile]
			job [job_relfile2name reports_histo- $bamfile] -optional 1 -deps {
				$dep1 $dep2 $dep1.$indexext
			} -targets {$target} -vars {bamroot bamfile} -code {
				analysisinfo_write $dep $target histo_tool genomecomb histo_version [version genomecomb]
				set tempfile [filetemp $target]
				set tempregionfile [filetemp $target]
				cg regcollapse $dep2 > $tempregionfile
				cg bam_histo $tempregionfile $dep {1 5 10 20 50 100 200 500 1000} > $tempfile
				file delete $tempregionfile
				file rename -force -- $tempfile $target
			}
		}
	}
	if {[inlist $reports predictgender]} {
		set target $sampledir/reports/report_predictgender-$sample.tsv
		set varfile [jobglob $sampledir/var-*[file_rootname $resultbamfile].tsv]
		set indexfile $resultbamfile.[indexext $resultbamfile]
		job predictgender-[file_rootname $resultbamfile] -optional 1 -deps {$resultbamfile $indexfile ($varfile)} -vars {resultbamfile dbdir sampledir} -targets {$target} -code {
			analysisinfo_write $dep $target predictgender_tool genomecomb predictgender_version [version genomecomb]
			cg predictgender -dbdir $dbdir $resultbamfile $target
		}
	}
	set fastqfiles [bsort [jobglob $sampledir/fastq/*.fastq.gz $sampledir/fastq/*.fastq $sampledir/fastq/*.fq.gz $sampledir/fastq/*.fq]]
	if {$paired} {
		set fastqfiles_fw [list_unmerge $fastqfiles 1 fastqfiles_rev]
	} else {
		set fastqfiles_fw $fastqfiles
		set fastqfiles_rev {}
	}
	if {[inlist $reports fastqc] && [llength $fastqfiles]} {
		foreach deps [list $fastqfiles_fw $fastqfiles_rev] dir {fw rev} {
			if {![llength $deps]} continue
			set target $sampledir/reports/fastqc_$dir-$sample.fastqc
			job reports_fastqc-$dir-$sample -time 2:00:00 -deps $deps -targets {
				$target $sampledir/reports/fastqc_$dir-$sample.fastqc/fastqc_data.txt
			} -code {
				analysisinfo_write $dep $target fastqc_version [version fastqc]
				file mkdir $target.temp
				set gzcat [gzcat [lindex $deps 0]]
				set foundreads 0
				foreach fastq $deps {
					if {[fastq_size $fastq 1]} {set foundreads 1 ; break}
				}
				if {!$foundreads} {
					file rename -force -- $target.temp $target
					file_write $target2 ""
				} else {
					exec -ignorestderr {*}$gzcat {*}$deps | fastqc -o $target.temp stdin
					exec unzip -o $target.temp/stdin_fastqc.zip -d $target.temp
					file rename -force -- {*}[glob $target.temp/stdin_fastqc/*] $target.temp
					file delete $target.temp/stdin_fastqc
					file delete $target.temp/stdin_fastqc.zip
					file delete $target.temp/stdin_fastqc.html
					if {[file exists $target]} {file delete -force $target}
					file rename -force -- $target.temp $target
				}
			}
		}
	}
	if {[inlist $reports fastqstats] && [llength $fastqfiles]} {
		foreach deps [list $fastqfiles_fw $fastqfiles_rev] dir {fw rev} {
			if {![llength $deps]} continue
			set target $sampledir/reports/report_fastq_$dir-$sample.tsv
			set target2 $sampledir/reports/fastq_stats_$dir-$sample.txt
			set target3 $sampledir/reports/fastx_$dir-$sample.tsv
			job reports_fastq-stats-$dir-$sample -time 5:00:00 -deps $deps -targets {$target $target2 $target3} -vars {sample dir} -code {
				analysisinfo_write $dep $target2 fastq_stats_version [version fastq-stats]
				analysisinfo_write $dep $target3 fastq_stats_version [version fastq-stats]
				set foundreads 0
				foreach fastq $deps {
					if {[fastq_size $fastq 1]} {set foundreads 1 ; break}
				}
				if {!$foundreads} {
					file_write $target1 ""
					file_write $target2 ""
					file_write $target3 ""
				} else {
					set gzcat [gzcat [lindex $deps 0]]
					if {[llength $deps] >= 1000} {
						set o [open [list | fastq-stats -x $target3 > $target2] w]
						foreach file $deps {
							exec {*}$gzcat $file >@ $o
						}
						close $o
					} else {
						exec -ignorestderr {*}$gzcat {*}$deps | fastq-stats -x $target3 > $target2
					}
					analysisinfo_write $dep $target fastq_stats_version [version fastq-stats]
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
					file rename -force -- $target.temp $target
				}
			}
		}
	}
	if {[inlist $reports singlecell]} {
		set deps [list]
		lappend deps [jobgzfile $sampledir/reports/report_fastq_fw-*$sample.tsv]
		lappend deps [jobgzfile $sampledir/umis_per_cell_raw-*$sample.tsv $sampledir/umis_per_cell_raw*.tsv]
		lappend deps [jobgzfile $sampledir/sc_gene_counts_raw-*$sample.tsv]
		lappend deps [jobgzfile $sampledir/sc_gene_counts_filtered-*$sample.tsv sc_gene_counts_filtered-*.tsv]
		lappend deps [jobgzfile $sampledir/sc_cellinfo_filtered-*$sample.tsv sc_cellinfo_filtered-*.tsv]
		# set deps "([join $deps ") ("])"
		job reports_singlecell-$sample -deps $deps -targets {
			$sampledir/reports/report_singlecell-$sample.tsv
			$sampledir/reports/report_singlecell-cellgrouping-$sample.tsv
			$sampledir/reports/singlecell-pseudobulk_isoforms-$sample.tsv
		} -vars {
			sampledir
		} -code {
			reports_singlecell $sampledir
		}
	}
	if {[inlist $reports vars]} {
		set refcodingfile [gzfile $dbdir/extra/reg_*_refcoding.tsv]
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
				analysisinfo_write $dep $target report_vars_tool genomecomb report_vars_version [version genomecomb]
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
				analysisinfo_write $dep $target report_covered_tool genomecomb report_covered_version [version genomecomb]
				set f [open $target.temp w]
				puts $f [join {sample source parameter value} \t]
				if {[llength [cg select -h $dep]] != 0} {
					set temp [split [exec cg covered $dep] \n]
					foreach line [lrange $temp 1 end] {
						foreach {chr cov} $line break
						puts $f $sample\tgenomecomb\tcovered_$chr\t$cov
					}
				}
				close $f
				file rename -force -- $target.temp $target
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

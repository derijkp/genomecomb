proc cg_meth_nanopolish_freqs {dep target {callthreshold 2.5}} {
	analysisinfo_write $dep $target meth_nanopolish_cg_version [version genomecomb]
	catch {close $o} ; catch {close $f}
	set f [gzopen $dep]
	set header [tsv_open $f]
	if {$header ne "chromosome strand start end read_name log_lik_ratio log_lik_methylated log_lik_unmethylated num_calling_strands num_motifs sequence"} {
		error "wrong header in file $dep"
	}
	set o [wgzopen $target]
	puts $o [join {chromosome start end num_motifs_in_group called_sites called_sites_methylated methylated_frequency group_sequence} \t]
	set curloc {}
	set num_reads 0
	set called_sites 0
	set called_sites_methylated 0
	while 1 {
		set read [gets $f line]
		set line [split $line \t]
		set loc [list_sub $line {0 2 3}]
		if {$loc ne $curloc} {
			if {$curloc ne ""} {
				incr end
				if {$called_sites > 0} {
					set methylated_frequency [expr {1.0*$called_sites_methylated/$called_sites}]
				} else {
					set methylated_frequency 0.0
				}
				puts $o [join [list $chromosome $begin $end \
					$num_motifs $called_sites $called_sites_methylated \
					[format %.3f $methylated_frequency] $sequence] \t]
			}
			set num_reads 0
			set called_sites 0
			set called_sites_methylated 0
			set curloc $loc
		}
		if {$read == -1} break
		foreach {chromosome strand begin end read_name log_lik_ratio log_lik_methylated log_lik_unmethylated num_calling_strands num_motifs sequence} $line break
	        # Skip ambiguous call
		if {abs($log_lik_ratio) < $callthreshold} continue
		incr numreads
		incr called_sites $num_motifs
		if {$log_lik_ratio > 0} {
			incr called_sites_methylated $num_motifs
		}
	}
	gzclose $o
	gzclose $f
}

proc meth_nanopolish_job {args} {
	# putslog [list meth_nanopolish_job {*}$args]
	global appdir
	upvar job_logdir job_logdir
	set cmdline "[list cd [pwd]] \; [list cg meth_nanopolish {*}$args]"
	set refseq {}
	set skips {}
	set resultfile {}
	set basecaller {}
	set callthreshold 2.5
	cg_options meth_nanopolish args {
		-callthreshold {
			set callthreshold $value
		}
		-refseq {
			set refseq $value
		}
		-skip {
			lappend skips -skip $value
		}
		-opts {
			set opts $value
		}
	} {fast5dir fastqdir bamfile resultfile} 3 4
	set bamfile [file_absolute $bamfile]
	set bamtail [file tail $bamfile]
	set fast5dir [file_absolute $fast5dir]
	set fastqdir [file_absolute $fastqdir]
	set refseq [refseq $refseq]
	if {$resultfile eq ""} {
		set root nanopolish-[file_rootname $bamtail]
		set resultfile [file dir $bamfile]/meth-$root.tsv.zst
	} else {
		set resultfile [file_absolute $resultfile]
	}
	file mkdir $resultfile.temp
	set destdir [file dir $resultfile]
	set tail [file tail $resultfile]
	if {[regexp ^meth- $tail]} {
		set smethfile [file dir $resultfile]/s$tail
	} else {
		set smethfile [file dir $resultfile]/smeth-$tail
	}
	job_logfile $destdir/meth_nanopolish_$bamtail $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools gatk picard java gnusort8 zst os]
	# start
	set keeppwd [pwd]
	set fastqfiles [gzfiles $fastqdir/*.fastq $fastqdir/*.fq]
	if {[file exists $fastqdir/info_basecaller.txt]} {
		set basecaller [file_read $fastqdir/info_basecaller.txt]
	}
	set seqmethfiles {}
	foreach fastqfile $fastqfiles {
		set fast5file $fast5dir/[file root [file tail [gzroot $fastqfile]]].fast5
		if {![jobfileexists $fast5file]} {
			puts stderr "Could not find fast5 for $fastqfile, skipping"
			continue
		}
		set root [file root [file tail [gzroot $fastqfile]]]
		set target $resultfile.temp/$root.smeth.tsv.zst
		lappend seqmethfiles $target
		job seqmeth_nanopolish-[file tail $target] {*}$skips -deps {
			$fastqfile $fast5file $refseq $bamfile
		} -targets {
			$target
		} -skip {$smethfile} -skip {$resultfile} -vars {
			fastqfile fast5file refseq bamfile basecaller
		} -code {
			analysisinfo_write $dep $target basecaller_version $basecaller meth_caller nanopolish meth_caller_version [version nanopolish]
			set tempdir [tempdir]
			foreach file [glob -nocomplain $tempdir/*] {file delete $file}
			mklink $fastqfile $tempdir/[file tail $fastqfile]
			mklink $fast5file $tempdir/[file tail $fast5file]
			catch_exec nanopolish index -d $tempdir $tempdir/[file tail $fastqfile]
			set error [catch {
				exec nanopolish call-methylation -r $tempdir/[file tail $fastqfile] \
					-b $bamfile -g $refseq | cg zst --compressionlevel 1 > $target.temp.zst
			} msg opt]
			if {$error} {
				if {$::errorCode ne "NONE"} {
					dict unset opt -level
					set errorInfo "$msg\n    while executing\nnanopolish call-methylation -r $tempdir/[file tail $fastqfile] -b $bamfile -g $refseq | cg zst --compressionlevel 1 > $target.temp.zst"
					return -code $error -errorcode $::errorCode -errorinfo $errorInfo $msg
				} else {
					puts stderr $msg
				}
			}
			result_rename $target.temp.zst $target
		}
	}
	set target $smethfile
	set root [file root [file tail [gzroot $smethfile]]]
	job meth_nanopolish_smethfinal-$root {*}$skips -deps $seqmethfiles -targets {
		$smethfile
	} -skip {$target} -code {
		analysisinfo_write $dep $target smeth_nanopolish_cg_version [version genomecomb]
		cg cat -c 0 {*}$deps | cg select -s - | cg zst > $target.temp
		file rename -- $target.temp $target
	}

	set root [file root [file tail [gzroot $resultfile]]]
	set dep $smethfile
	set target $resultfile
	job meth_nanopolish_methfinal-$root {*}$skips -deps {
		$dep
	} -targets {
		$target
	} -vars {
		callthreshold
	} -code {
		set tempresult [filetemp_ext $target]
		cg_meth_nanopolish_freqs $dep $tempresult $callthreshold
		result_rename $tempresult $target
	}
	return [list $resultfile $smethfile]
}

proc cg_meth_nanopolish {args} {
	set args [job_init {*}$args]
	meth_nanopolish_job {*}$args
	job_wait
}

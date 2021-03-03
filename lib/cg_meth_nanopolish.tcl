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
	if {[file tail $target] in ".gz .bgz"} {
		cg maketabix $target
	}
}

proc fastqs_mergename {fastqfiles} {
	set pbase [file_root [file tail [lindex $fastqfiles 0]]]
	if {[llength $fastqfiles] == 1} {
		return $pbase
	} else {
		set last [file_root [file tail [lindex $fastqfiles end]]]
		if {$last ne $pbase} {
			set pos 0
			string_foreach c1 $pbase c2 $last {
				if {$c1 ne $c2} break
				incr pos
			}
			set size [string length $pbase]
			set lastlen [string length $last]
			# -2 for cnnector (__), -30 for extensions and prefix
			set minpos [expr {$lastlen - (255-2-30-$size)}]
			if {$minpos > $pos} {set pos $minpos}
			append pbase __[string range $last $pos end]
		}
		return $pbase
	}
}

proc meth_nanopolish_distrfast5 {fast5dir fastqdir bamfile resultfile refseq skips basecaller callthreshold threads maxfastqdistr meth-compression nomatchingfast5error opts} {
	# putslog [list meth_nanopolish_job {*}$args]
	global appdir
	upvar job_logdir job_logdir
	# ignore nr threads, does not scale well enough with threads
	set threads 1
	file mkdir $resultfile.temp
	set destdir [file dir $resultfile]
	set tail [file tail $resultfile]
	if {[regexp ^meth- $tail]} {
		set smethfile [file dir $resultfile]/s$tail
	} else {
		set smethfile [file dir $resultfile]/smeth-$tail
	}
	# start
	set fastqfiles [gzfiles $fastqdir/*.fastq $fastqdir/*.fq]
	if {[file exists $smethfile] && ![jobtargetexists $smethfile $fastqfiles]} {
		putslog "$smethfile older than one of fastqfiles (renaming to .old)"
		file rename -force $smethfile $smethfile.old
	}
	if {[file exists $resultfile] && ![jobtargetexists $resultfile $fastqfiles]} {
		putslog "$resultfile older than one of fastqfiles (renaming to .old)"
		file rename -force $resultfile $resultfile.old
	}
	if {[file exists $fastqdir/info_basecaller.txt]} {
		set basecaller [file_read $fastqdir/info_basecaller.txt]
	}
	set seqmethfiles {}
	set processlist {}
	if {[isint $maxfastqdistr]} {
		set len [llength $fastqfiles]
		set perbatch [expr {($len + $maxfastqdistr -1)/$maxfastqdistr}]
		set pos 0
		for {set pos 0} {$pos < $len} {incr pos $perbatch} {
			lappend processlist [lrange $fastqfiles $pos [expr {$pos + $perbatch - 1}]]
		}
	} else {
		foreach {fastq} $fastqfiles {
			lappend processlist [list $fastq]
		}
	}
	set bamfileindex [index_file $bamfile]
	bam_index_job $bamfile
	if {[info exists ::env(CG_FAST_SHAREDSTORAGE)]} {
		set cachedir $::env(CG_FAST_SHAREDSTORAGE)/methcache
		set bamcache $cachedir/[file tail $bamfile]
		set bamcacheindex [index_file $bamcache]
		putslog "meth_nanopolish using cached bam: $bamcache"
		job [job_relfile2name meth_nanopolish_makecache- $bamfile] {*}$skips -skip {$smethfile} -skip {$resultfile} -deps {
			$bamfile
		} -targets {
			$bamcache
		} -vars {
			cachedir bamfileindex bamcacheindex
		} -code {
			catch {file mkdir $cachedir}
			set tempfile [filetemp $bamcacheindex]
			if {[catch {file link $dep} link]} {
				file copy -force $bamfileindex $tempfile
			} else {
				file copy -force $link $tempfile
			}
			file rename -force $tempfile $bamcacheindex
			exec touch -h -d [clock format [file mtime $bamfileindex]] $bamcacheindex
			set tempfile [filetemp $target]
			if {[catch {file link $dep} link]} {
				file copy -force $dep $tempfile
			} else {
				file copy -force $link $tempfile
			}
			file rename -force $tempfile $target
			exec touch -h -d [clock format [file mtime $dep]] $target
		}
	} else {
		set bamcache $bamfile
	}
	foreach fastqfiles $processlist {
		set usefastqfiles {}
		set fast5files {}
		foreach fastqfile $fastqfiles {
			set fast5file $fast5dir/[file root [file tail [gzroot $fastqfile]]].fast5
			if {![jobfileexists $fast5file]} {
				if {$nomatchingfast5error} {
					error "Could not find fast5 for $fastqfile"
				} else {
					puts stderr "Could not find fast5 for $fastqfile, skipping"
					continue
				}
			}
			lappend fast5files $fast5file
			lappend usefastqfiles $fastqfile
		}
		set root [fastqs_mergename $usefastqfiles]
		set target $resultfile.temp/$root.smeth.tsv.zst
		lappend seqmethfiles $target
		set deps [list {*}$usefastqfiles {*}$fast5files $refseq $bamcache ($bamfileindex)]
		job [job_relfile2name seqmeth_nanopolish- $target] {*}$skips -cores $threads -deps $deps -targets {
			$target
		} -skip {$smethfile} -skip {$resultfile} -vars {
			usefastqfiles fast5files refseq bamfile bamcache basecaller threads bamfileindex opts
		} -code {
			analysisinfo_write $dep $target basecaller_version $basecaller meth_caller nanopolish meth_caller_version [version nanopolish]
			set tempdir [tempdir]
			if {[llength $fast5files] == 1 && [file exists /dev/shm]} {
				# check if we can use ramdisk
				set fastqfile [lindex $usefastqfiles 0]
				set fast5file [lindex $fast5files 0]
				set size [expr {[file size $fast5file]+1.6*[file size $fastqfile]}]
				# tempdir will only be changed if enough space was free on /dev/shm
				catch {
					set tempdir [tempramdir $size]
				}
				file mkdir $tempdir
			}
			foreach file [glob -nocomplain $tempdir/*] {file delete $file}
			exec cp -fL {*}$fast5files $tempdir
			if {[llength $usefastqfiles] == 1} {
				set fastqfile [lindex $usefastqfiles 0]
				exec cp -fL $fastqfile $tempdir
				set fastqfile $tempdir/[file tail $fastqfile]
			} else {
				set fastqfile [tempfile].fastq.gz
				exec cg zcat {*}$usefastqfiles | gzip --fast > $fastqfile 2>@ stderr
			}
			set nanopolishprog [exec which nanopolish]
			if {[file exists $nanopolishprog.plugins]} {
				if {[get ::env(HDF5_PLUGIN_PATH) ""] ne ""} {
					set ::env(HDF5_PLUGIN_PATH) $nanopolishprog.plugins:$::env(HDF5_PLUGIN_PATH)
				} else {
					set ::env(HDF5_PLUGIN_PATH) $nanopolishprog.plugins
				}
			}
			catch_exec nanopolish index -d $tempdir $tempdir/[file tail $fastqfile]
			set error [catch {
				exec nanopolish call-methylation \
					{*}$opts \
					-t $threads -r $fastqfile \
					-b $bamcache -g $refseq | cg zst --compressionlevel 1 > $target.temp.zst 2>@ stderr
			} msg opt]
			if {$error} {
				if {$::errorCode ne "NONE"} {
					dict unset opt -level
					set errorInfo "$msg\n    while executing\nnanopolish call-methylation -r $tempdir/[file tail $fastqfile] -b $bamcache -g $refseq | cg zst --compressionlevel 1 > $target.temp.zst"
					return -code $error -errorcode $::errorCode -errorinfo $::errorInfo $msg
				} else {
					puts stderr $msg
				}
			}
			result_rename $target.temp.zst $target
			file delete -force $tempdir
		}
	}
	set target $smethfile
	set root [file root [file tail [gzroot $smethfile]]]
	job meth_nanopolish_smethfinal-$root {*}$skips -deps $seqmethfiles -targets {
		$smethfile
	} -vars {
		bamfile bamcache bamcacheindex meth-compression
	} -skip {$target} -code {
		analysisinfo_write $dep $target smeth_nanopolish_cg_version [version genomecomb]
		cg cat -c 0 {*}$deps | cg select -s - | cg ${meth-compression} > $target.temp
		file rename -- $target.temp $target
		if {${meth-compression} in "gz bgz"} {
			cg maketabix $target
		}
		if {$bamcache ne $bamfile} {
			file delete -force $bamcache
			file delete -force $bamcacheindex
		}
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
		file delete -force $target.temp
	}
	return [list $resultfile $smethfile]
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
	set callthreshold [get ::specialopt(-meth_nanopolish-callthreshold) 2.5]
	set threads 1
	set distrmethod fast5
	set distrreg 5000000
	set maxfastqdistr {}
	set nomatchingfast5error 1
	set meth-compression [get ::specialopt(-meth_nanopolish-compression) zst]
	set opts {}
	set callerroot nanopolish
	cg_options meth_nanopolish args {
		-callthreshold {
			set callthreshold $value
		}
		-preset {
			if {![inlist {cpg gpc dam dcm} $value]} {
				error "-preset $value unknown, must be one of: cpg gpc dam dcm"
			}
			lappend opts --methylation=$value
			append callerroot _$value
		}
		-refseq {set maxfastqdistr {}
			set refseq $value
		}
		-skip {
			lappend skips -skip $value
		}
		-threads {
			set threads $value
		}
		-distrmethod {
			set distrmethod $value
		}
		-distrreg {
			set distrreg $value
		}
		-maxfastqdistr {
			set maxfastqdistr $value
		}
		-meth-compression {
			set meth-compression $value
		}
		-nomatchingfast5error {
			set nomatchingfast5error [tru $value]
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
		set root $callerroot-[file_rootname $bamtail]
		set resultfile [file dir $bamfile]/meth-$root.tsv.${meth-compression}
	} else {
		set resultfile [file_absolute $resultfile]
	}
	file mkdir $resultfile.temp
	set fastqfile $resultfile.temp/merged-fastq.fastq.gz
	set destdir [file dir $resultfile]
	job_logfile $destdir/meth_nanopolish_[file root $bamtail] $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools gatk picard java gnusort8 ${meth-compression} os]
	# start
	set keeppwd [pwd]
	set fastqfiles [gzfiles $fastqdir/*.fastq $fastqdir/*.fq]
	if {[file exists $fastqdir/info_basecaller.txt]} {
		set basecaller [file_read $fastqdir/info_basecaller.txt]
	}
	set regions [list_remove [distrreg_regs $distrreg $refseq] unaligned]
	set bamfileindex [index_file $bamfile]
	bam_index_job $bamfile

	# distrmethod fast5
	# -------------------
	if {$distrmethod eq "fast5"} {
		return [meth_nanopolish_distrfast5 $fast5dir $fastqdir $bamfile $resultfile $refseq $skips $basecaller $callthreshold $threads $maxfastqdistr ${meth-compression} $nomatchingfast5error $opts]
	}

	# distrmethod regions
	# -------------------
	set tail [file tail $resultfile]
	if {[regexp ^meth- $tail]} {
		set smethfile [file dir $resultfile]/s$tail
	} else {
		set smethfile [file dir $resultfile]/smeth-$tail
	}
	set deps [list {*}$fastqfiles $fast5dir]
	job [job_relfile2name seqmeth_nanopolish_index- $bamfile] {*}$skips -cores $threads -deps $deps -targets {
		$fastqfile $fastqfile.index $fastqfile.index.fai $fastqfile.index.gzi $fastqfile.index.readdb
	} -skip {$smethfile} -skip {$resultfile} -vars {
		fastqfiles fastqfile fast5dir
	} -code {
		if {[llength $fastqfiles] == 1} {
			mklink [lindex $fastqfiles 0] $fastqfile
		} else {
			exec cg zcat {*}$fastqfiles | gzip --fast > $fastqfile
		}
		catch_exec nanopolish index -d $fast5dir $fastqfile
	}
	foreach region $regions {
		set target $resultfile.temp/$region.smeth.tsv.zst
		lappend seqmethfiles $target
		job [job_relfile2name seqmeth_nanopolish- $target] {*}$skips -cores $threads -mem 4G -deps {
			$fastqfile $fast5dir $refseq $bamfile ($bamfileindex)
		} -targets {
			$target
		} -skip {$smethfile} -skip {$resultfile} -vars {
			region fastqfile fast5dir refseq bamfile basecaller threads bamfileindex threads opts
		} -code {
			analysisinfo_write $dep $target basecaller_version $basecaller meth_caller nanopolish meth_caller_version [version nanopolish]
			set tempdir [tempdir]
			foreach file [glob -nocomplain $tempdir/*] {file delete $file}
			foreach {chr begin end} [split $region -] break
			set nregion ${chr}:${begin}-${end}
			set error [catch {
				exec nanopolish call-methylation -t $threads \
					-r $fastqfile -b $bamfile -g $refseq \
					-w $nregion \
					{*}$opts \
					 | cg zst --compressionlevel 1 > $target.temp.zst
			} msg opt]
			if {$error} {
				if {$::errorCode ne "NONE"} {
					dict unset opt -level
					set errorInfo "$msg\n    while executing\nnanopolish call-methylation -r $tempdir/[file tail $fastqfile] -b $bamfile -g $refseq | cg zst --compressionlevel 1 > $target.temp.zst"
					return -code $error -errorcode $::errorCode -errorinfo $::errorInfo $msg
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
		cg cat -c 0 {*}$deps | cg select -s - | cg ${meth-compression} > $target.temp
		file rename -- $target.temp $target
		if {${meth-compression} in "gz bgz"} {
			cg maketabix $target
		}
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
		file delete -force $target.temp
	}
	return [list $resultfile $smethfile]
}

proc cg_meth_nanopolish {args} {
	set args [job_init {*}$args]
	meth_nanopolish_job {*}$args
	job_wait
}

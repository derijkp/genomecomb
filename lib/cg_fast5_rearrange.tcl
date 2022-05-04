proc fast5_rearrange_job {args} {
	cg_options fast5_rearrange args {
	} {fastqdir fast5src fast5dst} 3 3 {
		rearrange (e.g. demultiplex) fast5s into matching fastq files (can be nested directories)
	}
	mkdir $fast5dst
	set base [join [lrange [file split $fast5dst] end-2 end] _]
	job_logfile $fast5dst/rearrange_fast5_$base [file_absolute $fast5dst] "cd [pwd] ; [list rearrange_fast5 $fastqdir $fast5src $fast5dst]"
	set nanopolishprog [exec which nanopolish]
	if {[file exists $nanopolishprog.plugins]} {
		if {[get ::env(HDF5_PLUGIN_PATH) ""] ne ""} {
			set ::env(HDF5_PLUGIN_PATH) $nanopolishprog.plugins:$::env(HDF5_PLUGIN_PATH)
		} else {
			set ::env(HDF5_PLUGIN_PATH) $nanopolishprog.plugins
		}
	}
	set fastqs [lsort -dict [split [exec find $fastqdir -name *.fastq.gz] \n]]
	set fast5s [lsort -dict [split [exec find $fast5src -name *.fast5] \n]]
	if {[string index $fastqdir end] ne "/"} {append fastqdir /}
	if {[string index $fast5src end] ne "/"} {append fast5src /}
	#
	# fast5 to slow5
	set workdir $fast5dst/tmp
	mkdir $workdir/src
	set srcs {}
	set num 0
	foreach fast5 $fast5s {
		if {![regsub ^$fast5src $fast5 {} efast5]} {
			error "fast5 $fast5 not in fast5dir $fast5dir?"
		}
		set dir [file dir $efast5]
		set tail [file tail $efast5]
		set target $workdir/src/$dir/[file root [gzroot $tail]].slow5
		lappend srcs $target
		job fast5_rearrange_2slow5_$target -deps {
			$fast5
		} -targets {
			$target
		} -vars {
			fast5
		} -code {
			mkdir [file dir $target]
			set nanopolishprog [exec which nanopolish]
			if {[file exists $nanopolishprog.plugins]} {
				if {[get ::env(HDF5_PLUGIN_PATH) ""] ne ""} {
					set ::env(HDF5_PLUGIN_PATH) $nanopolishprog.plugins:$::env(HDF5_PLUGIN_PATH)
				} else {
					set ::env(HDF5_PLUGIN_PATH) $nanopolishprog.plugins
				}
			}
			exec slow5tools f2s $fast5 -o $target.temp.slow5 >@ stdout 2>@ stderr
			file rename $target.temp.slow5 $target
			# exec slow5tools f2s $fast5 --to slow5 | zstd-mt -c -1 > $target 2>@ stderr
		}
	}
	# foreach file $srcs {if {![file exists $file]} {puts $file}}	
	# precalc ids
	mkdir $workdir/dst
	mkdir $workdir/ids
	set slow5targets {}
	set idfiles {}
	set num 0
	set tododeps {}
	set todotargets {}
	set batch 0
	set len [llength $fastqs]
	if {[string index $fastqdir end] ne "/"} {append fastqdir /}
	foreach fastq $fastqs {
		if {![regsub ^$fastqdir $fastq {} efastq]} {
			error "fastq $fastq not in fastqdir $fastqdir?"
		}
		set dir [file dir $efastq]
		set tail [file tail $efastq]
		set targetslow5 $workdir/dstslow/$dir/[file root [gzroot $tail]].slow5
		lappend slow5targets $targetslow5
		set root [file root [file tail $fastq]]
		set target $workdir/ids/[incr num]_$root.ids
		lappend idfiles $target
		lappend tododeps $fastq
		lappend todotargets $target
		incr batch
		if {$batch >= 100 || $num == $len} {
			job fast5_rearrange_2ids_$target -deps $tododeps -targets $todotargets -code {
				foreach fastq $deps target $targets {
					set idlines [exec cg fastq2tsv $fastq | cg select -f id -sh /dev/null]
					set list [regexp -all -inline {parent_read_id=([^ ]+)} $idlines]
					if {![llength $list]} {
						set list [regexp -all -inline {@([^ ]+)} $idlines]
					}
					list_unmerge $list 1 list
					file_write $target.temp [join $list \n]\n
					file rename -force $target.temp $target
				}
			}
			set batch 0
			set tododeps {}
			set todotargets {}
		}
	}
	# foreach file $idfiles {if {![file exists $file]} {puts $file}}	
	# create slow5targets
	job fast5_rearrange_main_$fast5dst \
	-deps [list {*}$srcs {*}$idfiles] \
	-targets [list {*}$slow5targets $workdir/dstslow/unidentified.slow5] \
	-vars {fastqdir fastqs srcs idfiles slow5targets workdir} \
	-code {
		# get slow5 header (based on first fast5)
		set slow5 [lindex $srcs 0]
		set slow5header {}
		set f [open $slow5]
		while {[gets $f line] != -1} {
			if {[string index $line 0] ni {# @}} break
			lappend slow5header $line
		}
		catch {close $f}
		set slow5header [join $slow5header \n]
	
		# go over fastqs and open (temp) matching slow5s
		mkdir $workdir/dstslow
		if {[string index $fastqdir end] ne "/"} {append fastqdir /}
		set len [llength $fastqs]
		set openfiles {}
		set slow5s {}
		set num 0
		unset -nocomplain a
		foreach fastq $fastqs idfile $idfiles {
			puts "[incr num]/$len $fastq"
			if {![regsub ^$fastqdir $fastq {} efastq]} {
				error "fastq $fastq not in fastqdir $fastqdir?"
			}
			set dir [file dir $efastq]
			set tail [file tail $efastq]
			set targetslow5 $workdir/dstslow/$dir/[file root [gzroot $tail]].slow5
			lappend slow5s $targetslow5
			mkdir [file dir $targetslow5]
			set o [open $targetslow5 w]
			lappend openfiles $o
			puts $o $slow5header
			set list [file_read $idfile]
			foreach id $list {
				set a($id) $o
			}
		}
		# go through all fast5s and distribute data to opened slow5s (matching fastqs)
		set ounid [open $workdir/dstslow/unidentified.slow5 w]
		puts $ounid $slow5header
		lappend slow5targets $workdir/dstslow/unidentified.slow5
		set len [llength $srcs]
		set num 0
		foreach slow5 $srcs {
			puts "[incr num]/$len $slow5"
			set f [open $slow5]
			while {[gets $f line] != -1} {
				if {[string index $line 0] ni {# @}} break
			}
			while 1 {
				set id [lindex [split $line \t] 0]
				if {![info exists a($id)]} {
					puts $ounid $line
				} else {
					puts $a($id) $line
				}
				if {[gets $f line] == -1} break
			}
			close $f
		}
		close $ounid
		foreach o $openfiles {
			close $o
		}
		unset a
	}
	# go over all generated slow5s and convert to fast5
	set len [llength $slow5targets]
	set num 0
	set resultfast5s {}
	mkdir $workdir/dst
	foreach slow5 [list {*}$slow5targets $workdir/dstslow/unidentified.slow5]  {
		puts "[incr num]/$len $slow5"
		if {![regsub ^$workdir/dstslow/ $slow5 {} eslow5]} {
			error "slow5 $slow5 not in workdir $workdir?"
		}
		set target $workdir/dst/[file root $eslow5].fast5
		lappend resultfast5s $target
		job fast5_rearrange_2fast5_$target -deps {
			$slow5
		} -targets {
			$target
		} -vars {
			slow5 eslow5
		} -code {
			set tempfile {}
			set nanopolishprog [exec which nanopolish]
			if {[file exists $nanopolishprog.plugins]} {
				if {[get ::env(HDF5_PLUGIN_PATH) ""] ne ""} {
					set ::env(HDF5_PLUGIN_PATH) $nanopolishprog.plugins:$::env(HDF5_PLUGIN_PATH)
				} else {
					set ::env(HDF5_PLUGIN_PATH) $nanopolishprog.plugins
				}
			}
			file delete -force $target.temp
			mkdir $target.temp
			exec slow5tools s2f -o $target.temp/result.fast5 $slow5 >@ stdout 2>@ stderr
			file rename -force $target.temp/result.fast5 $target
			file delete $target.temp
			# file delete $slow5
			if {$tempfile ne ""} {file delete $tempfile}
		}
	}
	job fast5_rearrange_mv -deps $resultfast5s -targets {
		$fast5dst/ready
	} -vars {
		slow5
	} -code {
		foreach dir [glob $fast5dst/tmp/dst/*] {
			exec cp -ral $dir $fast5dst/[file tail $dir]
		}
		if {[file exists $fast5dst/tmp/dst/unidentified.fast5]} {
			exec cp -ral $fast5dst/tmp/dst/unidentified.fast5 $fast5dst/unidentified.fast5
		}
		file_write $fast5dst/ready [timestamp]
		# file delete $fast5dst/temp2
	}
}

proc cg_fast5_rearrange {args} {
	set args [job_init {*}$args]
	fast5_rearrange_job {*}$args
	job_wait
}
proc refseq_minimap2_job {refseq {preset {}}} {
	if {$preset eq ""} {
		set preset map-ont
	}
	upvar job_logdir job_logdir
	job_logfile [file dir $refseq]/refseq_minimap2_[file tail $refseq]_$preset [file dir $refseq] "cg refseq_minimap2 $refseq $preset"
	set minimap2refseq $refseq.minimap2.$preset
	if {[file exists $minimap2refseq]} {return $minimap2refseq}
	set tail [file tail $refseq]
	if {[jobtargetexists [list $minimap2refseq] $refseq]} {	
		return $minimap2refseq
	}
	job [job_relfile2name minimap2_2refseq- $refseq] -deps {$refseq} -targets {$minimap2refseq} -vars {preset} -code {
		set size 0
		if {![file exists $dep.fai]} {
			catch_exec samtools faidx $dep
		}
		set f [open $dep.fai]
		while {[gets $f line] != -1} {
			incr size [lindex [split $line \t] 1]
		}
		close $f
		incr size 1000
		if {$size < 10000000000} {set size 10000000000}
		map_minimap2_presets $preset mpreset refpreset extraopts
		if {$preset eq "ontshort"} {
			set temp [catch_exec minimap2 -I $size -x map-ont -k 5 -w 1 -d $target.temp $dep]
		} elseif {$preset eq "short"} {
			set temp [catch_exec minimap2 -I $size -x sr -k 6 -w 2 -d $target.temp $dep]
		} else {
			set temp [catch_exec minimap2 -I $size -x $mpreset -d $target.temp $dep]
		}
		if {[regexp {loaded/built the index for 0 target sequence\(s\)} $temp]} {
			error "could not properly index $dep: contains no sequences"
		}
		file rename -- $target.temp $target
	}
	return $minimap2refseq
}

proc cg_refseq_minimap2 args {
	set args [job_init {*}$args]
	set return [refseq_minimap2_job {*}$args]
	job_wait
	return $return
}

proc refseq_minimap2 {refseq preset} {
	upvar job_logdir job_logdir
	if {$preset eq ""} {set preset map-ont}
	set refseq [file_absolute $refseq]
	set minimap2refseq $refseq.minimap2.$preset
	if {![jobfileexists $minimap2refseq]} {
		error "The minimap2 version for preset $preset of the refseq does not exist (should be at $minimap2refseq)
You can create it using:
cg refseq_minimap2 \'$refseq\' $preset"
	}
	return $minimap2refseq
}

proc map_mem_minimap2 {mem threads preset deps} {
	if {$mem eq ""} {
		set refseq [lindex $deps 0]
		if {[file exists $refseq.minimap2.$preset]} {
			# scale according to size index file
			set size [file size $refseq.minimap2.$preset]
			set mem [expr {round(2.5*$size)}]
			# but require minimum 6G
			if {$mem < 6442450944} {set mem 6442450944}
		} else {
			if {[regexp splice $preset]} {
				set mem 20G
			} else {
				set mem 10G
			}
		}
	}
	return $mem
}

# presets
# ont : Oxford Nanopore genomic reads (map-ont) (default)
# pb : PacBio genomic reads (map-pb)
# pacbio : PacBio genomic reads (map-pb)
# hifi : PacBio hifi genomic reads (map-hifi)
# asm20 : PacBio CCS genomic reads
# sr : short genomic paired-end reads
# splice : spliced long reads (strand unknown)
# splice : noisy Nanopore Direct RNA-seq
# splicehq : Final PacBio Iso-seq or traditional cDNA (splice:hq)
# asm5 : intra-species asm-to-asm alignment
# avaob : PacBio read overlap (ava-pb)
# avaont : Nanopore read overlap (ava-ont)
# any preset minimap2 version accepts can also be given
#
# custom presets
# ontshort : ONT optimized to find (very) short matches
# short : short read optimized to find (very) short matches
# splicesmall : minimap2 splice  preset with parameters set to detect small exons (can lead to errors in others)
# splicesens : minimap2 splice  preset with parameters set to be more sensitive
# splicesrsens : minimap2 splice:sr  preset with parameters set to be more sensitive

proc map_minimap2_presets {value mpresetVar refpresetVar extraoptsVar} {
	upvar $mpresetVar mpreset
	upvar $refpresetVar refpreset
	upvar $extraoptsVar extraopts
	set mpreset $value
	set refpreset $value
	if {$value eq "splicehq"} {
		set mpreset splice:hq
		set refpreset splice:hq
	} elseif {$value in "pb pacbio"} {
		set mpreset map-pb
		set refpreset map-pb
		set platform PACBIO
		lappend extraopts -L
	} elseif {$value in "hifi"} {
		set mpreset map-hifi
		set refpreset map-hifi
		set platform PACBIO
		lappend extraopts -L
	} elseif {$value in "ont"} {
		set mpreset map-ont
		set refpreset map-ont
		lappend extraopts -L
	} elseif {$value in "ontshort"} {
		# have to keep this and change just before using because it needs a different index
		set mpreset map-ont
		set refpreset ontshort
		lappend extraopts -n 1 -m 1 -k 5 -w 1 -s 20
	} elseif {$value in "short"} {
		# have to keep this and change just before using because it needs a different index
		set mpreset sr
		set refpreset short
		lappend extraopts -n 1 -m 5 -k 6 -w 2 -s 20 -r 50 --no-long-join --secondary=yes -N 50
	} elseif {$value in "avapb"} {
		set mpreset ava-pb
		set refpreset ava-pb
		set platform PACBIO
	} elseif {$value in "avaont"} {
		set mpreset ava-ont
		set refpreset ava-ont
		set platform PACBIO
	} elseif {$value eq "splicesmall"} {
		set mpreset splice
		set refpreset splice
		lappend extraopts -B3 -O3,6
	} elseif {$value eq "splicesens"} {
		set mpreset splice
		set refpreset splice
		lappend extraopts -N50 -p0.1 -A2 -B4 -O4,24 -E2,1
	} elseif {$value eq "splicesrsens"} {
		set mpreset splice:sr
		set refpreset splice:sr
		lappend extraopts -N50 -p0.1 -A2 -B4 -O4,24 -E2,1
	} else {
		set mpreset $value
		set refpreset $value
	}
	
}

proc cg_map_minimap2 {args} {
	if {[info exists ::cgextraopts(minimap2)]} {set extraopts $::cgextraopts(minimap2)} else {set extraopts {}}
	set paired 0
	set keepargs $args
	set preset {}
	set readgroupdata {}
	set threads 2
	set fixmate 1
	set aliformat bam
	set ali_keepcomments {}
	set platform ONT
	cg_options map_minimap2 args {
		-paired - -p {
			set paired $value
		}
		-x - -preset {
			set preset $value
			map_minimap2_presets $preset mpreset refpreset extraopts
		}
		-readgroupdata {
			set readgroupdata $value
		}
		-fixmate {
			set fixmate $value
		}
		-threads - -t {
			set threads $value
		}
		-extraopts {
			lappend extraopts {*}$value
		}
		-ali_keepcomments {
			set ali_keepcomments $value
		}
		-nohardclips {
			if {[true $value]} {
				lappend extraopts -Y
			}
		}
	} {result refseq sample fastqfile1} 4 ... {
		align reads in fastq files to a reference genome using minimap2
	}
	if {$ali_keepcomments eq ""} {
		if {[file extension $fastqfile1] in ".bam .ubam"} {set ali_keepcomments 1} else {set ali_keepcomments 0}
	}
	if {$ali_keepcomments} {
		lappend extraopts -y
	}
	if {$preset eq ""} {
		if {$paired} {
			set preset sr
			set mpreset sr
			set refpreset sr
			set platform illumina
		} else {
			set preset map-ont
			set mpreset map-ont
			set refpreset map-ont
			set platform ONT
		}
	}
	set files [list $fastqfile1 {*}$args]
	set result [file_absolute $result]
	set refseq [refseq $refseq]
	#
	set rg [sam_readgroup $readgroupdata $sample RG PL $platform]
	if {$rg ne ""} {
		lappend extraopts -R $rg
	}
	set minimap2refseq [refseq_minimap2 $refseq $preset]
	set outpipe [convert_pipe -.sam $result -endpipe 1 -refseq $refseq]
	analysisinfo_write $fastqfile1 $result sample [file tail $sample] aligner minimap2 aligner_version [version minimap2] aligner_preset $preset reference [file2refname $minimap2refseq] aligner_paired $paired
	if {!$paired} {
		putslog "making $result"
		if {[catch {
			exec minimap2 -a -x $mpreset -t $threads --MD \
				{*}$extraopts \
				$minimap2refseq {*}$files {*}$outpipe
		} msg]} {
			if {[regexp ERROR: $msg] || $::errorCode ne "NONE"} {
				puts stderr $msg
				error $msg
			}
		}
		# puts stderr $msg
		# puts stderr "previous is message, not error"
	} else {
		if {$fixmate} {
			set fixmate "| samtools fixmate -m -O sam - -"
		}
		if {[expr {[llength $files]%2}]} {
			error "minimap2 needs even number of files for paired analysis"
		}
		putslog "making $result"
		if {[catch {
			exec minimap2 -a -x $mpreset -t $threads --MD \
				{*}$extraopts \
				$minimap2refseq {*}$files {*}$fixmate {*}$outpipe
		} msg]} {
			if {[regexp ERROR: $msg] || $::errorCode ne "NONE"} {
				puts stderr $msg
				error $msg
			}
		}
		# puts stderr $msg
		# puts stderr "previous is message, not error"
	}
}

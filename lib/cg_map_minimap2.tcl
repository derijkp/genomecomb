proc refseq_minimap2_job {refseq preset} {
	upvar job_logdir job_logdir
	set minimap2refseq $refseq.minimap2.$preset
	if {[file exists $minimap2refseq]} {return $minimap2refseq}
	set tail [file tail $refseq]
	if {[jobtargetexists [list $minimap2refseq] $refseq]} return
	job [job_relfile2name minimap2_2refseq- $refseq] -deps {$refseq} -targets {$minimap2refseq} -vars {preset} -code {
		set temp [catch_exec minimap2 -x $preset -d $target.temp $dep]
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
	set refseq [file_absolute $refseq]
	set minimap2refseq $refseq.minimap2.$preset
	if {![file exists $minimap2refseq]} {
		error "The minimap2 version for preset $preset of the refseq does not exist (should be at $minimap2refseq)
You can create it using:
cg refseq_minimap2 \'$refseq\' $preset"
	}
	return $minimap2refseq
}

proc map_mem_minimap2 {mem threads} {
	if {$mem eq ""} {set mem 10G}
	return $mem
}

# presets
# map-pb : PacBio genomic reads
# map-ont : Oxford Nanopore genomic reads
# asm20 : PacBio CCS genomic reads
# sr : short genomic paired-end reads
# splice : spliced long reads (strand unknown)
# splice : noisy Nanopore Direct RNA-seq
# splice:hq : Final PacBio Iso-seq or traditional cDNA
# asm5 : intra-species asm-to-asm alignment
# va-pb : PacBio read overlap
# va-ont : Nanopore read overlap

proc cg_map_minimap2 {args} {
	set paired 1
	set keepargs $args
	set preset {}
	set readgroupdata {}
	set threads 2
	set mem 10G
	set fixmate 1
	set aliformat bam
	cg_options map_minimap2 args {
		-paired - -p {
			set paired $value
		}
		-x - -preset - -p {
			if {$value eq "splicehq"} {set value splice:hq}
			if {$value eq "splicesmall"} {
				set value splice
				lappend extraopts -B3 -O3,6
			}
			set preset $value
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
		-mem {
			set mem $value
		}
	} {result refseq sample fastqfile1} 4 5 {
		align reads in fastq files to a reference genome using minimap2
	}
	if {$preset eq ""} {
		if {$paired} {set preset sr} else {set preset map-ont}
	}
	set files [list $fastqfile1 {*}$args]
	set result [file_absolute $result]
	set refseq [refseq $refseq]
	#
	set readgroupdata [map_readgroupdata $readgroupdata $sample]
	set minimap2refseq [refseq_minimap2 $refseq $preset]
	set outpipe [convert_pipe -.sam $result -endpipe 1 -refseq $refseq]
	analysisinfo_write $fastqfile1 $result aligner minimap2 aligner_version [version minimap2] reference [file2refname $minimap2refseq] aligner_paired $paired
	if {!$paired} {
		putslog "making $result"
		set rg {}
		foreach {key value} $readgroupdata {
			lappend rg "$key:$value"
		}
		if {[catch {
			exec minimap2 -a -x $preset -t $threads --MD \
				-R @RG\\tID:$sample\\t[join $rg \\t] \
				$minimap2refseq {*}$files {*}$outpipe
		} msg]} {
			if {[regexp ERROR: $msg]} {
				puts stderr $msg
				error $msg
			}
		}
		puts stderr $msg
	} else {
		if {$fixmate} {
			set fixmate "| samtools fixmate -m -O sam - -"
		}
		if {[expr {[llength $files]%2}]} {
			error "minimap2 needs even number of files for paired analysis"
		}
		putslog "making $result"
		set rg {}
		foreach {key value} $readgroupdata {
			lappend rg "$key:$value"
		}
		if {[catch {
			exec minimap2 -a -x $preset -t $threads --MD \
				-R @RG\\tID:$sample\\t[join $rg \\t] \
				$minimap2refseq {*}$files {*}$fixmate {*}$outpipe
		} msg]} {
			if {[regexp ERROR: $msg]} {
				puts stderr $msg
				error $msg
			}
		}
		puts stderr $msg
	}
}

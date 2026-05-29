proc map_mem_bwa {mem threads preset deps} {
	upvar job_logdir job_logdir
# todo
	if {$mem eq ""} {
		set refseq [lindex $deps 0]
		set bwarefseq [refseq_bwa_job $refseq $preset]
		if {[file exists $bwarefseq.bwt]} {
			# scale according to size index file
			set size [file size $bwarefseq.bwt]
			set mem [expr {round(2*$size) + $threads*100000000}]
			if {$preset eq "short"} {set mem [expr {2*$mem}]}
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

proc refseq_bwa_job {refseq {preset {}}} {
	upvar job_logdir job_logdir
	set refseq [file_absolute $refseq]
	set bwarefseq $refseq.bwa/[file tail $refseq]
	set bwarefseqfa [file root $bwarefseq].fa
	if {[file exists $bwarefseqfa]} {return $bwarefseqfa}
	set tail [file tail $refseq]
	set targets [list $bwarefseq $bwarefseq.fai]
	foreach ext {amb ann bwt pac sa} {
		lappend targets $bwarefseq.$ext
	}
	if {[file extension $refseq] ne ".fa"} {
		set fatargets [list $bwarefseqfa $bwarefseqfa.fai]
		foreach ext {amb ann bwt pac sa} {
			lappend fatargets $bwarefseqfa.$ext
		}
	} else {
		set fatargets {}
	}
	if {[jobtargetexists $targets $refseq] && [jobtargetexists $fatargets $refseq]} {return $bwarefseqfa}
	# set dep $refseq ; set target $bwarefseq
	job [job_relfile2name bwa2refseq- $refseq] -deps {$refseq} -targets $targets -code {
		file delete -force $target.temp
		file mkdir $target.temp
		mklink $dep $target.temp/[file tail $dep]
		if {![file exists $dep.fai]} {
			catch_exec samtools faidx $dep
		}
		mklink $dep.fai $target.temp/[file tail $dep].fai
		exec bwa index $target.temp/[file tail $dep] 2>@ stderr
		set targetdir [file dir $target]
		foreach file [glob $target.temp/*] {
			file delete $targetdir/[file tail $file]
			file rename -- $file $targetdir/[file tail $file]
		}
		file delete -force $target $target.fai
		mklink $dep $target
		mklink $dep.fai $target.fai
		file delete $target.temp
	}
	if {[llength $fatargets]} {
		# with .fa as extension
		job [job_relfile2name bwa2refseq_fa- $refseq] -deps $targets -targets $fatargets -vars {
			bwarefseq bwarefseqfa refseq fatargets
		} -code {
			foreach file $deps nfile $fatargets {
				if {$file eq "$refseq.bwa"} continue
				mklink $file $nfile
			}
			mklink $bwarefseq $bwarefseqfa
		}
	}
	return $bwarefseqfa
}

proc cg_refseq_bwa args {
	set args [job_init {*}$args]
	set return [refseq_bwa_job {*}$args]
	job_wait
	return $return
}

proc refseq_bwa {refseq {preset {}}} {
	upvar job_logdir job_logdir
	set refseq [file_absolute $refseq]
	set bwarefseq $refseq.bwa/[file tail $refseq]
	set bwarefseqfa [file root $bwarefseq].fa
	if {![jobfileexists $bwarefseqfa]} {
		error "The bwa version of the refseq does not exist (should be at $bwarefseqfa)
You can create it using:
cg refseq_bwa \'$refseq\'"
	}
	return $bwarefseqfa
}

proc cg_map_bwa {args} {
	set paired 1
	set readgroupdata {}
	set threads 2
	set fixmate 1
	set extraopts {}
	set ali_keepcomments {}
	cg_options map_bwa args {
		-paired - -p {
			set paired $value
		}
		-x - -preset - -p {
			if {$value eq "short"} {
				lappend extraopts -k 12 -T 18 -a
			} elseif {$value ne ""} {
				error "unknown bwa preset $value"
			}
		}		
		-readgroupdata {
			set readgroupdata $value
		}
		-fixmate {
			set fixmate $value
		}
		-ali_keepcomments {
			# not used yet
			set ali_keepcomments $value
		}
		-nohardclips {
			if {[true $value]} {
				lappend extraopts -Y
			}
		}
		-threads - -t {
			set threads $value
		}
		-extraopts {
			lappend extraopts {*}$value
		}
	} {result refseq sample fastqfile1 fastqfile2} 4 5 {
		align reads in fastq files to a reference genome using bwa-mem
	}
	set refseq [refseq $refseq]
	dbdir [file dir $refseq]
	set bwarefseq [refseq_bwa $refseq]
	set outpipe [convert_pipe -.sam $result -endpipe 1 -refseq $refseq]
	putslog "making $result"
	analysisinfo_write $fastqfile1 $result sample [file tail $sample] aligner bwa aligner_version [version bwa] reference [file2refname $bwarefseq] aligner_paired $paired
	set rg [sam_readgroup $readgroupdata $sample]
	if {$rg ne ""} {
		lappend extraopts -R $rg
	}
	if {!$paired} {
		catch_exec bwa mem -t $threads \
			{*}$extraopts \
			$bwarefseq $fastqfile1 {*}$outpipe 2>@ stderr
	} else {
		if {$fixmate} {
			set fixmate "| samtools fixmate -m -O sam - -"
		}
		if {![info exists fastqfile2]} {
			error "bwa needs 2 files for paired analysis"
		}
		catch_exec bwa mem -t $threads -M \
			{*}$extraopts \
			$bwarefseq $fastqfile1 $fastqfile2 {*}$fixmate {*}$outpipe 2>@ stderr
	}
}

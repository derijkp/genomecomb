proc validate_var_medaka {refseq distrreg datatype} {
	# seperate because command is not the same as cmd
	if {[catch {exec which medaka_variant}]} {
		error "command \"medaka_variant\" not available, make sure medaka is installed, e.g. using \"cg install medaka\""
	}
}

proc var_medaka_tools {} {
	return {medaka}
}

proc var_medaka_distrreg {distrreg} {
	# per chromosome takes too much memory, split further
	# we do keep 0, as we suppose it would only be used for small organisms
	if {$distrreg in {chr chromosome 1 schr schromosome}} {
		return s10000000
	} else {
		return $distrreg
	}
}

proc version_medaka {} {
	set version ?
	catch {exec medaka_variant -h} c
	regexp {medaka ([^\n]+)} $c temp version
	return $version
}

proc var_medaka_job {args} {
	# putslog [list var_medaka_job {*}$args]
	global appdir
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg var_medaka {*}$args]
	set pre ""
	set split 1
	set deps {}
	set region {}
	set threads 2
	set cleanup 1
	set mincoverage 8
	set mingenoqual 25
	set resultfiles 0
	set rootname {}
	set skips {}
	set opts {}
	set index 1
	set resultfile {}
	set batchsize 40
	set mem 5G
	set time 3:00:00
	cg_options var_medaka args {
		-L - -deps {
			lappend deps [file_absolute $value]
		}
		-region {
			set region $value
		}
		-pre {
			set pre $value
		}
		-split {
			set split $value
		}
		-threads - -t {
			if {$value > 2} {
				putslog "-threads $value reduced to 2 because more do not help for medaka"
				set threads 2
			} else {
				set threads $value
			}
		}
		-cleanup {
			set cleanup $value
		}
		-mincoverage {
			set mincoverage $value
		}
		-resultfiles {
			set resultfiles $value
		}
		-rootname {
			set rootname $value
		}
		-datatype {
			# not actually used
		}
		-index {
			set index [true $value]
		}
		-opts {
			set opts $value
		}
		-skip {
			lappend skips -skip $value
		}
		-mem {
			set mem $value
		}
		-time {
			set time $value
		}
	} {bamfile refseq resultfile} 2 3
	foreach {key value} [specialopts -medaka] {
		switch $key {
			-b {set batchsize $value}
			default {lappend opts $key $value}
		}
	}
	lappend opts -b $batchsize
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	if {$resultfile eq ""} {
		set resultfile [file dir $bamfile]/${pre}var-medaka-[file_rootname $bamfile].tsv.zst
	} else {
		set resultfile [file_absolute $resultfile]
	}
	set resulttail [file tail $resultfile]
	set destdir [file dir $resultfile]
	if {$rootname eq ""} {
		set root [file_rootname $resultfile]
	} else {
		set root $rootname
	}
	# resultfiles
	set varfile $resultfile
	set vcffile [file root [gzroot $varfile]].vcf.gz
	set sregfile $destdir/${pre}sreg-$root.tsv.zst
	set resultlist [list $varfile $sregfile {} $vcffile]
	set medakatargets [list $varfile $vcffile]
	if {$resultfiles} {
		return $resultlist
	}
	lappend skips -skip $resultlist
	# logfile
	job_logfile $destdir/var_medaka_$resulttail $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools gatk picard java gnusort8 zst os]
	# start
	## Produce medaka SNP calls
	set dep $bamfile
	set bamindex $bamfile.[indexext $bamfile]
	set deps [list $bamfile $refseq $bamindex {*}$deps]
#putsvars deps medakatargets vcffile region refseq root varfile split opts region outbam mincoverage index
#error stop
	job [job_relfile2name medaka- $varfile] {*}$skips -mem $mem -time $time -cores $threads \
	-deps $deps \
	-targets $medakatargets -vars {
		vcffile region refseq root varfile split opts region outbam mincoverage index threads
	} -code {
		analysisinfo_write $dep $varfile analysis $root sample $root varcaller medaka varcaller_version [version medaka] varcaller_cg_version [version genomecomb] varcaller_region $region
		set regions [samregions $region $refseq]
		# set tempdir [gzroot $vcffile].temp
		set tempdir [filetemp $vcffile]
		catch {file delete -force $tempdir}
		mkdir $tempdir
		set tempvcf $tempdir/result.vcf.gz
		if {![llength $regions]} {set regions {{}}}
		set todo {}
		foreach region $regions {
			mkdir $tempdir/$region
			set runopts $opts
			if {$region ne ""} {
				lappend runopts -r $region
			}
			set tempfile [tempfile].vcf
			if {[catch {
				catch_exec medaka_variant {*}$runopts \
					-l -p \
					-t $threads \
					-i $dep \
					-f $refseq \
					-o $tempdir/$region
			} msg]} {
				error $msg
			}
			lappend todo $tempdir/$region/round_1.vcf
		}
		if {[llength $regions] > 1} {
			exec cg vcfcat {*}$todo | bgzip > $tempvcf
			file rename $tempvcf $vcffile
		} else {
			compress [lindex $todo 0] $vcffile 1
		}
		cg vcf2tsv -split $split $vcffile $varfile
		cg_zindex $varfile
		file delete -force $tempdir
	}
	# make sreg
	job [job_relfile2name medaka-sreg- $varfile] {*}$skips \
	-deps $deps \
	-targets {
		$sregfile
	} -vars {
		bamfile varfile refseq sregfile region mincoverage refseq
	} -code {
		set compress [compresspipe $sregfile 1]
		set temptarget [filetemp $sregfile]
		set opts {}
		if {$region ne ""} {
			lappend opts -region $region
		}
		exec cg regextract -stack 1 {*}$opts -refseq $refseq -min $mincoverage $bamfile {*}$compress > $temptarget
		file rename -force -- $temptarget $sregfile
		cg_zindex $sregfile
	}
	return $resultlist
}

proc cg_var_medaka {args} {
	set args [job_init {*}$args]
	var_medaka_job {*}$args
	job_wait
}

proc var_longshot_tools {} {
	return {longshot}
}

proc version_longshot {} {
	set version ?
	catch {exec longshot -V} version
	set version [string trim $version]
	regsub {Longshot *} $version {} version
	return $version
}

proc var_longshot_job {args} {
	# putslog [list var_longshot_job {*}$args]
	global appdir
	upvar job_logdir job_logdir
	set cmdline "[list cd [pwd]] \; [list cg var_longshot {*}$args]"
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
	set tech ont
	set opts {}
	cg_options var_longshot args {
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
			putslog "-threads $value ignored because longshot does not support threads for now"
			set threads $value
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
		-tech {
			if {$value ni "ont pacbio"} {error "-tech $value not supported, must be one of: ont pacbio"}
			set tech $value
		}
		-opts {
			set opts $value
		}
		-skip {
			lappend skips -skip $value
		}
	} {bamfile refseq}
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	set destdir [file dir $bamfile]
	set bamtail [file tail $bamfile]
	if {$rootname eq ""} {
		set root longshot-[file_rootname $bamtail]
	} else {
		set root $rootname
	}
	# resultfiles
	set varfile ${pre}var-$root.tsv.zst
	set sregfile ${pre}sreg-$root.tsv.zst
	set vcffile [file root [gzroot $varfile]].vcf.gz
	set resultlist [list $destdir/$varfile $destdir/$sregfile {} $destdir/$vcffile]
	if {$resultfiles} {
		return $resultlist
	}
	# logfile
	job_logfile $destdir/var_longshot_$bamtail $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools gatk picard java gnusort8 zst os]
	# start
	## Produce longshot SNP calls
	set keeppwd [pwd]
	cd $destdir
	set dep $bamtail
	set bamtailindex $bamtail.[indexext $bamtail]
	set deps [list $bamtail $refseq $bamtailindex {*}$deps]
	job longshot-[file tail $varfile] {*}$skips -mem 15G -deps $deps -targets {
		$varfile $vcffile
	} -vars {
		vcffile region refseq root varfile split tech opts region
	} -code {
		if {$tech eq "ont"} {
			lappend opts --strand_bias_pvalue_cutoff 0.01
		}
		if {$region ne ""} {
			lappend opts --region [samregion $region $refseq]
		}
		analysisinfo_write $dep $varfile sample $root varcaller longshot varcaller_version [version longshot] varcaller_cg_version [version genomecomb] varcaller_region $region
		catch_exec longshot {*}$opts \
			--bam $dep \
			--ref $refseq \
			--out [gzroot $vcffile].temp
		exec gzip [gzroot $vcffile].temp
		file rename -force -- [gzroot $vcffile].temp.gz $vcffile
		cg vcf2tsv -split $split $vcffile $varfile
		cg_zindex $varfile
	}
	# make sreg
	job longshot-sreg-[file tail $varfile] {*}$skips -deps {
		$bamfile
	} -targets {
		$sregfile
	} -vars {
		bamfile varfile refseq sregfile region mincoverage
	} -code {
		set compress [compresspipe $sregfile 1]
		set temptarget [filetemp $sregfile]
		set opts {}
		if {$region ne ""} {
			lappend opts -region $region
		}
		exec cg regextract -stack 1 {*}$opts -min $mincoverage $bamfile {*}$compress > $temptarget
		file rename -force -- $temptarget $sregfile
		cg_zindex $sregfile
	}
	cd $keeppwd
	return $resultlist
}

proc cg_var_longshot {args} {
	set args [job_init {*}$args]
	var_longshot_job {*}$args
	job_wait
}
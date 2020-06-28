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

proc longshot_empty_vcf {vcffile} {
	set o [open $vcffile w]
	puts $o [deindent {
		##fileformat=VCFv4.2
		##source=Longshot v0.4.0
		##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth of reads passing MAPQ filter">
		##INFO=<ID=AC,Number=R,Type=Integer,Description="Number of Observations of Each Allele">
		##INFO=<ID=AM,Number=1,Type=Integer,Description="Number of Ambiguous Allele Observations">
		##INFO=<ID=MC,Number=1,Type=Integer,Description="Minimum Error Correction (MEC) for this single variant">
		##INFO=<ID=MF,Number=1,Type=Float,Description="Minimum Error Correction (MEC) Fraction for this variant.">
		##INFO=<ID=MB,Number=1,Type=Float,Description="Minimum Error Correction (MEC) Fraction for this variant's haplotype block.">
		##INFO=<ID=AQ,Number=1,Type=Float,Description="Mean Allele Quality value (PHRED-scaled).">
		##INFO=<ID=GM,Number=1,Type=Integer,Description="Phased genotype matches unphased genotype (boolean).">
		##INFO=<ID=DA,Number=1,Type=Integer,Description="Total Depth of reads at any MAPQ (but passing samtools filter 0xF00).">
		##INFO=<ID=MQ10,Number=1,Type=Float,Description="Fraction of reads (passing 0xF00) with MAPQ>=10.">
		##INFO=<ID=MQ20,Number=1,Type=Float,Description="Fraction of reads (passing 0xF00) with MAPQ>=20.">
		##INFO=<ID=MQ30,Number=1,Type=Float,Description="Fraction of reads (passing 0xF00) with MAPQ>=30.">
		##INFO=<ID=MQ40,Number=1,Type=Float,Description="Fraction of reads (passing 0xF00) with MAPQ>=40.">
		##INFO=<ID=MQ50,Number=1,Type=Float,Description="Fraction of reads (passing 0xF00) with MAPQ>=50.">
		##INFO=<ID=PH,Number=G,Type=Integer,Description="PHRED-scaled Probabilities of Phased Genotypes">
		##INFO=<ID=SC,Number=1,Type=String,Description="Reference Sequence in 21-bp window around variant.">
		##FILTER=<ID=dn,Description="In a dense cluster of variants">
		##FILTER=<ID=dp,Description="Exceeds maximum depth">
		##FILTER=<ID=sb,Description="Allelic strand bias">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
		##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase Set">
		##FORMAT=<ID=UG,Number=1,Type=String,Description="Unphased Genotype (pre-haplotype-assembly)">
		##FORMAT=<ID=UQ,Number=1,Type=Float,Description="Unphased Genotype Quality (pre-haplotype-assembly)">
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
	}]
	close $o
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
	set maxcov {}
	set mingenoqual 25
	set resultfiles 0
	set rootname {}
	set skips {}
	set tech ont
	set opts {}
	set hap_bam 0
	set index 1
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
		-maxcoverage {
			set maxcov $value
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
		-hap_bam {
			set hap_bam [true $value]
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
	set longshottargets [list $destdir/$varfile $destdir/$vcffile]
	if {$hap_bam} {
		set hap_base $destdir/${pre}map-$root
		lappend resultlist $hap_base.hap1.bam $hap_base.hap2.bam $hap_base.unassigned.bam
		lappend longshottargets $hap_base.hap1.bam $hap_base.hap2.bam $hap_base.unassigned.bam
	}
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
	job longshot-[file_part $varfile end-1 end] {*}$skips -mem 8G -deps $deps -targets $longshottargets -vars {
		vcffile region refseq root varfile split tech opts region hap_bam hap_base maxcov mincoverage index
	} -code {
		if {$tech eq "ont"} {
			lappend opts --strand_bias_pvalue_cutoff 0.01
		}
		if {$maxcov ne ""} {
			lappend opts --max_cov $maxcov
		}
		analysisinfo_write $dep $varfile sample $root varcaller longshot varcaller_version [version longshot] varcaller_cg_version [version genomecomb] varcaller_region $region
		set regions [samregions $region $refseq]
		if {[llength $regions] > 1} {
			set todo {}
			foreach region $regions {
				set runopts $opts
				set tempfile [tempfile].vcf
				if {$hap_bam} {
					lappend runopts --hap_bam_prefix $tempfile
				}
				if {[catch {
					catch_exec longshot {*}$runopts --region $region -F \
						--min_cov $mincoverage \
						--bam $dep \
						--ref $refseq \
						--out $tempfile
				} msg]} {
					if {[regexp "^error: Chromosome name for region is not in BAM file" $msg]} {
						putslog "longshot warning: Chromosome name ($region) for region is not in BAM file, writing empty"
						longshot_empty_vcf $tempfile
						exec samtools view -H -b $dep > $tempfile.hap1.bam
						exec samtools view -H -b $dep > $tempfile.hap2.bam
						exec samtools view -H -b $dep > $tempfile.unassigned.bam
					} else {
						error $msg
					}
				}
				if {![file exists $tempfile] && [regexp {0 potential variants identified.} $msg]} {
					longshot_empty_vcf $tempfile
					exec samtools view -H -b $dep > $tempfile.hap1.bam
					exec samtools view -H -b $dep > $tempfile.hap2.bam
					exec samtools view -H -b $dep > $tempfile.unassigned.bam
				}
				lappend todo $tempfile
			}
			exec cg vcfcat {*}$todo > [gzroot $vcffile].temp
			if {$hap_bam} {
				foreach part {hap1 hap2 unassigned} {
					set files {}
					foreach file $todo {
						lappend files $file.$part.bam
					}
					exec cg sam_catmerge -sort nosort $hap_base.$part.bam {*}$files
				}
			}
		} else {
			if {$hap_bam} {
				lappend opts --hap_bam_prefix $hap_base
			}
			if {[llength $regions]} {lappend opts --region [lindex $regions 0]}
			set tempfile [gzroot $vcffile].temp
			if {[catch {
				catch_exec longshot {*}$opts -F \
					--bam $dep \
					--ref $refseq \
					--out $tempfile
			} msg]} {
				if {[regexp "^error: Chromosome name for region is not in BAM file" $msg]} {
					putslog "longshot warning: Chromosome name ($region) for region is not in BAM file, writing empty"
					longshot_empty_vcf $tempfile
				} else {
					error $msg
				}
			}
			if {![file exists $tempfile] && [regexp {0 potential variants identified.} $msg]} {
				longshot_empty_vcf $tempfile
			}
			if {$hap_bam} {
				foreach file [list $hap_base.hap1.bam $hap_base.hap2.bam $hap_base.unassigned.bam] {
					if {![file exists $file]} {
						exec samtools view -H -b $dep > $file
					}
				}
				if {$index} {
					foreach file [glob [file root $dep]*.bam] {
						exec samtools index $file
					}
				}
			}
		}
		exec gzip [gzroot $vcffile].temp
		file rename -force -- [gzroot $vcffile].temp.gz $vcffile
		if {[file size $vcffile] == 0} {
			file_write $varfile ""
		} else {
			cg vcf2tsv -split $split $vcffile $varfile
			cg_zindex $varfile
		}
	}
	# make sreg
	job longshot-sreg-[file_part $varfile end-1 end] {*}$skips -deps {
		$bamfile
	} -targets {
		$sregfile
	} -vars {
		bamfile varfile refseq sregfile region mincoverage refseq maxcov
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
	cd $keeppwd
	return $resultlist
}

proc cg_var_longshot {args} {
	set args [job_init {*}$args]
	var_longshot_job {*}$args
	job_wait
}

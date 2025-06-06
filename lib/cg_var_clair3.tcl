proc var_clair3_tools {} {
	return {clair3}
}

proc validate_var_clair3 {refseq distrreg datatype} {
	# seperate because command is not clair3
	if {[catch {exec which run_clair3.sh}]} {
		error "command \"run_clair3.sh\" not available, try installing clair3 using e.g. \"cg install clair3\""
	}
}

proc version_clair3 {} {
	set version ?
	catch {exec run_clair3.sh -v} version
	regexp {Clair3 ([^\n]+)} $version temp version
	set version [lindex $version end]
}

proc clair3_dir {} {
	set exe [file_resolve [exec which run_clair3.sh]]
	set clair3dir [file dir $exe]
}

proc clair3_replacebam {finalbam oribam} {
	file lstat $oribam a
	set time $a(mtime)
	catch {
		file lstat $oribam.[indexext $oribam] a
		set btime $a(mtime)
	}
        file rename -force $oribam $oribam.old
	catch {file rename -force $oribam.[indexext $oribam] $oribam.[indexext $oribam].old}
        mklink $finalbam $oribam
        exec touch -h -d [clock format $time] $oribam
	if {[info exists btime]} {
	        mklink $finalbam.[indexext $finalbam] $oribam.[indexext $oribam]
	        exec touch -h -d [clock format $btime] $oribam.[indexext $oribam]
	}
	file delete $oribam.old
	file delete $oribam.[indexext $oribam].old
}

proc clair3_empty_vcf {vcffile} {
	set o [wgzopen $vcffile]
	puts $o [deindent {
		##fileformat=VCFv4.2
		##source=clair3 v0.4.0
		##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth of reads passing MAPQ filter">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
		##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
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
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
		##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase Set">
		##FORMAT=<ID=UG,Number=1,Type=String,Description="Unphased Genotype (pre-haplotype-assembly)">
		##FORMAT=<ID=UQ,Number=1,Type=Float,Description="Unphased Genotype Quality (pre-haplotype-assembly)">
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
	}]
	gzclose $o
}

proc var_clair3_find_model {sampledir} {
	set clairscript [follow_links [exec which run_clair3.sh]]
	set clairdir [file dir $clairscript]
	set bcmodel {}
	set fastq [gzfile $sampledir/ubam/*.bam]
	if {[file exists $fastq]} {
		set f [open "| samtools view $fastq"]
		set line [gets $f]
		catch {close $f}
		regexp {RG:Z:([^ \t]+)} $line temp bcmodel
		regsub {^[^_]+_} $bcmodel {} bcmodel
	} else {
		set fastq [gzfile $sampledir/fastq/*.f*q]
		if {[file exists $fastq]} {
			set f [gzopen $fastq]
			set line [gets $f]
			gzclose $f
			regexp {basecall_model_version_id=([^ ]+)} $line temp bcmodel
		}
	}
	if {$bcmodel ne ""} {
		if {[file exists $clairdir/models/$bcmodel]} {
			set model $clairdir/models/$bcmodel
			return $model
		}
		regsub {^[dr]na_} $bcmodel {} bcmodel
		set bcmodel2 {}
		if {[regsub {_barcode[0-9]+$} $bcmodel {} bcmodel2]} {
			set bcmodel2 [string_change $bcmodel2 {. {} @ _}]
			set model $clairdir/models/$bcmodel2
			if {[file exists $model]} {
				return $model
			}
			regsub {_[^_]+$} $bcmodel2 {} bcmodel2
			set model $model
			if {[file exists $clairdir/models/$bcmodel2]} {
				return $model
			}
		}
		if {[regsub {(@v[0-9.]+)[^@]*$} $bcmodel {\1} bcmodel2]} {
			set bcmodel2 [string_change $bcmodel2 {. {} @ _}]
			set model $clairdir/models/$bcmodel2
			if {[file exists $model]} {
				return $model
			}
			set model [gzfile $clairdir/models/$bcmodel2*]
			if {[file exists $model]} {
				return $model
			}
		}
		set bcmodel2 [string_change $bcmodel {. {}}]
		if {[regexp {^([^_]+)_(.+)_([^@]+)@([^_]+)} $bcmodel2 temp pore machine mode version]} {
			set model [gzfile $clairdir/models/${pore}_*_${mode}_$version]
			if {[file exists $model]} {
				return $model
			}
			set model [gzfile $clairdir/models/${pore}_*_${mode}_*]
			if {[file exists $model]} {
				return $model
			}
		}
	}
	set runinfofile [gzfile $sampledir/runinfo*.tsv]
	if {[file exists $runinfofile]} {
		set guppyversion [lindex [cg select -f "basecaller_version" $runinfofile] end]
	} else {
		if {![regexp {_v([0-9]\.[0-9]+\.[0-9]+)} [file tail $sampledir] temp guppyversion]} {
			set guppyversion ?
		}
	}
	if {[string index $guppyversion 0] in "5 6"} {
		set usemodel r941_prom_sup_g5014
	} elseif {[string index $guppyversion 0] in "3 4"} {
		set usemodel r941_prom_hac_g360+g422
	} elseif {[string index $guppyversion 0] in "2"} {
		set usemodel r941_prom_hac_g238
	} else {
		puts stderr "warning: -model for clair3 not given, could not be found in the fastq/ubam, and no runinfo file; using default r941_prom_sup_g5014"
		set usemodel r941_prom_sup_g5014
	}
	return $usemodel
}

proc var_clair3_job {args} {
	# putslog [list var_clair3_job {*}$args]
	global appdir
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg var_clair3 {*}$args]
	set preset {}
	set platform {}
	set model {}
	set pre ""
	set split 1
	set deps {}
	set region {}
	set phasing 1
	set threads 2
	set mincoverage 8
	set mingenoqual 6
	set resultfiles 0
	set rootname {}
	set skips {}
	set opts {}
	set index 1
	set resultfile {}
	set cleanup 1
	set mem {}
	set time {}
	cg_options var_clair3 args {
		-preset {
			set preset $value
		}
		-platform {
			if {$value ni "ont hifi ilmn"} {
				error "unsupported platform $value, must be one of: ont hifi ilmn"
			}
			set platform $value
		}
		-model {
			set model $value
		}
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
			set threads $value
		}
		-mincoverage {
			set mincoverage $value
		}
		-mingenoqual {
			set mingenoqual $value
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
		-phasing {
			set phasing [true $phasing]
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
		-cleanup {
			# not actually used here, but must be an option
			set cleanup $value
		}
		-mem {
			set mem $value
		}
		-time {
			set time $value
		}
	} {bamfile refseq resultfile} 2 3
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	if {$time eq ""} {set time ${threads}:00:00}
	# if {$mem eq ""} {set mem [expr {2*$threads}]G}
	if {$mem eq ""} {set mem 24G}
	# if {$phasing} {lappend opts	--enable_phasing --use_longphase_for_final_output_phasing}
	if {$preset ne ""} {
		switch $preset {
			ont {
				set platform ont
				set tech ont
			}
			hifi - pacbio {
				set platform hifi
				set tech pacbio
				set model hifi
			}
			ilmn {
				set platform ilmn
				set tech {}
				set model ilmn
			}
			default {
				error "preset $preset not supported by var_clair3"
			}
		}
	}
	if {$platform eq ""} {
		puts stderr "warning: -platform for clair3 not given; using default ont"
		set platform ont
	}
	if {$resultfile eq ""} {
		if {$rootname eq ""} {
			set resultfile [file dir $bamfile]/${pre}var-clair3-[file_rootname $bamfile].tsv.zst
		} else {
			set resultfile [file dir $bamfile]/${pre}var-clair3-$rootname.tsv.zst
		}
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
	set varallfile $destdir/${pre}varall-$root.gvcf.gz
	set sregfile $destdir/${pre}sreg-$root.tsv.zst
	set resultlist [list $varfile $sregfile $varallfile $vcffile]
	if {$resultfiles} {
		return $resultlist
	}
	lappend skips -skip $resultlist
	# logfile
	job_logfile $destdir/var_clair3_$resulttail $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools gatk picard java gnusort8 zst os]
	# start
	## Produce clair3 SNP calls
	set dep $bamfile
	set bamindex $bamfile.[indexext $bamfile]
	set deps [list $bamfile $refseq $bamindex {*}$deps]
#putsvars deps vcffile region refseq root varfile split tech opts region index
#error stop
	job clair3-[file_rootname $varfile] {*}$skips -mem $mem -time $time -cores $threads \
	-deps $deps -targets {
		$varfile $vcffile $varallfile
	} -vars {
		vcffile varfile varallfile root refseq
		region opts index threads
		mincoverage mingenoqual split platform model phasing
	} -code {
		if {$model eq ""} {
			set sampledir [file dir $dep]
			set usemodel [var_clair3_find_model $sampledir]
		} else {
			set usemodel $model
		}
		if {[file exists $usemodel]} {
			set usemodel [file_absolute $usemodel]
		} else {
			set clairscript [follow_links [exec which run_clair3.sh]]
			set clairdir [file dir $clairscript]
			if {![file exists $clairdir/models/$usemodel]} {
				error "model does not exists: $clairdir/models/$usemodel"
			}
			set usemodel $clairdir/models/$usemodel
		}
		analysisinfo_write $dep $varfile \
			analysis $root sample $root \
			varcaller clair3 varcaller_version [version clair3] \
			varcaller_cg_version [version genomecomb] varcaller_region $region \
			varcaller_platform $platform varcaller_model $usemodel \
			varcaller_mincoverage $mincoverage varcaller_mingenoqual $mingenoqual
		set regions [samregions $region $refseq]
		set tempvcfdir [gzroot $vcffile].temp
		mkdir $tempvcfdir
		if {$region eq "unmapped"} {
			clair3_empty_vcf $tempvcfdir/temp.vcf
			# add unmapped reads
			set tempoutbam [tempfile].unmapped.bam
			catch_exec samtools view -h -b -f 4 $dep > $tempoutbam
			catch_exec samtools index $tempoutbam
			file rename -force -- $tempoutbam $outbam
			file rename -force -- $tempoutbam.[indexext $tempoutbam] $outbam.[indexext $outbam]
		} else {
			if {[llength $regions]} {
				set tempbed [tempfile].bed
				distrreg_reg2bed $tempbed $regions $refseq
				lappend opts --bed_fn=$tempbed
			}
			putslog "Running clair3"
			if {![file exists $usemodel]} {
				error "model not found: $usemodel"
			}
			set result [catch_exec run_clair3.sh {*}$opts \
				--include_all_ctgs \
				--threads $threads \
				--platform=$platform \
				--model_path=$usemodel \
				--bam_fn=$dep \
				--ref_fn=$refseq \
				--output=$tempvcfdir \
				--gvcf \
			]
			if {[regexp {[Ee]rror} $result]} {
				error $result
			} else {
				puts stderr $result
			}
		}
		set gvcf [gzfile $tempvcfdir/merge_output.gvcf]
		if {[file exists $gvcf]} {
			putslog "Sorting gvcf ($gvcf.s.gz)"
			cg sortvcf -threads $threads $gvcf $gvcf.s.gz
			file rename -force -- $gvcf.s.gz $varallfile
		} else {
			empty_gvcf $varallfile [list $root]
		}
		if {[file exists $tempvcfdir/merge_output.vcf.gz]} {
			if {$phasing} {
				# check if empty -> longphase gives incorrect output on empty vcf file
				set ft [gzopen $tempvcfdir/merge_output.vcf.gz]
				set header [tsv_open $ft]
				if {[gets $ft line] == -1} {set empty 1} else {set empty 0}
				gzclose $ft
				if {$empty} {
					clair3_empty_vcf $vcffile
				} else {
					if {$platform eq "ont"} {
						set seqplatformopt --ont
					} else {
						set seqplatformopt --pb
					}
					putslog "Running longphase ($vcffile)"
					set result [catch_exec [clair3_dir]/bin/longphase phase \
						-s $tempvcfdir/merge_output.vcf.gz \
						-b $dep \
						-r $refseq \
						-t $threads -o $tempvcfdir/phased-merge_output \
						$seqplatformopt]
					putslog $result
					cg sortvcf -threads $threads $tempvcfdir/phased-merge_output.vcf $vcffile.s.gz
					file rename -force -- $vcffile.s.gz $vcffile
				}
			} else {
				cg sortvcf -threads $threads $tempvcfdir/merge_output.vcf.gz $vcffile.s.gz
				file rename -force -- $vcffile.s.gz $vcffile
			}
		} else {
			clair3_empty_vcf $vcffile
		}
		putslog "Making varfile ($varfile)"
		set tempfile [filetemp $varfile]
		set fields {chromosome begin end type ref alt quality alleleSeq1 alleleSeq2}
		lappend fields [subst {sequenced=if(\$genoqual < $mingenoqual || \$coverage < $mincoverage,"u","v")}]
		lappend fields [subst {zyg=if(\$genoqual < $mingenoqual || \$coverage < $mincoverage,"u",\$zyg)}]
		lappend fields *
		exec cg vcf2tsv -split $split -meta [list refseq [file tail $refseq]] -removefields {
			name filter
		} $vcffile | cg select -f $fields | cg zst > $tempfile
		file rename -force -- $tempfile $varfile
		cg_zindex $varfile
		catch {file delete -force $tempvcfdir}
	}
	# make sreg
	job clair3-sreg-[file_rootname $varfile] {*}$skips \
	-deps {
		$varallfile
	} -targets {
		$sregfile
	} -vars {
		bamfile varfile refseq sregfile region refseq
		mincoverage mingenoqual
	} -code {
		set temp [filetemp $target]
		if {![file size $dep]} {
			file_write $temp ""
		} else {
			exec cg vcf2tsv $dep \
				| cg select -q [subst {
					\$genoqual >= $mingenoqual && \$coverage >= $mincoverage && \$type ne "ins"
				}] -f {chromosome begin end} -s - \
				| cg regjoin {*}[compresspipe $target] > $temp
		}
		file rename -force -- $temp $target
		cg_zindex $target
	}
	return $resultlist
}

proc cg_var_clair3 {args} {
	set args [job_init {*}$args]
	var_clair3_job {*}$args
	job_wait
}

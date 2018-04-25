proc sreg_gatkh_job {job varallfile resultfile {mincoverage 5} {mingenoqual 12}} {
	upvar job_logdir job_logdir
	job $job -deps {$varallfile} -targets {$resultfile} -vars {mincoverage mingenoqual} -code {
		set temp [filetemp $target]
		file_write $temp "# regions selected from [gzroot $dep]: \$genoqual >= $mingenoqual && \$coverage >= $mincoverage\n"
		exec cg vcf2tsv $dep \
			| cg select -q [subst {
				\$genoqual >= $mingenoqual && \$coverage >= $mincoverage && \$type ne "ins"
			}] -f {chromosome begin end} \
			| cg regjoin {*}[compresspipe $target] > $temp
		file rename $temp $target
		if {[file extension $target] eq ".lz4"} {
			exec lz4index $target
		}
	}
}

proc var_gatkh_job {args} {
	global appdir
	upvar job_logdir job_logdir
	set pre ""
	set opts {}
	set split 1
	set deps {}
	set regionfile {}
	set threads 2
	set cleanup 1
	set regmincoverage 3
	set ERC BP_RESOLUTION
	set distrchr 0
	set mincoverage 5
	set mingenoqual 12
	cg_options var_gatkh args {
		-L - -deps {
			lappend deps [file_absolute $value]
		}
		-regionfile {
			set regionfile $value
		}
		-regmincoverage {
			set regmincoverage $value
		}
		-pre {
			set pre $value
		}
		-split {
			set split $value
		}
		-threads - -t {
			puts "-threads $value ignored because gatkh does not support threads for now"
			set threads $value
		}
		-cleanup {
			set cleanup $value
		}
		-ERC - -emitRefConfidence {
			if {$value ni {BP_RESOLUTION GVCF NONE}} {error "option $value not supported for -ERC, must be one of: BP_RESOLUTION, GVCF, NONE"}
			set ERC $value
		}
		-distrchr {
			set distrchr $value
		}
		-mincoverage {
			set mincoverage $value
		}
		-mingenoqual {
			set mingenoqual $value
		}
		default {
			lappend opts $key $value
		}
	} {bamfile refseq}
	set bamfile [file_absolute $bamfile]
	set refseq [file_absolute $refseq]
	set destdir [file dir $bamfile]
	if {$regionfile ne ""} {
		set regionfile [file_absolute $regionfile]
	} else {
		set regionfile [bam2reg_job -mincoverage $regmincoverage $bamfile]
	}
	lappend deps $regionfile
	# logfile
	set cmdline [list cg var_gatkh]
	foreach option {
		split deps bed pre regionfile regmincoverage
	} {
		if {[info exists $option]} {
			lappend cmdline -$option [get $option]
		}
	}
	lappend cmdline {*}$opts $bamfile $refseq
	job_logfile $destdir/var_gatkh_[file tail $bamfile] $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools gatk picard java gnusort8 lz4 os]
	set file [file tail $bamfile]
	set root [file_rootname $file]
	if {$root eq ""} {set root [file root $file]}
	# start
	## Produce gatkh SNP calls
	set keeppwd [pwd]
	cd $destdir
	set gatkrefseq [gatk_refseq_job $refseq]
	set dep $file
	set resultgvcf ${pre}varall-gatkh-$root.gvcf.gz
	set resultname ${pre}varall-gatkh-$root.gvcf
	set deps [list $file $gatkrefseq $file.bai {*}$deps]
	if {!$distrchr} {
		job $resultname -mem 5G -cores $threads -deps $deps -targets {
			$resultgvcf $resultgvcf.tbi ${pre}varall-gatkh-$root.gvcf.analysisinfo
		} -vars {
			opts regionfile gatkrefseq refseq threads root ERC resultgvcf
		} -code {
			analysisinfo_write $dep $resultgvcf sample gatkh-$root varcaller gatkh varcaller_version [version gatk] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
			if {$regionfile ne ""} {
				set bedfile [tempbed $regionfile $refseq]
				lappend opts -L $bedfile
			}
			gatkexec {-XX:ParallelGCThreads=1 -d64 -Xms512m -Xmx4g} HaplotypeCaller \
				{*}$opts -R $gatkrefseq \
				-I $dep \
				-O $resultgvcf.temp.gz \
				--annotate-with-num-discovered-alleles \
				-ERC $ERC \
				-G StandardAnnotation \
				-G StandardHCAnnotation \
				-G AS_StandardAnnotation \
				2>@ stderr >@ stdout
			file rename -force $resultgvcf.temp.gz $resultgvcf
			file rename -force $resultgvcf.temp.gz.tbi $resultgvcf.tbi
			# file delete $resultgvcf.temp
		}
	} else {
		set indexdir ${pre}varall-gatkh-$root.gvcf.index
		file mkdir $indexdir
		set chromosomes [exec cut -f 1 $refseq.fai]
		set basename [gzroot $resultgvcf]
		set regfiles {}
		foreach chromosome $chromosomes {
			lappend regfiles $indexdir/$basename-$chromosome.bed
		}
		job $resultname-distrchr-beds -deps {
			$regionfile
		} -skip {
			$resultgvcf
		} -targets $regfiles -vars {
			regionfile chromosomes appdir basename indexdir
		} -code {
			cg select -f {chromosome begin end} $regionfile | $appdir/bin/distr2chr $indexdir/$basename-
			foreach chromosome $chromosomes {
				if {[file exists $indexdir/$basename-$chromosome]} {
					file rename -force $indexdir/$basename-$chromosome $indexdir/$basename-$chromosome.bed
				} else {
					file_write $indexdir/$basename-$chromosome.bed {}
				}
			}
			file delete $indexdir/$basename-chromosome
		}
		set todo {}
		foreach chromosome $chromosomes {
			lappend todo $indexdir/$basename.$chromosome.gvcf
			job $resultname-$chromosome -mem 4G \
			-deps [list {*}$deps $indexdir/$basename-$chromosome.bed] \
			-targets {
				$indexdir/$basename.$chromosome.gvcf
			} -skip {
				$resultgvcf
			} -vars {
				opts regionfile gatkrefseq refseq threads basename ERC chromosome indexdir
			} -code {
				set bedfile $indexdir/$basename-$chromosome.bed
				if {![file size $bedfile]} {
					file_write $indexdir/$basename.$chromosome.gvcf {}
					return
				}
				gatkexec {-XX:ParallelGCThreads=1 -d64 -Xms512m -Xmx4g} HaplotypeCaller \
					{*}$opts -R $gatkrefseq \
					-I $dep \
					-O $indexdir/$basename.$chromosome.gvcf.temp \
					-L $bedfile \
					--annotate-with-num-discovered-alleles \
					-ERC $ERC \
					-G StandardAnnotation \
					-G StandardHCAnnotation \
					-G AS_StandardAnnotation \
					2>@ stderr >@ stdout
				file rename -force $indexdir/$basename.$chromosome.gvcf.temp $indexdir/$basename.$chromosome.gvcf
				file delete $indexdir/$basename.$chromosome.gvcf.temp.idx
			}
		}
		job $resultname-cat -deps $todo -targets {
			$resultgvcf
		} -vars {
			todo regfiles basename regionfile indexdir resultgvcf
		} -rmtargets [list {*}$todo $regfiles] -code {
			analysisinfo_write $dep $resultgvcf sample gatkh-$basename varcaller gatkh varcaller_version [version gatk] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
			set temp [filetemp $resultgvcf].gz
			file delete [file root $temp]
			set o [open "| bgzip > $temp" w]
			set header 1
			foreach file $todo {
				if {![file size $file]} continue
				set f [open $file]
				if {$header} {
					fcopy $f $o
					set header 0
				} else {
					while {[gets $f line] != -1} {
						if {[string index $line 0] eq {#}} continue
						puts $o $line
						break
					}
					fcopy $f $o
				}
				close $f
			}
			close $o
			file rename $temp $resultgvcf
			file delete {*}$todo
			file delete {*}$regfiles
			if {![llength [glob -nocomplain $indexdir/*]]} {catch {file delete $indexdir}}
		}
		job index-$resultname -deps {$resultgvcf} -targets {$resultgvcf.tbi} -vars {
			todo regfiles basename regionfile indexdir resultgvcf
		} -code {
			gatkexec {-XX:ParallelGCThreads=1 -d64 -Xms512m -Xmx4g} IndexFeatureFile \
				-F $resultgvcf \
				2>@ stderr >@ stdout
		}
	}
	job ${pre}gvcf2tsv-$root -deps {
		$resultgvcf $gatkrefseq
	} -targets {
		${pre}uvar-gatkh-$root.tsv ${pre}uvar-gatkh-$root.tsv.analysisinfo
	} -vars {
		sample split pre root gatkrefseq resultgvcf mincoverage mingenoqual
	} -skip {
		${pre}var-gatkh-$root.tsv.lz4 ${pre}var-gatkh-$root.tsv.analysisinfo
	} -code {
		analysisinfo_write $dep $target varcaller_mincoverage $mincoverage varcaller_mingenoqual $mingenoqual varcaller_cg_version [version genomecomb]
		gatkexec {-XX:ParallelGCThreads=1 -d64 -Xms512m -Xmx4g} GenotypeGVCFs \
			-R $gatkrefseq \
			-V $resultgvcf \
			-O ${pre}uvar-gatkh-$root.temp.vcf \
			-G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation
		catch {file delete ${pre}uvar-gatkh-$root.temp.vcf.idx}
		file rename -force ${pre}uvar-gatkh-$root.temp.vcf ${pre}uvar-gatkh-$root.vcf
		set fields {chromosome begin end type ref alt quality alleleSeq1 alleleSeq2}
		lappend fields [subst {sequenced=if(\$genoqual < $mingenoqual || \$coverage < $mincoverage,"u","v")}]
		lappend fields [subst {zyg=if(\$genoqual < $mingenoqual || \$coverage < $mincoverage,"u",\$zyg)}]
		lappend fields *
		exec cg vcf2tsv -split $split -removefields {
			name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR
		} ${pre}uvar-gatkh-$root.vcf | cg select -f $fields > ${pre}uvar-gatkh-$root.tsv.temp
		file rename -force ${pre}uvar-gatkh-$root.tsv.temp ${pre}uvar-gatkh-$root.tsv
	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job ${pre}uvar-gatkh-$root.tsv ${pre}var-gatkh-$root.tsv.lz4
	# make sreg
	sreg_gatkh_job ${pre}sreg-gatkh-$root ${pre}varall-gatkh-$root.gvcf ${pre}sreg-gatkh-$root.tsv.lz4 $mincoverage $mingenoqual
	# cleanup
	if {$cleanup} {
		set cleanupfiles [list \
			${pre}var-gatkh-$root.vcf \
			${pre}uvar-gatkh-$root.tsv ${pre}uvar-gatkh-$root.tsv.index
		]
		set cleanupdeps [list ${pre}var-gatkh-$root.tsv $resultgvcf]
		cleanup_job clean_${pre}var-gatkh-$root $cleanupfiles $cleanupdeps
	}
	cd $keeppwd
	return [file join $destdir ${pre}var-gatkh-$root.tsv]
}

proc cg_var_gatkh {args} {
	set args [job_init {*}$args]
	var_gatkh_job {*}$args
	job_wait
}

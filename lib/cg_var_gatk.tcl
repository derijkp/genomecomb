proc gatk_refseq_job refseq {
	upvar job_logdir job_logdir
	set nrefseq [file root $refseq].fa
	if {![file exists $nrefseq] && $refseq ne $nrefseq} {
		mklink $refseq $nrefseq
	}
	job gatkrefseq_faidx-[file tail $nrefseq] -deps $nrefseq -targets {$nrefseq.fai} -code {
		exec samtools faidx $dep
	}
	set dict [file root $nrefseq].dict
	job gatkrefseq-[file tail $nrefseq] -deps $nrefseq -targets {$dict} -vars {nrefseq} -code {
		file delete $target.temp
		picard CreateSequenceDictionary R= $nrefseq O= $target.temp 2>@ stderr >@ stdout
		file rename -force $target.temp $target
	}
	return $nrefseq
}

proc annotvar_clusters_job {file resultfile} {
	upvar job_logdir job_logdir
	set root [join [lrange [split [file root $file] -] 1 end] -]
	job annotvar-clusters-$root -deps $file -targets reg_cluster-$root.tsv.lz4 -code {
		cg clusterregions < $dep > $target.temp
		cg lz4 $target.temp
		file rename -force $target.temp.lz4 $target
		cg lz4index $target
	}
	job annotvar-annotclusters-$root -deps {$file reg_cluster-$root.tsv.lz4} -targets {$resultfile} -code {
		cg annotate $dep $target {*}[list_remove [lrange $deps 1 end] {}]
	}
}

proc sreg_gatk_job {job varallfile resultfile} {
	upvar job_logdir job_logdir
	job $job -deps {$varallfile} -targets {$resultfile} -code {
		set temp [filetemp $target]
		set temp2 [filetemp $target]
		cg select -q {$quality >= 30 && $totalcoverage >= 5 && $type ne "ins"} -f {chromosome begin end} $dep $temp
		file_write $temp2 "# regions selected from [gzroot $dep]: \$quality >= 30 && \$totalcoverage >= 5\n"
		cg regjoin $temp >> $temp2
		if {[file extension $target] eq ".lz4"} {
			cg_lz4 -keep 0 -i 1 -o $target $temp2
		} else {
			file rename $temp2 $target
		}
		file delete $temp
	}
}

proc var_gatk_job {args} {
	upvar job_logdir job_logdir
	set pre ""
	set opts {}
	set split 0
	set deps {}
	set regionfile {}
	cg_options var_sam args {
		-L {
			lappend deps $value
		}
		-bed {
			set regionfile $value
			lappend deps $value
		}
		-pre {
			set pre $value
		}
		-split {
			set split $value
		}
		default {
			lappend opts $key $value
		}
	} {bamfile refseq}
	set gatk [gatk]
	## Produce gatk SNP calls
	set dir [file dir $bamfile]
	set keeppwd [pwd]
	cd $dir
	set file [file tail $bamfile]
	set root [join [lrange [split [file root $file] -] 1 end] -]
	set gatkrefseq [gatk_refseq_job $refseq]
	set deps [list $file $gatkrefseq $file.bai {*}$deps]
	job ${pre}varall-gatk-$root -deps $deps \
	-targets ${pre}varall-gatk-$root.vcf -skip ${pre}varall-gatk-$root.tsv -vars {gatk opts regionfile gatkrefseq refseq} -code {
		if {$regionfile ne ""} {
			set bedfile [tempbed $regionfile $refseq]
			lappend opts -L $bedfile
		}
		exec [gatkjava] -d64 -Xms512m -Xmx4g -jar $gatk -T UnifiedGenotyper \
			{*}$opts -R $dep2 -I $dep -o $target.temp \
			-stand_call_conf 10.0 -dcov 1000 \
			--annotateNDA \
			-glm SNP --output_mode EMIT_ALL_CONFIDENT_SITES 2>@ stderr >@ stdout
		file rename -force $target.temp $target
		catch {file delete $target.temp.idx}
		# file delete $target.temp
	}
	job ${pre}varall-gatk2sft-$root -deps [list ${pre}varall-gatk-$root.vcf] -targets ${pre}varall-gatk-$root.tsv.lz4 -vars {sample split} -code {
		cg vcf2tsv -split $split -removefields {name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR} $dep $target.temp.lz4
		file rename -force $target.temp.lz4 $target
	}
	# lz4_job ${pre}varall-gatk-$root.tsv -i 1
	lz4index_job ${pre}varall-gatk-$root.tsv.lz4
	# predict deletions separately, because gatk will not predict snps in a region where a deletion
	# was predicted in the varall
	job ${pre}delvar-gatk-$root -deps $deps \
	-targets ${pre}delvar-gatk-$root.vcf -skip {${pre}delvar-gatk-$root.tsv} \
	-skip {${pre}var-gatk-$root.tsv} -vars {gatk opts regionfile gatkrefseq refseq} -code {
		if {$regionfile ne ""} {
			set bedfile [tempbed $regionfile $refseq]
			lappend opts -L $bedfile
		}
		exec [gatkjava] -d64 -Xms512m -Xmx4g -jar $gatk -T UnifiedGenotyper \
			{*}$opts -R $dep2 -I $dep -o $target.temp \
			-stand_call_conf 10.0 -dcov 1000 \
			--annotateNDA \
			-glm INDEL 2>@ stderr >@ stdout
		file rename -force $target.temp $target
		catch {file delete $target.temp.idx}
		# file delete $target.temp
	}
	job ${pre}delvar-gatk2sft-$root -deps [list ${pre}delvar-gatk-$root.vcf] \
	-targets ${pre}delvar-gatk-$root.tsv \
	-skip {${pre}var-gatk-$root.tsv} -vars {sample split} -code {
		cg vcf2tsv -split $split -removefields {name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR} $dep $target.temp
		cg select -q {$alt ne "." && $alleleSeq1 ne "." &&$quality >= 10 && $totalcoverage > 4} \
		-f {
			chromosome begin end type ref alt quality alleleSeq1 alleleSeq2 
			{sequenced=if($quality < 30 || $totalcoverage < 5,"u",if($zyg eq "r","r","v"))}
			{zyg=if($quality < 30 || $totalcoverage < 5,"u",$zyg)}
			*
		} $target.temp $target.temp2
		file rename -force $target.temp2 $target
		file delete $target.temp
	}
	job ${pre}uvar-gatk-$root -deps {${pre}varall-gatk-$root.tsv ${pre}delvar-gatk-$root.tsv} \
	-targets ${pre}uvar-gatk-$root.tsv \
	-skip {${pre}var-gatk-$root.tsv} -code {
		cg select -q {$alt ne "." && $alleleSeq1 ne "." &&$quality >= 10 && $totalcoverage > 4} \
		-f {
			chromosome begin end type ref alt quality alleleSeq1 alleleSeq2 
			{sequenced=if($quality < 30 || $totalcoverage < 5,"u",if($zyg eq "r","r","v"))}
			{zyg=if($quality < 30 || $totalcoverage < 5,"u",$zyg)}
			*
		} $dep $target.temp
		cg cat $target.temp $dep2 > $target.temp2
		cg select -s - $target.temp2 $target.temp3
		file rename -force $target.temp3 $target
		file delete $target.temp
		file delete $target.temp2
	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job ${pre}uvar-gatk-$root.tsv ${pre}var-gatk-$root.tsv.lz4
	sreg_gatk_job ${pre}sreg-gatk-$root ${pre}varall-gatk-$root.tsv ${pre}sreg-gatk-$root.tsv.lz4
	## filter SNPs (according to seqanswers exome guide)
	# java -d64 -Xms512m -Xmx4g -jar $gatk -R $reference -T VariantFiltration -B:variant,VCF snp.vcf.recalibrated -o $outprefix.snp.filtered.vcf --clusterWindowSize 10 --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "DP < 5 " --filterName "LowCoverage" --filterExpression "QUAL < 30.0 " --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0 " --filterName "LowQual" --filterExpression "QD < 1.5 " --filterName "LowQD" --filterExpression "SB > -10.0 " --filterName "StrandBias"
	# cleanup
	job clean_${pre}var-gatk-$root -deps {${pre}var-gatk-$root.tsv} -vars {pre root} -targets {} \
	-rmtargets {${pre}uvar-gatk-$root.tsv ${pre}uvar-gatk-$root.tsv.index ${pre}varall-gatk-$root.vcf ${pre}delvar-gatk-$root.vcf ${pre}delvar-gatk-$root.tsv} -code {
		catch {file delete ${pre}uvar-gatk-$root.tsv}
		catch {file delete -force ${pre}uvar-gatk-$root.tsv.index}
		catch {file delete ${pre}varall-gatk-$root.vcf}
		catch {file delete ${pre}delvar-gatk-$root.vcf}
		catch {file delete ${pre}delvar-gatk-$root.tsv}
	}
	cd $keeppwd
	return [file join $dir ${pre}var-gatk-$root.tsv]
}

proc cg_var_gatk {args} {
	set args [job_init {*}$args]
	var_gatk_job {*}$args
	job_wait
}

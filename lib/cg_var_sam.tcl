proc sreg_sam_job {job varallfile resultfile} {
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
			file rename -force $temp2 $target
		}
		file delete $temp
	}
}

proc var_sam_job {args} {
	upvar job_logdir job_logdir
	set pre ""
	set opts {}
	set split 0
	set BQ 0
	cg_options var_sam args {
		-l {
			lappend deps $value
		}
		-bed {
			lappend opts -l $value
			lappend deps $value
		}
		-pre {
			set pre $value
		}
		-split {
			set split $value
		}
		-BQ {
			set BQ $value
		}
		default {
			lappend opts $key $value
		}
	} {bamfile refseq}
	set dir [file_absolute [file dir $bamfile]]
	set keeppwd [pwd]
	cd $dir
	set file [file tail $bamfile]
	set root [join [lrange [split [file root $file] -] 1 end] -]
	# make sure reference sequence is indexed
	job ${pre}var_sam_faidx -deps $refseq -targets {$refseq.fai} -code {
		exec samtools faidx $dep
	}
	set deps [list $file $refseq $refseq.fai {*}$deps]
	job ${pre}varall-sam-$root -deps $deps -targets {${pre}varall-sam-$root.vcf} \
		-vars {refseq opts BQ} -skip ${pre}varall-sam-$root.tsv -code {
		if {[catch {version samtools 1}]} {
			exec samtools mpileup -uDS -Q $BQ -f $refseq {*}$opts $dep 2>@ stderr | bcftools view -cg - > $target.temp 2>@ stderr
		} else {
			# bcftools -v for variant only
			# -t DP: Number of high-quality bases (per sample)
			# -t SP: Phred-scaled strand bias P-value
			exec samtools mpileup --uncompressed -t DP,SP --min-BQ $BQ --fasta-ref $refseq {*}$opts $dep 2>@ stderr | bcftools call -c - > $target.temp 2>@ stderr
		}
		file rename -force $target.temp $target
	}
	job ${pre}varall-sam2sft-$root -deps ${pre}varall-sam-$root.vcf -targets ${pre}varall-sam-$root.tsv.lz4 -vars split -code {
		cg vcf2tsv -split $split -removefields {name filter AN AC AF AA INDEL G3 HWE CLR UGT CGT PCHI2 QCHI2 PR} $dep $target.temp.lz4
		file rename -force $target.temp.lz4 $target
	}
	# lz4_job ${pre}varall-sam-$root.tsv -i 1
	lz4index_job ${pre}varall-sam-$root.tsv.lz4
	job ${pre}var-sam-$root -deps ${pre}varall-sam-$root.tsv -targets {${pre}uvar-sam-$root.tsv} \
	-skip {${pre}var-sam-$root.tsv} \
	-code {
		cg select -q {
				$alt ne "." && $alleleSeq1 ne "." && $quality >= 10 && $totalcoverage > 4
				&& $zyg != "r"
			} \
			-f {
				chromosome begin end type ref alt quality alleleSeq1 alleleSeq2
				{sequenced=if($quality < 30 || $totalcoverage < 5,"u",if($zyg eq "r","r","v"))}
				{zyg=if($quality < 30 || $totalcoverage < 5,"u",$zyg)}
				*
			} \
			$dep $target.temp
		file rename -force $target.temp $target
	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job ${pre}uvar-sam-$root.tsv ${pre}var-sam-$root.tsv.lz4
	# find regions
	sreg_sam_job ${pre}sreg-sam-$root ${pre}varall-sam-$root.tsv ${pre}sreg-sam-$root.tsv.lz4
	# cleanup
	job clean_${pre}var-sam-$root -deps {${pre}var-sam-$root.tsv ${pre}varall-sam-$root.tsv} -vars {pre root} -targets {} \
	-rmtargets {${pre}uvar-sam-$root.tsv ${pre}uvar-sam-$root.tsv.index ${pre}varall-sam-$root.vcf ${pre}varall-sam-$root.vcf.idx} -code {
		catch {file delete ${pre}uvar-sam-$root.tsv}
		catch {file delete -force ${pre}uvar-sam-$root.tsv.index}
		catch {file delete ${pre}varall-sam-$root.vcf}
		catch {file delete ${pre}varall-sam-$root.vcf.idx}
	}
	cd $keeppwd
	return [file join $dir var-sam-$root.tsv]
}

proc cg_var_sam {args} {
	set args [job_init {*}$args]
	var_sam_job {*}$args
	job_wait
}

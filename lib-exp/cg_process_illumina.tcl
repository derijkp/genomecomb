proc picard {} {
	global picard
	if {![info exists picard]} {
		if {[info exists ::env(PICARD)]} {
			set picard $::env(PICARD)
		} else {
 			set picard $::externdir/picard
		}
	}
	return $picard
}

proc gatk {} {
	global gatk
	if {![info exists gatk]} {
		if {[info exists ::env(GATK)]} {
			set gatk $::env(GATK)
		} else {
 			set gatk $::externdir/gatk/GenomeAnalysisTK.jar
		}
	}
	return $gatk
}

proc fastq_clipadapters {files targets {adapterfile {}}} {
	if {$adapterfile eq ""} {
		set adapterfile $::externdir/adaptors.fa
	}
	# clip primers, quality
	set out {}
	foreach target $targets {
		file mkdir [file dir $target]
		lappend out -o $target.temp
	}
	exec fastq-mcf {*}$out $adapterfile {*}$files 2>@ stderr
	foreach target $targets {
		file rename $target.temp $target
	}
}

proc fastq_clipadapters_job {files {adapterfile {}}} {
	upvar job_logdir job_logdir
	set targets {}
	foreach file $files {
		set file [file normalize [gzroot $file]]
		set root [file root $file]
		lappend targets [file dir $root].clipped/[file tail $root].clipped.fastq
	}
	job clip-[file dir [file dir $root]] -deps $files -targets $targets -code {
		fastq_clipadapters $deps $targets
	}
	return $targets
}

proc bowtie2refseq_job {refseq} {
	upvar job_logdir job_logdir
	set bowtie2refseq $refseq.bowtie2/[file tail $refseq]
	job bowtie2refseq-[file tail $refseq] -deps $refseq -targets {$refseq.bowtie2 $bowtie2refseq.1.bt2} \
	-vars {refseq} -code {
		file mkdir $refseq.bowtie2.temp
		mklink $refseq $refseq.bowtie2.temp/[file tail $refseq]
		exec bowtie2-build $refseq $refseq.bowtie2.temp/[file tail $refseq]
		file rename $refseq.bowtie2.temp $refseq.bowtie2
	}
	return $bowtie2refseq
}

proc map_bowtie2_job {refseq files sample {readgroupdata {}} {pre {}}} {
	array set a [list PL illumina LB solexa-123 PU $sample SM $sample]
	if {$readgroupdata ne ""} {
		array set a $readgroupdata
	}
	set result ${pre}map-bowtie2-$sample
	set readgroupdata [array get a]
	upvar job_logdir job_logdir
	set bowtie2refseq [bowtie2refseq_job $refseq]
	job bowtie2-$sample -deps [list $bowtie2refseq {*}$files] -targets $result.sam -vars {bowtie2refseq readgroupdata sample} -skip $result.bam -code {
		puts "making $target"
		list_shift deps
		set rg {}
		foreach {key value} $readgroupdata {
			lappend rg --rg "$key:$value"
		}
		set files1 {}
		set files2 {}
		foreach {file1 file2} $deps {
			lappend files1 $file1
			lappend files2 $file2
		}
		exec bowtie2 -p 2 --sensitive -x $bowtie2refseq -1 [join $files1 ,] -2 [join $files2 ,] \
		--rg-id "$sample" {*}$rg \
		-S $target.temp >@ stdout 2>@ stderr
		file rename $target.temp $target
	}
	job bowtie2_bam -deps $result.sam -targets $result.bam -vars {result} -code {
		puts "making $target"
		catch {exec samtools view -S -h -b -o $result.ubam $result.sam >@ stdout 2>@ stderr}
		catch {exec samtools sort $result.ubam $result.temp >@ stdout 2>@ stderr}
		file rename $result.temp.bam $result.bam
		file delete $result.ubam
		file delete $result.sam
	}
	job bowtie2_index -deps $result.bam -targets $result.bam.bai -code {
		exec samtools index $dep >@ stdout 2>@ stderr
		puts "making $target"
	}
#	job bowtie2_coverage -deps $result.bam -targets {coverage coverage/coverage-$sample-.FINISHED} -vars {sample} -code {
#		cg bam2coverage $dep coverage-bowtie2-$sample/coverage-bowtie2-$sample
#	}
#	job bowtie2_coverage -deps $result.bam -targets sreg-$sample.tsv -vars {sample} -code {
#		cg regextract -above 1 7 $dep > $target.temp
#		file rename $target.temp $target
#	}
}

proc bwarefseq_job {refseq} {
	upvar job_logdir job_logdir
	set bwarefseq $refseq.bwa/[file tail $refseq]
	job bowtie2refseq-[file tail $refseq] -deps $refseq -targets {$refseq.bwa} -code {
		file delete -force $target.temp
		file mkdir $target.temp
		mklink $dep $target.temp/[file tail $dep]
		exec bwa index $target.temp/[file tail $dep] 2>@ stderr
		file rename $target.temp $target
	}
	return $bwarefseq
}

proc map_bwa_job {refseq files sample {readgroupdata {}} {pre {}}} {
	upvar job_logdir job_logdir
	array set a [list PL illumina LB solexa-123 PU $sample SM $sample]
	if {$readgroupdata ne ""} {
		array set a $readgroupdata
	}
	set result ${pre}map-bwa-$sample
	set readgroupdata [array get a]
	set bwarefseq [bwarefseq_job $refseq]
	job bwa -deps [list $bwarefseq {*}$files] -targets $result.sam -vars {readgroupdata sample} -skip $result.bam -code {
		puts "making $target"
		set bwarefseq [list_shift deps]
		set rg {}
		foreach {key value} $readgroupdata {
			lappend rg "$key:$value"
		}
		set files1 {}
		set files2 {}
		foreach {file1 file2} $deps {
			lappend files1 $file1
			lappend files2 $file2
		}
		file delete bwa1.fastq bwa2.fastq
		if {[llength $files1] > 1} {
			exec cat {*}$files1 > bwa1.fastq
			exec cat {*}$files2 > bwa2.fastq
		} else {
			mklink $file1 bwa1.fastq
			mklink $file2 bwa2.fastq
		}
		exec bwa mem -t 2 -a -M -R @RG\tID:$sample\t[join $rg \t] $bwarefseq bwa1.fastq bwa2.fastq > $target.temp 2>@ stderr
		file rename $target.temp $target
		file delete bwa1.fastq bwa2.fastq
	}
	job bwa_bam -deps $result.sam -targets $result.bam -vars {result} -code {
		puts "making $target"
		catch {exec samtools view -S -h -b -o $result.ubam $result.sam >@ stdout 2>@ stderr}
		catch {exec samtools sort $result.ubam $result.temp >@ stdout 2>@ stderr}
		file rename $result.temp.bam $result.bam
		file delete $result.ubam
		file delete $result.sam
	}
	job bwa_index -deps $result.bam -targets $result.bam.bai -code {
		exec samtools index $dep >@ stdout 2>@ stderr
		puts "making $target"
	}
#	job bwa_coverage -deps $result.bam -targets {coverage coverage/coverage-$sample-.FINISHED} -vars {sample} -code {
#		cg bam2coverage $dep coverage-bwa-$sample/coverage-bwa-$sample
#	}
#	job bwa_coverage -deps $result.bam -targets sreg-$sample.tsv -vars {sample} -code {
#		cg regextract -above 1 7 $dep > $target.temp
#		file rename $target.temp $target
#	}
}

proc gatk_refseq_job refseq {
	upvar job_logdir job_logdir
	set nrefseq [file root $refseq].fa
	file delete $nrefseq
	exec ln -s $refseq $nrefseq
	file mkdir $refseq.gatk
	set picard [picard]
	job bam_clean_sam_faidx-[file tail $nrefseq] -deps $nrefseq -targets {$nrefseq.fai} -code {
		exec samtools faidx $dep
	}
	set dict [file root $nrefseq].dict
	job gatkrefseq-[file tail $nrefseq] -deps $nrefseq -targets {$dict} -vars {nrefseq picard} -code {
		file delete $target.temp
		exec java -jar $picard/CreateSequenceDictionary.jar R= $nrefseq O= $target.temp 2>@ stderr > stdout
		file rename $target.temp $target
	}
	return $nrefseq
}

proc bam_clean_job {bamfile refseq sample args} {
	set gatk [gatk]
	set picard [picard]
	upvar job_logdir job_logdir
	array set opt $args
	set removeduplicates 1
	set realign 1
	foreach {key value} $args {
		switch -- $key {
			-removeduplicates {set removeduplicates $value}
			-realign {set realign $value}
			default {error "bam_clean_job: unknown option $key"}
		}
	}
	set pre [lindex [split $bamfile -] 0]
	set dir [file dir $bamfile]
	set file [file tail $bamfile]
	set root [join [lrange [split [file root $file] -] 1 end] -]
	# make gatk refseq
	set gatkrefseq [gatk_refseq_job $refseq]
	set dict [file root $gatkrefseq].dict
	# clean bamfile: sort using picard, mark duplicates, and add readgroup info (needed for gatk)
	if {![file exists $dir/$pre-ds$root.bam] && ![file exists $dir/$pre-rs$root.bam]} {
		job bamsort-$root -deps {$bamfile} -targets {$dir/$pre-s$root.bam} \
		-vars {removeduplicates sample picard} -skip {$dir/$pre-rds$root.bam} -code {
			file delete $target.temp
			exec java -jar $picard/SortSam.jar	I=$dep	O=$target.temp	SO=coordinate 2>@ stderr > stdout
			file rename $target.temp $target
		}
	}
	set root s$root
	if {$removeduplicates} {
		job bamremdup-$root -deps {$dir/$pre-$root.bam} -targets {$dir/$pre-d$root.bam} \
		-vars {sample picard} -skip {$dir/$pre-rd$root.bam} -code {
			puts "removing duplicates"
			exec java -jar $picard/MarkDuplicates.jar	I=$dep	O=$target.temp METRICS_FILE=$target.dupmetrics 2>@ stderr > stdout
			file rename $target.temp $target
	#		set src $target.temp2
	#		# exec java -jar $picard/AddOrReplaceReadGroups.jar	I=$src	O=$target.temp3	RGID=$sample	RGLB=solexa-123	RGPL=illumina	RGPU=$sample RGSM=$sample 2>@ stderr > stdout
	#		# file delete $target.temp $target.temp2
	#		# file rename $target.temp3 $target
	#		file rename $src $target
	#		file delete $target.temp $target.temp2
		}
		set root d$root
		# index cleaned result
		job bamrs_index-$root -deps $dir/$pre-$root.bam -targets $dir/$pre-$root.bam.bai \
		-code {
			exec samtools index $dep >@ stdout 2>@ stderr
			puts "making $target"
		}
	}
	if {$realign} {
		# realign around indels
		job bamrealign-$root -deps {$dir/$pre-$root.bam $dict} -targets {$dir/$pre-r$root.bam} \
		-vars {gatkrefseq gatk pre} -code {
			exec java -jar $gatk -T RealignerTargetCreator -R $gatkrefseq -I $dep -o $target.intervals 2>@ stderr >@ stdout
			exec java -jar $gatk -T IndelRealigner -R $gatkrefseq -targetIntervals $target.intervals -I $dep -o $target.temp 2>@ stderr >@ stdout
			file rename $target.temp $target
			file delete $dep $dep.bai
		}
		set root r$root
		job bamrealign_index-$root -deps $dir/$pre-$root.bam -targets $dir/$pre-$root.bam.bai -code {
			exec samtools index $dep >@ stdout 2>@ stderr
			puts "making $target"
		}
	}
	return $dir/$pre-$root.bam
}

proc annotvar_clusters_job {file resultfile} {
	upvar job_logdir job_logdir
	set root [join [lrange [split [file root $file] -] 1 end] -]
	job annotvar-clusters-$root -deps $file -targets reg_cluster-$root.tsv -code {
		cg clusterregions < $dep > $target.temp
		file rename -force $target.temp $target
	}
	job annotvar-annotclusters-$root -deps {$file reg_cluster-$root.tsv} -targets {$resultfile} -code {
		cg annotate $dep $target {*}[list_remove [lrange $deps 1 end] {}]
	}
}

proc sreg_sam_job {job varallfile resultfile} {
	upvar job_logdir job_logdir
	job $job -deps {$varallfile} -targets {$resultfile} -code {
		cg select -q {$quality >= 20 && $totalcoverage >= 5 && $type ne "ins"} -f {chromosome begin end} $dep $target.temp
		file_write $target.temp2 "# regions selected from $dep: \$quality >= 20 && \$totalcoverage >= 5\n"
		cg regjoin $target.temp >> $target.temp2
		file rename $target.temp2 $target
		file delete $target.temp
	}
}

proc var_sam_job {bamfile refseq {pre {}}} {
	upvar job_logdir job_logdir
	set dir [file normalize [file dir $bamfile]]
	set keeppwd [pwd]
	cd $dir
	set file [file tail $bamfile]
	set root [join [lrange [split [file root $file] -] 1 end] -]
	# make sure reference sequence is indexed
	job ${pre}var_sam_faidx -deps $refseq -targets {$refseq.fai} -code {
		exec samtools faidx $dep
	}
	job ${pre}varall-sam-$root -deps {$file $refseq.fai} -targets {${pre}varall-sam-$root.vcf} -vars {refseq} -skip varall-sam-$root.tsv -code {
		# bcftools -v for variant only
		exec samtools mpileup -uDS -f $refseq $dep 2>@ stderr | bcftools view -cg - > $target.temp 2>@ stderr
		file rename $target.temp $target
	}
	job ${pre}varall-sam_2sft-$root -deps ${pre}varall-sam-$root.vcf -targets ${pre}varall-sam-$root.tsv -code {
		cg vcf2sft $dep $target.temp
		cg select -s - $target.temp $target.temp2
		file delete $target.temp
		file rename $target.temp2 $target
	}
	job ${pre}var-sam-$root -deps ${pre}varall-sam-$root.tsv -targets {${pre}var-sam-$root.tsv} -code {
		cg select -s - -q {$alt ne "." && $alleleSeq1 ne "." &&$quality >= 5 && $totalcoverage > 3} \
			-f {chromosome begin end type ref alt name quality filter alleleSeq1 alleleSeq2 {sequenced=if($quality < 30 || $totalcoverage < 5,"u",if($zyg eq "r","r","v"))} *} \
			$dep $target.temp
		file rename $target.temp $target
	}
	sreg_sam_job ${pre}sreg-sam-$root ${pre}varall-sam-$root.tsv ${pre}sreg-sam-$root.tsv
#	job ${pre}var-sam-$root -deps $file -targets {var-sam-$root.vcf} -vars {refseq} -code {
#		# bcftools -v for variant only
#		exec samtools mpileup -uDS -f $refseq $dep 2>@ stderr | bcftools view -vcg - > $target.temp 2>@ stderr
#		file rename $target.temp $target
#	}
#	job ${pre}var-sam2sft-$root -deps var-sam-$root.vcf -targets var-sam-$root.tsv -vars {refseq} -code {
#		cg vcf2sft $dep $target.temp
#		file rename $target.temp $target
#	}
	cd $keeppwd
	return [file join $dir var-sam-$root.tsv]
}

proc sreg_gatk_job {job varallfile resultfile} {
	upvar job_logdir job_logdir
	job $job -deps {$varallfile} -targets {$resultfile} -code {
		cg select -q {$quality >= 30 && $totalcoverage >= 5 && $type ne "ins"} -f {chromosome begin end} $dep $target.temp
		file_write $target.temp2 "# regions selected from $dep: \$quality >= 30 && \$totalcoverage >= 5\n"
		cg regjoin $target.temp >> $target.temp2
		file rename $target.temp2 $target
		file delete $target.temp
	}
}

proc var_gatk_job {bamfile refseq {pre {}}} {
	upvar job_logdir job_logdir
	set gatk [gatk]
	## Produce gatk SNP calls
	set dir [file dir $bamfile]
	set keeppwd [pwd]
	cd $dir
	set file [file tail $bamfile]
	set root [join [lrange [split [file root $file] -] 1 end] -]
	set gatkrefseq [gatk_refseq_job $refseq]
	job ${pre}varall-gatk-$root -deps [list $file $gatkrefseq] \
	-targets ${pre}varall-gatk-$root.vcf -skip ${pre}varall-gatk-$root.tsv -vars {gatk} -code {
		exec java -d64 -Xms512m -Xmx4g -jar $gatk -T UnifiedGenotyper -R $dep2 -I $dep -o $target.temp \
			-stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1000 \
			-glm BOTH --output_mode EMIT_ALL_CONFIDENT_SITES 2>@ stderr
		file rename $target.temp $target
		catch {file rename $target.temp.idx $target.idx}
		# file delete $target.temp
	}
	job ${pre}varall-gatk2sft-$root -deps [list ${pre}varall-gatk-$root.vcf] -targets ${pre}varall-gatk-$root.tsv -vars {sample} -code {
		cg vcf2sft $dep $target.temp
		file rename $target.temp $target
	}
	job ${pre}uvar-gatk-$root -deps ${pre}varall-gatk-$root.tsv -targets ${pre}uvar-gatk-$root.tsv -skip {${pre}var-gatk-$root.tsv} -code {
		cg select -s - -q {$alt ne "." && $alleleSeq1 ne "." &&$quality >= 10 && $totalcoverage > 3} \
			-f {chromosome begin end type ref alt name quality filter alleleSeq1 alleleSeq2 {sequenced=if($quality < 30 || $totalcoverage < 5,"u",if($zyg eq "r","r","v"))} *} \
			$dep $target.temp
		file rename $target.temp $target
	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job ${pre}uvar-gatk-$root.tsv ${pre}var-gatk-$root.tsv
	sreg_gatk_job ${pre}sreg-gatk-$root ${pre}varall-gatk-$root.tsv ${pre}sreg-gatk-$root.tsv
	## filter SNPs (according to seqanswers exome guide)
	# java -d64 -Xms512m -Xmx4g -jar $gatk -R $reference -T VariantFiltration -B:variant,VCF snp.vcf.recalibrated -o $outprefix.snp.filtered.vcf --clusterWindowSize 10 --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "DP < 5 " --filterName "LowCoverage" --filterExpression "QUAL < 30.0 " --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0 " --filterName "LowQual" --filterExpression "QD < 1.5 " --filterName "LowQD" --filterExpression "SB > -10.0 " --filterName "StrandBias"
	cd $keeppwd
	return [file join $dir ${pre}var-gatk-$root.tsv]
}

proc multicompar_job {experiment dbdir todo} {
	upvar job_logdir job_logdir
	file mkdir compar
	if {[catch {cg select -n compar/compar-$experiment.tsv} done]} {set done {}}
	set stilltodo {}
	set deps {}
	foreach sample $todo {
		set name [lindex [split $sample -] end]
		if {![inlist $done $sample]} {
			lappend stilltodo $name/var-$sample.tsv
			lappend deps $name/sreg-$sample.tsv
		}
	}
	if {[llength $stilltodo]} {
		file delete compar/compar-$experiment.tsv.temp
		if {[file exists compar/compar-$experiment.tsv]} {
			file rename compar/compar-$experiment.tsv compar/compar-$experiment.tsv.temp
		}
		job multicompar-$experiment -deps $stilltodo -targets compar/compar-$experiment.tsv -vars stilltodo -code {
			# should maybe better recheck todo here
			cg multicompar $target.temp {*}$stilltodo
			cg multicompar_reannot $target.temp
			file rename $target.temp $target
		}
	}
	job annotcompar-$experiment -deps compar/compar-$experiment.tsv \
	-targets compar/annot_compar-$experiment.tsv -vars dbdir -code {
		cg annotate $dep $target.temp $dbdir
		file rename $target.temp $target
	}
	job indexannotcompar-$experiment \
	-deps compar/annot_compar-$experiment.tsv \
	-targets compar/annot_compar-$experiment.tsv.index/info.tsv -vars dbdir -code {
		cg index $dep
	}
	if {[catch {cg select -n compar/sreg-$experiment.tsv} done]} {set done {}}
	set stilltodo {}
	foreach sample $todo {
		set name [lindex [split $sample -] end]
		if {![inlist $done $sample]} {
			lappend stilltodo $name/sreg-$sample.tsv
		}
	}
	if {[llength $stilltodo]} {
		file delete compar/sreg-$experiment.tsv.temp
		if {[file exists compar/sreg-$experiment.tsv]} {
			file rename compar/sreg-$experiment.tsv compar/sreg-$experiment.tsv.temp
		}
		job sreg-$experiment -deps $stilltodo -targets compar/sreg-$experiment.tsv -vars stilltodo -code {
			cg multireg $target.temp {*}$stilltodo
			file rename $target.temp $target
		}
	}
}

proc process_illumina {destdir dbdir} {
	set refseq [glob $dbdir/genome_*.ifas]
	set destdir [file normalize $destdir]
	set samples {}
	foreach dir [dirglob $destdir */fastq] {
		lappend samples [file dir $dir]
	}
	set samples [lsort -dict $samples]
	set keeppwd [pwd]
	cd $destdir
	job_logdir $destdir/log_jobs
	foreach sample $samples {
		puts $sample
		set dir $destdir/$sample
		catch {file mkdir $dir}
		puts $dir
		cd $dir
		job_logdir $dir/log_jobs
		# job_logdir $dir/log_jobs
		set files [glob -nocomplain fastq/*.fastq.gz fastq/*.fastq]
		if {![llength $files]} continue
		# quality and adapter clipping
		set files [fastq_clipadapters_job $files]
		#
#		# map using bowtie2
#		map_bowtie2_job $refseq $sample $files
#		# clean bamfile (mark duplicates, realign)
#		bam_clean_job map-bowtie2-$sample.bam $refseq $sample
#		# samtools variant calling on map-bowtie2
#		#var_sam_job map-bowtie2-$sample.bam $refseq
#		# samtools variant calling on map-rdsbowtie2
#		var_sam_job map-rdsbowtie2-$sample.bam $refseq
#		# gatk variant calling on map-rdsbowtie2
#		var_gatk_job map-rdsbowtie2-$sample.bam $refseq
		#
		# map using bwa
		map_bwa_job $refseq $files $sample
		# clean bamfile (mark duplicates, realign)
		set cleanedbam [bam_clean_job map-bwa-$sample.bam $refseq $sample -removeduplicates 1 -realign 0]
		# samtools variant calling on map-rdsbwa
		var_sam_job $cleanedbam $refseq
		# gatk variant calling on map-rdsbwa
		var_gatk_job $cleanedbam $refseq
	}
	job_logdir $destdir/log_jobs
	cd $destdir
	set todo {}
	foreach sample $samples {
		lappend todo sam-rdsbowtie2-$sample gatk-rdsbowtie2-$sample
	}
	multicompar_job $experiment $dbdir $todo
	cd $keeppwd

}

proc cg_process_illumina {args} {
	set args [job_init {*}$args]
	if {[llength $args] < 2} {
		errorformat process_illumina
		exit 1
	}
	process_illumina {*}$args
	job_wait
}

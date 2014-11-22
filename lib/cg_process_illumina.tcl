proc bam2covstats_job {bamfile regionfile} {
	upvar job_logdir job_logdir
	set bamfile [file_absolute $bamfile]
	set dir [file dir $bamfile]
	set file [file tail $bamfile]
	set root [join [lrange [split [file root $file] -] 1 end] -]
#	job bam2coverage-$root -deps $bamfile -targets {$dir/coverage-$root $dir/coverage-$root/coverage-$root.FINISHED} -vars {root} -code {
#		cg bam2coverage $dep $target/coverage-$root
#	}
	job make_histo-$root -deps {$bamfile $bamfile.bai $regionfile} -targets $dir/$root.histo -vars {regionfile} -code {
		set tempfile [file_tempwrite $target]
		cg bam_histo $regionfile $dep {1 5 10 20 50 100 200 500 1000} > $tempfile
		file rename -force $tempfile $target
	}
}

proc cg_bcl2fastq {rundir outdir {rtr 6} {dtr 6} {ptr 6} {wtr 6} } {
	#-r, --loading-threads Number of threads used for loading BCL data.
	#-d, --demultiplexing-threads Number of threads used for demultiplexing.
	#-p, --processing-threads  Number of threads used for processing demultiplexed data.
	#-w, --writing-threads Number of threads used for writing FASTQ data. This must not be higher than number of samples.
	exec nohup /complgen/bin/bcl2fastq --create-fastq-for-index-reads -r $rtr -d $dtr -p $ptr -w $wtr --runfolder-dir $rundir --output-dir $outdir  
}

proc cg_process_conv_illnextseq {illsrc destdir} {
	set illsrc [file_absolute $illsrc]
	set destdir [file_absolute $destdir]
	file mkdir $destdir
	set keeppwd [pwd]
	cd $destdir
	# copy files from illumina dir
	set files [glob $illsrc/*_R*.fastq*]
	foreach file $files {
		set sample [file tail $file]
		regsub {_[^_]+_[^_]+_[^_]+_[^_]+\.fastq.*} $sample {} sample
		regsub -all -- - $sample _ sample 
		if {$sample != "Undetermined"} {
			file mkdir $destdir/$sample/ori
			file mkdir $destdir/$sample/ori/fastq
			file mkdir $destdir/$sample
			file mkdir $destdir/$sample/fastq
			exec cp -al $file $destdir/$sample/ori/fastq
			cplinked $destdir/$sample/ori/fastq $destdir/$sample/fastq
		}
	}
	cd $keeppwd
}

proc searchpath {envvar args} {
	set name [lindex $args 0]
	if {[info exists ::env($envvar)]} {
		if {![file exists $::env($envvar)]} {
			error "$name not found at $::env($envvar) (from env var $envvar)"
		}
		return $::env($envvar)
	} else {
		set dirlist [split [get ::env(PATH) ""] :]
		list_addnew dirlist $::externdir
		foreach pattern $args {
			foreach dir $dirlist {
				if {![catch {glob $dir/$pattern} dirs]} {
					return [lindex [lsort -dict $dirs] 0]
				}
			}
		}
		error "$name not found in PATH"
	}
}

proc picard {} {
	global picard
	if {![info exists picard]} {
		set picard [searchpath PICARD picard picard*]
	}
	return $picard
}

proc gatk {} {
	global gatk
	if {![info exists gatk]} {
		set gatk [searchpath GATK gatk GenomeAnalysisTK*]/GenomeAnalysisTK.jar
	}
	return $gatk
}

proc fastq_clipadapters {files targets args} {
	set adapterfile {}
	set paired 1
	foreach {key value} $args {
		if {$key eq "-adapterfile"} {
			set adapterfile $value
		} elseif {$key eq "-paired"} {
			set paired $value
		} else {
			lappend opts $key $value
		}
	}
	if {$adapterfile eq ""} {
		set adapterfile $::externdir/adaptors.fa
	}
	# clip primers, quality
	set temptargets {}
	if {[llength $files] == 1 || !$paired} {
		foreach {f1} $files {t1} $targets {
			set tempout1 [file_tempwrite $t1]
			exec fastq-mcf -a -o $tempout1 $adapterfile $f1 2>@ stderr
			lappend temptargets $tempout1
		}
	} else {
		foreach {f1 f2} $files {t1 t2} $targets {
			set tempout1 [file_tempwrite $t1]
			set tempout2 [file_tempwrite $t2]
			exec fastq-mcf -a -o $tempout1 -o $tempout2 $adapterfile $f1 $f2 2>@ stderr
			lappend temptargets $tempout1 $tempout2
		}
	}
	foreach target $targets temptarget $temptargets {
		file rename -force $temptarget $target
	}
}

proc fastq_clipadapters_job {files args} {
	upvar job_logdir job_logdir
	set targets {}
	set skips {}
	set adapterfile {}
	set paired 1
	set files [ssort -natural $files]
	foreach {key value} $args {
		if {$key eq "-adapterfile"} {
			set adapterfile $value
		} elseif {$key eq "-paired"} {
			set paired $value
		} elseif {$key eq "-skips"} {
			set skips $value
		} else {
			lappend opts $key $value
		}
	}
	foreach file $files {
		set file [file_absolute [gzroot $file]]
		set root [file root $file]
		file mkdir [file dir $root].clipped
		lappend targets [file dir $root].clipped/[file tail $root].clipped.fastq
	}
	job clip-[file dir [file dir $root]] -deps $files -targets $targets \
	-vars {adapterfile paired} {*}$skips -code {
		fastq_clipadapters $deps $targets -adapterfile $adapterfile -paired $paired
	}
	return $targets
}

proc gatkworkaround_tsv2bed_job {file refseq} {
	upvar job_logdir job_logdir
	job tsv2bed-[file tail $file] -deps {$file $refseq.index} -targets [file root $file].bed -code {
		set f [open $dep2]
		while {![eof $f]} {
			set line [split [gets $f] \t]
			set chr [lindex $line 0]
			set chr [chr_clip $chr]
			set maxa($chr) [lindex $line 1 1]
		}
		close $f
		set f [gzopen $dep]
		set header [tsv_open $f]
		set poss [tsv_basicfields $header 3]
		set temptarget [file_tempwrite $target]
		set o [open $temptarget w]
		while {![eof $f]} {
			set line [split [gets $f] \t]
			set line [list_sub $line $poss]
			foreach {chr begin end} $line break
			set cchr [chr_clip $chr]
			if {[info exists maxa($cchr)]} {
				if {$end > $maxa($cchr)} {set end $maxa($cchr)}
			}
			if {$end == $begin} continue
			puts $o $chr\t$begin\t$end
		}
		close $o
		close $f
		file rename -force $temptarget $target
	}
	return [file root $file].bed
}

proc bowtie2refseq_job {refseq} {
	upvar job_logdir job_logdir
	set bowtie2refseq $refseq.bowtie2/[file tail $refseq]
	job bowtie2refseq-[file tail $refseq] -deps $refseq -targets {$refseq.bowtie2 $bowtie2refseq.1.bt2} \
	-vars {refseq} -code {
		file mkdir $refseq.bowtie2.temp
		mklink $refseq $refseq.bowtie2.temp/[file tail $refseq]
		exec bowtie2-build $refseq $refseq.bowtie2.temp/[file tail $refseq]
		file rename -force $refseq.bowtie2.temp $refseq.bowtie2
	}
	return $bowtie2refseq
}

proc map_bowtie2_job {args} {
	oargs map_bowtie2_job {refseq files sample 
		{paired 1}
		{readgroupdata {}}
		{pre {}}
		{skips {}}
	} $args
	array set a [list PL illumina LB solexa-123 PU $sample SM $sample]
	if {$readgroupdata ne ""} {
		array set a $readgroupdata
	}
	set result ${pre}map-bowtie2-$sample
	set readgroupdata [array get a]
	upvar job_logdir job_logdir
	set bowtie2refseq [bowtie2refseq_job $refseq]
	job bowtie2-$sample -deps [list $bowtie2refseq {*}$files] -targets $result.sam \
	-vars {paired bowtie2refseq readgroupdata sample} \
	-skip $result.bam {*}$skips -code {
		puts "making $target"
		list_shift deps
		set rg {}
		foreach {key value} $readgroupdata {
			lappend rg --rg "$key:$value"
		}
		set temptarget [file_tempwrite $target]
		if {$paired} {
			set files1 {}
			set files2 {}
			foreach {file1 file2} $deps {
				lappend files1 $file1
				lappend files2 $file2
			}
			exec bowtie2 -p 2 --sensitive -x $bowtie2refseq -1 [join $files1 ,] -2 [join $files2 ,] \
			--rg-id "$sample" {*}$rg \
			-S $temptarget >@ stdout 2>@ stderr
		} else {
			exec bowtie2 -p 2 --sensitive -x $bowtie2refseq -U [join $$deps ,] \
			--rg-id "$sample" {*}$rg \
			-S $temptarget >@ stdout 2>@ stderr
		}
		file rename -force $temptarget $target
	}
	job bowtie2_bam-$sample -deps $result.sam -targets $result.bam -vars {result} {*}$skips -code {
		puts "making $target"
		catch {exec samtools view -S -h -b -o $result.ubam $result.sam >@ stdout 2>@ stderr}
		catch {exec samtools sort $result.ubam $result.temp >@ stdout 2>@ stderr}
		file rename -force $result.temp.bam $result.bam
		file delete $result.ubam
		file delete $result.sam
	}
	job bowtie2_index-$sample -deps $result.bam -targets $result.bam.bai {*}$skips -code {
		exec samtools index $dep >@ stdout 2>@ stderr
		puts "making $target"
	}
}


proc bam2reg_job {bamfile {mincoverage 5}} {
	upvar job_logdir job_logdir
	set bamfile [file_absolute $bamfile]
	set pre [lindex [split $bamfile -] 0]
	set dir [file dir $bamfile]
	set file [file tail $bamfile]
	set root [join [lrange [split [file root $file] -] 1 end] -]
#	job bam2coverage-$root -deps $bamfile -targets {$dir/coverage-$root $dir/coverage-$root/coverage-$root.FINISHED} -vars {root} -code {
#		cg bam2coverage $dep $target/coverage-$root
#	}
	job cov$mincoverage-$root -deps $bamfile -targets $dir/sreg-cov$mincoverage-$root.tsv -vars {mincoverage} -code {
		set temptarget [file_tempwrite $target]
		cg regextract -above 1 [expr {$mincoverage-1}] $dep > $temptarget
		file rename -force $temptarget $target
	}
	return $dir/sreg-cov$mincoverage-$root.tsv
}

proc bwarefseq_job {refseq} {
	upvar job_logdir job_logdir
	set bwarefseq $refseq.bwa/[file tail $refseq]
	if {[file exists $bwarefseq]} {return $bwarefseq}
	job bwa2refseq-[file tail $refseq] -deps $refseq -targets {$refseq.bwa} -code {
		file delete -force $target.temp
		file mkdir $target.temp
		mklink $dep $target.temp/[file tail $dep]
		exec bwa index $target.temp/[file tail $dep] 2>@ stderr
		file rename -force $target.temp $target
	}
	return $bwarefseq
}

proc map_bwa_job {args} {
	oargs map_bwa_job {refseq files sample
		{paired 1}
		{readgroupdata {}}
		{pre {}}
		{skips {}}
	} $args
	upvar job_logdir job_logdir
	array set a [list PL illumina LB solexa-123 PU $sample SM $sample]
	if {$readgroupdata ne ""} {
		array set a $readgroupdata
	}
	set result ${pre}map-bwa-$sample
	set readgroupdata [array get a]
	set bwarefseq [bwarefseq_job $refseq]
	job bwa-$sample -deps [list $bwarefseq {*}$files] -targets $result.sam -vars {readgroupdata sample paired} \
	-skip $result.bam {*}$skips -code {
		puts "making $target"
		set bwarefseq [list_shift deps]
		set rg {}
		foreach {key value} $readgroupdata {
			lappend rg "$key:$value"
		}
		set tempdir [scratchdir]
		if {!$paired} {
			if {[llength $deps] > 1} {
				exec cat {*}$deps > $tempdir/bwa1.fastq
			} else {
				mklink $deps $tempdir/bwa1.fastq
			}
			exec bwa mem -t 2 -a  -R @RG\tID:$sample\t[join $rg \t] $bwarefseq $tempdir/bwa1.fastq > $target.temp 2>@ stderr
		} else {
			set files1 {}
			set files2 {}
			foreach {file1 file2} $deps {
				lappend files1 $file1
				lappend files2 $file2
			}
			if {[llength $files1] > 1} {
				exec cat {*}$files1 > $tempdir/bwa1.fastq
				exec cat {*}$files2 > $tempdir/bwa2.fastq
			} else {
				mklink $file1 $tempdir/bwa1.fastq
				mklink $file2 $tempdir/bwa2.fastq
			}
			exec bwa mem -t 2 -a -M -R @RG\tID:$sample\t[join $rg \t] $bwarefseq $tempdir/bwa1.fastq $tempdir/bwa2.fastq > $target.temp 2>@ stderr
		}
		file rename -force $target.temp $target
		file delete bwa1.fastq bwa2.fastq
	}
	job bwa2bam-$sample -deps $result.sam -targets $result.bam {*}$skips -vars {result} -code {
		puts "making $target"
		catch {exec samtools view -S -h -b -o $result.ubam $result.sam >@ stdout 2>@ stderr}
		catch {exec samtools sort $result.ubam $result.temp >@ stdout 2>@ stderr}
		file rename -force $result.temp.bam $result.bam
		file delete $result.ubam
		file delete $result.sam
	}
	job bwa_index-$sample -deps $result.bam -targets $result.bam.bai {*}$skips -code {
		exec samtools index $dep >@ stdout 2>@ stderr
		puts "making $target"
	}
#	job bwa_coverage-$sample -deps $result.bam -targets {coverage coverage/coverage-$sample-.FINISHED} -vars {sample} -code {
#		cg bam2coverage $dep coverage-bwa-$sample/coverage-bwa-$sample
#	}
#	job bwa_coverage-$sample -deps $result.bam -targets sreg-$sample.tsv -vars {sample} -code {
#		cg regextract -above 1 7 $dep > $target.temp
#		file rename -force $target.temp $target
#	}
}

proc gatk_refseq_job refseq {
	upvar job_logdir job_logdir
	set nrefseq [file root $refseq].fa
	if {![file exists $nrefseq] && $refseq ne $nrefseq} {
		mklink $refseq $nrefseq
	}
	set picard [picard]
	job gatkrefseq_faidx-[file tail $nrefseq] -deps $nrefseq -targets {$nrefseq.fai} -code {
		exec samtools faidx $dep
	}
	set dict [file root $nrefseq].dict
	job gatkrefseq-[file tail $nrefseq] -deps $nrefseq -targets {$dict} -vars {nrefseq picard} -code {
		file delete $target.temp
		exec java -jar $picard/CreateSequenceDictionary.jar R= $nrefseq O= $target.temp 2>@ stderr >@ stdout
		file rename -force $target.temp $target
	}
	return $nrefseq
}

proc bam_clean_job {args} {
	oargs bam_clean_job {bamfile refseq sample
		{removeduplicates 1}
		{realign 1}
		{realignopts {}}
		{realigndeps {}}
		{clipamplicons {}}
		{cleanup 1}
		{bed {}}
		args
	} $args
	if {$bed ne ""} {
		lappend realigndeps $bed
	}
	if {[llength $args]} {
		array set opt $args
	}
	set gatk [gatk]
	set picard [picard]
	upvar job_logdir job_logdir
	if {$realign eq "srma"} {
		foreach value $realigndeps {
			lappend realignopts RANGES=$value
		}
	} else {
		foreach value $realigndeps {
			lappend realignopts -L $value
		}
	}
	set dir [file dir $bamfile]
	set file [file tail $bamfile]
	set pre [lindex [split $file -] 0]
	set root [join [lrange [split [file root $file] -] 1 end] -]
	# make gatk refseq
	set gatkrefseq [gatk_refseq_job $refseq]
	set dict [file root $gatkrefseq].dict
	# precalc skips and cleanuplist
	set skips {}
	set cleanuplist {}
	lappend cleanuplist $dir/$pre-$root.bam $dir/$pre-$root.bam.bai
	set temproot s$root
	lappend skips -skip [list $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.bai]
	if {$removeduplicates} {
		lappend cleanuplist $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.bai
		set temproot d$temproot
		lappend skips -skip [list $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.bai]
	}
	if {$realign ne "0"} {
		lappend cleanuplist $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.bai
		set temproot r$temproot
		lappend skips -skip [list $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.bai]
	}
	if {$clipamplicons ne ""} {
		lappend cleanuplist $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.bai
		set temproot c$temproot
		lappend skips -skip [list $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.bai]
	}
	# start jobs
	# sort using picard
	job bamsort-$root -deps {$bamfile} -targets {$dir/$pre-s$root.bam} \
	-vars {removeduplicates sample picard} {*}$skips -code {
		file delete $target.temp
		exec java -jar $picard/SortSam.jar	I=$dep	O=$target.temp	SO=coordinate 2>@ stderr >@ stdout
		file rename -force $target.temp $target
	#	# exec java -jar $picard/AddOrReplaceReadGroups.jar	I=$src	O=$target.temp3	RGID=$sample	RGLB=solexa-123	RGPL=illumina	RGPU=$sample RGSM=$sample 2>@ stderr >@ stdout
	#	# file delete $target.temp $target.temp2
	#	# file rename -force $target.temp3 $target
	}
	set root s$root
	if {$removeduplicates} {
		list_pop skips 0; list_pop skips 0;
		job bamremdup-$root -deps {$dir/$pre-$root.bam} -targets {$dir/$pre-d$root.bam} \
		-vars {sample picard} {*}$skips -code {
			puts "removing duplicates"
			exec java -jar $picard/MarkDuplicates.jar	I=$dep	O=$target.temp METRICS_FILE=$target.dupmetrics 2>@ stderr >@ stdout
			file rename -force $target.temp $target
		}
		set root d$root
	}
	# index intermediate result
	job bamindex-$pre-$root -deps $dir/$pre-$root.bam -targets $dir/$pre-$root.bam.bai {*}$skips -code {
		exec samtools index $dep >@ stdout 2>@ stderr
		puts "making $target"
	}
	if {$realign ne "0"} {
		list_pop skips 0; list_pop skips 0;
		# realign around indels
		set deps [list $dir/$pre-$root.bam $dir/$pre-$root.bam.bai $dict $gatkrefseq {*}$realigndeps]
		if {$realign eq "srma"} {
			set srma [srma]
			job bamrealign-$root -deps $deps -targets {$dir/$pre-r$root.bam} \
			-vars {gatkrefseq srma pre realignopts} {*}$skips -code {
				exec java -jar $srma I=$dep O=$target.temp R=$gatkrefseq {*}$realignopts 2>@ stderr >@ stdout
				catch {file rename -force $target.temp.bai $target.bai}
				catch {file delete $target.intervals}
				file rename -force $target.temp $target
			}
		} else {
			job bamrealign-$root -deps $deps -targets {$dir/$pre-r$root.bam} \
			-vars {gatkrefseq gatk pre realignopts} {*}$skips -code {
				exec java -jar $gatk -T RealignerTargetCreator -R $gatkrefseq -I $dep -o $target.intervals {*}$realignopts 2>@ stderr >@ stdout
				exec java -jar $gatk -T IndelRealigner -R $gatkrefseq -targetIntervals $target.intervals -I $dep -o $target.temp 2>@ stderr >@ stdout
				catch {file rename -force $target.temp.bai $target.bai}
				catch {file delete $target.intervals}
				file rename -force $target.temp $target
			}
		}
		set root r$root
		job bamrealign_index-$root -deps $dir/$pre-$root.bam -targets $dir/$pre-$root.bam.bai -code {
			exec samtools index $dep >@ stdout 2>@ stderr
			puts "making $target"
		}
	}
	if {$clipamplicons ne ""} {
		list_pop skips 0; list_pop skips 0;
		job bamclean_clipamplicons-$root -deps {$dir/$pre-$root.bam $clipamplicons} -targets {$dir/$pre-c$root.bam} {*}$skips -code {
			cg sam_clipamplicons $dep2 $dep $target.temp
			file rename $target.temp $target
		}
		set root c$root
		job bamclean_clipamplicons_index-$root -deps $dir/$pre-$root.bam -targets $dir/$pre-$root.bam.bai -code {
			exec samtools index $dep >@ stdout 2>@ stderr
			puts "making $target"
		}
	}
	if {$cleanup} {
		cleanup_job bamclean_remtemp-$root $cleanuplist [list $dir/$pre-$root.bam]
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
		cg select -q {$quality >= 30 && $totalcoverage >= 5 && $type ne "ins"} -f {chromosome begin end} $dep $target.temp
		file_write $target.temp2 "# regions selected from $dep: \$quality >= 30 && \$totalcoverage >= 5\n"
		cg regjoin $target.temp >> $target.temp2
		file rename -force $target.temp2 $target
		file delete $target.temp
	}
}

proc var_sam_job {bamfile refseq args} {
	upvar job_logdir job_logdir
	set pre ""
	set opts {}
	set split 0
	foreach {key value} $args {
		if {$key eq "-l"} {lappend deps $value}
		if {$key eq "-bed"} {
			lappend opts -l $value
			lappend deps $value
		} elseif {$key eq "-pre"} {
			set pre $value
		} elseif {$key eq "-split"} {
			set split $value
		} else {
			lappend opts $key $value
		}
	}
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
		-vars {refseq opts} -skip ${pre}varall-sam-$root.tsv -code {
		# bcftools -v for variant only
		exec samtools mpileup -uDS -f $refseq {*}$opts $dep 2>@ stderr | bcftools view -cg - > $target.temp 2>@ stderr
		file rename -force $target.temp $target
	}
	job ${pre}varall-sam2sft-$root -deps ${pre}varall-sam-$root.vcf -targets ${pre}varall-sam-$root.tsv -vars split -code {
		cg vcf2tsv -split $split $dep $target.temp
		file rename -force $target.temp $target
	}
	razip_job ${pre}varall-sam-$root.tsv
	job ${pre}var-sam-$root -deps ${pre}varall-sam-$root.tsv -targets {${pre}uvar-sam-$root.tsv} \
	-skip {${pre}var-sam-$root.tsv} \
	-code {
		cg select -q {
				$alt ne "." && $alleleSeq1 ne "." && $quality >= 10 && $totalcoverage > 4
				&& $zyg != "r"
			} \
			-f {
				chromosome begin end type ref alt name quality filter alleleSeq1 alleleSeq2
				{sequenced=if($quality < 30 || $totalcoverage < 5,"u",if($zyg eq "r","r","v"))}
				{zyg=if($quality < 30 || $totalcoverage < 5,"u",$zyg)}
				*
			} \
			$dep $target.temp
		file rename -force $target.temp $target
	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job ${pre}uvar-sam-$root.tsv ${pre}var-sam-$root.tsv
	# find regions
	sreg_sam_job ${pre}sreg-sam-$root ${pre}varall-sam-$root.tsv ${pre}sreg-sam-$root.tsv
	# cleanup
	job clean_${pre}var-sam-$root -deps {${pre}var-sam-$root.tsv ${pre}varall-sam-$root.tsv} -vars {pre root} -targets {} \
	-rmtargets {${pre}uvar-sam-$root.tsv ${pre}varall-sam-$root.vcf ${pre}varall-sam-$root.vcf.idx} -code {
		catch {file delete ${pre}uvar-sam-$root.tsv}
		catch {file delete ${pre}varall-sam-$root.vcf}
		catch {file delete ${pre}varall-sam-$root.vcf.idx}
	}
	cd $keeppwd
	return [file join $dir var-sam-$root.tsv]
}

proc sreg_gatk_job {job varallfile resultfile} {
	upvar job_logdir job_logdir
	job $job -deps {$varallfile} -targets {$resultfile} -code {
		cg select -q {$quality >= 30 && $totalcoverage >= 5 && $type ne "ins"} -f {chromosome begin end} $dep $target.temp
		file_write $target.temp2 "# regions selected from $dep: \$quality >= 30 && \$totalcoverage >= 5\n"
		cg regjoin $target.temp >> $target.temp2
		file rename -force $target.temp2 $target
		file delete $target.temp
	}
}

proc var_gatk_job {bamfile refseq args} {
	upvar job_logdir job_logdir
	set pre ""
	set opts {}
	set split 0
	foreach {key value} $args {
		if {$key eq "-L"} {lappend deps $value}
		if {$key eq "-bed"} {
			lappend opts -L $value
			lappend deps $value
		} elseif {$key eq "-pre"} {
			set pre $value
		} elseif {$key eq "-split"} {
			set split $value
		} else {
			lappend opts $key $value
		}
	}
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
	-targets ${pre}varall-gatk-$root.vcf -skip ${pre}varall-gatk-$root.tsv -vars {gatk opts} -code {
		exec java -d64 -Xms512m -Xmx4g -jar $gatk -T UnifiedGenotyper \
			{*}$opts -R $dep2 -I $dep -o $target.temp \
			-stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1000 \
			--annotateNDA \
			-glm SNP --output_mode EMIT_ALL_CONFIDENT_SITES 2>@ stderr >@ stdout
		file rename -force $target.temp $target
		catch {file rename -force $target.temp.idx $target.idx}
		# file delete $target.temp
	}
	job ${pre}varall-gatk2sft-$root -deps [list ${pre}varall-gatk-$root.vcf] -targets ${pre}varall-gatk-$root.tsv -vars {sample split} -code {
		cg vcf2tsv -split $split $dep $target.temp
		file rename -force $target.temp $target
	}
	razip_job ${pre}varall-gatk-$root.tsv
	# predict deletions separately, because gatk will not predict snps in a region where a deletion
	# was predicted in the varall
	job ${pre}delvar-gatk-$root -deps $deps \
	-targets ${pre}delvar-gatk-$root.vcf -skip ${pre}delvar-gatk-$root.tsv -skip ${pre}var-gatk-$root.tsv -vars {gatk opts} -code {
		exec java -d64 -Xms512m -Xmx4g -jar $gatk -T UnifiedGenotyper \
			{*}$opts -R $dep2 -I $dep -o $target.temp \
			-stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1000 \
			--annotateNDA \
			-glm INDEL 2>@ stderr >@ stdout
		file rename -force $target.temp $target
		catch {file rename -force $target.temp.idx $target.idx}
		# file delete $target.temp
	}
	job ${pre}delvar-gatk2sft-$root -deps [list ${pre}delvar-gatk-$root.vcf] -targets ${pre}delvar-gatk-$root.tsv -vars {sample split} -code {
		cg vcf2tsv -split $split $dep $target.temp
		cg select -q {$alt ne "." && $alleleSeq1 ne "." &&$quality >= 10 && $totalcoverage > 4} \
		-f {
			chromosome begin end type ref alt name quality filter alleleSeq1 alleleSeq2 
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
			chromosome begin end type ref alt name quality filter alleleSeq1 alleleSeq2 
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
	annotvar_clusters_job ${pre}uvar-gatk-$root.tsv ${pre}var-gatk-$root.tsv
	sreg_gatk_job ${pre}sreg-gatk-$root ${pre}varall-gatk-$root.tsv ${pre}sreg-gatk-$root.tsv
	## filter SNPs (according to seqanswers exome guide)
	# java -d64 -Xms512m -Xmx4g -jar $gatk -R $reference -T VariantFiltration -B:variant,VCF snp.vcf.recalibrated -o $outprefix.snp.filtered.vcf --clusterWindowSize 10 --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "DP < 5 " --filterName "LowCoverage" --filterExpression "QUAL < 30.0 " --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0 " --filterName "LowQual" --filterExpression "QD < 1.5 " --filterName "LowQD" --filterExpression "SB > -10.0 " --filterName "StrandBias"
	# cleanup
	job clean_${pre}var-gatk-$root -deps {${pre}var-gatk-$root.tsv} -vars {pre root} -targets {} \
	-rmtargets {${pre}uvar-gatk-$root.tsv ${pre}varall-gatk-$root.vcf ${pre}varall-gatk-$root.vcf.idx  ${pre}delvar-gatk-$root.vcf ${pre}delvar-gatk-$root.tsv} -code {
		catch {file delete ${pre}uvar-gatk-$root.tsv}
		catch {file delete ${pre}varall-gatk-$root.vcf}
		catch {file delete ${pre}varall-gatk-$root.vcf.idx}
		catch {file delete ${pre}delvar-gatk-$root.vcf}
		catch {file delete ${pre}delvar-gatk-$root.tsv}
	}
	cd $keeppwd
	return [file join $dir ${pre}var-gatk-$root.tsv]
}

proc multicompar_job {experiment dbdir todo args} {
	upvar job_logdir job_logdir
	set skipincomplete 1
	set split 0
	set dbfiles {}
	set addtargets 0
	foreach {key value} $args {
		if {$key eq "-skipincomplete"} {
			set skipincomplete $value
		} elseif {$key eq "-split"} {
			set split $value
		} elseif {$key eq "-dbfiles"} {
			set dbfiles $value
		} elseif {$key eq "-targetsfile"} {
			set addtargets 1
			set targetsfile $value
		} else {
			lappend opts $key $value
		}
	}
	file mkdir compar
	if {[catch {cg select -n compar/compar-$experiment.tsv} done]} {set done {}}
	set stilltodo {}
	set deps {}
	foreach sample $todo {
		set name [lindex [split $sample -] end]
		if {![inlist $done $sample]} {
			lappend stilltodo $name/var-$sample.tsv
			lappend deps \($name/sreg-$sample.tsv\) \($name/varall-$sample.tsv\)
			lappend deps \($name/coverage/coverage-*.bcol\) \($name/coverage/refScore-*.bcol\)
			lappend deps \($name/coverage/coverage-*.tsv\)
			lappend deps \($name/reg_refcons-$sample.tsv\) \($name/reg_nocall-$sample.tsv\) \($name/reg_cluster-$sample.tsv\)
		}
	}
	if {$addtargets} {
		if {[catch {cg select -n compar/compar-$experiment.tsv} header]} {set header {}}
		if {![llength $stilltodo] && [inlist $header [lindex [split $targetsfile -] end]]} {
			set addtargets 0
		} else {
			lappend deps $targetsfile
		}
	}
	if {[llength $stilltodo] || $addtargets} {
		file delete compar/compar-$experiment.tsv.temp
		if {[file exists compar/compar-$experiment.tsv]} {
			file rename -force compar/compar-$experiment.tsv compar/compar-$experiment.tsv.temp
		}
		job multicompar-$experiment -deps [list_concat $stilltodo $deps] -targets compar/compar-$experiment.tsv \
		-vars {stilltodo skipincomplete split addtargets targetsfile} -code {
			# should maybe better recheck todo here
			if {$addtargets} {
				cg multicompar -split $split -targetsfile $targetsfile $target.temp {*}$stilltodo
			} else {
				cg multicompar -split $split $target.temp {*}$stilltodo
			}
			if {$skipincomplete} {
				cg multicompar_reannot -paged 100 $target.temp skipincomplete
			} else {
				cg multicompar_reannot -paged 100 $target.temp
			}
			file rename -force $target.temp $target
		}
	}
	job annotcompar-$experiment -deps [list compar/compar-$experiment.tsv {*}$dbfiles] \
	-targets compar/annot_compar-$experiment.tsv -vars {dbdir dbfiles} -code {
		cg annotate $dep $target.temp $dbdir {*}$dbfiles
		file rename -force $target.temp $target
	}
	job indexannotcompar-$experiment \
	-deps compar/annot_compar-$experiment.tsv \
	-targets compar/annot_compar-$experiment.tsv.index/info.tsv -vars dbdir -code {
		cg index -colinfo $dep
	}
	if {[catch {cg select -n compar/sreg-$experiment.tsv} done]} {set done {}}
	set stilltodo {}
	foreach sample $todo {
		set name [lindex [split $sample -] end]
		if {![inlist $done $sample]} {
			lappend stilltodo \($name/sreg-$sample.tsv\)
		}
	}
	if {[llength $stilltodo]} {
		file delete compar/sreg-$experiment.tsv.temp
		if {[file exists compar/sreg-$experiment.tsv]} {
			file rename -force compar/sreg-$experiment.tsv compar/sreg-$experiment.tsv.temp
		}
		job sreg-$experiment -deps $stilltodo -targets compar/sreg-$experiment.tsv -vars stilltodo -code {
			cg multireg $target.temp {*}[list_remove $deps {}]
			file rename -force $target.temp $target
		}
	}
}

proc process_illumina {args} {
	set dbdir {}
	set dbfiles {}
	set realign 1
	set paired 1
	set adapterfile {}
	set conv_nextseq 0
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-realign {
				set realign $value
			}
			-dbdir {
				set dbdir $value
			}
			-split {
				set split $value
			}
			-dbfile {
				lappend dbfiles $value
			}
			-paired {
				set paired $value
			}
			-adapterfile {
				set adapterfile $value
			}
			-conv_nextseq {
				set conv_nextseq $value
			}
			default break
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	set len [llength $args]
	if {$len == 1} {
		set destdir [lindex $args 0]
	} elseif {$len == 2} {
		foreach {destdir dbdir} $args break
	} else {
		errorformat process_illumina
		exit 1
	}
	lappend dbfiles [glob $dbdir/extra/*dbnsfp*.tsv]
	lappend dbfiles [glob $dbdir/extra/var_*_evs.tsv]
	set destdir [file_absolute $destdir]
	# check projectinfo
	projectinfo $destdir dbdir split
	# start
	##in case of nextseq500 data - generate fastqs & distribute data
	if {$conv_nextseq} {
		set rundir [glob $destdir/*NS500*]
		cg_bcl2fastq $rundir fastq 4 4 4 4
		cg_process_conv_illnextseq	fastq $destdir
	}
	set refseq [glob $dbdir/genome_*.ifas]
	set samples {}
	set experiment [file tail $destdir]
	foreach dir [dirglob $destdir */fastq] {
		lappend samples [file dir $dir]
	}
	set samples [ssort -natural $samples]
	set keeppwd [pwd]
	cd $destdir
	job_logdir $destdir/log_jobs
	set todo {}
	foreach sample $samples {
		puts $sample
		set dir $destdir/$sample
		catch {file mkdir $dir}
		puts $dir
		cd $dir
		job_logdir $dir/log_jobs
		# convert existing vcfs
		set files [gzfiles var-*.vcf]
		foreach file $files {
			set target [file root [gzroot $file]].tsv
			if {![file exists $target]} {
				job vcf2tsv-$file -deps $file -targets $target -vars split -code {
					cg vcf2tsv -split $split $dep $target.temp
					file rename -force $target.temp $target
				}
				lappend todo [string range $target 4 end-4]
			}
		}
		# add existing var files to todo
		set files [gzfiles var-*.tsv]
		foreach file $files {
			set target [file root [gzroot $file]].tsv
			lappend todo [string range $target 4 end-4]
		}
		# job_logdir $dir/log_jobs
		set files [ssort -natural [glob -nocomplain fastq/*.fastq.gz fastq/*.fastq fastq/*.fq.gz fastq/*.fq]]
		if {[llength $files]} {
			# quality and adapter clipping
			set files [fastq_clipadapters_job $files -adapterfile $adapterfile -paired $paired]
			#
			# map using bwa
			map_bwa_job $refseq $files $sample $paired
		}
		# extract regions with coverage >= 5
		set cov5reg [bam2reg_job map-bwa-$sample.bam 5]
		set cov5bed [gatkworkaround_tsv2bed_job $cov5reg $refseq]
		# clean bamfile (mark duplicates, realign)
		set cleanedbam [bam_clean_job map-bwa-$sample.bam $refseq $sample -removeduplicates 1 -realign $realign -bed $cov5bed]
		# samtools variant calling on map-rdsbwa
		var_sam_job $cleanedbam $refseq -bed $cov5bed -split $split
		lappend todo sam-rdsbwa-$sample
		# gatk variant calling on map-rdsbwa
		var_gatk_job $cleanedbam $refseq -bed $cov5bed -split $split
		lappend todo gatk-rdsbwa-$sample
	}
	job_logdir $destdir/log_jobs
	cd $destdir
	set todo [list_remdup $todo]
	multicompar_job $experiment $dbdir $todo -skipincomplete 1 -split $split -dbfiles $dbfiles
	cd $keeppwd

}

proc cg_process_illumina {args} {
	set args [job_init {*}$args]
	if {[llength $args] < 1} {
		errorformat process_illumina
		exit 1
	}
	process_illumina {*}$args
	job_wait
}

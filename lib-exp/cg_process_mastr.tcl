proc tsv2bed {file bedfile chromname args} {
	if {$chromname eq ""} {
		cg select -f "$args" $file $bedfile.temp
	} else {
		cg select -f "\{chrom=\"$chromname\"\} $args" $file $bedfile.temp
	}
	set f [open $bedfile.temp]
	set o [open $bedfile w]
	set header [gets $f]
	puts $o #$header
	fcopy $f $o
	close $f; close $o
	file delete $bedfile.temp
}

proc makeminigenome {dbdir name ampliconsfile namefield {adaptorseq TGGAGAACAGTGACGATCGCAAGACTCGGCAGCATCTCCA}} {
	# using the real adaptorseq leads to false positives near the adaptor
	# better to have a bit lower coverage in these regions
	regsub -all {.} $adaptorseq N adaptorseq
	# sort and collapse regions
	set dir [file dir $ampliconsfile]
	set tail [file tail $ampliconsfile]
	cg select -s {chromosome begin end} $ampliconsfile $dir/s$tail
	cg regcollapse $dir/s$tail > $dir/reg-$name.tsv
	# the resulting reg file is used to make minigenome
	# The reg file should contain chromosome,begin,end and $namefield
	# data of the mapping is stored to reg-$name.map, for later remapping to genomic coordinates
	cg genome_seq -n 1 --namefield $namefield -m $dir/reg-$name.map -c $adaptorseq -cn $name -e $adaptorseq $dir/reg-$name.tsv $dbdir > $dir/seq-$name.fa
	# make bed files
	tsv2bed $dir/reg-$name.tsv $dir/reg-$name.bed {} chromosome begin end $namefield
	tsv2bed $dir/reg-$name.map $dir/reg-mini_$name.bed $name begin end name
	catch {
		cg select -f {chromosome begin=$primer1_end end=$primer2_begin name} $ampliconsfile $dir/inner_$tail.temp
		cg select -s - $dir/inner_$tail.temp $dir/inner_$tail.temp2
		file rename -force $dir/inner_$tail.temp2 $dir/inner_$tail
		file delete $dir/inner_$tail.temp
	}
}

proc cg_process_conv_illmastr {illsrc destdir} {
	set illsrc [file normalize $illsrc]
	set destdir [file normalize $destdir]
	file mkdir $destdir
	set keeppwd [pwd]
	cd $destdir
	# copy files from illumina dir
	if {[file exists $illsrc/Data/Intensities/BaseCalls]} {
		set illsrc $illsrc/Data/Intensities/BaseCalls
	} elseif {[file exists $illsrc/BaseCalls]} {
		set illsrc $illsrc/BaseCalls
	}
	set files [glob $illsrc/*.fastq*]
	foreach file $files {
		set sample [lindex [split [file tail $file] _] 0]
		file mkdir $destdir/$sample
		file mkdir $destdir/$sample/fastq
		exec cp -al $file $destdir/$sample/fastq
	}
	set alidir [lindex [lsort -dict [glob $illsrc/Alignment*]] end]
	set files [glob $alidir/*.bam $alidir/*.bam.bai]
	foreach file $files {
		set tail [file tail $file]
		foreach {sample barcode} [split $tail _.] break
		file mkdir $destdir/$sample
		if {[file ext $file] eq ".bam"} {set ext bam} else {set ext bam.bai}
		exec cp -al $file $destdir/$sample/map-ill-$sample.$ext
	}
	set files [glob $alidir/*.vcf]
	foreach file $files {
		set tail [file tail $file]
		foreach {sample barcode} [split $tail _.] break
		if {$sample eq "tempvariants"} continue
		file mkdir $destdir/$sample
		exec cp -al $file $destdir/$sample/var-ill-ill-$sample.vcf
	}
	cd $keeppwd
}

proc mastr_refseq_job {mastrdir dbdir} {
	set keeppwd [pwd]
	cd $mastrdir
	set mastrname [file root [file tail $mastrdir]]
	set genome [file tail $dbdir]
	job_logdir log_jobs
	set refseq seq-$mastrname.fa
	set mapfile reg-$mastrname.map
	job makeminigenome-$mastrname -deps amplicons-$mastrname.tsv \
	-targets {$refseq reg-$mastrname.bed reg-$mastrname.map reg-$mastrname.tsv} \
	-vars {dbdir mastrname} -code {
		puts stderr "makeminigenome $dbdir $mastrname $dep name"
		makeminigenome $dbdir $mastrname $dep name
	}
#	set refseq [lindex [glob $dbdir/genome_*.ifas] 0]
	# index refseq for bowtie2
	bowtie2refseq_job $refseq
	# index refseq for gatk
	gatk_refseq_job $refseq
	# index refseq for bwa
	bwarefseq_job $refseq
	cd $keeppwd
	return [list $mastrname $mastrdir/seq-$mastrname.fa $mastrdir/reg-$mastrname.map]
}

proc process_mastr_job {mastrdir destdir dbdir} {
	#
#	# make minigenome
	set mastrdir [file normalize $mastrdir]
	set destdir [file normalize $destdir]
	set dbdir [file normalize $dbdir]
	# make sure mastrdir contains everything needed
	foreach {mastrname refseq mapfile} [mastr_refseq_job $mastrdir $dbdir] break
	# start mastr analysis
	set keeppwd [pwd]
	cd $destdir
	catch {mklink $mastrdir [file tail $mastrdir]}
	set experiment [file tail $destdir]
	# which samples are there
	job_logdir $destdir/log_jobs
	set samples {}
	foreach file [dirglob $destdir */fastq] {
		lappend samples [file dir $file]
	}
	set samples [lsort -dict $samples]
	foreach sample $samples {
		puts $sample
		set name ${sample}
		set dir $destdir/$name
		catch {file mkdir $dir}
		puts $dir
		cd $dir
		job_logdir $dir/log_jobs
		# convert illumina
		job vcf2sft-ill-$name -deps {var-ill-ill-$name.vcf} -targets var-ill-ill-$name.tsv -code {
			cg vcf2sft $dep $target
		}
		# do own alignment
		set files [glob -nocomplain fastq/*.fastq.gz fastq/*.fastq]
#		# clip primers, quality
#		set adapterfile ../primers.fa
#		set outfiles {}
#		set out {}
#		foreach file $files {
#			set outfile [file root $file].clipped.fastq
#			lappend out -o $outfile
#			lappend outfiles $outfile
#		}
#		exec fastq-mcf -t 0 -k 0 -x 0 {*}$out $adapterfile {*}$files 2>@ stderr

#		set files {}
#		set tempfiles {}
#		foreach file $gzfiles {
#			if {[file extension $file] eq ".gz"} {
#				gunzip $file $file.temp
#				lappend files $file.temp
#				lappend tempfiles $file.temp
#			} else {
#				lappend files $file
#			}
#		}
		if {![llength $files]} continue
		#
		# map using bowtie2
		map_bowtie2_job $refseq $files $name {PL illumina LB solexa-123} reg_
		# clean bamfile (do not mark duplicates, realign)
		bam_clean_job reg_map-bowtie2-$name.bam $refseq $sample -removeduplicates 0
		# samtools variant calling on map-rclbowtie2
		var_sam_job reg_map-rclbowtie2-$name.bam $refseq reg_
		job remapsam-varall-$name -deps {reg_varall-sam-rclbowtie2-$name.tsv $mapfile} -targets varall-sam-rclbowtie2-$name.tsv -code {
			cg remap $dep1 $dep2 $target
		}
		job remapsam-var-$name -deps {reg_var-sam-rclbowtie2-$name.tsv $mapfile} -targets var-sam-rclbowtie2-$name.tsv -code {
			cg remap $dep1 $dep2 $target
		}
		sreg_sam_job sreg-sam-rclbowtie2-$name varall-sam-rclbowtie2-$name.tsv sreg-sam-rclbowtie2-$name.tsv
		# gatk variant calling on map-rclbowtie2
		var_gatk_job reg_map-rclbowtie2-$name.bam $refseq reg_
		job remapgatk-varall-$name -deps {reg_varall-gatk-rclbowtie2-$name.tsv $mapfile} -targets varall-gatk-rclbowtie2-$name.tsv -code {
			cg remap $dep1 $dep2 $target
		}
		job remapgatk-var-$name -deps {reg_var-gatk-rclbowtie2-$name.tsv $mapfile} -targets var-gatk-rclbowtie2-$name.tsv -code {
			cg remap $dep1 $dep2 $target
		}
		sreg_gatk_job sreg-gatk-rclbowtie2-$name varall-gatk-rclbowtie2-$name.tsv sreg-gatk-rclbowtie2-$name.tsv
	}
	job_logdir $destdir/log_jobs
	cd $destdir
	set todo {}
	foreach sample $samples {
		lappend todo sam-rclbowtie2-$sample gatk-rclbowtie2-$sample
	}
	multicompar_job $experiment $dbdir $todo
	cd $keeppwd
}

proc cg_process_mastr {args} {
	set args [job_init {*}$args]
	if {[llength $args] < 2} {
		puts "Wrong number of arguments"
		errorformat process_mastr
		exit 1
	}
	process_mastr_job {*}$args
	job_wait
}

if 0 {
cd ~/projects/mastr.data
export GATK=/home/peter/bio/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar
export PICARD=/home/peter/bio/picard-tools-1.87
cg process_conv_illmastr illuminaseqdir mastrprojectdir
cg process_mastr -d 2 mastrdesigndir mastrprojectdir /complgen/refseq/hg19
}

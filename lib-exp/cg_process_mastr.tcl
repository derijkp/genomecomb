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
	set alidir [lindex [ssort -natural [glob $illsrc/Alignment*]] end]
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

proc mastr_refseq_job {mastrdir dbdir useminigenome} {
	set keeppwd [pwd]
	cd $mastrdir
	set mastrname [file root [file tail $mastrdir]]
	set genome [file tail $dbdir]
	job_logdir log_jobs
	set refseq seq-$mastrname.fa
	set mapfile reg-$mastrname.map
	job makeminigenome-$mastrname -deps amplicons-$mastrname.tsv \
	-targets {$refseq reg-$mastrname.bed $mapfile reg-$mastrname.tsv} \
	-vars {dbdir mastrname} -code {
		puts stderr "makeminigenome $dbdir $mastrname $dep name"
		makeminigenome $dbdir $mastrname $dep name
	}
	if {!$useminigenome} {
		set refseq [glob $dbdir/genome_*.ifas]
	}
	# index refseq for bowtie2
	bowtie2refseq_job $refseq
	# index refseq for gatk
	gatk_refseq_job $refseq
	# index refseq for bwa
	bwarefseq_job $refseq
	cd $keeppwd
	return [list $mastrname $refseq $mastrdir/reg-$mastrname.map]
}

proc process_mastr_job {mastrdir destdir dbdir {useminigenome 0} {aligner bwa}} {
	#
#	# make minigenome
	set mastrdir [file normalize $mastrdir]
	set destdir [file normalize $destdir]
	set dbdir [file normalize $dbdir]
	if {$useminigenome} {set pre reg_} else {set pre {}}
	# make sure mastrdir contains everything needed
	foreach {mastrname refseq mapfile} [mastr_refseq_job $mastrdir $dbdir $useminigenome] break
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
	set samples [ssort -natural $samples]
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
		if {![llength $files]} continue
		#
		# map using $aligner
		map_${aligner}_job $refseq $files $name {PL illumina LB solexa-123} $pre
		# clean bamfile (do not mark duplicates, realign)
		set cleanbam [bam_clean_job ${pre}map-${aligner}-$name.bam $refseq $sample -removeduplicates 0]
		# samtools variant calling on map-rs${aligner}
		var_sam_job $cleanbam $refseq $pre
		if {$useminigenome} {
			job remapsam-varall-$name -deps {reg_varall-sam-rs${aligner}-$name.tsv $mapfile} -targets varall-sam-rs${aligner}-$name.tsv -code {
				cg remap $dep1 $dep2 $target
			}
			job remapsam-var-$name -deps {reg_var-sam-rs${aligner}-$name.tsv $mapfile} -targets var-sam-rs${aligner}-$name.tsv -code {
				cg remap $dep1 $dep2 $target
			}
		}
		sreg_sam_job sreg-sam-rs${aligner}-$name varall-sam-rs${aligner}-$name.tsv sreg-sam-rs${aligner}-$name.tsv
		job_razip varall-sam-rs${aligner}-$name.tsv
		if {$useminigenome} {
			# gatk variant calling on map-rs${aligner}
			var_gatk_job $cleanbam $refseq $pre
			job remapgatk-varall-$name -deps {reg_varall-gatk-rs${aligner}-$name.tsv $mapfile} -targets varall-gatk-rs${aligner}-$name.tsv -code {
				cg remap $dep1 $dep2 $target
			}
			job remapgatk-var-$name -deps {reg_var-gatk-rs${aligner}-$name.tsv $mapfile} -targets var-gatk-rs${aligner}-$name.tsv -code {
				cg remap $dep1 $dep2 $target
			}
		} else {
			# gatk variant calling on map-rs${aligner}
			var_gatk_job $cleanbam $refseq $pre -L $mastrdir/reg-$mastrname.bed
		}
		sreg_gatk_job sreg-gatk-rs${aligner}-$name varall-gatk-rs${aligner}-$name.tsv sreg-gatk-rs${aligner}-$name.tsv
		job_razip varall-gatk-rs${aligner}-$name.tsv
	}
	job_logdir $destdir/log_jobs
	cd $destdir
	set todo {}
	foreach sample $samples {
		lappend todo sam-rs${aligner}-$sample gatk-rs${aligner}-$sample
	}
	multicompar_job $experiment $dbdir $todo
	cd $keeppwd
}

proc cg_process_mastr {args} {
	set args [job_init {*}$args]
	set useminigenome 0
	set aligner bwa
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-m - --minigenome {
				set useminigenome $value
			}
			-a - --aligner {
				set aligner $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {[llength $args] < 3} {
		puts "Wrong number of arguments"
		errorformat process_mastr
		exit 1
	}
	foreach {mastrdir destdir dbdir} $args break
	process_mastr_job $mastrdir $destdir $dbdir $useminigenome $aligner
	job_wait
}

if 0 {
cd ~/projects/mastr.data
export GATK=/home/peter/bio/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar
export PICARD=/home/peter/bio/picard-tools-1.87
cg process_conv_illmastr illuminaseqdir mastrprojectdir
cg process_mastr -d 2 mastrdesigndir mastrprojectdir /complgen/refseq/hg19
}

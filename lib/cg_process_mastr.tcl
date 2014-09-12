proc tsv2bed {file bedfile args} {
	if {[llength $args]} {
		set chromname [list_shift args]
		if {$chromname eq ""} {
			file_write $bedfile.temp \#[join $args \t]\n
			cg select -sh /dev/null -f "$args" $file >> $bedfile.temp
		} else {
			file_write $bedfile.temp \#chrom\t[join $args \t]\n
			cg select -sh /dev/null -f "\{chrom=\"$chromname\"\} $args" $file >> $bedfile.temp
		}
	} else {
		set f [gzopen $file]
		set header [tsv_open $f]
		close $f
		set poss [tsv_basicfields $header 3]
		set fields [list_sub $header $poss]
		file_write $bedfile.temp \#[join $fields \t]\n
		cg select -sh /dev/null -f "$fields" $file >> $bedfile.temp		
	}
	file rename -force $bedfile.temp $bedfile
}

proc tsv2bed_job {file} {
	upvar job_logdir job_logdir
	job tsv2bed-[file tail $file] -deps $file -targets [file root $file].bed -code {
		tsv2bed $dep $target
	}
	return [file root $file].bed
}

proc make_alternative_compar_job {experiment} {
	upvar job_logdir job_logdir
	job altcompar-$experiment -deps compar/annot_compar-$experiment.tsv -targets [list compar/annot_compar_gatk-${experiment}.tsv compar/annot_compar_gatk-${experiment}_long.tsv] -code {
		set target1 [lindex $targets 0]
		set target2 [lindex $targets 1]
		##remove samtools analysis & move specific annotfields forwards		
		cg select -rf {*-sam-*} $dep $target1.temp1
		set cfields [cg select -h $target1.temp1]
		set fields [list_common {chromosome begin end type ref alt amplicons dbnsfp_SIFT_score dbnsfp_Polyphen2_HDIV_score dbnsfp_Polyphen2_HDIV_pred dbnsfp_Polyphen2_HVAR_score dbnsfp_Polyphen2_HVAR_pred snp138_name 1000gCEU refGene_impact refGene_gene refGene_descr dbnsfp_MutationTaster_score dbnsfp_MutationTaster_pred} $cfields]
		list_addnew fields  {*}$cfields
		cg select -f $fields $target1.temp1 $target1.temp2
		file delete $target1.temp1
		file rename -force $target1.temp2 $target1
		##depivot compar file
		cg long $target1 $target2.temp
		file rename -force $target2.temp $target2
	}
}

proc makeminigenome {dbdir name ampliconsfile namefield {adaptorseq TGGAGAACAGTGACGATCGCAAGACTCGGCAGCATCTCCA}} {
	# using the real adaptorseq leads to false positives near the adaptor
	# better to have a bit lower coverage in these regions
	regsub -all {.} $adaptorseq N adaptorseq
	# sort and collapse regions
	set dir [file dir $ampliconsfile]
	set tail [file tail $ampliconsfile]
	set header [cg select -h $ampliconsfile]
	if {[inlist $header upprobelen] && [inlist $header downprobelen]} {
		# clipped files
		cg select -f {chromosome {begin=$begin+$upprobelen} {end=$end - $downprobelen} name {outer_begin=$begin} {outer_end=$end} *} $ampliconsfile $dir/s$tail.temp
	} elseif {[inlist $header primer1_end] && [inlist $header primer2_begin]} {
		cg select -f {chromosome begin=$primer1_end end=$primer2_begin name outer_begin=$begin outer_end=$end *} $ampliconsfile $dir/s$tail.temp
	} else {
		set fields {chromosome begin end name}
		if {![inlist $header outer_begin]} {lappend fields outer_begin=\$begin} else {lappend fields outer_begin}
		if {![inlist $header outer_end]} {lappend fields outer_end=\$end} else {lappend fields outer_end}
		cg select -f $fields $ampliconsfile $dir/s$tail.temp
	}
	cg select -s {chromosome begin end} $dir/s$tail.temp $dir/s$tail.temp2
	file rename $dir/s$tail.temp2 $dir/s$tail
	file delete $dir/s$tail.temp
	set ampliconsfile $dir/s$tail
	cg select -f {chromosome begin=$outer_begin end=$outer_end name} $ampliconsfile $dir/reg-$name.tsv.temp
	cg select -s {chromosome begin end} $dir/reg-$name.tsv.temp $dir/reg-$name.tsv.temp2
	cg regcollapse $dir/reg-$name.tsv.temp2 > $dir/reg-$name.tsv
	file delete $dir/reg-$name.tsv.temp $dir/reg-$name.tsv.temp2
	# the resulting reg file is used to make minigenome
	# The reg file should contain chromosome,begin,end and $namefield
	# data of the mapping is stored to reg-$name.map, for later remapping to genomic coordinates
	cg genome_seq -n 1 --namefield $namefield -m $dir/reg-$name.map -c $adaptorseq -cn $name -e $adaptorseq $dir/reg-$name.tsv $dbdir > $dir/seq-$name.fa
	# make bed files
	tsv2bed $dir/reg-$name.tsv $dir/reg-$name.bed {} chromosome begin end $namefield
	tsv2bed $dir/reg-$name.map $dir/reg-mini_$name.bed $name begin end name
	set header [cg select -h $ampliconsfile]
	cg select -f {chromosome begin end name} $ampliconsfile $dir/inner_$tail.temp
	cg select -s - $dir/inner_$tail.temp $dir/inner_$tail.temp2
	file delete $dir/inner_$tail.temp
	file rename -force $dir/inner_$tail.temp2 $dir/inner_$tail
	cg regcollapse $dir/inner_$tail > $dir/reg-inner-$name.tsv
	tsv2bed $dir/reg-inner-$name.tsv $dir/reg-inner-$name.bed {} chromosome begin end $namefield
}

proc cg_process_conv_illmastr {illsrc destdir} {
	set illsrc [file_absolute $illsrc]
	set destdir [file_absolute $destdir]
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
		set sample [file tail $file]
		regsub {_[^_]+_[^_]+_[^_]+_[^_]+\.fastq.*} $sample {} sample
		regsub -all -- - $sample _ sample 
		file mkdir $destdir/$sample/ori
		file mkdir $destdir/$sample/ori/fastq
		file mkdir $destdir/$sample
		file mkdir $destdir/$sample/fastq
		exec cp -al $file $destdir/$sample/ori/fastq
		cplinked $destdir/$sample/ori/fastq $destdir/$sample/fastq
	}
	#set alidir [lindex [ssort -natural [glob $illsrc/Alignment*]] end]
	#set files [glob $alidir/*.bam $alidir/*.bam.bai]
	#foreach file $files {
		#set tail [file tail $file]
		#foreach {sample barcode} [split $tail _.] break
		#regsub -all -- - $sample _ sample
		#file mkdir $destdir/$sample
		#file mkdir $destdir/$sample/ori
		#if {[file ext $file] eq ".bam"} {set ext bam} else {set ext bam.bai}
		#exec cp -al $file $destdir/$sample/ori/map-ill-$sample.$ext
		#cplinked $destdir/$sample/ori/map-ill-$sample.$ext $destdir/$sample/map-ill-$sample.$ext
	#}
	#set files [glob $alidir/*.vcf]
	#foreach file $files {
		#set tail [file tail $file]
		#foreach {sample barcode} [split $tail _.] break
		#regsub -all -- - $sample _ sample
		#if {$sample eq "tempvariants"} continue
		#file mkdir $destdir/$sample
		#file mkdir $destdir/$sample/ori
		#exec cp -al $file $destdir/$sample/ori/var-ill-ill-$sample.vcf
		#cplinked $destdir/$sample/ori/var-ill-ill-$sample.vcf $destdir/$sample/var-ill-ill-$sample.vcf
	#}
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
	-targets {
		$refseq reg-$mastrname.bed reg-$mastrname.tsv
		inner_amplicons-$mastrname.tsv reg-inner-$mastrname.tsv reg-inner-$mastrname.bed
		reg-mini_$mastrname.bed $mapfile samplicons-$mastrname.tsv
	} \
	-vars {dbdir mastrname} -code {
		puts stderr "makeminigenome $dbdir $mastrname $dep name"
		makeminigenome $dbdir $mastrname $dep name
	}
	if {!$useminigenome} {
		set refseq [glob $dbdir/genome_*.ifas]
	}
	# index refseq for bowtie2
	# bowtie2refseq_job $refseq
	# index refseq for gatk
	gatk_refseq_job $refseq
	# index refseq for bwa
	bwarefseq_job $refseq
	cd $keeppwd
	return [list $mastrname $refseq $mastrdir/reg-$mastrname.map]
}

proc process_mastr_job {args} {
	oargs process_mastr_job {mastrdir destdir dbdir
		{useminigenome 0}
		{aligner bwa}
		{split 0}
		{paired 1}
		{cleanup 1}
		{clipamplicons {}}
	} $args
	set mastrdir [file_absolute $mastrdir]
	set destdir [file_absolute $destdir]
	set dbdir [file_absolute $dbdir]
	if {$useminigenome} {set pre reg_} else {set pre {}}
	# make sure mastrdir contains everything needed
	foreach {mastrname refseq mapfile} [mastr_refseq_job $mastrdir $dbdir $useminigenome] break
	# start mastr analysis
	set keeppwd [pwd]
	cd $destdir
	catch {mklink $mastrdir [file tail $mastrdir]}
	set experiment [file tail $destdir]
	#set additional annotation files 
	set dbfiles {}
	lappend dbfiles [glob *.mastr/reg-inner-*.tsv]
	lappend dbfiles [glob $dbdir/extra/*dbnsfp*.tsv]
	# which samples are there
	job_logdir $destdir/log_jobs
	set samples {}
	foreach file [dirglob $destdir */fastq] {
		lappend samples [file dir $file]
	}
	set samples [ssort -natural $samples]
	set todo {}
	if {$clipamplicons ne ""} {
		set resultbamprefix crs
	} else {
		set resultbamprefix rs
	}

	foreach sample $samples {
		puts $sample
		set name ${sample}
		set dir $destdir/$name
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
		# do own alignment
		set files [glob -nocomplain fastq/*.fastq.gz fastq/*.fastq]
		if {![llength $files]} continue
		# do not do any of preliminaries if end product is already there
		set bamfile ${pre}map-${aligner}-$name.bam
		set resultbamfile ${pre}map-$resultbamprefix${aligner}-$name.bam
		# quality and adapter clipping
		set files [fastq_clipadapters_job $files -skips [list -skip $bamfile -skip $resultbamfile]]
		#
		# map using $aligner
		map_${aligner}_job $refseq $files $name $paired -readgroupdata {PL illumina LB solexa-123} -pre $pre -skips [list -skip $resultbamfile]
		if {$cleanup} {
			# clean clipped files
			lappend files [file dir [lindex $files 0]]
			cleanup_job bamclean_clipped-$name $files [list $resultbamfile] [list $bamfile]
		}
		# clean bamfile (do not mark duplicates, realign)
		set cleanbam [bam_clean_job $bamfile $refseq $sample \
			-removeduplicates 0 \
			-clipamplicons $clipamplicons \
			-bed $mastrdir/reg-inner-$mastrname.bed \
			-cleanup $cleanup]
		# coverage statistics
		bam2covstats_job $cleanbam $mastrdir/reg-inner-$mastrname.tsv
		# samtools variant calling on map-rs${aligner}
		if {$useminigenome} {
			var_sam_job $cleanbam $refseq -pre $pre -split $split
			job remapsam-varall-$name -deps {reg_varall-sam-rs${aligner}-$name.tsv $mapfile} -targets varall-sam-rs${aligner}-$name.tsv -code {
				cg remap $dep1 $dep2 $target
			}
			razip_job varall-sam-rs${aligner}-$name.tsv
			job remapsam-var-$name -deps {reg_var-sam-rs${aligner}-$name.tsv $mapfile} -targets var-sam-rs${aligner}-$name.tsv -code {
				cg remap $dep1 $dep2 $target
			}
			sreg_sam_job sreg-sam-rs${aligner}-$name varall-sam-rs${aligner}-$name.tsv sreg-sam-rs${aligner}-$name.tsv
		} else {
			var_sam_job $cleanbam $refseq -pre $pre -bed $mastrdir/reg-inner-$mastrname.bed -split $split
		}
		lappend todo sam-$resultbamprefix${aligner}-$sample
		if {$useminigenome} {
			# gatk variant calling on map-rs${aligner}
			var_gatk_job $cleanbam $refseq -pre $pre -dt NONE -split $split
			job remapgatk-varall-$name -deps {reg_varall-gatk-rs${aligner}-$name.tsv $mapfile} -targets varall-gatk-rs${aligner}-$name.tsv -code {
				cg remap $dep1 $dep2 $target
			}
			razip_job varall-gatk-rs${aligner}-$name.tsv
			job remapgatk-var-$name -deps {reg_var-gatk-rs${aligner}-$name.tsv $mapfile} -targets var-gatk-rs${aligner}-$name.tsv -code {
				cg remap $dep1 $dep2 $target
			}
			sreg_gatk_job sreg-gatk-rs${aligner}-$name varall-gatk-rs${aligner}-$name.tsv sreg-gatk-rs${aligner}-$name.tsv
		} else {
			# gatk variant calling on map-rs${aligner}
			var_gatk_job $cleanbam $refseq -pre $pre -dt NONE -bed $mastrdir/reg-inner-$mastrname.bed -split $split
		}
		lappend todo gatk-$resultbamprefix${aligner}-$sample
	}
	job_logdir $destdir/log_jobs
	cd $destdir
	set todo [list_remdup $todo]
	multicompar_job $experiment $dbdir $todo -split $split -dbfiles $dbfiles
	make_alternative_compar_job $experiment 
	cd $keeppwd
}

proc cg_process_mastr {args} {
	set args [job_init {*}$args]
	set useminigenome 0
	set aligner bwa
	set cleanup 1
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-m - --minigenome {
				set useminigenome $value
			}
			-a - --aligner {
				set aligner $value
			}
			-c - --cleanup {
				set cleanup $value
			}
			-clipamplicons {
				set clipamplicons $value
			}
			-split {
				set split $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	set len [llength $args]
	if {$len == 2} {
		foreach {mastrdir destdir} $args break
	} elseif {$len == 3} {
		foreach {mastrdir destdir dbdir} $args break
	} else {
		puts stderr "Wrong number of arguments"
		errorformat process_mastr
		exit 1
	}
	set mastrname [file root [file tail $mastrdir]]
	if {![info exists clipamplicons]} {set clipamplicons $mastrdir/samplicons-$mastrname.tsv}
	# check projectinfo
	set mastrdir [file_absolute $mastrdir]
	set dbdir [file_absolute $dbdir]
	projectinfo $destdir dbdir mastrdir {split 1}
	process_mastr_job $mastrdir $destdir $dbdir -useminigenome $useminigenome \
		-aligner $aligner -split $split -cleanup $cleanup \
		-clipamplicons $clipamplicons
	job_wait
}

proc cg_process_mastrdesign {args} {
	set args [job_init {*}$args]
	# useminigenome is only important for which reference sequence is returned by mastr_refseq_job
	set useminigenome 0
	if {[llength $args] < 2} {
		puts "Wrong number of arguments, should be cg process_mastrdesign mastrdir dbdir"
		exit 1
	}
	foreach {mastrdir dbdir} $args break
	mastr_refseq_job $mastrdir $dbdir $useminigenome
	job_wait
}

if 0 {
cd ~/projects/mastr.data
export GATK=/home/peter/bio/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar
export PICARD=/home/peter/bio/picard-tools-1.87
cg process_conv_illmastr illuminaseqdir mastrprojectdir
cg process_mastr -d 2 mastrdesigndir mastrprojectdir /complgen/refseq/hg19
}

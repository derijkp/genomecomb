#####
### Auxiliary procedures
#####

proc tophat_job {sample files libtype outdir bowtie_index} {
	upvar job_logdir job_logdir
	set targets {}
	set files [ssort -natural $files]
	set target map-tophat-${sample}.bam
	job tophat-$sample -deps $files -targets {$target} -vars {libtype bowtie_index outdir}  -code {
		exec tophat  -p 8 -o $outdir --library-type $libtype $bowtie_index $dep1 $dep2 >@ stdout 2>@ stderr
		exec ln -s accepted_hits.bam $target
		}
}

proc bam_sort_rseq_job {bam} {
	upvar job_logdir job_logdir
	set pre [lindex [split $bam -] 0]
	set root [join [lrange [split [file root $bam] -] 1 end] -]
	set dir [file dir $bam]
	job bamsort-$root -deps {$bam} -targets {$dir/$pre-s$root.bam $dir/$pre-sn$root.bam} -vars {dir pre root} -code {
		file delete $target1.temp.bam
		file delete $target2.temp.bam
		samtools_sort $dep $target1.temp
		samtools_sort -n $dep $target2.temp
		file rename -force $target1.temp.bam $target1
		file rename -force $target2.temp.bam $target2
	}	
}

proc bam_index_job {bam} { 
	upvar job_logdir job_logdir
	set pre [lindex [split $bam -] 0]
	set root [join [lrange [split [file root $bam] -] 1 end] -]
	set dir [file dir $bam]
	job bamindex-$pre-$root -deps [list $dir/$pre-$root.bam] -targets [list $dir/$pre-$root.bam.bai] -code {
		exec samtools index $dep >@ stdout 2>@ stderr
		puts "making $target"
	}
}

proc htseqcount_job {bam gff order stranded mode} {
	upvar job_logdir job_logdir
	set pre count
	set root [join [lrange [split [file root $bam] -] 1 end] -]
	set dir [file dir $bam]
	job htseqcount-$bam -deps {$bam} -targets [list $dir/$pre-$root.tsv] -vars {pre root gff order stranded mode} -code {
		file delete $target.temp
		#use -q quiet option? or redirect
		exec echo "id\t$pre-$root" > $target.temp
		exec htseq-count --quiet --format=bam --order=$order --stranded=$stranded --mode=$mode $dep $gff >> $target.temp
		file rename -force $target.temp $target
	}
}

proc make_count_table_job {dir experiment samples} {
	upvar job_logdir job_logdir
	job make_count_table -deps $samples -targets {$dir/counts-$experiment.tsv}  -code {
		file delete $target.temp1
		file delete $target.temp2
		exec paste {*}$deps > $target.temp1
		##keep 1st id column, remove empty ids, remove __no_feature __ambiguous __too_low_aQual __not_aligned __alignment_not_unique
		cg select -q {$id != "" && regexp($id,"^\[^_\]+") } -f {id *-*} $target.temp1 $target.temp2
		file rename -force $target.temp2 $target
		file delete $target.temp1
	}
}
		
proc generate_fastqc {files outdir} {
	upvar job_logdir job_logdir
	catch {file mkdir $outdir}
	job generate_fastqc-$outdir -deps $files -vars {outdir}  -code {
		exec fastqc -q -o $outdir {*}$deps 2>@ stderr >@ stdout
	}
}


#####
### Main RNA seq workflow
### 1) generate fastqc (optional)
### 2) clip adapters (optional)
### 3) generate fastqc of clipped fastq  (optional)
### 4) map with tophat (bowtie2)
### 5) sort bam
### 6) index bam file
### 7) get count table using htseqcount per sample
### 8) create project count table
#####

proc process_rnaseq_job {destdir libtype bowtie_index gff fastqc adapterfile paired} {
	set destdir [file_absolute $destdir]
	set keeppwd [pwd]
	cd $destdir
	set experiment [file tail $destdir]
	if {$adapterfile == ""} {set clip 0} else {set clip 1}
	# set variables for htseqcount
		#info on htseq stranded option
		# fr-unstranded -> no ; fr-firststrand -> reverse; fr-secondstrand -> yes
		# -s <yes/no/reverse>, --stranded=<yes/no/reverse>-
  		# 		whether the data is from a strand-specific assay (default: yes)
		#		For stranded=no, a read is considered overlapping with a feature regardless of whether it is mapped to the same or 
		#		the opposite strand as the feature. For stranded=yes and single-end reads, the read has to be mapped to the same strand as the feature. 
		#		For paired-end reads, the first read has to be on the same strand and the second read on the opposite strand. For stranded=reverse, these rules are reversed.
		#
		#		Important: The default for strandedness is yes. 
		#			If your RNA-Seq data has not been made with a strand-specific protocol, this causes half of the reads to be lost. 
		#			Hence, make sure to set the option --stranded=no unless you have strand-specific data!
	switch -- $libtype  {
		fr-unstranded
			{set stranded no}
		fr-firststrand
			{set stranded reverse}
		fr-secondstrand
			{set stranded yes}
		default 
			{puts "$libtype :unknown library type"
			 exit 1}
	}
	set mode union
		#set order pos -- htseqcount gives error when file to big if sorted by pos
	set order name
	
	# which samples are there
	job_logdir $destdir/log_jobs
	set samples {}
	foreach file [dirglob $destdir */fastq] {
		lappend samples [file dir $file]
	}
	set samples [ssort -natural $samples]
	set todo {}
	foreach sample $samples {
		puts $sample
		set name ${sample}
		set dir $destdir/$name
		catch {file mkdir $dir}
		puts $dir
		cd $dir
		job_logdir $dir/log_jobs
		# do alignment
		set files [glob -nocomplain fastq/*.fastq.gz fastq/*.fastq fastq/*.fastq.bz2 fastq/*.slx.gz]
		if {![llength $files]} continue
		if $clip { 
			#clip_adapters
			set files_clipped [fastq_clipadapters_job $files -adapterfile $adapterfile -paired $paired]
			#tophat
			tophat_job $sample $files_clipped $libtype $dir $bowtie_index 
		} else {
			tophat_job $sample $files $libtype $dir $bowtie_index
		}
		if $fastqc {
			#fastqc
			generate_fastqc $files $dir/fastqc
			if $clip {
				#fastqc clipped files 
				generate_fastqc $files_clipped $dir/fastqc
			}
		}

		#sort by coordinate and by name -- latter is needed for htseqcount
		bam_sort_rseq_job map-tophat-$sample.bam
		bam_index_job map-stophat-$sample.bam
		htseqcount_job map-sntophat-$sample.bam $gff $order $stranded $mode
		lappend todo $dir/count-sntophat-$sample.tsv 
	}
	job_logdir $destdir/log_jobs
	cd $destdir
	set todo [list_remdup $todo]
	make_count_table_job $destdir $experiment $todo
	cd $keeppwd
}

proc cg_process_rnaseq {args} {
	set args [job_init {*}$args]
	#set default values for optional args
	set adapterfile {}
	set fastqc 0
	set paired 1
	#parse options
	cg_options process_rnaseq args {
		-adapterfile {
			set adapterfile $value
		}
		-fastqc {
			set fastqc $value
		}
		-paired {
			set paired $value
		}
	} {destdir libtype bowtie_index gff} 4 4
	# logfile
	set cmdline [list cg process_rnaseq]
	foreach option {
		adapterfile fastqc paired
	} {
		if {[info exists $option]} {
			lappend cmdline -$option [get $option]
		}
	}
	lappend cmdline $destdir $libtype $bowtie_index $gff
	job_logfile $destdir/process_rnaseq_[file tail $destdir] $destdir $cmdline \
		{*}[versions dbdir fastqc bowtie2 samtools gatk picard java gnusort8 lz4 os]
	#get required arguments
	process_rnaseq_job $destdir $libtype $bowtie_index $gff $fastqc $adapterfile $paired
	job_wait
}



if {0} {
#process_rnaseq -d sge NG7522 fr-firststrand /complgen3/refseq/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome /complgen3/refseq/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
#process_rnaseq -d sge -fastqc 1 -adapterfile /complgen3/project-work/rnaseq/adaptors.fa test_pipeline fr-secondstrand /complgen3/project-work/rnaseq/refseq/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome /complgen3/project-work/rnaseq/refseq/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
}
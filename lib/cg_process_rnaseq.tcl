#####
### Auxiliary procedures
#####

proc tophat_job {sample files libtype outdir bowtie_index} {
	upvar job_logdir job_logdir
	set targets {}
	set files [ssort -natural $files]
	set target map-tophat-${sample}.bam
	job tophat-$sample -deps $files -targets $target -vars {libtype bowtie_index outdir}  -code {
		exec tophat  -p 8 -o $outdir --library-type $libtype $bowtie_index $dep1 $dep2 >@ stdout 2>@ stderr
		exec ln -s accepted_hits.bam $target
		}
}

proc bam_sort_job {bam} {
	upvar job_logdir job_logdir
	set pre [lindex [split $bam -] 0]
	set root [join [lrange [split [file root $bam] -] 1 end] -]
	set dir [file dir $bam]
	job bamsort-$root -deps {$bam} -targets {$dir/$pre-s$root.bam $dir/$pre-sn$root.bam} -vars {dir pre root} -code {
		file delete $target1.temp.bam
		file delete $target2.temp.bam
		exec samtools sort $dep $target1.temp 2>@ stderr >@ stdout
		exec samtools sort -n $dep $target2.temp 2>@ stderr >@ stdout
		file rename -force $target1.temp.bam $target1
		file rename -force $target2.temp.bam $target2
	}	
}

proc bam_index_job {bam} { 
	upvar job_logdir job_logdir
	set pre [lindex [split $bam -] 0]
	set root [join [lrange [split [file root $bam] -] 1 end] -]
	set dir [file dir $bam]
	job bamindex-$pre-$root -deps $dir/$pre-$root.bam -targets $dir/$pre-$root.bam.bai -code {
		exec samtools index $dep >@ stdout 2>@ stderr
		puts "making $target"
	}
}

proc htseqcount_job {bam gff} {
	upvar job_logdir job_logdir
	set pre count
	set root [join [lrange [split [file root $bam] -] 1 end] -]
	set dir [file dir $bam]
	job htseqcount-$bam -deps $bam -targets $dir/$pre-$root.tsv -vars {gff} -code {
		file delete $target.temp
		#exec htseq-count --format=bam --order=name --stranded=yes --mode=union $dep $gff > $target.temp
		exec htseq-count --format=bam --order=name --stranded=reverse --mode=union $dep $gff > $target.temp
		file rename -force $target.temp $target
	}
}

#####
### Main RNA seq workflow
### 1) map with tophat (bowtie2)
### 2) sort bam
### 3) index bam file
### 4) get count table using htseqcount
#####

proc process_rnaseq_job {destdir libtype bowtie_index gff} {
	set destdir [file_absolute $destdir]
	set keeppwd [pwd]
	cd $destdir
	set experiment [file tail $destdir]
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
		# do own alignment
		set files [glob -nocomplain fastq/*.fastq.gz fastq/*.fastq fastq/*.fastq.bz2]
		if {![llength $files]} continue
		# tophat
		tophat_job $sample $files $libtype $dir $bowtie_index
		#sort by coordinate and by name -- latter is needed for htseqcount
		bam_sort_job map-tophat-$sample.bam
		bam_index_job map-stophat-$sample.bam
		htseqcount_job map-sntophat-$sample.bam $gff
	}
	job_logdir $destdir/log_jobs
	cd $destdir
	cd $keeppwd
}

proc process_rnaseq {args} {
	set args [job_init {*}$args]
	set pos 0
	set args [lrange $args $pos end]
	if {[llength $args] < 4} {
		puts "Wrong number of arguments"
		errorformat process_rnaseq
		exit 1
	}
	foreach {destdir libtype bowtie_index gff} $args break
	process_rnaseq_job $destdir $libtype $bowtie_index $gff
	job_wait
}


process_rnaseq -d sge NG7522 fr-firststrand /complgen3/refseq/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome /complgen3/refseq/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf




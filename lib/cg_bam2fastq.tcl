proc cg_bam2fastq {args} {
	set pos 0
	set method picard
	set fastqfile2 {}
	cg_options bam2fastq args {
		-m - -method {
			set method $value
		}
	} {bamfile fastqfile1 fastqfile2} 2 3
	set compress 0
	if {[file extension $fastqfile1] eq ".gz"} {
		set fastqfile1 [file root $fastqfile1]
		set compress 1
	}
	if {[file extension $fastqfile2] eq ".gz"} {
		set fastqfile2 [file root $fastqfile2]
		set compress 1
	}
	set fastqfile1 [file_absolute $fastqfile1]
	set fastqfile2 [file_absolute $fastqfile2]
	set destdir [file dir $fastqfile1]
	set tempbam [file root [tempfile]].bam
	# Aligning the generated fastq files may give problems/biases if the bam is sorted on position
	# Sorting based on name should avoid this
	putslog "Sorting bam file on name"
	samtools_sort -n $bamfile $tempbam
	if {$method eq "picard"} {
		putslog "Using picard to convert bam to fastq"
		if {$fastqfile2 ne ""} {
			set picard [findpicard]
			if {[catch {
				exec samtools view -hf 0x2 $tempbam | java -jar $picard/SamToFastq.jar I=/dev/stdin F=$fastqfile1.temp F2=$fastqfile2.temp VALIDATION_STRINGENCY=SILENT
			} msg] && ![regexp "done. Elapsed time:" $msg]} {
				error $msg
			}
		} else {
			picard SamToFastq I=$tempbam F=$fastqfile1 VALIDATION_STRINGENCY=SILENT >@ stdout
		}
		putslog $msg
	} elseif {$method eq "sam"} {
		putslog "Using samtools to convert bam to fastq (mate pairs are not always kept synced!)"
		# just a try; this does not always keep mate pairs synced, so not very reliable
		exec samtools view -uf64 $tempbam | samtools bam2fq - > $fastqfile1.temp
		exec samtools view -uf128 $tempbam | samtools bam2fq - > $fastqfile2.temp
	} else {
		error "unknown method \"$method\", must be picard or sam"
	}
	if {$compress} {
		exec gzip $fastqfile1.temp
		exec gzip $fastqfile2.temp
		file rename -force $fastqfile1.temp.gz $fastqfile1.gz
		file rename -force $fastqfile2.temp.gz $fastqfile2.gz
	} else {
		file rename -force $fastqfile1.temp $fastqfile1
		file rename -force $fastqfile2.temp $fastqfile2
	}
	file delete $tempbam
}

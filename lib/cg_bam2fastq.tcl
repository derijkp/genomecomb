proc cg_bam2fastq {args} {
	set pos 0
	set method biobambam
	set namesort 1
	set fastqfile2 {}
	set singlefile [scratchfile]
	set unmatchedfile [scratchfile]
	set unmatchedfile2 [scratchfile]
	cg_options bam2fastq args {
		-m - -method {
			set method $value
		}
		-namesort {
			set namesort $value
		}
		-s - -single {
			set singlefile $value
		}
		-u - -unmatched {
			set unmatchedfile $value
		}
		-u2 - -unmatched2 {
			set unmatchedfile2 $value
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
	# Aligning the generated fastq files may give problems/biases if the bam is sorted on position
	# Sorting based on name should avoid this
	if {$namesort} {
		putslog "Sorting bam file on name"
		set tempbam [file root [scratchfile]].bam
		bam_sort -sort name $bamfile $tempbam
	} else {
		set tempbam $bamfile
	}
	if {$method eq "biobambam"} {
		putslog "Using biobambam to convert bam to fastq"
		if {$fastqfile2 ne ""} {
			biobambam bamtofastq filename=$tempbam F=$fastqfile1.temp F2=$fastqfile2.temp S=$singlefile O=$unmatchedfile O2=$unmatchedfile2 collate=1 exclude=SECONDARY T=[scratchfile] gz=$compress
		} else {
			biobambam bamtofastq filename=$tempbam F=$fastqfile1.temp S=$singlefile O=$unmatchedfile collate=0 exclude=SECONDARY T=[scratchfile] gz=$compress
		}
	} elseif {$method eq "picard"} {
		putslog "Using picard to convert bam to fastq"
		if {$fastqfile2 ne ""} {
			set picard [findpicard]
			if {[catch {
				exec samtools view -hf 0x2 $tempbam | java -jar $picard/SamToFastq.jar I=/dev/stdin F=$fastqfile1.temp F2=$fastqfile2.temp VALIDATION_STRINGENCY=SILENT
			} msg] && ![regexp "done. Elapsed time:" $msg]} {
				error $msg
			}
		} else {
			picard SamToFastq I=$tempbam F=$fastqfile1.temp VALIDATION_STRINGENCY=SILENT >@ stdout
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
	if {$compress && $method ne "biobambam"} {
		exec gzip $fastqfile1.temp
		file rename -force $fastqfile1.temp.gz $fastqfile1.gz
		if {$fastqfile2 ne ""} {
			exec gzip $fastqfile2.temp
			file rename -force $fastqfile2.temp.gz $fastqfile2.gz
		}
	} elseif {$compress} {
		file rename -force $fastqfile1.temp $fastqfile1.gz
		if {$fastqfile2 ne ""} {
			file rename -force $fastqfile2.temp $fastqfile2.gz
		}
	} else {
		file rename -force $fastqfile1.temp $fastqfile1
		if {$fastqfile2 ne ""} {
			file rename -force $fastqfile2.temp $fastqfile2
		}
	}
	if {$tempbam ne $bamfile} {file delete $tempbam}
}

proc cg_bam2fastq {args} {
	set pos 0
	set method sam
	set sortmethod collate
	set namesort 1
	set fastqfile2 {}
	set singlefile [scratchfile]
	set unmatchedfile [scratchfile]
	set unmatchedfile2 [scratchfile]
	set threads 1
	cg_options bam2fastq args {
		-m - -method {
			set method $value
		}
		-sortmethod {
			set sortmethod $value
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
		-threads {
			set threads $value
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
	if {$sortmethod eq "collate" && $method ni "sam samtools"} {
		error "sortmethod collate only supported for samtools"
	}
	# Aligning the generated fastq files may give problems/biases if the bam is sorted on position
	# Sorting based on name should avoid this
	if {$namesort && !($method in "sam samtools" && $sortmethod eq "collate")} {
		putslog "Sorting bam file on name"
		set tempbam [file root [scratchfile]].bam
		cg_bam_sort -threads $threads -sort name -method $sortmethod $bamfile $tempbam
	} else {
		set tempbam $bamfile
	}
	set tempfastq1 $fastqfile1.temp[gzext $fastqfile1]
	set tempfastq2 $fastqfile2.temp[gzext $fastqfile2]
	if {$method eq "biobambam"} {
		putslog "Using biobambam to convert bam to fastq"
		if {$fastqfile2 ne ""} {
			biobambam bamtofastq filename=$tempbam F=$tempfastq1 F2=$tempfastq2 S=$singlefile O=$unmatchedfile O2=$unmatchedfile2 collate=1 exclude=SECONDARY T=[scratchfile] gz=$compress
		} else {
			biobambam bamtofastq filename=$tempbam F=$tempfastq1 S=$singlefile O=$unmatchedfile collate=0 exclude=SECONDARY T=[scratchfile] gz=$compress
		}
	} elseif {$method eq "picard"} {
		putslog "Using picard to convert bam to fastq"
		if {$fastqfile2 ne ""} {
			set picard [findpicard]
			if {[catch {
				exec samtools view --no-PG -hf 0x2 $tempbam | java -jar $picard/SamToFastq.jar I=/dev/stdin F=$tempfastq1 F2=$tempfastq2 VALIDATION_STRINGENCY=SILENT
			} msg] && ![regexp "done. Elapsed time:" $msg]} {
				error $msg
			}
		} else {
			picard SamToFastq I=$tempbam F=$tempfastq1 VALIDATION_STRINGENCY=SILENT >@ stdout
		}
		putslog $msg
	} elseif {$method in "sam samtools" && $sortmethod eq "collate"} {
		putslog "Using samtools to convert bam to fastq"
		catch_exec samtools collate -u -O $tempbam | samtools fastq -1 $tempfastq1 -2 $tempfastq2 -0 /dev/null -s $singlefile -n -N -F 0x900 -
	} elseif {$method in "sam samtools"} {
		putslog "Using samtools to convert bam to fastq"
		catch_exec samtools fastq -1 $tempfastq1 -2 $tempfastq2 -0 /dev/null -s $singlefile -n -N -F 0x900 $tempbam
	} else {
		error "unknown method \"$method\", must be picard or sam"
	}
	if {$compress && $method ni "biobambam"} {
		catch {file delete $tempfastq1.gz}
		cg_gzip $tempfastq1
		file rename -force -- $tempfastq1.gz $fastqfile1.gz
		if {$fastqfile2 ne ""} {
			catch {file delete $tempfastq2.gz}
			cg_gzip $tempfastq2
			file rename -force -- $tempfastq2.gz $fastqfile2.gz
		}
	} elseif {$compress} {
		file rename -force -- $tempfastq1 $fastqfile1.gz
		if {$fastqfile2 ne ""} {
			file rename -force -- $tempfastq2 $fastqfile2.gz
		}
	} else {
		file rename -force -- $tempfastq1 $fastqfile1
		if {$fastqfile2 ne ""} {
			file rename -force -- $tempfastq2 $fastqfile2
		}
	}
	if {$tempbam ne $bamfile} {file delete $tempbam}
}

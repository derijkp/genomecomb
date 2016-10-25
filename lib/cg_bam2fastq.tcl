proc cg_bam2fastq {args} {
	set pos 0
	set method picard
	while 1 {
		set key [lindex $args $pos]
		switch -- $key {
			-m {
				incr pos
				set method [lindex $args $pos]
			}
			-- {
				incr pos
				break
			}
			default {
				if {[string index $key 0] eq "-"} {error "unknown option \"$key\""}
				break
			}
		}
		incr pos
	}
	set args [lrange $args $pos end]
	if {[llength $args] < 2 || [llength $args] > 3} {
		errorformat bam2fastq
		exit 1
	}
	set fastqfile2 {}
	foreach {bamfile fastqfile1 fastqfile2} $args break
	set destdir [file dir $fastqfile1]
	set tempbam $destdir/[file tail $bamfile].temp_namesorted.bam
	# Aligning the generated fastq files may give problems/biases if the bam is sorted on position
	# Sorting based on name should avoid this
	exec samtools sort -n $bamfile [file root $tempbam] >@ stdout 2>@ stderr
#	if {[catch {
#		picard SortSam	I=$bamfile	O=$tempbam	SO=queryname VALIDATION_STRINGENCY=SILENT 2>@ stderr >@ stdout
#	} msg]} {
#		if {![regexp "done. Elapsed time:" $msg]} {
#			error $msg
#		}
#	}
	if {$method eq "picard"} {
		if {[catch {
			if {$fastqfile2 ne ""} {
				set picard [findpicard]
				exec samtools view -hf 0x2 Sample.sorted.bam chr21 | java -jar $picard/SamToFastq.jar I=$tempbam F=$fastqfile1 F2=$fastqfile2 VALIDATION_STRINGENCY=SILENT >@ stdout
			} else {
				picard SamToFastq I=$tempbam F=$fastqfile1 VALIDATION_STRINGENCY=SILENT >@ stdout
			}
		} msg]} {
			if {![regexp "done. Elapsed time:" $msg]} {
				error $msg
			}
		}
		puts stderr $msg
	} elseif {$method eq "sam"} {
		# just a try; this does not always keep mate pairs synced, so not very reliable
		exec samtools view -uf64 x.bam | samtools bam2fq - > $fastqfile1.temp
		exec samtools view -uf128 x.bam | samtools bam2fq - > $fastqfile2.temp
		file rename $fastqfile1.temp $fastqfile1
		file rename $fastqfile2.temp $fastqfile2
	} else {
		error "unknown method \"$method\", must be picard or sam"
	}
	file delete $tempbam
}

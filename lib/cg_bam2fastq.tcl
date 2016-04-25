proc cg_bam2fastq {bamfile fastqfile1 {fastqfile2 {}}} {
	set destdir [file dir $fastqfile1]
	set tempbam $destdir/[file tail $bamfile].temp_namesorted.bam
	# Aligning the generated fastq files may give problems/biases if the bam is sorted on position
	# Sorting based on name should avoid this
	exec samtools sort -n $bamfile [file root $tempbam] >@ stdout 2>@ stderr
#	if {[catch {
#		exec java -jar [picard]/SortSam.jar	I=$bamfile	O=$tempbam	SO=queryname VALIDATION_STRINGENCY=SILENT 2>@ stderr >@ stdout
#	} msg]} {
#		if {![regexp "done. Elapsed time:" $msg]} {
#			error $msg
#		}
#	}
	if {[catch {
		if {$fastqfile2 ne ""} {
			exec java -jar [picard]/SamToFastq.jar I=$tempbam F=$fastqfile1 F2=$fastqfile2 VALIDATION_STRINGENCY=SILENT >@ stdout
		} else {
			exec java -jar [picard]/SamToFastq.jar I=$tempbam F=$fastqfile1 VALIDATION_STRINGENCY=SILENT >@ stdout
		}
	} msg]} {
		if {![regexp "done. Elapsed time:" $msg]} {
			error $msg
		}
	}
	puts stderr $msg
	file delete $tempbam
}

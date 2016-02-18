proc cg_bam2fastq {bamfile fastqfile1 {fastqfile2 {}}} {
	set dir [file dir $destprefix]
	if {$fastqfile2 ne ""} {
		exec java -jar [picard]/SamToFastq.jar I=$bamfile F=$fastqfile1 F2=$fastqfile2 >@ stdout 2>@ stderr
	} else {
		exec java -jar [picard]/SamToFastq.jar I=$bamfile F=$fastqfile1 >@ stdout 2>@ stderr
	}
}

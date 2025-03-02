proc fastq_size {fastq maxsize} {
	incr maxsize
	if {[file ext $fastq] in ".cram .bam"} {
		catch {
			cg bam2tsv $fastq | head -$maxsize | wc -l
		} out
	} else {
		catch {
			cg fastq2tsv $fastq | head -$maxsize | wc -l
		} out
	}
	set size [expr {[lindex [split $out \n] 0]-1}]
}

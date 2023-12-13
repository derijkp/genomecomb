proc fastq_size {fastq maxsize} {
	incr maxsize
	catch {
		cg fastq2tsv $fastq | head -$maxsize | wc -l
	} out
	set size [expr {[lindex [split $out \n] 0]-1}]
}

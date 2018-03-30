proc cg_versions {args} {
	# extra commands used, but not in def list: {wget gzip gunzip zcat cat paste tail wc cp ln bigWigToBedGraph find grep chmod tar}
	# give no version: bgzip bzcat razip
	cg_options versions args {
	} {} 0 ... {
		returns the (current) versions of the given programs as a tsv file
	}
	if {![llength $args]} {
		set args {genomecomb dbdir fastqc fastq-stats fastq-mcf bwa bowtie2 samtools gatk biobambam picard plink primer3 java R gnusort8 tabix lz4 os}
	}
	puts "item\tversion"
	foreach item $args {
		puts $item\t[version $item]
	}
}

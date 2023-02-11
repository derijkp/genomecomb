proc cg_gatk_index {args} {
	cg_options makepvt args {
	} {file} 1 ... {
		index files using gatk IndexFeatureFile
	}
	if {[catch {version gatk 4.1.4}]} {
		set opt -F
	} else {
		set opt -I
	}
	foreach file [list $file {*}$args] {
		if {[file extension $file] eq ".gz" && [file extension [file root $file]] in ".vcf .gvcf"} {
			exec tabix -p vcf $file
		} else {
			gatkexec [list -XX:ParallelGCThreads=1 -Xms1g -Xmx1g] IndexFeatureFile $opt $file
		}
	}
}

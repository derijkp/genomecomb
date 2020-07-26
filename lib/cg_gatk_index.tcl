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
		if {[string range $file end-6 end] eq ".vcf.gz"} {
			exec tabix -p vcf $file
		} else {
			gatkexec [list -XX:ParallelGCThreads=1 -d64 -Xms1g -Xmx1g] IndexFeatureFile $opt $file
		}
	}
}

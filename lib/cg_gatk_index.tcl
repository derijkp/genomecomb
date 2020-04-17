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
	gatkexec [list -XX:ParallelGCThreads=1 -d64 -Xms1g -Xmx1g] IndexFeatureFile $opt $file
	foreach file $args {
		gatkexec [list -XX:ParallelGCThreads=1 -d64 -Xms1g -Xmx1g] IndexFeatureFile $opt $file
	}
}

proc cg_gatk_index {args} {
	cg_options makepvt args {
	} {file} 1 ... {
		index files using gatk IndexFeatureFile
	}
	gatkexec [list -XX:ParallelGCThreads=1 -d64 -Xms1g -Xmx1g] IndexFeatureFile -F $file
	foreach file $args {
		gatkexec [list -XX:ParallelGCThreads=1 -d64 -Xms1g -Xmx1g] IndexFeatureFile -F $file
	}
}

proc compresscmd_rz {{threads 1} {compressionlevel {}} {blocksize {}}} {
	list razip -c
}

proc compress_rz {file {destfile {}} {index 1} {keep 1} {threads 1} {compressionlevel {}} {blocksize {}} args} {
	# putsvars file destfile index keep threads compressionlevel blocksize
	set cmd [compresscmd_rz $threads $compressionlevel $blocksize]
	compress_template $file $destfile rz $cmd $index $keep
}

proc cg_razip args {
	set args [job_init {*}$args]
	cg_compress_job -method rz {*}$args
	job_wait
}

proc cg_rz args {
	set args [job_init {*}$args]
	cg_compress_job -method rz {*}$args
	job_wait
}

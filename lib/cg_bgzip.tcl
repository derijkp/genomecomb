proc compresscmd_bgz {{threads {}} {compressionlevel {}} {blocksize {}}} {
	set threads [compressionthreads $threads]
	list bgzip -@ $threads -c
}

proc compress_bgz {file {destfile {}} {index 1} {keep 1} {threads {}} {compressionlevel {}} {blocksize {}} args} {
	# putsvars file destfile index keep threads compressionlevel blocksize
	set cmd [compresscmd_bgz $threads $compressionlevel $blocksize]
	compress_template $file $destfile bgz $cmd $index $keep
}

proc cg_bgzip args {
	cg_compress_job -method bgz {*}$args
}

proc cg_bgz args {
	cg_compress_job -method bgz {*}$args
}

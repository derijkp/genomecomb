proc compresscmd_gz {{threads 1} {compressionlevel {}} {blocksize {}}} {
	set compressionlevel [compressionlevel $compressionlevel 6 1 9]
	list gzip -$compressionlevel -c
}

proc compress_gz {file {destfile {}} {index 1} {keep 1} {threads 1} {compressionlevel {}} {blocksize {}} args} {
	# putsvars file destfile index keep threads compressionlevel blocksize
	set cmd [compresscmd_gz $threads $compressionlevel $blocksize]
	compress_template $file $destfile gz $cmd $index $keep
}

proc cg_gzip args {
	cg_compress_job -method gz {*}$args
}

proc cg_gz args {
	cg_compress_job -method gz {*}$args
}

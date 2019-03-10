proc compresscmd_bz2 {{threads 1} {compressionlevel {}} {blocksize 5}} {
	list bzip2 -q -$compressionlevel -c
}

proc compress_bz2 {file {destfile {}} {index 1} {keep 1} {threads 1} {compressionlevel {}} {blocksize 5} args} {
	set cmd [compresscmd_bz2 $threads $compressionlevel $blocksize]
	compress_template $file $destfile bz2 $cmd $index $keep
}

proc cg_bzip2 args {
	set args [job_init {*}$args]
	cg_compress_job -method bz2 {*}$args
	job_wait
}

proc cg_bz2 args {
	set args [job_init {*}$args]
	cg_compress_job -method bz2 {*}$args
	job_wait
}

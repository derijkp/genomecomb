proc compresscmd_zst {{threads {}} {compressionlevel {}} {blocksize {}}} {
	set threads [compressionthreads $threads]
	set compressionlevel [compressionlevel $compressionlevel 8 1 22]
	if {$blocksize eq ""} {set blocksize 512}
	set blocksize [expr {$blocksize/1024.0}]
	list zstd-mt -k -q -$compressionlevel -b $blocksize -T $threads -c
}

proc decompresscmd_zst {{threads {}}} {
	set threads [compressionthreads $threads]
	list zstd-mt -k -q -d -T $threads -c
}

proc compress_zst {file {destfile {}} {index 1} {keep 1} {threads {}} {compressionlevel {}} {blocksize {}} args} {
	# putsvars file destfile index keep threads compressionlevel blocksize
	set cmd [compresscmd_zst $threads $compressionlevel $blocksize]
	compress_template $file $destfile zst $cmd $index $keep
}

proc cg_zst args {
	set args [job_init {*}$args]
	cg_compress_job -method zst {*}$args
	job_wait
}

proc zstindex {file} {
	exec zstdindex $file
}

proc cg_zstindex {args} {
	foreach file $args {
		exec zstdindex $file
	}
}

proc cg_zstcat {args} {
	if {![llength $args]} {
		exec zstd-mt -T 1 -k -q -d -c <@ stdin >@ stdout
	} else {
		foreach file $args {
			exec zstd-mt -T 1 -k -q -d -c $file >@ stdout
		}
	}
}

proc cg_zstra {args} {
	set start 0
	cg_options zstra args {
	} {filename start size} 1 3
	if {![info exists size]} {
		exec zstdra $filename $start >@ stdout
	} else {
		exec zstdra $filename $start $size >@ stdout
	}
}

proc cg_zstless {args} {
	if {![llength $args]} {
		set f [open "| zstd-mt -T 1 -k -q -d -c | less" w]
		fconfigure stdin -translation binary
		fconfigure $f -translation binary
		fcopy stdin $f
		close $f
	} else {
		foreach file $args {
			set f [open "| zstd-mt -T 1 -k -q -d -c $file | less" w]
			close $f
		}
	}
}

proc zstversion {{minversion {}}} {
	global zstversion
	if {![info exists zstversion]} {
		catch {exec zstd-mt -V} temp
		regexp { v([0-9.]+)} $temp temp zstversion
		if {[string index $zstversion 0] eq "v"} {set zstversion [string range $zstversion  1 end]}
	}
	if {$minversion ne ""} {
		if {[lindex [bsort [list $minversion $zstversion]] 0] ne "$minversion"} {
			error "zst version ($zstversion) smaller than $minversion"
		}
	}
	return $zstversion
}

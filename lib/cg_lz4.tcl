proc index_lz4 {file} {
	exec lz4index $file
}

proc compresscmd_lz4 {{threads 1} {compressionlevel {}} {blocksize 5}} {
	set compressionlevel [compressionlevel $compressionlevel 9 1 9]
	if {$blocksize eq ""} {
		set blocksize 5
	} elseif {$blocksize < 4} {
		error "blocksize $blocksize not supported"
	} elseif {$blocksize >= 16048} {
		set blocksize 7
	} elseif {$blocksize >= 1024} {
		set blocksize 7
	} elseif {$blocksize >= 1024} {
		set blocksize 6
	} elseif {$blocksize >= 256} {
		set blocksize 5
	} elseif {$blocksize >= 64} {
		set blocksize 4
	}
	list lz4 -q -$compressionlevel -B$blocksize -c
}

proc compress_lz4 {file {destfile {}} {index 1} {keep 1} {threads 1} {compressionlevel {}} {blocksize {}} args} {
	# putsvars file destfile index keep threads compressionlevel blocksize
	set cmd [compresscmd_lz4 $threads $compressionlevel $blocksize]
	compress_template $file $destfile lz4 $cmd $index $keep
}

proc cg_lz4 args {
	set args [job_init {*}$args]
	cg_compress_job -method lz4 {*}$args
	job_wait
}

proc cg_lz4index {args} {
	foreach file $args {
		exec lz4index $file
	}
}

proc cg_lz4cat {args} {
	if {![llength $args]} {
		exec lz4 -q -d -c <@ stdin >@ stdout
	} else {
		foreach file $args {
			exec lz4 -q -d -c $file >@ stdout
		}
	}
}

proc cg_lz4ra {args} {
	set start 0
	cg_options lz4ra args {
	} {filename start size} 1 3
	if {![info exists size]} {
		exec lz4ra $filename $start >@ stdout
	} else {
		exec lz4ra $filename $start $size >@ stdout
	}
}

proc cg_lz4less {args} {
	if {![llength $args]} {
		set f [open "| lz4 -d -c | less" w]
		fconfigure stdin -translation binary
		fconfigure $f -translation binary
		fcopy stdin $f
		close $f
	} else {
		foreach file $args {
			set f [open "| lz4 -d -c $file | less" w]
			close $f
		}
	}
}

proc lz4version {{minversion {}}} {
	global lz4version
	if {![info exists lz4version]} {
		catch {exec lz4 -h} temp
		regexp { 32-bits ([^, \n\t]+)[, \n\t]} $temp temp lz4version
		if {[string index $lz4version 0] eq "v"} {set lz4version [string range $lz4version  1 end]}
	}
	if {$minversion ne ""} {
		if {[lindex [ssort -natural [list $minversion $lz4version]] 0] ne "$minversion"} {
			error "lz4 version ($lz4version) smaller than $minversion"
		}
	}
	return $lz4version
}

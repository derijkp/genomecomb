proc compressblocksize_lz4 {} {
	return 5
}

proc index_lz4 {file} {
	exec lz4index $file
}

proc compresscmd_lz4 {{threads 1} {compressionlevel {}} {blocksize 5}} {
	list lz4c -q -$compressionlevel -B$blocksize -c
}

proc compress_lz4 {file {destfile {}} {index 1} {keep 1} {threads 1} {compressionlevel {}} {blocksize 5} args} {
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
		exec lz4c -q -d -c <@ stdin >@ stdout
	} else {
		foreach file $args {
			exec lz4c -q -d -c $file >@ stdout
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
		set f [open "| lz4c -d -c | less" w]
		fconfigure stdin -translation binary
		fconfigure $f -translation binary
		fcopy stdin $f
		close $f
	} else {
		foreach file $args {
			set f [open "| lz4c -d -c $file | less" w]
			close $f
		}
	}
}

proc lz4version {{minversion {}}} {
	global lz4version
	if {![info exists lz4version]} {
		catch {exec lz4c -h} temp
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

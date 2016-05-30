proc cg_unzip args {
	foreach file $args {
		putslog "Uncompressing $file"
		set ext [file extension $file]
		if {[inlist {.lz4} $ext]} {
			exec lz4c -d -c $file > [file root $file].temp
			file rename [file root $file].temp [file root $file]
		} elseif {[inlist {.bz2} $ext]} {
			exec bunzip2 $file
		} else {
			gunzip $file
		}
	}
}

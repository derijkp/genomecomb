proc cg_lz4 args {
	set pos 0
	set keep 0
	set compressionlevel 9
	foreach {key value} $args {
		switch -- $key {
			-k {
				set keep $value
			}
			-c {
				set compressionlevel $value
			}
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	foreach file $args {
		set ext [file extension $file]
		switch $ext {
			.gz {
				putslog "lz4 $file"
				set result [file root $file].lz4
				exec gunzip -d -c $file > $result.temp2
				exec lz4c -$compressionlevel -c $result.temp2 > $result.temp
				file delete $result.temp2
				file rename -force $result.temp $result
				if {!$keep} {file delete $file}
			}
			.rz {
				putslog "lz4 $file"
				set result [file root $file].lz4
				exec razip -d -c $file > $result.temp2
				exec lz4c -$compressionlevel -c $result.temp2 > $result.temp
				file delete $result.temp2
				file rename -force $result.temp $result
				if {!$keep} {file delete $file}
				putslog "$file already lz4"
			}
			.lz4 {
				putslog "$file already lz4"
			}
			.bz2 {
				putslog "lz4c $file"
				set result [file root $file].lz4
				exec bzcat $file > $result.temp2
				exec lz4c -$compressionlevel -c $result.temp2 > $result.temp
				file delete $result.temp2
				file rename -force $result.temp $result
				if {!$keep} {file delete $file}
			}
			default {
				putslog "lz4c $file"
				exec lz4c -$compressionlevel -c $file > $file.lz4.temp
				file rename -force $file.lz4.temp $file.lz4
				if {!$keep} {file delete $file}
			}
		}
	}
}

proc cg_lz4index {args} {
	foreach file $args {
		exec lz4index $file
	}
}

proc cg_lz4cat {args} {
	if {![llength $args]} {
		exec lz4c -d -c <@ stdin >@ stdout
	} else {
		foreach file $args {
			exec lz4c -d -c $file >@ stdout
		}
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

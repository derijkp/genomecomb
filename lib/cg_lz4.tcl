proc cg_lz4 args {
	set keep 0
	set compressionlevel 9
	set blocksize 5
	set index 0
	set outputfile {}
	cg_options lz4 args {
		-k - --keep {
			set keep $value
		}
		-c - --compressionlevel {
			set compressionlevel $value
		}
		-b - --blocksize {
			set blocksize $value
		}
		-i - --index {
			set index $value
		}
		-o - --outputfile {
			set outputfile $value
		}
	}
	if {$outputfile ne "" && [llength $args] > 1} {
		error "option -o can only be used for compressing one file"
	}
	foreach file $args {
		set ext [file extension $file]
		if {$outputfile eq ""} {
			set result [file root $file].lz4
		} else {
			set result $outputfile
		}
		switch $ext {
			.gz {
				putslog "lz4 $file"
				exec gunzip -d -c $file > $result.temp2
				exec lz4c -$compressionlevel -B$blocksize -c $result.temp2 > $result.temp
			}
			.rz {
				putslog "lz4 $file"
				exec razip -d -c $file > $result.temp2
				exec lz4c -$compressionlevel -B$blocksize -c $result.temp2 > $result.temp
			}
			.lz4 {
				putslog "$file already lz4"
				continue
			}
			.bz2 {
				putslog "lz4c $file"
				exec bzcat $file > $result.temp2
				exec lz4c -$compressionlevel -B$blocksize -c $result.temp2 > $result.temp
			}
			default {
				if {$outputfile eq ""} {
					set result $file.lz4
				} else {
					set result $outputfile
				}
				putslog "lz4c $file"
				exec lz4c -$compressionlevel -B$blocksize -c $file > $result.temp
			}
		}
		if {$index} {exec lz4index $result.temp}
		file delete $result.temp2
		if {$index} {file rename -force $result.temp.lz4i [file root $result].lz4i}
		file rename -force $result.temp $result
		if {!$keep} {file delete $file}
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

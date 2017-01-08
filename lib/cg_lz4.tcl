proc cg_lz4 args {
	set keep {}
	set compressionlevel 9
	set blocksize 5
	set index 0
	set outputfile {}
	cg_options lz4 args {
		-k - -keep - --keep {
			set keep $value
		}
		-c - -compressionlevel - --compressionlevel {
			set compressionlevel $value
		}
		-b - -blocksize - --blocksize {
			set blocksize $value
		}
		-i - -index - --index {
			set index $value
		}
		-o - -outputfile - --outputfile {
			set outputfile $value
			if {$keep eq ""} {set keep 1}
		}
	}
	if {$keep eq ""} {set keep 0}
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
		set temp [filetemp $result]
		set temp2 [filetemp $result]
		switch $ext {
			.gz {
				putslog "lz4 $file"
				set temp2 [filetemp $result]
				exec gunzip -d -c $file > $temp2
				exec lz4c -$compressionlevel -B$blocksize -c $temp2 > $temp
			}
			.rz {
				putslog "lz4 $file"
				exec razip -d -c $file > $temp2
				exec lz4c -$compressionlevel -B$blocksize -c $temp2 > $temp
			}
			.lz4 {
				putslog "$file already lz4"
				continue
			}
			.bz2 {
				putslog "lz4c $file"
				exec bzcat $file > $temp2
				exec lz4c -$compressionlevel -B$blocksize -c $temp2 > $temp
			}
			.lz4i {
				putslog "not compressin lz4 index file $file"
			}
			default {
				if {$outputfile eq ""} {
					set result $file.lz4
				} else {
					set result $outputfile
				}
				putslog "lz4c $file"
				exec lz4c -$compressionlevel -B$blocksize -c $file > $temp
			}
		}
		if {$index} {exec lz4index $temp}
		file delete $temp2
		if {$index} {file rename -force $temp.lz4i [file root $result].lz4i}
		file rename -force $temp $result
		if {!$keep} {file delete $file}
	}
}

proc cg_lz4index {args} {
	foreach file $args {
		exec lz4index $file
		if {[file extension $file] eq ".lz4"} {file rename $file.lz4i [file root $file].lz4i}
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

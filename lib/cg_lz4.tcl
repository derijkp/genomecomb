proc cg_lz4 args {
	set keep {}
	set compressionlevel 9
	set blocksize 5
	set index 0
	set outputfile {}
	cg_options lz4 args {
		-k - -keep {
			set keep $value
		}
		-c - -compressionlevel {
			set compressionlevel $value
		}
		-b - -blocksize {
			set blocksize $value
		}
		-i - -index {
			set index $value
		}
		-o - -outputfile {
			set outputfile $value
			if {$keep eq ""} {set keep 1}
		}
	}
	if {$keep eq ""} {set keep 0}
	if {$outputfile ne "" && [llength $args] > 1} {
		error "option -o can only be used for compressing one file"
	}
	if {![llength $args]} {
		exec lz4c -q -$compressionlevel -B$blocksize -c <@ stdin >@ stdout 2>@ stderr
		return
	}
	foreach file $args {
		set ext [file extension $file]
		if {$outputfile eq ""} {
			set result [file root $file].lz4
		} else {
			set result $outputfile
		}
		set temp [filetemp $result]
		switch $ext {
			.gz {
				putslog "lz4 $file"
				set temp2 [filetemp $result]
				exec gunzip -d -c $file > $temp2
				exec lz4c -q -$compressionlevel -B$blocksize -c $temp2 > $temp
			}
			.rz {
				putslog "lz4 $file"
				set temp2 [filetemp $result]
				exec razip -d -c $file > $temp2
				exec lz4c -q -$compressionlevel -B$blocksize -c $temp2 > $temp
			}
			.lz4 {
				putslog "$file already lz4"
				continue
			}
			.bz2 {
				putslog "lz4c $file"
				set temp2 [filetemp $result]
				exec bzcat $file > $temp2
				exec lz4c -q -$compressionlevel -B$blocksize -c $temp2 > $temp
			}
			.lz4i {
				putslog "not compressing lz4 index file $file"
			}
			default {
				if {$outputfile eq ""} {
					set result $file.lz4
				} else {
					set result $outputfile
				}
				putslog "lz4c $file"
				exec lz4c -q -$compressionlevel -B$blocksize -c $file > $temp
			}
		}
		if {$index} {exec lz4index $temp}
		if {$index} {file rename -force $temp.lz4i $result.lz4i}
		if {[info exists temp2]} {file delete $temp2}
		file rename -force $temp $result
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

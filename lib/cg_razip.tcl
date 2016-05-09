proc cg_razip args {
	set pos 0
	set keep 0
	foreach {key value} $args {
		switch -- $key {
			-k {
				set keep 1
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
				putslog "razip $file"
				set result [file root $file].rz
				exec gunzip -d -c $file > $result.temp2
				exec razip -c $result.temp2 > $result.temp
				file delete $result.temp2
				file rename -force $result.temp $result
				if {!$keep} {file delete $file}
			}
			.rz {
				putslog "$file already razip"
			}
			.lz4 {
				putslog "razip $file"
				set result [file root $file].rz
				exec lz4c -d $file > $result.temp2
				exec razip -c $result.temp2 > $result.temp
				file delete $result.temp2
				file rename -force $result.temp $result
				if {!$keep} {file delete $file}
			}
			.bz2 {
				putslog "razip $file"
				set result [file root $file].rz
				exec bzcat $file > $result.temp2
				exec razip -c $result.temp2 > $result.temp
				file delete $result.temp2
				file rename -force $result.temp $result
				if {!$keep} {file delete $file}
			}
			default {
				putslog "razip $file"
				exec razip -c $file > $file.rz.temp
				file rename -force $file.rz.temp $file.rz
				if {!$keep} {file delete $file}
			}
		}
	}
}

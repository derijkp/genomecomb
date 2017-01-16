proc cg_bgzip args {
	set keep 0
	cg_options bgzip args {
		-k {
			set keep 1
		}
	}
	foreach file $args {
		set ext [file extension $file]
		switch $ext {
			.gz {
				set error [catch {exec tabix $file} errormsg]
				if {[regexp {was bgzip used to compress} $errormsg]} {
					putslog "bgzip $file"
					exec gunzip -d -c $file > $file.temp2
					exec bgzip -c $file > $file.temp
					file delete $file.temp2
					if {$keep} {file rename -force $file $file.old}
					file rename -force $file.temp $file
				}
			}
			.rz {
				putslog "bgzip $file"
				set result [file root $file].gz
				exec razip -d -c $file > $result.temp2
				exec bgzip -c $result.temp2 > $result.temp
				file delete $result.temp2
				file rename -force $result.temp $result
				if {!$keep} {file delete $file}
			}
			.lz4 {
				putslog "bgzip $file"
				set result [file root $file].gz
				exec lz4c -q -d -c $file > $result.temp2
				exec bgzip -c $result.temp2 > $result.temp
				file delete $result.temp2
				file rename -force $result.temp $result
				if {!$keep} {file delete $file}
			}
			.bz2 {
				putslog "bgzip $file"
				set result [file root $file].gz
				exec bzcat $file > $result.temp2
				exec bgzip -c $result.temp2 > $result.temp
				file delete $result.temp2
				file rename -force $result.temp $result
				if {!$keep} {file delete $file}
			}
			default {
				putslog "bgzip $file"
				exec bgzip -c $file > $file.gz.temp
				file rename -force $file.gz.temp $file.gz
				if {!$keep} {file delete $file}
			}
		}
	}
}

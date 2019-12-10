proc cg_unzip args {
	set keep 0
	set force 0
	cg_options unzip args {
		-k - -keep {
			set keep $value
		}
		-f - -force {
			set force $value
		}
	} {} 1 ...
	foreach file $args {
		putslog "Uncompressing $file"
		set ext [file extension $file]
		set target [file root $file]
		set tempfile [filetemp $target]
		if {!$force && [file exists $target]} {error "uncompressed file $target already exists"}
		exec {*}[gzcat $file] $file > $tempfile
		file rename -force -- $tempfile $target
		if {!$keep} {file delete $file}
	}
}

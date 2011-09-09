#proc checksumdir {basedir o dir} {
#	set files [dirglob $basedir [file join $dir *]]
#	
#	foreach file $files {
#		if {[file isdir $file]} {
#			checksumdir $basedir $o $file
#		} else {
#			puts "Checking $file"
##			set sha256 [lindex [exec sha256sum $file] 0]
#set sha256 test
#			puts $o "$sha256 *$file"
#		}
#	}
#}

proc cg_checksum {args} {
	if {[llength $args] < 1} {
		error "format is: cg checksum dir\n or: cg checksum manifestfile manifestfile ..."
		exit 1
	}
	if {[llength $args] == 1} {
		set basedir [lindex $args 0]
		if {[file isdir $basedir]} {
			set files [glob $basedir/manifest.all $basedir/*/manifest.all $basedir/*/*/manifest.all]
			if {[llength $files] == 0} {
				set files [exec find $basedir -iname manifest.all]
			}
		} else {
			set files $args
		}
	} else {
		set files $args
	}
	set numok 0
	set numfailed 0
	foreach file $files {
		set dir [file dir $file]
		cd $dir
		if {![file exists checksum.txt]} {
			puts "Checking $file"
			set error [catch {exec sha256sum -c manifest.all > checksum.txt.temp} result]
			file rename checksum.txt.temp checksum.txt
		}
		if {![catch {exec grep FAILED checksum.txt}]} {
			incr numfailed
			puts "FAILED $file"
		} else {
			incr numok
			puts "OK $file"
		}
	}
	puts "ok: $numok"
	puts "failed: $numfailed"
}

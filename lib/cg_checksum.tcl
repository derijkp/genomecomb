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
	set check 1; set outdir {}
	set pos 0
	while 1 {
		set key [lindex $args $pos]
		switch -- $key {
			-n {
				set check 0
				incr pos
			}
			-o {
				incr pos
				set outdir [file normalize [lindex $args $pos]]
				incr pos
			}
			-- break
			default {
				break
			}
		}
	}
	set args [lrange $args $pos end]
	if {[llength $args] < 1} {
		error "format is: cg checksum dir\n or: cg ?-n? checksum manifestfile manifestfile ...\n (-c option for dryrun: check existing checksum files only)"
		exit 1
	}
	set files {}
	foreach file $args {
		set file [file normalize $file]
		if {[file isdir $file]} {
			set temp [glob $file/manifest.all $file/*/manifest.all $file/*/*/manifest.all]
			if {[llength $files] == 0} {
				set files [split [exec find $file -iname manifest.all] \n]
			}
		} else {
			lappend files $file
		}
	}
	set numok 0
	set numfailed 0
	foreach file $files {
		set dir [file dir $file]
		cd $dir
		if {$outdir eq ""} {
			set checksumfile checksum.txt
		} else {
			set checksumfile [file join $outdir [file tail $dir]-checksum.txt]
		}
		if {$check && ![file exists $checksumfile]} {
			puts "Checking $file"
			set error [catch {exec sha256sum -c manifest.all > $checksumfile.temp} result]
			if {[regexp {no properly formatted SHA256 checksum} $result]} {
				set error [catch {exec md5sum -c manifest.all > $checksumfile.temp} result]
			}
			if {![file size $checksumfile.temp]} {
				file_write $checksumfile.temp "FAILED: $result\n"
			}
			file rename -force $checksumfile.temp $checksumfile
			puts "$checksumfile written"
		}
		# grep returns error if nothing is found
		set line {}
		catch {
			set f [open $dir/checksum.txt]
			set line [gets $f]
			close $f
		}
		if {![catch {exec grep FAILED checksum.txt}] || [catch {exec grep OK checksum.txt}]} {
			incr numfailed
			puts "FAILED $dir/checksum.txt ($line)"
		} else {
			incr numok
			puts "OK $dir/checksum.txt ($line)"
		}
	}
	puts "ok: $numok"
	puts "failed: $numfailed"
}

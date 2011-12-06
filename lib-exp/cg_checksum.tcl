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
	if {[lindex $args 0] eq "-n"} {
		set check 0
		set args [lrange $args 1 end]
	} else {
		set check 1
	}
	if {[llength $args] < 1} {
		error "format is: cg checksum dir\n or: cg ?-n? checksum manifestfile manifestfile ...\n (-c option for dryrun: check existing checksum files only)"
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
		if {$check && ![file exists checksum.txt]} {
			puts "Checking $file"
			set error [catch {exec sha256sum -c manifest.all > checksum.txt.temp} result]
			if {[regexp {no properly formatted SHA256 checksum} $result]} {
				set error [catch {exec md5sum -c manifest.all > checksum.txt.temp} result]
			}
			if {![file size checksum.txt.temp]} {
				file_write checksum.txt.temp "FAILED: $result\n"
			}
			file rename -force checksum.txt.temp checksum.txt
			puts "$dir/checksum.txt written"
		}
		# grep returns error if nothing is found
		set line {}
		catch {
			set f [open $dir/checksum.txt]
			set line [gets $f]
			close $f
		}
		if {![catch {exec grep FAILED checksum.txt}] && [catch {exec grep OK checksum.txt}]} {
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

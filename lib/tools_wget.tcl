proc wgetfile {url {resultfile {}} {force 0}} {
	if {$resultfile eq ""} {
		set resultfile [file tail $url]
	}
	if {!$force && [file exists $resultfile]} {return $resultfile}
	file delete -force $resultfile.temp
	set tail [file tail $url]
	set webcache [get ::env(webcache)]
	if {$webcache ne "" && [file exists $webcache/$tail]} {
		putslog "Getting from webcache: $tail"
		if {[catch {hardlink $webcache/$tail $resultfile.temp}]} {
			file copy $webcache/$tail $resultfile.temp
		}
	} else {
		if {[catch {
			exec wget -c --tries=45 -O $resultfile.temp $url 2>@ stderr
		} errmsg]} {
			if {[file size $resultfile.temp] == 0} {
				file delete $resultfile.temp
			}
			return {}
		}
		if {![file exists $resultfile.temp]} {
			return {}
		}
		if {[regexp "No such file" $errmsg]} {
			file delete $resultfile.temp
			return {}
		}
		if {$webcache ne ""} {
			file copy $resultfile.temp $webcache/$tail
		}
	}
	file rename -force $resultfile.temp $resultfile
	return $resultfile
}

proc wgetfiles {url resultdir} {
	if {[catch {
		exec wget -c --tries=45 -P $resultdir $url 2>@ stderr
	} errmsg]} {
		puts $errmsg
		return {}
	}
}

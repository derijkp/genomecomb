proc wgetfile {url {resultfile {}} {force 0}} {
	if {$resultfile eq ""} {
		set resultfile [file tail $url]
	}
	if {!$force && [file exists $resultfile]} {return $resultfile}
	file delete -force $resultfile.temp
	set tail [file tail $url]
	if {[info exists ::env(webcache)]} {
		set webcache $::env(webcache)
	} else {
		set webcache [tempdir]/webcache
	}
	file mkdir $webcache
	regsub -all {[:/]} $url _ temp
	set webcachename $webcache/$temp
	if {$webcache ne "" && [file exists $webcachename]} {
		putslog "Getting from webcache: $tail"
		if {[catch {hardlink $webcachename $resultfile.temp}]} {
			file copy $webcachename $resultfile.temp
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
		if {$webcache ne "" && [file exists $webcache]} {
			if {[catch {hardlink $resultfile.temp $webcachename.temp}]} {
				file copy $resultfile.temp $webcachename.temp
			}
			file rename $webcachename.temp $webcachename
		}
	}
	file rename -force $resultfile.temp $resultfile
	return $resultfile
}

proc wgetfiles {url resultdir {force 0}} {
	if {!$force && [file exists $resultdir]} {return $resultdir}
	set tail [file tail $resultdir]
	if {[info exists ::env(webcache)]} {
		set webcache $::env(webcache)
	} else {
		set webcache [tempdir]/webcache
	}
	file mkdir $webcache
	set webcachename $webcache/$tail
	if {$webcache ne "" && [file exists $webcachename]} {
		putslog "Getting from webcache: $tail"
		file delete -force $resultdir.temp
		if {[catch {hardlink $webcachename $resultdir.temp}]} {
			file copy $webcachename $resultdir.temp
		}
	} else {
		if {[catch {
			exec wget -c --tries=45 -P $resultdir.temp $url 2>@ stderr
		} errmsg]} {
			return {}
		}
		if {[catch {glob $resultdir.temp/*}]} {
			return {}
		}
		if {[regexp "No such file" $errmsg]} {
			file delete $resultdir.temp
			return {}
		}
		if {$webcache ne ""} {
			if {[catch {hardlink $resultdir.temp $webcachename.temp}]} {
				file copy $resultdir.temp $webcachename.temp
			}
			file rename $webcachename.temp $webcachename
		}
	}
	file rename -force $resultdir.temp $resultdir
	return $resultdir
}


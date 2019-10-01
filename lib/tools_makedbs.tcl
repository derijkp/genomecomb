proc mkdbs_write_info {target category description args} {
	array set a $args
	if {[info exists a(build)]} {
		set build $a(build)
	} else {
		upvar build build
	}
	if {[info exists a(dbname)]} {
		set dbname $a(dbname)
	} else {
		set dbname [file root [file tail $target]]
		regsub ^reg_${build}_ $dbname {} dbname
	}
	set source [get a(source) ""]
	if {[info exists a(version)]} {
		set version $a(version)
	} else {
		set version [timestamp]
		catch {
			regexp {version\t([^\n]+)\n} [file_read [gzroot $source].info] temp version
		}
	}
	set result [subst [deindent {
		= $dbname (build $build) =
		
		== Download info ==
		dbname	$dbname
		version	$version
		license	free
		source	$source
		time	[timestamp]
		
		== Description ==
	}]]
	append result \n[deindent $description]
	append result "\n\n== Category ==\n$category"
	file_write [gzroot $target].info $result
}

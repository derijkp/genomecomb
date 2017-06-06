proc adapterfile {{adapterfile {}}} {
	if {$adapterfile eq ""} {
		return $::externdir/adaptors.fa
	} elseif {![file exists $adapterfile]} {
		error "adapterfile $adapterfile does not exists"
	} else {
		return [file_absolute $adapterfile]
	}
}

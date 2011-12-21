proc query_tool {w} {
	Classy::Entry $w
	return [varsubst w {
		$w configure -command [list %W query] -textvariable [privatevar %W tdata(query)]
	}]
}

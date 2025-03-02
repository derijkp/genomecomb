proc reghisto {regfile} {
	global cache
	unset -nocomplain a
	set v [opencgifile $regfile header]
	set num 0
	while 1 {
		incr num
		if {![expr $num%100000]} {putsprogress $num}
		set vline [cggets $v]
		if {![llength $vline]} break
		foreach {vchr vstart vend} $vline break
		set len [expr {$vend-$vstart}]
		if {![info exists a($vchr,$len)]} {
			set a($vchr,$len) 1
		} else {
			incr a($vchr,$len)
		}
	}
	close $v
	unset -nocomplain tota
	set o stdout
	set names [bsort -sortchromosome [array names a]]
	foreach name $names {
		foreach {chr len} [split $name ,] break
		puts $o $chr\t$len\t$a($name)
		if {![info exists tota($len)]} {
			set tota($len) 1
		} else {
			incr tota($len) $a($name)
		}
	}
	set names [bsort [array names tota]]
	set tot 0
	foreach len $names {
		puts $o total\t$len\t$tota($len)
		incr tot $tota($len)
	}
	puts $o all\tall\t$tot
	puts $o all\tallnonull\t[expr {$tot-[get tota(0) 0]}]
}

proc cg_reghisto {args} {
	global scriptname action
	if {[llength $args] != 1} {
		error "format is: $scriptname $action region_file\n"
	}
	foreach {region_file} $args break
	reghisto $region_file
}

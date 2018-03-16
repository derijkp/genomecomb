proc findbam {dir sample} {
	set samples {}
	set temp $sample
	while 1 {
		lappend samples $temp
		if {![regsub -- {^[^-]+-} $temp {} temp]} break
	}
	set esample $temp
	set patterns {}
	foreach s $samples {
		lappend patterns [file dir $dir]/samples/$esample/map-$s.bam
	}
	foreach s $samples {
		lappend patterns $dir/map-$s.bam
	}
	foreach s $samples {
		lappend patterns [file dir $dir]/samples/*-$esample/map-$s.bam
	}
	foreach s $samples {
		lappend patterns [file dir $dir]/map-$s.bam
	}
	lappend patterns $dir/map-*-$esample.bam
	lappend patterns [file dir $dir]/samples/$esample/map-*$esample.bam
	lindex [glob -nocomplain {*}$patterns] 0
}

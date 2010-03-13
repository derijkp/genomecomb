package require Tk

proc redraw {w pwidth pheight} {
	set width [winfo width $w]
	set height [winfo height $w]
	set scalex [expr {double($width)/$pwidth}]
	set scaley [expr {double($height)/$pheight}]
	$w scale all 0 0 $scalex $scaley
	bind $w <Configure> [list redraw $w $width $height]
}

set drawnum 0

proc draw {list {start 0} {step 0} {highlights {}} {maxmax 0} {highlightnames {}}} {
	global drawnum
	incr drawnum
	set w .d$drawnum.c
	toplevel .d$drawnum
	canvas $w -bd 0 -highlightthickness 0
	pack $w -fill both -expand yes
	tkwait visibility $w
	catch {$w delete all}
	set len [llength $list]
	if {!$len} {error "list of length 0"}
	set width [winfo width $w]
	set height [winfo height $w]
	set uheight [expr {$height - 20}]
	set max [lmath_max $list]
	set min [lmath_min $list]
	if {$maxmax && ($max > $maxmax)} {set max $maxmax}
	if {$step == 0} {
		set step [expr {$len/10}]
		if {$step == 0} {set step 1}
	}
	set scalex [expr {double($width)/$len}]
	set scaley [expr {-double($uheight)/($max)}]
	set lh [expr {round(20/-$scaley)}]
	foreach pos $highlights name $highlightnames {
		set pos [expr {$pos-$start}]
		$w create line $pos 0 $pos $max -fill orange
		$w create text $pos [expr {round($lh*0.4)}] -anchor w -text [expr {$start+$pos}]
		$w create text $pos [lindex $list $pos] -anchor sw -text $name
	}
	for {set pos 0} {$pos < $len} {incr pos $step} {
		$w create line $pos $max $pos [expr {round(-$lh*0.3)}] -fill gray
		$w create text $pos [expr {round(-$lh*0.4)}] -anchor n -text [expr {$start+$pos}]
	}
	foreach pos {50 100 500 1000 2000} {
		$w create line 0 $pos $len $pos
		$w create text 1 $pos -anchor sw -text $pos
	}
	$w create text 1 $max -anchor nw -text $max
	set coords [list_merge [list_fill $len 0 1] $list]
	$w create line $coords
	$w scale all 0 0 $scalex $scaley
#	$w move all 0 [expr {$uheight - $min * $scaley}]
	$w move all 0 $uheight
	bind $w <Configure> [list redraw $w $width $height]
	return $w
}


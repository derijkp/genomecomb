# taken from the wiki http://wiki.tcl.tk/666

proc hls2rgb {h l s} {
    # h, l and s are floats between 0.0 and 1.0, ditto for r, g and b
    # h = 0   => red
    # h = 1/3 => green
    # h = 2/3 => blue 

    set h6 [expr {($h-floor($h))*6}]
    set r [expr {  $h6 <= 3 ? 2-$h6
                            : $h6-4}]
    set g [expr {  $h6 <= 2 ? $h6
                            : $h6 <= 5 ? 4-$h6
                            : $h6-6}]
    set b [expr {  $h6 <= 1 ? -$h6
                            : $h6 <= 4 ? $h6-2
                            : 6-$h6}]
    set r [expr {$r < 0.0 ? 0.0 : $r > 1.0 ? 1.0 : double($r)}]
    set g [expr {$g < 0.0 ? 0.0 : $g > 1.0 ? 1.0 : double($g)}]
    set b [expr {$b < 0.0 ? 0.0 : $b > 1.0 ? 1.0 : double($b)}]

    set r [expr {(($r-1)*$s+1)*$l}]
    set g [expr {(($g-1)*$s+1)*$l}]
    set b [expr {(($b-1)*$s+1)*$l}]
    return [list $r $g $b]
}

proc hls2tk {h l s} {
    set rgb [hls2rgb $h $l $s]
    foreach c $rgb {
       set intc [expr {int($c * 256)}]
       if {$intc == 256} { set intc 255 }
       set c1 [format %1X $intc]
       if {[string length $c1] == 1} {set c1 "0$c1"}
       append init $c1
    }
    return #$init
}

#
# Same as distinctLabels2, but it returns a list
# of colors rather than drawing buttons
#
proc distinctColors {n} {
    set nn 1
    set hue_increment .15
    set s 1.0 ;# non-variable saturation

    set lum_steps [expr $n * $hue_increment]
    set int_lum_steps [expr int($lum_steps)]
    if {$lum_steps > $int_lum_steps} { ;# round up
        set lum_steps [expr $int_lum_steps + 1]
    }
    set lum_increment [expr .7 / $lum_steps]

    for {set l 1.0} {$l > 0.3} {set l [expr {$l - $lum_increment}]} {
        for {set h 0.0} {$h < 1.0} {set h [expr {$h + $hue_increment}]} {
            lappend rc [hls2tk $h $l $s]
            incr nn
            if {$nn > $n} { return $rc }
        }
    }
    return $rc
}

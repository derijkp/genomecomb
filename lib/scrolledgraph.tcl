#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require ClassyTk
package require Extral
package require rbc
namespace import ::rbc::*

bind all <MouseWheel> [list ::tk::MouseWheel %W %X %Y %D 0]
bind all <Shift-MouseWheel> [list ::tk::MouseWheel %W %X %Y %D 1]
if {[tk windowingsystem] eq "x11"} {
    # Support for mousewheels on Linux/Unix commonly comes through
    # mapping the wheel to the extended buttons.
    bind all <4> [list ::tk::MouseWheel %W %X %Y 120]
    bind all <5> [list ::tk::MouseWheel %W %X %Y -120]
    bind all <Shift-4> [list ::tk::MouseWheel %W %X %Y 120 1]
    bind all <Shift-5> [list ::tk::MouseWheel %W %X %Y -120 1]
}

Widget subclass scrolledgraph

scrolledgraph method init {args} {
	super init frame
	graph $object.g
	bind $object.g <2> "$object paste; break"
	bind $object.g <Control-1> "[list $object point %x %y]; break"
	scrollbar $object.scy -orient vertical
	scrollbar $object.scx -orient horizontal
	frame $object.b
	grid $object.b -sticky nwse -row 0 -columnspan 2
	grid $object.g $object.scy -row 1 -sticky nwse
	grid $object.scx -row 2 -sticky nwse
	grid columnconfigure $object 0 -weight 1
	grid rowconfigure $object 1 -weight 1
	# link graph --> scrollbars
	$object.g axis configure x -scrollcommand [list $object.scx set] -min 0 -max 1
	$object.g axis configure x -scrollcommand [list $object xset] -min 0 -max 1
	$object.g axis configure y -scrollcommand [list $object yset] -min 0 -max 1
	# link scrollbars --> graph
	$object.scx configure -command [list $object xview]
	$object.scy configure -command [list $object yview]
	$object.g grid configure -hide no
	Rbc_ZoomStack $object.g
	set style [$object gradientstyle blue]
	bind $object.g <MouseWheel> {
		if {%D < 0} {
			[winfo parent %W] zoomy 1.05
		} else {
			[winfo parent %W] zoomy 0.95
		}
		break
	}
	
	bind $object.g <Shift-MouseWheel> {
		if {%D < 0} {
			[winfo parent %W] zoomx 1.05
		} else {
			[winfo parent %W] zoomx 0.95
		}
		break
	}
	button $object.b.print -text "Print" -command [list $object print]
	pack $object.b.print -side left
	foreach type {xmin xmax ymin ymax} {
		destroy $object.b.$type
		Classy::NumEntry $object.b.$type -width 10 -label $type -orient horizontal \
			-textvariable [privatevar $object region($type)] -command [list Classy::todo $object redraw]
		pack $object.b.$type -side left
	}
	foreach type {xlog ylog} {
		checkbutton $object.b.$type  -text $type -onvalue 1 -offvalue 0 \
			-variable [privatevar $object region($type)] -command [list Classy::todo $object redraw]
		pack $object.b.$type -side left
	}
#	Classy::ProgressWidget $object.progress
#	grid $object.progress -sticky nwse -columnspan 2
#	Classy::Progress display $object.progress
	vector create ::$object.ends.x
	vector create ::$object.ends.y
	::$object.ends.x set {0 1}
	::$object.ends.y set {0 0}
	$object clear
	$object.g element create legend -xdata ::$object.ends.x -ydata ::$object.ends.y -symbol none -linewidth 0
	bindtags $object.g [list_concat [list $object.g] [list_remove [bindtags $object.g] $object.g]]
	Classy::todo $object start
}

scrolledgraph method destroy {} {
	$object clear
}

scrolledgraph chainoptions {$object.g}

scrolledgraph addoption -datafile {DataFile dataFile plot} {}

catch {PushZoom}
rename PushZoom ori.PushZoom

proc PushZoom {graph} {
	ori.PushZoom $graph
	[winfo parent $graph] _fillregion
	[winfo parent $graph] _setvars
}

catch {PopZoom}
rename PopZoom ori.PopZoom

proc PopZoom {graph} {
	ori.PopZoom $graph
	[winfo parent $graph] _fillregion
	[winfo parent $graph] _setvars
}

scrolledgraph method paste {} {
	private $object region
	set region(xrange) [::tk::GetSelection $object PRIMARY]
	Classy::todo $object redraw
}

scrolledgraph method start {} {
	private $object region
	Classy::todo $object redraw
}

scrolledgraph method clear {} {
	private $object xmin xmax ymin ymax data
	foreach el [list_remove [$object.g element names] legend] {
		$object.g element	delete $el
	}
	foreach v [vector names $object.*] {
		vector destroy ::$v
	}
if 0 {
	foreach v [info commands $object.*] {
		rename ::$v {}
	}
}
	unset -nocomplain data
	lassign {0 0 0 0} xmin xmax ymin ymax
	set data(entries) {}
	::$object.ends.x set {0 1}
	::$object.ends.y set {0 0}
	$object.scx set 0 1
	event generate $object.g <Configure>
}

scrolledgraph method gradientstyle {basecolor {min 0} {max 100} {null 0} {symbol plus}} {
	private $object colors
	if {![llength [$object.g pen names $basecolor-0]]} {
		$object.g pen create $basecolor-0 -fill $basecolor -outline $basecolor -outlinewidth 0 -pixels 0 -symbol $symbol
	}
	if {$basecolor eq "gray"} {
		set num 1
		foreach {gnum} {
			90 80 70 60 50 40 30 20 10 0
		} {
			set color gray-$num
			set colorcode gray$gnum
			if {![llength [$object.g pen names $color]]} {
				$object.g pen create $color -fill $colorcode -outline $colorcode -outlinewidth 1 -pixels 2 -symbol $symbol
			} else {
				# $object.g pen configure $color -fill $colorcode -outline $colorcode -outlinewidth 1 -pixels 2 -symbol $symbol
			}
			incr num
		}
	} else {
		set num 1
		foreach {factor} {
			0.4 0.3 0.2 0.1 0.0 -0.1 -0.2 -0.3 -0.4 -0.5
		} {
			set color $basecolor-$num
			set colorcode [gradient $basecolor $factor]
			if {![llength [$object.g pen names $color]]} {
				$object.g pen create $color -fill $colorcode -outline $colorcode -outlinewidth 1 -pixels 2 -symbol $symbol
			} else {
				# $object.g pen configure $color -fill $colorcode -outline $colorcode -outlinewidth 1 -pixels 2 -symbol $symbol
			}
			incr num
		}
	}
	set style {}
	set step [expr {($max-$min)/10.0}]
	set num 1
	set cur $min
	for {set num 1} {$num <= 10} {incr num} {
		set next [expr {$cur+$step}]
		lappend style [list $basecolor-$num $cur $next]
		set cur $next
		incr num
	}
	lset style end end [expr {$max+1}]
	lappend style [list $basecolor-0 $null $null]
	return $style
}

scrolledgraph method _configureevent {} {
	private $object data region
	event generate $object.g <Configure>
	set region(status) ""
	set region(cancel) 0
}

scrolledgraph method getlabel {win num} {
putsvars win num
	private $object labels
	if {[isint $num]} {
		return [lindex $labels $num]
	} else {
		return ""
	}
}

scrolledgraph method add {table {xtitle {}} {ytitle {}}} {
	private $object xv yv wv data colors labels region
	set showregion 0
	# we will add to existing (data(entries) already there)
	set vnum [llength $data(entries)]
	set header [list_shift table]
	set table [ssort -natural -index 0 $table]
	set xs [list_subindex $table 0]
	vector create ::$object.x
	$object.g axis configure y -title $ytitle
	if {[catch {::$object.x set $xs}]} {
		set labels $xs
		set xs [list_fill [llength $xs] 0 1]
		::$object.x set $xs
		$object.g axis configure x -command [list $object getlabel] -title $xtitle
		# -subdivisions 0 -stepsize 1
	} else {
		unset -nocomplain labels
		$object.g axis configure x -command {} -title $xtitle
	}
#	set amin [min [lmath_min $xs] [::$object.ends.x index 0]]
#	set amax [max [lmath_max $xs] [::$object.ends.x index 1]]
	set amin [lmath_min [list_lremove $xs {Inf -Inf NaN "" "?" "-"}]]
	set amax [lmath_max [list_lremove $xs {Inf -Inf NaN "" "?" "-"}]]
	if {[info exists labels]} {
		set amin [expr {$amin - 1}]
		set amax [expr {$amax + 1}]
	}
	if {$region(xmax) == Inf} {
		set region(xmax) [::$object.ends.x index 1]
	}
	::$object.ends.x set [list $amin $amax]
	set pos 0
	set symbols {cross circle square diamond plus splus scross triangle}
	# set dcolors {blue red gray orange yellow green violet}
	set dcolors [distinctColors [expr {[llength $header]-1}]]
	set amin {}
	set amax {}
	foreach field [lrange $header 1 end] {
		incr pos
		set name $field
		foreach field {y w i} {
			vector create ::$object.$vnum.$field
		}
		set ys [list_subindex $table $pos]
		::$object.$vnum.y set $ys
		set amin [min [lmath_min [list_lremove $ys {Inf -Inf NaN "" "?" "-"}]] $amin]
		set amax [max [lmath_max [list_lremove $ys {Inf -Inf NaN "" "?" "-"}]] $amax]
		set basecolor [lindex $dcolors $vnum]
		if {$basecolor eq ""} {set basecolor blue}
		lappend data(entries) $name
		set data($name,vnum) $vnum
		set data($name,lstart) 0
		set data($name,lend) 0
		set data($name,header) $header
		set xv ::$object.x
		set yv ::$object.$vnum.y
		set wv ::$object.$vnum.w
		$object.g element create e$name -label $name -xdata $xv -ydata $yv -symbol none
		# $object.g element configure e$name -weights ::$wv
		set style [$object gradientstyle $basecolor $amin $amax 0]
		set color $basecolor
		# $object.g element configure e$name -linewidth 1 -fill $color -outlinewidth 1 -pixels 1 -symbol circle
		$object.g element configure e$name -linewidth 1 -color $basecolor -outline $basecolor -outlinewidth 1 -pixels 2 -symbol plus \
			-styles $style
		if {$showregion} {
			$object.g element configure e$name -linewidth 10 -trace decreasing
		}
		set data($name,color) $color
		$object.g legend bind $name <1> [list $object elconf $name]
		incr vnum
	}
	if {![isdouble $amin]} {set amin 0}
	if {![isdouble $amax]} {set amax 0}
	::$object.ends.y set [list $amin $amax]
	set region(xmin) [::$object.ends.x index 0]
	set region(xmax) [::$object.ends.x index 1]
	set region(ymin) [::$object.ends.y index 0]
	set region(ymax) [::$object.ends.y index 1]
	Classy::todo $object redraw
}

scrolledgraph method addscatter {table xcol ycol datacols} {
	private $object xv yv wv data colors labels region
	# we will add to existing (data(entries) already there)
	set vnum [llength $data(entries)]
	set header [list_shift table]
	set xpos [lsearch $header $xcol]
	set ypos [lsearch $header $ycol]
	# set table [ssort -natural -index $xpos $table]
	# x
	set xs [list_subindex $table $xpos]
	vector create ::$object.x
	if {[catch {::$object.x set $xs}]} {
		error "x column may only contain numbers"
	}
	set amin [lmath_min [list_lremove $xs {Inf -Inf NaN "" "?" "-"}]]
	set amax [lmath_max [list_lremove $xs {Inf -Inf NaN "" "?" "-"}]]
	if {![isdouble $amin]} {set amin 0}
	if {![isdouble $amax]} {set amax 0}
	::$object.ends.x set [list $amin $amax]
	$object.g axis configure x -title $xcol
	# y
	set ys [list_subindex $table $ypos]
	vector create ::$object.y
	if {[catch {::$object.y set $ys}]} {
		error "y column may only contain numbers"
	}
	set amin [lmath_min [list_lremove $ys {Inf -Inf NaN "" "?" "-"}]]
	set amax [lmath_max [list_lremove $ys {Inf -Inf NaN "" "?" "-"}]]
	if {![isdouble $amin]} {set amin 0}
	if {![isdouble $amax]} {set amax 0}
	::$object.ends.y set [list $amin $amax]
	$object.g axis configure y -title $ycol
	set symbols {cross circle square diamond plus splus scross triangle}
	set pos 0
	if {[llength $datacols] <= 7} {
		set dcolors {blue red gray orange yellow green violet}
	} else {
		set dcolors [distinctColors [llength $datacols]]
	}
	set amin {}
	set amax {}
	if {[llength $datacols]} {
		foreach field $datacols {
			set datapos [lsearch $header $field]
			incr pos
			set name $field
			vector create ::$object.$vnum.w
			set ws [list_subindex $table $datapos]
			set min [lmath_min [list_lremove $ws {Inf -Inf NaN "" "?" "-"}]]
			set max [lmath_max [list_lremove $ws {Inf -Inf NaN "" "?" "-"}]]
			::$object.$vnum.w set $ws
			set basecolor [lindex $dcolors $vnum]
			if {$basecolor eq ""} {set basecolor blue}
			lappend data(entries) $name
			set data($name,vnum) $vnum
			set data($name,lstart) 0
			set data($name,lend) 0
			set data($name,header) $header
			set xv ::$object.x
			set yv ::$object.y
			set wv ::$object.$vnum.w
			$object.g element create e$name -xdata $xv -ydata $yv -weight $wv -symbol none
			# $object.g element configure e$name -weights ::$wv
			set style [$object gradientstyle $basecolor $min $max 0]
			set color $basecolor
			# $object.g element configure e$name -linewidth 1 -fill $color -outlinewidth 1 -pixels 1 -symbol circle
			$object.g element configure e$name -linewidth 0 -color $basecolor -outline $basecolor -outlinewidth 1 -pixels 2 -symbol plus \
				-styles $style
			set data($name,color) $color
			$object.g legend bind $name <1> [list $object elconf $name]
			incr vnum
		}
	} else {
		incr pos
		set name x
		set data($name,vnum) $vnum
		set data($name,lstart) 0
		set data($name,lend) 0
		set data($name,header) $header
		set xv ::$object.x
		set yv ::$object.y
		set wv ::$object.$vnum.w
		$object.g element create e$name -xdata $xv -ydata $yv -symbol none
		set color black
		$object.g element configure e$name -linewidth 0 -color $color -outline $color -outlinewidth 1 -pixels 2 -symbol plus
		set data($name,color) $color
		$object.g legend bind $name <1> [list $object elconf $name]
		incr vnum
	}
	set region(xmin) [::$object.ends.x index 0]
	set region(xmax) [::$object.ends.x index 1]
	set region(ymin) [::$object.ends.y index 0]
	set region(ymax) [::$object.ends.y index 1]
	set region(xlog) [get region(xlog) 0]
	set region(ylog) [get region(ylog) 0]
	Classy::todo $object redraw
}

proc ::tk::MouseWheel {wFired X Y D {shifted 0}} {
    # Set event to check based on call
    set evt "<[expr {$shifted?{Shift-}:{}}]MouseWheel>"
    # do not double-fire in case the class already has a binding
    if {[bind [winfo class $wFired] $evt] ne ""} { return }
    # obtain the window the mouse is over
    set w [winfo containing $X $Y]
    # if we are outside the app, try and scroll the focus widget
    if {![winfo exists $w]} { catch {set w [focus]} }
    if {[winfo exists $w]} {
        if {[bind $w $evt] ne ""} {
            # Awkward ... this widget has a MouseWheel binding, but to
            # trigger successfully in it, we must give it focus.
            catch {focus} old
            if {$w ne $old} { focus $w }
            event generate $w $evt -rootx $X -rooty $Y -delta $D
            if {$w ne $old} { focus $old }
            return
        }
        # aqua and x11/win32 have different delta handling
        if {[tk windowingsystem] ne "aqua"} {
            set delta [expr {- ($D / 30)}]
        } else {
            set delta [expr {- ($D)}]
        }
        # scrollbars have different call conventions
        if {[string match "*Scrollbar" [winfo class $w]]} {
            catch {tk::ScrollByUnits $w \
                       [string index [$w cget -orient] 0] $delta}
        } else {
            set cmd [list $w [expr {$shifted ? "xview" : "yview"}] \
                         scroll $delta units]
            # Walking up to find the proper widget handles cases like
            # embedded widgets in a canvas
            while {[catch $cmd] && [winfo toplevel $w] ne $w} {
                set w [winfo parent $w]
            }
        }
    }
}

scrolledgraph method redraw {args} {
puts "----------redraw $object----------"
	private $object region
	set w $object.g
	Classy::canceltodo $object redraw
	# $object _xrange
	if {$region(xmax) <= $region(xmin)} {set region(xmax) [expr {$region(xmin)+1}]}
	if {$region(ymax) <= $region(ymin)} {set region(ymax) [expr {$region(ymin)+1}]}
	$w axis configure x -min $region(xmin) -max $region(xmax) -logscale $region(xlog)
	$w axis configure y -min $region(ymin) -max $region(ymax) -logscale $region(ylog)
	if {$region(xlog)} {
		$object.scx configure -command {}
		$object.g axis configure x -scrollcommand {}
		grid forget $object.scx
	} else {
		$object.scx configure -command [list $object xview]
		$object.g axis configure x -scrollcommand [list $object xset]
		grid $object.scx -row 2 -sticky nwse
	}
	if {$region(ylog)} {
		$object.scy configure -command {}
		$object.g axis configure y -scrollcommand {}
		grid forget $object.scy
	} else {
		$object.scy configure -command [list $object yview]
		$object.g axis configure y -scrollcommand [list $object yset]
		grid $object.scy -row 1 -column 2 -sticky nwse
	}
	Classy::todo $object _configureevent
	# Classy::todo $object reload
}

scrolledgraph method _setvars {} {
	private $object region
	set w $object.g
	set xmin [expr {round([$w axis cget x -min])}]
	set xmax [expr {round([$w axis cget x -max])}]
	if {[isdouble $xmin] && [isdouble $xmax]} {
		set region(xmin) $xmin
		set region(xmax) $xmax
	}
	set ymin [$w axis cget y -min]
	set ymax [$w axis cget y -max]
	catch {set ymin [format %.0f $ymin]}
	catch {set ymax [format %.0f $ymax]}
	set region(ymin) $ymin
	set region(ymax) $ymax
	Classy::todo $object _configureevent
}

scrolledgraph method _fillregion {args} {
	private $object region
	set w $object.g
	foreach type {xmin xmax ymin ymax} {
		set val [$w axis cget [string index $type 0] -[string range $type 1 end]]
		if {[isdouble $val]} {
			set val [format %.1f $val]
		}
		set region($type) $val
	}
}

scrolledgraph method zoomx {{factor 2}} {
	private $object xmin xmax
	set w $object.g
	set min [$w axis cget x -min]
	set max [$w axis cget x -max]
	if {$min eq ""} {set min $xmin}
	if {$max eq ""} {set max $xmax}
	set newwidth [expr {$factor * ($max-$min)}]
	set extra [expr {($newwidth - ($max-$min))/2}]
	set min [expr {$min-$extra}]
	set max [expr {$max+$extra}]
#	if {$min < $xmin} {
#		set min $xmin
#		set max [expr {$min+$newwidth}]
#	}
#	if {$max > $xmax} {
#		set max $xmax
#	}
	$w axis configure x -min $min
	$w axis configure x -max $max
	set range(xmin) $min
	set range(xmax) $max
	Classy::todo $object _configureevent
}

scrolledgraph method zoomy {{factor 2}} {
	private $object ymin ymax range
	set w $object.g
	set min [$w axis cget y -min]
	set max [$w axis cget y -max]
	if {$min eq ""} {set min $ymin}
	if {$max eq ""} {set max $ymax}
	set newheight [expr {$factor * ($max-$min)}]
	set extra [expr {($newheight - ($max-$min))/2}]
	set min [expr {$min-$extra}]
	set max [expr {$max+$extra}]
#	if {$min < $ymin} {
#		set min $ymin
#		set max [expr {$min+$newheight}]
#	}
#	if {$max > $ymax} {
#		set max $ymax
#	}
	$w axis configure y -min $min
	$w axis configure y -max $max
	set range(ymin) $min
	set range(ymax) $max
	Classy::todo $object _configureevent
}

proc gradient {rgb factor {window .}} {
	foreach {r g b} [winfo rgb $window $rgb] {break}
	### Figure out color depth and number of bytes to use in
	### the final result.
	if {($r > 255) || ($g > 255) || ($b > 255)} {
	    set max 65535
	    set len 4
	} else {
	    set max 255
	    set len 2
	}
	### Compute new red value by incrementing the existing
	### value by a value that gets it closer to either 0 (black)
	### or $max (white)
	set range [expr {$factor >= 0.0 ? $max - $r : $r}]
	set increment [expr {int($range * $factor)}]
	incr r $increment
	### Compute a new green value in a similar fashion
	set range [expr {$factor >= 0.0 ? $max - $g : $g}]
	set increment [expr {int($range * $factor)}]
	incr g $increment
	### Compute a new blue value in a similar fashion
	set range [expr {$factor >= 0.0 ? $max - $b : $b}]
	set increment [expr {int($range * $factor)}]
	incr b $increment
	### Format the new rgb string
	set rgb \
	    [format "#%.${len}X%.${len}X%.${len}X" \
	         [expr {($r>$max)?$max:(($r<0)?0:$r)}] \
	         [expr {($g>$max)?$max:(($g<0)?0:$g)}] \
	         [expr {($b>$max)?$max:(($b<0)?0:$b)}]]
	### Return the new rgb string
	return $rgb
}

scrolledgraph method shift {name {shift 1}} {
	private $object data
	if {![info exists data($name,shift)]} {
		set data($name,shift) 0
	}
	set size [expr {$shift-$data($name,shift)}]
	set yv [$object.g element cget $name -ydata]
	$yv expr {$yv + $size}
	set data($name,shift) $shift
}

scrolledgraph method reconf {args} {
	private $object data conf
	set name $conf(entry)
	set data($name,color) [get conf(color) gray]
	set style [$object gradientstyle $data($name,color)]
	set color [get colors($data($name,color)-0) $data($name,color)]
	$object.g element configure e$name -linewidth $conf(linewidth) -fill $color -outline $color -color $color \
		-outlinewidth 1 -pixels 2 -symbol $conf(symbol) \
		-styles $style
	Classy::todo $object _configureevent
}

scrolledgraph method confcurrent {args} {
	private $object data conf
	set name $conf(entry)
	set conf(symbol) [$object.g element cget e$name -symbol]
	set conf(linewidth) [$object.g element cget e$name -linewidth]
	set conf(color) $data($name,color)
	set conf(fields) $data($name,header)
	lappend conf(fields) ""
	set conf(Y) [lindex $data($name,elements) 1]
	set conf(W) [lindex $data($name,elements) 2]
}

scrolledgraph method elconf {name} {
return
	puts "conf $name"
	private $object conf
	catch {destroy $object.conf}
	Classy::Dialog $object.conf -title "Configure data"
	set conf(entry) $name
	set conf(symbols) {plus splus square circle diamond cross scross triangle arrow none}
	$object confcurrent
	$object.conf option select "Entry" [privatevar $object conf(entry)] [privatevar $object data(entries)] \
		-command "$object confcurrent"
	$object.conf option color "Color" [privatevar $object conf(color)] \
		-command "$object reconf"
	$object.conf option select "Symbol" [privatevar $object conf(symbol)] [privatevar $object conf(symbols)] \
		-command "$object reconf"
	$object.conf option numentry "Line width" [privatevar $object conf(linewidth)] \
		-command "$object reconf"
	$object.conf option numentry "Y shift" [privatevar $object conf(yshift)] \
		-command "$object shift \[getprivate $object conf(entry)\]"
	$object.conf option select "Y" [privatevar $object conf(Y)] [privatevar $object conf(fields)] \
		-command "$object changefields Y"
	$object.conf option select "Weight" [privatevar $object conf(W)] [privatevar $object conf(fields)] \
		-command "$object changefields W"
	$object.conf option string "Trans" [privatevar $object conf(trans)] \
		-command "$object changefields W"
	$object.conf option button "Remove" \
		"$object delelement \[getprivate $object conf(entry)\]"
}

scrolledgraph method changefields {what args} {
	private $object data conf
	set name $conf(entry)
	set header $data($name,header)
	lset data($name,elements) 1 $conf(Y)
	lset data($name,elements) 2 $conf(W)
	lset data($name,poss) 1 [lsearch $header $conf(Y)]
	lset data($name,poss) 2 [lsearch $header $conf(W)]
	set data($name,lstart) 0
	set data($name,lend) 0
	set data($name,trans) $conf(trans)
	Classy::todo $object redraw
}

scrolledgraph method delelement {name} {
	private $object data
	set vnum $data($name,vnum)
	foreach field {x y w region i} {
		catch {vector destroy ::$object.$vnum.$field}
	}
	$object.g element	delete $name
	set list [array names data $name,*]
	foreach el $list {
		unset data($el)
	}
	set data(entries) [list_remove $data(entries) $name]
}

scrolledgraph method reload {} {
return
}

scrolledgraph method xset {args} {
	private $object pxv
	if {[get pxv ""] ne $args} {
		$object.scx configure -command {}
		$object.scx set {*}$args
		$object.scx configure -command [list $object xview]
		set pxv $args
	}
}

scrolledgraph method xview {args} {
	private $object data region
	set w $object.g
	set amin [::$object.ends.x index 0]
	set amax [::$object.ends.x index 1]
	set xmin [expr {round([$w axis cget x -min])}]
	set xmax [expr {round([$w axis cget x -max])}]
	set cursize [expr {$xmax-$xmin}]
	if {$cursize <= 0} {set cursize 1000}
	set first [lindex $args 0]
	switch $first {
		"" {
			set cursx [expr {}]]
			return [$w axis view x]
		}
		moveto {
			set fraction [lindex $args 1]
			set xmin [expr {round($amin + $fraction*($amax-$amin+1))}]
		}
		scroll {
			set number [lindex $args 1]
			set what [lindex $args 2]
			if {$what eq "pages"} {
				set xmin [expr {round($xmin+$number*$cursize - ($cursize/10))}]
			} else {
				set xmin [expr {round($xmin+$number*($cursize/10))}]
			}
		}
	}
	if {$xmin < $amin} {set xmin $amin}
	if {$xmin > $amax} {set xmin $amax}
	set xmax [expr {$xmin+$cursize}]
	$w axis configure x -min $xmin -max $xmax
	$object _setvars
	Classy::todo $object reload
}

scrolledgraph method yset {args} {
	private $object pyv
	if {[get pyv ""] ne $args} {
		$object.scy configure -command {}
		$object.scy set {*}$args
		$object.scy configure -command [list $object yview]
		set pyv $args
	}
}

scrolledgraph method yview {args} {
	private $object data region
	set w $object.g
	set amin [::$object.ends.y index 0]
	set amax [::$object.ends.y index 1]
	set ymin [expr {round([$w axis cget y -min])}]
	set ymax [expr {round([$w axis cget y -max])}]
	set cursize [expr {$ymax-$ymin}]
	set first [lindex $args 0]
	switch $first {
		"" {
			set cursy [expr {}]]
			return [$w axis view y]
		}
		moveto {
			set fraction [lindex $args 1]
			set ymin [expr {$amax - $fraction*($amax-$amin+1)}]
		}
		scroll {
			set number [lindex $args 1]
			set what [lindex $args 2]
			if {$what eq "pages"} {
				set ymin [expr {$ymin - $number*$cursize + ($cursize/10)}]
			} else {
				set ymin [expr {$ymin - $number*($cursize/10)}]
			}
		}
	}
	if {$ymin < $amin} {set ymin $amin}
	if {$ymin > $amax} {set ymin $amax}
	set ymax [expr {$ymin+$cursize}]
	$w axis configure y -min $ymin -max $ymax
	$object _setvars
	Classy::todo $object reload
}

scrolledgraph method print {args} {
	global printdialog
	private $object region options
	if {![llength $args]} {
		destroy $object.ps
		set file [file root [file tail [get options(-datafile) plot]]]__[$object.g axis cget y -title]_vs_[$object.g axis cget x -title]__$region(xmin)-$region(xmax).ps
		regsub -all {[:\?\*/]} $file _ file
		set printdialog(file) [file dir [get options(-datafile) [pwd]]]/$file
		Classy::Dialog $object.ps -title "Print to postscript"
		$object.ps option file "Filename" printdialog(file)
		$object.ps add print "Print" "[list $object] print \$printdialog(file)" default
		$object.ps persistent remove print
	} else {
		foreach {file} $args break
		$object.g postscript configure -landscape yes -maxpect yes
		$object.g postscript output $file
	}
}

scrolledgraph method point {x y} {
	private $object data points ptable
	puts $x,$y
	$object.g element closest $x $y pd
	set name $pd(name)
	set index $pd(index)
	set ext [file extension $data($name,file)]
	if {[inlist {.gz .bgz} $ext]} {
		set f [open "| tabix [list $data($name,file)] $data($name,bgzregion)"]
	} else {
		set f [gzopen $data($name,file) $data($name,fpos)]
	}
	incr index
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		incr index -1
		if {!$index} break
	}
	catch {close $f}
	if {![winfo exists [get ptable]]} {
		unset -nocomplain points
		set ptable [tableedit]
	}
	if {![info exists points] || [llength $points] <= 1} {
		set header $data($name,header)
		list_unshift header sample
		set points [list $header]
	}
	list_unshift line $name
	lappend points $line
	$ptable.editor configure -variable [privatevar $object points]
	$ptable.editor autosize
}

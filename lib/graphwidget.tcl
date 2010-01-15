package require ClassyTk
package require Extral
package require rbc
namespace import ::rbc::*

Widget subclass graphwidget

graphwidget method init {args} {
	super init frame
	graph $object.g
	scrollbar $object.scy -orient vertical
	scrollbar $object.scx -orient horizontal
	frame $object.b
	grid $object.b -sticky nwse -row 0 -columnspan 2
	grid $object.g $object.scy -row 1 -sticky nwse
	grid $object.scx -row 2 -sticky nwse
	grid columnconfigure $object 0 -weight 1
	grid rowconfigure $object 1 -weight 1
	# link graph --> scrollbars
	$object.g axis configure x -scrollcommand [list $object.scx set]
	$object.g axis configure y -scrollcommand [list $object.scy set]
	# link scrollbars --> graph
	$object.scx configure -command [list $object.g axis view x]
	$object.scy configure -command [list $object.g axis view y]
	$object.g grid configure -hide no
	Rbc_ZoomStack $object.g
	foreach num {90 80 70 60 50 40 30} {
		$object.g pen create gray$num -fill gray$num -outline gray$num -outlinewidth 1 -pixels 2 -symbol plus
	}
	if 0 {
		foreach num {90 80 70 60 50 40 30} {
			$object.g pen configure gray$num -fill gray$num -outline gray$num -outlinewidth 1 -pixels 2 -symbol plus
		}
	}
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
	bind $object.g <Configure> [list $object _fillregion]
	button $object.b.open -text "Open" -command [list $object open]
	pack $object.b.open -side left
	button $object.b.clear -text "Clear" -command [list $object clear]
	pack $object.b.clear -side left
	foreach type {xmin xmax ymin ymax} {destroy $object.b.$type
		Classy::NumEntry $object.b.$type -width 10 -label $type -orient horizontal \
			-textvariable [privatevar $object region($type)] -command [list $object redraw]
		pack $object.b.$type -side left
	}
	$object clear
}

graphwidget chainoptions {$object.g}

graphwidget method clear {} {
	private $object xmin xmax ymin ymax data
	foreach el [$object.g element names] {
		$object.g element	delete $el
	}
	foreach v [vector names $object.*] {
		vector destroy ::$v
	}
	unset -nocomplain data
	lassign {0 0 0 0} xmin xmax ymin ymax
	set data(entries) {}
}

graphwidget method gradientstyle {basecolor} {
	private $object colors
	if {$basecolor eq "gray"} {
		foreach {num} {
			100 90 80 70 60 50 40 30 20 10 0
		} {
			set color gray-$num
			set colorcode gray$num
			set colors($color) $colorcode
			if {![llength [$object.g pen names $color]]} {
				$object.g pen create $color -fill $colorcode -outline $colorcode -outlinewidth 1 -pixels 2 -symbol plus
			} else {
				# $object.g pen configure $color -fill $colorcode -outline $colorcode -outlinewidth 1 -pixels 2 -symbol plus
			}
		}
	} else {
		foreach {num factor} {
			100 0.5 90 0.4 80 0.3 70 0.2 60 0.1 50 0.0 40 -0.1 30 -0.2 20 -0.3 10 -0.4 0 -0.5
		} {
			set color $basecolor-$num
			set colorcode [gradient $basecolor $factor]
			if {![llength [$object.g pen names $color]]} {
				$object.g pen create $color -fill $colorcode -outline $colorcode -outlinewidth 1 -pixels 2 -symbol plus
			} else {
				# $object.g pen configure $color -fill $colorcode -outline $colorcode -outlinewidth 1 -pixels 2 -symbol plus
			}
		}
	}
	set style {}
	lappend style [list $basecolor-80 0 10]
	lappend style [list $basecolor-60 11 20]
	lappend style [list $basecolor-40 21 30]
	lappend style [list $basecolor-0 31 10000000]
	return $style
}

graphwidget method open {{file {}} {region 0}} {
	private $object xmin xmax ymin ymax xv yv wv data colors
	global graphd
	if {$file eq ""} {set file [Classy::selectfile]}
	set name [file root [file tail $file]]
	set f [open $file]
	set header [tsv_open $f]
	set graphd(header) $header
	lappend graphd(header) {}
	set graphd(region) 0
	set graphd(do) 0
	set graphd(nan) -50
	foreach {graphd(xfield) graphd(yfield) graphd(wfield)} $graphd(header) break
	if {[lsearch $header begin] != -1} {
		set graphd(xfield) begin
		set graphd(region) 1
	} elseif {[lsearch $header start] != -1} {
		set graphd(xfield) start
		set graphd(region) 1
	}
	if {[lsearch $header end] != -1} {
		set graphd(yfield) end
		set graphd(region) 1
	}
	if {$graphd(region) == 1} {
		set graphd(wfield) ""
	}
	Classy::Dialog $object.open -title "Open file $file"
	$object.open option select "X" graphd(xfield) graphd(header)
	$object.open option select "Y" graphd(yfield) graphd(header)
	$object.open option select "Weight" graphd(wfield) graphd(header)
	$object.open option check "Regions" graphd(region) ""
	$object.open option entry "Nan" graphd(nan)
	$object.open add do Do {set graphd(do) 1} default
	tkwait window $object.open
	if {!$graphd(do)} return
	set region $graphd(region)
	set elements [list $graphd(xfield) $graphd(yfield) $graphd(wfield)]
	set poss {}
	foreach el $elements {
		lappend poss [lsearch $header $el]
	}
	set vnum [llength $data(entries)]
	if {$region} {
		set elements [list range y $graphd(wfield)]
	}
	foreach field [list_remove $elements {}] {
		vector create ::$object.$vnum.$field
	}
	set basecolor [lindex {gray blue yellow green orange red violet} $vnum]
	lappend data(entries) $name
	set data($name,header) $header
	set data($name,elements) $elements
	if {$region} {
		while {![eof $f]} {
			set line [split [gets $f] \t]
			if {![llength $line]} continue
			foreach {begin end} [list_sub $line $poss] break
			::$object.$vnum.region append $end
			::$object.$vnum.region append $begin
			::$object.$vnum.y append 0
			::$object.$vnum.y append 0
		}
		set graphd(xfield) region
		set graphd(yfield) y
	} else {
		while {![eof $f]} {
			set line [split [gets $f] \t]
			if {![llength $line]} continue
			foreach el $elements p $poss v [list_sub $line $poss] {
				if {$p == -1} continue
				if {[isdouble $v]} {
					::$object.$vnum.$el append $v
				} else {
					::$object.$vnum.$el append $graphd(nan)
				}
			}
		}
	}
	close $f
	set xv ::$object.$vnum.$graphd(xfield)
	set yv ::$object.$vnum.$graphd(yfield)
	if {$graphd(wfield) ne ""} {
		set wv ::$object.$vnum.$graphd(wfield)
	}
	if {[get ${xv}(min)] < $xmin} {
		set xmin [get ${xv}(min)]
	}
	if {[get ${xv}(max)] > $xmax} {
		set xmax [get ${xv}(max)]
	}
	if {[get ${yv}(min)] < $ymin} {
		set ymin [get ${yv}(min)]
	}
	if {[get ${yv}(max)] > $ymax} {
		set ymax [get ${yv}(max)]
	}
	$object.g element create $name -xdata ::$xv -ydata ::$yv -symbol none
	if {[llength $header] == 3} {
		$object.g element configure $name -weights ::$wv
	}
	set style [$object gradientstyle $basecolor]
	set color [get colors(${basecolor}0) $basecolor]
	$object.g element configure $name -linewidth 0 -fill $color -outlinewidth 1 -pixels 1 -symbol circle
	$object.g element configure $name -linewidth 0 -outline $color -outlinewidth 1 -pixels 2 -symbol plus \
		-styles $style
	if {$region} {
		$object.g element configure $name -linewidth 5 -trace decreasing
	}
	set data($name,color) $color
	event generate $object.g <Configure>
	$object.g legend bind $name <1> [list $object elconf $name]
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

graphwidget method redraw {args} {
	private $object region
	set w $object.g
	$w axis configure x -min $region(xmin) -max $region(xmax)
	$w axis configure y -min $region(ymin) -max $region(ymax)
	event generate $object.g <Configure>
}

graphwidget method _fillregion {args} {
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

graphwidget method zoomx {{factor 2}} {
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
	if {$min < $xmin} {
		set min $xmin
		set max [expr {$min+$newwidth}]
	}
	if {$max > $xmax} {
		set max $xmax
	}
	$w axis configure x -min $min
	$w axis configure x -max $max
	event generate $object.g <Configure>
}

graphwidget method zoomy {{factor 2}} {
	private $object ymin ymax
	set w $object.g
putsvars w ymin ymax factor
	set min [$w axis cget y -min]
	set max [$w axis cget y -max]
	if {$min eq ""} {set min $ymin}
	if {$max eq ""} {set max $ymax}
	set newheight [expr {$factor * ($max-$min)}]
	set extra [expr {($newheight - ($max-$min))/2}]
	set min [expr {$min-$extra}]
	set max [expr {$max+$extra}]
	if {$min < $ymin} {
		set min $ymin
		set max [expr {$min+$newheight}]
	}
	if {$max > $ymax} {
		set max $ymax
	}
	$w axis configure y -min $min
	$w axis configure y -max $max
	event generate $object.g <Configure>
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

graphwidget method shift {name {shift 1}} {
	private $object data
	if {![info exists data($name,shift)]} {
		set data($name,shift) 0
	}
	set size [expr {$shift-$data($name,shift)}]
	set yv [$object.g element cget $name -ydata]
	$yv expr {$yv + $size}
	set data($name,shift) $shift
}

graphwidget method reconf {args} {
	private $object data conf
	set name $conf(entry)
	set data($name,color) [get conf(color) gray]
	set style [$object gradientstyle $data($name,color)]
	set color [get colors($data($name,color)-0) $data($name,color)]
	$object.g element configure $name -linewidth $conf(linewidth) -fill $color -outline $color -color $color \
		-outlinewidth 1 -pixels 2 -symbol $conf(symbol) \
		-styles $style
	event generate $object.g <Configure>
}

graphwidget method confcurrent {args} {
	private $object data conf
	set name $conf(entry)
	set conf(symbol) [$object.g element cget $name -symbol]
	set conf(linewidth) [$object.g element cget $name -linewidth]
	set conf(color) $data($name,color)
}

graphwidget method elconf {name} {
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
}

if 0 {

lappend auto_path /home/peter/bin/tcl
lappend auto_path /home/peter/dev/completegenomics/lib

set object .g

	package require Tk
	graphwidget .g
	pack .g -fill both -expand yes

}

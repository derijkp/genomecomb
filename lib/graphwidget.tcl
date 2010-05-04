package require ClassyTk
package require Extral
package require rbc
namespace import ::rbc::*

Widget subclass graphwidget

graphwidget method init {args} {
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
	button $object.b.open -text "Open" -command [list $object opendialog]
	pack $object.b.open -side left
	button $object.b.clear -text "Clear" -command [list $object clear]
	pack $object.b.clear -side left
	button $object.b.print -text "Print" -command [list $object print]
	pack $object.b.print -side left
	button $object.b.check -text "Check" -command [list $object checkfile]
	pack $object.b.check -side left
	button $object.b.cmd -text "Cmd" -command {Classy::cmd}
	pack $object.b.cmd -side left
	Classy::Entry $object.b.xrange -width 20 -label xrange -orient horizontal \
		-textvariable [privatevar $object region(xrange)] -command [list Classy::todo $object redraw]
	bind $object.b.xrange <3> "$object paste; break"
	pack $object.b.xrange -side left
	foreach type {xrangeextra ymin ymax} {destroy $object.b.$type
		Classy::NumEntry $object.b.$type -width 10 -label $type -orient horizontal \
			-textvariable [privatevar $object region($type)] -command [list Classy::todo $object redraw]
		pack $object.b.$type -side left
	}
	Classy::ProgressWidget $object.progress
	grid $object.progress -sticky nwse -columnspan 2
	Classy::Progress display $object.progress
	$object.b.xrange configure -width 20
	$object.b.xrangeextra configure -width 5
	vector create ::$object.ends.x
	vector create ::$object.ends.y
	::$object.ends.x set {0 1}
	::$object.ends.y set {0 0}
	$object clear
	$object.g element create legend -xdata ::$object.ends.x -ydata ::$object.ends.y -symbol none -linewidth 0
	bindtags $object.g [list_concat [list $object.g] [list_remove [bindtags $object.g] $object.g]]
	Classy::todo $object start
}

graphwidget chainoptions {$object.g}

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

graphwidget method paste {} {
	private $object region
	set region(xrange) [::tk::GetSelection $object PRIMARY]
	Classy::todo $object redraw
}

graphwidget method start {} {
	private $object region
	array set region {xrange {1000 9000} xrangeextra 1000 ymin -2 ymax 600}
	Classy::todo $object redraw
}

graphwidget method clear {} {
	private $object xmin xmax ymin ymax data
	foreach el [list_remove [$object.g element names] legend] {
		$object.g element	delete $el
	}
	foreach v [vector names $object.*] {
		vector destroy ::$v
	}
	unset -nocomplain data
	lassign {0 0 0 0} xmin xmax ymin ymax
	set data(entries) {}
	::$object.ends.x set {0 1}
	::$object.ends.y set {0 0}
	$object.scx set 0 1
	event generate $object.g <Configure>
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
	lappend style [list $basecolor-80 -1000000000 10]
	lappend style [list $basecolor-60 11 20]
	lappend style [list $basecolor-40 21 30]
	lappend style [list $basecolor-0 31 1000000000]
	return $style
}

graphwidget method opendialog {{file {}}} {
	global graphd
	if {$file eq ""} {set file [Classy::selectfile]}
	set f [rzopen $file]
	set header [tsv_open $f]
	catch {close $f}
	set graphd(file) [file normalize $file]
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
	} elseif {[lsearch $header end1] != -1} {
		set graphd(xfield) end1
		set graphd(yfield) dist
		set graphd(region) 0
		set graphd(wfield) "weight1"
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
	$object.open option entry "Chromosome" graphd(chr)
	$object.open add do Do [list $object open $file ::graphd] default
}

graphwidget method _configureevent {} {
	private $object data region
	event generate $object.g <Configure>
	set region(status) ""
	set region(cancel) 0
}

graphwidget method open {file settingsVar} {
	private $object xv yv wv data colors
	upvar $settingsVar graphd
	set showregion $graphd(region)
	set vnum [llength $data(entries)]
	foreach field {x y w i} {
		vector create ::$object.$vnum.$field
	}
	set elements [list $graphd(xfield) $graphd(yfield) $graphd(wfield)]
	set indexname [rzroot $file].$graphd(xfield)_index
	set name [file root [file tail [rzroot $file]]]
	set f [rzopen $file]
	set header [tsv_open $f]
	set poss {}
	foreach el $elements {
		lappend poss [lsearch $header $el]
	}
	set xpos [lindex $poss 0]
	set ypos [lindex $poss 1]
	set basecolor [lindex {blue red gray orange yellow green violet} $vnum]
	lappend data(entries) $name
	set data($name,chr) $graphd(chr)
	set data($name,vnum) $vnum
	set data($name,lstart) 0
	set data($name,lend) 0
	set data($name,file) [file normalize $file]
	set data($name,header) $header
	set data($name,elements) $elements
	set data($name,region) $graphd(region)
	set data($name,nan) $graphd(nan)
	set data($name,poss) $poss
	# create index
	if {$graphd(region)} {
		$object loadregtype $file $f $name
	} else {
		set indexed 1
		set index ::$object.$vnum.i
		if {![file exists $indexname]} {
			tsv_index $file $graphd(xfield)
		}
		# read index
		set o [open $indexname]
		set data($name,step) [gets $o]
		set data($name,findex) [gets $o]
		set xmin [gets $o]
		set data($name,xmin) $xmin
		set data($name,xmax) [gets $o]
		set data($name,fx) [expr {$data($name,xmin)-$xmin%10000}]
		set temp [split [string trim [read $o]] \n]
		::$index set $temp
		close $o
		#
		set list [array names data *,xmin]
		set amin $data([lindex $list 0])
		foreach n $list {
			if {$data($n) < $amin} {set amin $data($n)}
		}
		set list [array names data *,xmax]
		set amax $data([lindex $list 0])
		foreach n [array names data *,xmax] {
			if {$data($n) > $amin} {set amax $data($n)}
		}
		::$object.ends.x set [list $amin $amax]
		::$object.ends.y set {0 0}
		set data($name,indexed) $indexed
	}
	catch {close $f}
	set xv ::$object.$vnum.x
	set yv ::$object.$vnum.y
	set wv ::$object.$vnum.w
	$object.g element create $name -xdata ::$xv -ydata ::$yv -symbol none
	$object.g element configure $name -weights ::$wv
	set style [$object gradientstyle $basecolor]
	set color [get colors(${basecolor}0) $basecolor]
	$object.g element configure $name -linewidth 0 -fill $color -outlinewidth 1 -pixels 1 -symbol circle
	$object.g element configure $name -linewidth 0 -outline $color -outlinewidth 1 -pixels 2 -symbol plus \
		-styles $style
	if {$showregion} {
		$object.g element configure $name -linewidth 10 -trace decreasing
	}
	set data($name,color) $color
	$object.g legend bind $name <1> [list $object elconf $name]
	Classy::todo $object reload
	Classy::todo $object _configureevent
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

graphwidget method _xrange {args} {
	private $object region
	set xmin ""
	set xmax ""
	if {![isdouble $region(xrangeextra)]} {set region(xrangeextra) 0}
	set temp [regexp -inline -all {[0-9.]+} $region(xrange)]
	foreach {xmin xmax} $temp break
	if {[isdouble $xmin] && ![isdouble $xmax]} {
		set xmax [expr {$xmin+1}]
	}
	if {[isdouble $xmin] && [isdouble $xmax] && [isdouble $region(xrangeextra)]} {
		set xmin [expr {$xmin - $region(xrangeextra)}]
		set xmax [expr {$xmax + $region(xrangeextra)}]
	}
	set region(xmin) [expr {round($xmin)}]
	set region(xmax) [expr {round($xmax)}]
	if {[llength $temp] > 2} {
		set y [lindex $temp 2]
		if {$y < 3000} {
			set region(ymin) -1
			set region(ymax) [expr {$y+500}]
		} else {
			set region(ymin) [expr {$y-500}]
			set region(ymax) [expr {$y+500}]
		}
	}
	set region(xrange) $xmin-$xmax
}

graphwidget method redraw {args} {
puts ----------redraw----------
	private $object region
	set w $object.g
	Classy::canceltodo $object redraw
	$object _xrange
	$w axis configure x -min $region(xmin) -max $region(xmax)
	$w axis configure y -min $region(ymin) -max $region(ymax)
	Classy::todo $object reload
}

graphwidget method _setvars {} {
	private $object region
	set w $object.g
	set xmin [expr {round([$w axis cget x -min])}]
	set xmax [expr {round([$w axis cget x -max])}]
puts _setvars($xmin-$xmax)
	if {![isdouble $region(xrangeextra)]} {set region(xrangeextra) 0}
	if {[isdouble $xmin] && [isdouble $xmax]} {
		set region(xmin) $xmin
		set region(xmax) $xmax
		set region(xrange) [expr {$xmin + $region(xrangeextra)}]-[expr {$xmax - $region(xrangeextra)}]
	}
	set ymin [$w axis cget y -min]
	set ymax [$w axis cget y -max]
	catch {set ymin [format %.0f $ymin]}
	catch {set ymax [format %.0f $ymax]}
	set region(ymin) $ymin
	set region(ymax) $ymax
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

graphwidget method zoomy {{factor 2}} {
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
putsvars name color
	$object.g element configure $name -linewidth $conf(linewidth) -fill $color -outline $color -color $color \
		-outlinewidth 1 -pixels 2 -symbol $conf(symbol) \
		-styles $style
	Classy::todo $object _configureevent
}

graphwidget method confcurrent {args} {
	private $object data conf
	set name $conf(entry)
	set conf(symbol) [$object.g element cget $name -symbol]
	set conf(linewidth) [$object.g element cget $name -linewidth]
	set conf(color) $data($name,color)
	set conf(fields) $data($name,header)
	lappend conf(fields) ""
	set conf(Y) [lindex $data($name,elements) 1]
	set conf(W) [lindex $data($name,elements) 2]
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
	$object.conf option select "Y" [privatevar $object conf(Y)] [privatevar $object conf(fields)] \
		-command "$object changefields Y"
	$object.conf option select "Weight" [privatevar $object conf(W)] [privatevar $object conf(fields)] \
		-command "$object changefields W"
	$object.conf option string "Trans" [privatevar $object conf(trans)] \
		-command "$object changefields W"
	$object.conf option button "Remove" \
		"$object delelement \[getprivate $object conf(entry)\]"
}

graphwidget method changefields {what args} {
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

graphwidget method delelement {name} {
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

graphwidget method loadregion {name} {
	private $object region data
	array set trans [get data($name,trans) ""]
	if {[get region(cancel) 0]} return
	set start $region(xmin)
	set end $region(xmax)
	set start [expr {round([$object.g axis cget x -min])}]
	set end [expr {round([$object.g axis cget x -max])}]
	if {($start >= $data($name,lstart)) && ($end <= $data($name,lend))} return
	incr start -5000
	incr end 5000
	set start [expr {round($start)-round($start)%10000}]
	if {$start < $data($name,findex)} {set start $data($name,findex)}
puts "load $name $start $end $data($name,lstart) $data($name,lend)"
	if {![info exists data($name,vnum)]} return
	set vnum $data($name,vnum)
	set poss $data($name,poss)
	set index ::$object.$vnum.i
	$object.progress configure -message "Loading $name"
	if {$data($name,region)} {
#		while {![eof $f]} {
#			set line [split [gets $f] \t]
#			if {![llength $line]} continue
#			foreach {begin end} [list_sub $line $poss] break
#			::$object.$vnum.region append $end
#			::$object.$vnum.region append $begin
#			::$object.$vnum.y append 0
#			::$object.$vnum.y append 0
#		}
#		set graphd(xfield) region
#		set graphd(yfield) y
	} else {
		foreach el {x y w} v [list $data($name,xmin) 0 0] {
			::$object.$vnum.$el set {}
		}
		update
		if {[get region(cancel) 0]} {
			$object.progress configure -message "Loading $name canceled"
			puts "Loading $name canceled"
			return
		}
		set xpos [lindex $poss 0]
		set data($name,lstart) 0
		set data($name,lend) 0
		set fpos [expr {round([::$index index [expr {($start-$data($name,findex))/10000}]])}]
		set data($name,fpos) $fpos
		set f [rzopen $data($name,file) $fpos]
		set pnext [expr {$start+200}]
		set tot [expr {$end-$start}]
		while {![eof $f]} {
			set line [split [gets $f] \t]
			if {![llength $line]} continue
			set x [lindex $line $xpos]
			if {$x > $pnext} {
				update
				$object.progress configure -message "Loading $name [format %.1f [expr {100*($x-$start)/$tot}]]% ([get region(cancel) 0])"
				puts "$x [format %.1f [expr {100*($x-$start)/$tot}]]%"
				if {[get region(cancel) 0]} {
					$object.progress configure -message "Loading $name canceled"
					puts "Loading $name canceled"
					return
				}
				incr pnext 2000
			}
			if {$x > $end} break
			foreach el {x y w} p $poss v [list_sub $line $poss] {
				if {$p == -1} continue
				if {[info exists trans($v)]} {set v $trans($v)}
				if {[isdouble $v]} {
					::$object.$vnum.$el append $v
				} else {
					::$object.$vnum.$el append $data($name,nan)
				}
			}
		}
		catch {close $f}
		set data($name,lstart) $start
		set data($name,lend) $end
		if {[get data($name,shift) 0] != 0} {
			set yv [$object.g element cget $name -ydata]
			$yv expr {$yv + $data($name,shift)}
		}
	}
	puts "Finished loading $name"
	$object.progress configure -message "Finished loading $name"
}

graphwidget method loadregtype {file f name} {
	private $object region data
	if {![info exists data($name,vnum)]} return
	set vnum $data($name,vnum)
	set poss $data($name,poss)
	set chr $data($name,chr)
	set chrpos [lsearch $data($name,header) chromosome]
	chrindexseek $file $f $chr
	list_pop poss
	lappend poss $chrpos
	set index ::$object.$vnum.i
	::$object.$vnum.x set {}
	::$object.$vnum.y set {}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		foreach {begin endw cchr} [list_sub $line $poss] break
		if {$cchr ne $chr} break
		::$object.$vnum.x append $endw
		::$object.$vnum.x append $begin
		::$object.$vnum.y append 1
		::$object.$vnum.y append 1
	}
	update
	$object.progress configure -message "Finished loading $name"
}

graphwidget method reload {} {
	private $object data region
	if {[get region(status) ""] eq "loading"} {
		set region(cancel) 1
		return
	}
	set region(status) "loading"
	unset -nocomplain region(cancel)
	set error [catch {
		set miny 0
		set maxy 0
		foreach name $data(entries) {
			set vnum $data($name,vnum)
			set index ::$object.$vnum.y
			if {$data($name,region)} continue
			$object loadregion $name
			if {[get region(cancel) 0]} {
				puts "Loading canceled"
				set region(cancel) 0
				set region(status) ""
				puts "Loading $name canceled"
				Classy::todo $object reload
				return
				break
			}
			if {[get region(cancel) 0]} break
			set y [vector expr max($index)]
			if {$y > $maxy} {set maxy $y}
			set y [vector expr min($index)]
			if {$y < $miny} {set miny $y}
			puts "$name loaded"
		}
		::$object.ends.y set [list $miny $maxy]
	} e]
	if {$error} {
		puts ERROR:$e
	}
	Classy::todo $object _configureevent
	set region(status) "finished"
}

graphwidget method xset {args} {
	private $object pxv
	if {[get pxv ""] ne $args} {
		$object.scx configure -command {}
		$object.scx set {*}$args
		$object.scx configure -command [list $object xview]
		set pxv $args
	}
}

graphwidget method xview {args} {
	private $object data region
	set w $object.g
	set amin [::$object.ends.x index 0]
	set amax [::$object.ends.x index 1]
	set xmin [expr {round([$w axis cget x -min])}]
	set xmax [expr {round([$w axis cget x -max])}]
	set cursize [expr {$xmax-$xmin}]
	set first [lindex $args 0]
	switch $first {
		"" {
			set cursx [expr {}]]
			return [$w axis view x]
		}
		moveto {
			set fraction [lindex $args 1]
			set xmin [expr {$amin + $fraction*($amax-$amin+1)}]
		}
		scroll {
			set number [lindex $args 1]
			set what [lindex $args 2]
			if {$what eq "pages"} {
				set xmin [expr {$xmin+$number*$cursize - ($cursize/10)}]
			} else {
				set xmin [expr {$xmin+$number*($cursize/10)}]
			}
		}
	}
	set xmax [expr {$xmin+$cursize}]
	$w axis configure x -min $xmin -max $xmax
	$object _setvars
	Classy::todo $object reload
}

graphwidget method yset {args} {
	private $object pyv
	if {[get pyv ""] ne $args} {
		$object.scy configure -command {}
		$object.scy set {*}$args
		$object.scy configure -command [list $object yview]
		set pyv $args
	}
}

graphwidget method yview {args} {
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

graphwidget method print {args} {
	global printdialog
	private $object region
	if {![llength $args]} {
		destroy $object.ps
		set file $region(xrange)
		Classy::Dialog $object.ps -title "Print to postscript"
		$object.ps option file "Filename" printdialog(file)
		$object.ps add print "Print" "[list $object] print \$printdialog(file)" default
		$object.ps persistent remove print
	}
	foreach {file} $args break
	$object.g postscript configure -landscape yes -maxpect yes
	$object.g postscript output $file
}

graphwidget method checkbrowse {args} {
	private $object checkw region
	set table $checkw.editor
	upvar #0 [$table cget -variable] checkdata
	set row [$table index active row]
	$table selection clear
	$table selection set $row,0 $row,[$table cget -cols]
	set header [lindex $checkdata 0]
	set line [lindex $checkdata $row]
	set poss [list_cor $header {patchstart pos size}]
	set cur [list_sub $line $poss]
	set region(xrange) $cur
	Classy::todo $object redraw
}

graphwidget method checkset {value} {
	private $object checkw region
	set table $checkw.editor
	upvar #0 [$table cget -variable] checkdata
	set row [$table index active row]
	set header [lindex $checkdata 0]
	set line [lindex $checkdata $row]
	set pos [lsearch $header check]
	switch $value {
		d {set value diff}
		m {set value mdiff}
		u {set value uncert}
		n {set value noind}
		b {set value both}
		a {set value artefact}
		s {set value same}
		p {set value pdip}
		2 {set value double}
		t {set value true}
		f {set value false}
		l {set value likely}
		default {error "unknown value"}
	}
	lset checkdata $row $pos $value
	incr row
 	$table see $row,1
	$table activate $row,1
}

graphwidget method checkfile {} {
	private $object checkw
	set file [Classy::selectfile]

	# set file /complgen/sv/svcompar_sv78_sv79-20.tsv
	# set file /complgen/sv/svcompar_GS103_GS102-9.tsv
	private $object checkw
	catch {destroy $checkw}
	set checkw [tableedit]
	set table $checkw.editor
	$table load $file
	upvar #0 [$table cget -variable] checkdata
	set header [lindex $checkdata 0]
	set chpos [lsearch $header check]
	set end1pos [lsearch $header pos]
	if {$chpos == -1} {
		set size 0
		set temp {}
		foreach line $checkdata {
			lappend temp [linsert $line 1 {}]
			set size [max $size [llength $line]]
		}
		incr size
		lset temp 0 1 check
		set checkdata $temp
		$table autosize
		$table configure -browsecommand "[list $object] checkbrowse" -cols [llength $header]
	} else {
		$table configure -browsecommand "[list $object] checkbrowse"
	}
	# sort
	set header [lindex $checkdata 0]
	set cor [list_cor $header {type quality problems size match}]
	set temp [list $header]
	set endtemp {}
	set matchpos [lindex $cor end]
	foreach line [lsort -index $matchpos [lrange $checkdata 1 end]] {
		foreach {type quality problems size match} [list_sub $line $cor] break
		set move 0
		if {($quality <= 3) || ($type eq "trans") \
			|| ([llength [list_common {msmall trfbad trfartefact} $problems]]) \
			|| (($type eq "ins") && ($size > 250))
		} {
			set move 1
		} elseif {$type eq "trans"} {
			set move 1
		} elseif {$match eq "db"} {
			set move 1
		}
		if {$move} {
			lappend endtemp $line
		} else {
			lappend temp $line
		}
	}
	lappend temp {*}$endtemp
	# set checkdata
	set checkdata $temp
	set header [lindex $checkdata 0]
	#
	$table configure -titlerows 1 -labels $header
	$table configure -labelcommand [list $table sort]
	foreach key {d m u n b a s p 2 t f l} {
		bind $table <KeyPress-$key> "[list $object] checkset %K; break"
	}
}

graphwidget method point {x y} {
	private $object data points ptable
	puts $x,$y
	$object.g element closest $x $y pd
	set name $pd(name)
	set index $pd(index)
	set f [rzopen $data($name,file) $data($name,fpos)]
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

if 0 {

	package require Tclx
	signal -restart error SIGINT
lappend auto_path /home/peter/bin/tcl
lappend auto_path /home/peter/dev/completegenomics/lib
cd /complgen/sv
set object .g
	package require Tk
	graphwidget .g
	pack .g -fill both -expand yes
set file sv79-20-pairs.tsv
set file /complgen/sv/GS103/GS103-9-paired.tsv.rz
.g opendialog $file


.g clear
set chr 20
.g open /complgen/sv/GS103/GS103-$chr-paired.tsv.rz ::graphd
.g open /complgen/sv/GS102/GS102-$chr-paired.tsv.rz ::graphd


set file sv79-20-pairs.tsv
set file sv78-20-pairs.tsv
catch {close $f}
catch {close $o}

$object xview scroll 1 pages

$object reload
$object.scx configure -command [list $object xview]
$object.scx configure -command {}

$object.g axis configure x -scrollcommand [list $object xset]
$object.g axis configure x -scrollcommand {}

}

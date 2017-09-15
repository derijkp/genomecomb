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
			-textvariable [privatevar $object graphsettings($type)] -command [list Classy::todo $object redraw]
		pack $object.b.$type -side left
	}
	button $object.b.conf  -text "Configure" -command [list Classy::todo $object graphconfigure]
	pack $object.b.conf -side left
	Classy::NumEntry $object.b.linew -width 2 -label lw -orient horizontal \
		-textvariable [privatevar $object graphsettings(linew)] -command [list Classy::todo $object rearrange]
	pack $object.b.linew -side left
	button $object.b.max  -text "Maxview" -command [subst {
		setprivate $object graphsettings(xmin) {}
		setprivate $object graphsettings(xmax) {}
		setprivate $object graphsettings(ymin) {}
		setprivate $object graphsettings(ymax) {}
		$object redraw
	}]
	pack $object.b.max -side left
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

if {[info commands ori.PushZoom] eq ""} {
	catch {PushZoom}
	rename PushZoom ori.PushZoom
}

proc PushZoom {graph} {
	ori.PushZoom $graph
	[winfo parent $graph] _fillregion
	[winfo parent $graph] _setvars
}

if {[info commands ori.PopZoom] eq ""} {
	catch {PopZoom}
	rename PopZoom ori.PopZoom
}

proc PopZoom {graph} {
	ori.PopZoom $graph
	[winfo parent $graph] _fillregion
	[winfo parent $graph] _setvars
}

scrolledgraph method paste {} {
	private $object graphsettings
	set graphsettings(xrange) [::tk::GetSelection $object PRIMARY]
	Classy::todo $object redraw
}

scrolledgraph method start {} {
	private $object graphsettings
	set graphsettings(xmin) {}
	set graphsettings(xmax) {}
	set graphsettings(ymin) {}
	set graphsettings(ymax) {}
	set graphsettings(xlog) 0
	set graphsettings(ylog) 0
	set graphsettings(linew) 0
	$object.g axis bind x <1> [list $object graphconfigure x]
	$object.g axis bind y <1> [list $object graphconfigure y]
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

scrolledgraph method _retable {} {
	private $object coldata header
	set len [[lindex $coldata([lindex $header 0]) 1] length]
	set lastpos [expr {$len-1}]
	set result [list_fill $len {}]
	foreach col $header {
		if {[lindex $coldata($col) 0] eq "labels"} {
			set temp [lindex $coldata($col) 2]
		} else {
			set v [lindex $coldata($col) 1]
			set temp [$v range 0 $lastpos]
		}
		set pos 0
		foreach v $temp {
			lset result $pos [list {*}[lindex $result $pos] $v]
			incr pos
		}
	}
	return $result
}

scrolledgraph method sort {sortcol} {
	private $object coldata header
	set type [lindex $coldata($sortcol) 0]
	set sortv [lindex $coldata($sortcol) 1]
	set lastpos [expr {[$sortv length]-1}]
	if {$type eq "labels"} {
		set list [lindex $coldata($sortcol) 2]
	} else {
		set list [$sortv range 0 $lastpos]
	}
	set slist [ssort -natural $list]
	if {$list eq $slist} return
	set neworder [list_cor $list $slist]
	set others {}
	foreach col $header {
		if {[lindex $coldata($col) 0] eq "labels"} {
			set list [lindex $coldata($col) 2]
			lset coldata($col) 2 [list_sub $list $neworder]
		} else {
			lappend others [lindex $coldata($col) 1]
		}
	}
	if {[llength $others]} {
		vector create $object.sortvector
		$object.sortvector set [list_cor $slist $list]
		$object.sortvector sort {*}$others
	}
}

scrolledgraph method gradientstyle {basecolor {wv {}} {null 0} {symbol plus}} {
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
				$object.g pen create $color -fill $colorcode -outline $colorcode -outlinewidth 1 -pixels $num -symbol $symbol
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
				$object.g pen create $color -fill $colorcode -outline $colorcode -outlinewidth 1 -pixels $num -symbol $symbol
			} else {
				# $object.g pen configure $color -fill $colorcode -outline $colorcode -outlinewidth 1 -pixels 2 -symbol $symbol
			}
			incr num
		}
	}
#	set style {}
#	set step [expr {($max-$min)/10.0}]
#	set num 1
#	set cur $min
#	for {set num 1} {$num <= 10} {incr num} {
#		set next [expr {$cur+$step}]
#		lappend style [list $basecolor-$num $cur $next]
#		set cur $next
#		incr num
#	}
	set style {}
	if {$wv ne ""} {
		$wv dup ::tv
		::tv sort
		set min $::tv(0)
		set max $::tv(end)
		set len [::tv length]
		set step [expr {round(($len)/10.0)}]
		set num 1
		set cur 0
		for {set num 1} {$num <= 10} {incr num} {
			set next [expr {$cur+$step}]
			if {$next >= $len} {set next [expr {$len-1}]}
			lappend style [list $basecolor-$num $::tv($cur) $::tv($next)]
			set cur $next
		}
		lset style end end [expr {$max+1}]
	} else {
		lappend style [list $basecolor -Inf Inf]
	}
	lappend style [list $basecolor-0 $null $null]
	# join $style \n
	return $style
}

scrolledgraph method _configureevent {} {
	private $object graphsettings
	event generate $object.g <Configure>
	set graphsettings(status) ""
	set graphsettings(cancel) 0
}

scrolledgraph method getlabel {win num} {
	private $object labels
	if {[isint $num]} {
		return [lindex $labels $num]
	} else {
		return ""
	}
}

scrolledgraph method add {table {xtitle {}} {ytitle {}}} {
#set ::table $table
#putsvars object
	private $object data colors labels graphsettings header coldata
	set showregion 0
	# we will add to existing (data(entries) already there)
	set header [list_shift table]
	unset -nocomplain coldata
	set vnum 0
	foreach colname $header {
		vector create ::$object.col$vnum
		set col [list_subindex $table $vnum]
		if {[catch {::$object.col$vnum set $col}]} {
			set tempcol [list_regsub -all {^(NaN||[?-])$} $col -Inf]
			if {[catch {::$object.col$vnum set $tempcol}]} {
				set coldata($colname) [list labels ::$object.col$vnum $col]
				::$object.col$vnum set [list_cor $col $col]
			} else {
				set coldata($colname) [list vectsub ::$object.col$vnum $col]
			}
		} else {
			set coldata($colname) [list vect ::$object.col$vnum {}]
		}
		incr vnum
	}
	set graphsettings(xmin) {}
	set graphsettings(xmax) {}
	set graphsettings(ymin) {}
	set graphsettings(ymax) {}
	Classy::todo $object rearrange
}

scrolledgraph method addscatter {table xcol ycol datacols} {
	private $object xv yv wv data colors labels graphsettings
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
	set graphsettings(xlog) [get graphsettings(xlog) 0]
	set graphsettings(ylog) [get graphsettings(ylog) 0]
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

scrolledgraph method rearrange_clear {args} {
	foreach el [list_remove [$object.g element names] legend] {
		$object.g element delete $el
	}
	foreach el [list_remove [$object.g marker names] legend] {
		$object.g marker delete $el
	}
}

# change display according to settings (x, y, ...)
scrolledgraph method rearrange {args} {
	private $object data colors labels graphsettings header coldata xcol ycols fcols wcols scols
	if {[get graphsettings(type)] in {roc pr}} {
		$object rearrange_roc
		return
	}
	set showregion 0
	set linew [get graphsettings(linew) 0]
	$object rearrange_clear
	set numcols {}
	foreach field $header {
		if {[lindex $coldata($field) 0] ne "labels"} {
			lappend numcols $field
		}
	}
	set xpos [lsearch $header [get xcol ""]]
	if {$xpos == -1} {
		if {[llength $numcols] > 1} {
			set xcol [lindex $numcols 0]
		}
	}
	set xpos [lsearch $header $xcol]
	set yposs [list_remove [list_cor $header [get ycols ""]] -1 $xpos]
	if {![llength $yposs]} {
		set ycols [list_remove $numcols $xcol]
	} else {
		set ycols [list_sub $header $yposs]
	}
	set fposs [list_remove [list_cor $header [get fcols ""]] -1 $xpos]
	if {![llength $fposs]} {
		set fcols [list_remove $header $xcol {*}$ycols {*}$numcols]
	} else {
		set fcols [list_sub $header $fposs]
	}
	set wposs [list_remove [list_cor $header [get wcols ""]] -1]
	set wcols [list_sub $header $wposs]
	set sposs [list_remove [list_cor $header [get scols ""]] -1]
	set scols [list_sub $header $sposs]
	# sort according to x or scols
	set others {}
	if {[llength $scols]} {
		$object sort $scols
	} else {
		$object sort $xcol
	}
	# create axis
	set xs [lindex $coldata($xcol) 1]
	$object.g axis configure y -title [join $ycols ,]
	if {[lindex $coldata($xcol) 0] eq "labels"} {
		set labels [lindex $coldata($xcol) end]
		$object.g axis configure x -command [list $object getlabel] -title $xcol
	} else {
		unset -nocomplain labels
		$object.g axis configure x -command {} -title $xcol
	}
	set xmin [get ${xs}(min)]
	set xmax [get ${xs}(max)]
#	if {[info exists labels]} {
#		set amin [expr {$xmin - 1}]
#		set amax [expr {$xmax + 1}]
#	}
	::$object.ends.x set [list $xmin $xmax]
	if {$graphsettings(xmax) == Inf} {
		set graphsettings(xmax) [::$object.ends.x index 1]
	}
	set vnum 0
	set pos 0
	set ymin {}
	set ymax {}
	set todo {}
	if {[llength $wcols] == 0} {
		foreach ycol $ycols {
			lappend todo $ycol {}
		}
	} elseif {[llength $ycols] == 1} {
		set ycol [lindex $ycols 0]
		foreach wcol $wcols {
			lappend todo $ycol $wcol
		}
	} elseif {[llength $wcols] == 1} {
		set wcol [lindex $wcols 0]
		foreach ycol $ycols {
			lappend todo $ycol $wcol
		}
	} else {
		foreach ycol $ycols wcol $wcols {
			lappend todo $ycol $wcol
		}
	}
	if {$fcols eq ""} {
		set fcolvals {{}}
		set dcolors [distinctColors [llength $todo]]
	} else {
		set fcol [lindex $fcols 0]
		if {[lindex $coldata($fcol) 0] eq "labels"} {
			set temp [lindex $coldata($fcol) end]
		} else {
			set temp [[lindex $coldata($fcol) 1] range 0 end]
		}
		foreach fcol [lrange $fcols 1 end] {
			set newtemp {}
			if {[lindex $coldata($fcol) 0] eq "labels"} {
				set vals [lindex $coldata($fcol) end]
			} else {
				set vals [[lindex $coldata($fcol) 1] range 0 end]
			}
			foreach f $temp v $vals {
				lappend f $v
				lappend newtemp $f
			}
			set temp $newtemp
		}
		set allfcolvals $temp
		set fcolvals [list_remdup $allfcolvals]
		set dcolors [distinctColors [expr {[llength $todo]*[llength $fcolvals]}]]
	}
	set symbols {cross circle square diamond plus splus scross triangle}
	foreach {field wcol} $todo {
		foreach fcolval $fcolvals {
			incr pos
			set name $field
			if {[llength $wcol]} {
				append name -$wcol
			}
			if {[llength $fcols]} {
				append name -$fcolval
			}
			set ys [lindex $coldata($field) 1]
			if {[llength $fcols]} {
				# split out xs, ys and ws for this value of fill
				set poss [list_find $allfcolvals $fcolval]
				set uxs [vector create ::$object.fxs$vnum]
				$uxs set [list_sub [$xs range 0 end] $poss]
				set uys [vector create ::$object.fys$vnum]
				$uys set [list_sub [$ys range 0 end] $poss]
				if {[llength $wcol]} {
					set wv [lindex $coldata($wcol) 1]
					set uwv [vector create ::$object.fwv$vnum]
					$uwv set [list_sub [$wv range 0 end] $poss]
				}
				set elementlabel $fcolval
			} else {
				set uxs $xs
				set uys $ys
				if {[llength $wcol]} {
					set wv [lindex $coldata($wcol) 1]
				}
				set elementlabel $field
			}			
			set curmin [get ${uys}(min)]
			set ymin [min $ymin $curmin]
			set curmax [get ${uys}(max)]
			set ymax [max $ymax $curmax]
			set basecolor [lindex $dcolors $vnum]
			if {$basecolor eq ""} {set basecolor blue}
			lappend data(entries) $name
			set data($name,vnum) $vnum
			set data($name,lstart) 0
			set data($name,lend) 0
			set data($name,header) $header
			$object.g element create e$name -label $elementlabel -xdata $uxs -ydata $uys -symbol none
			set color $basecolor
			# $object.g element configure e$name -linewidth $linew -fill $color -outlinewidth 1 -pixels 1 -symbol circle
			set symbol plus
			$object.g element configure e$name -linewidth $linew -color $basecolor -outline $basecolor -outlinewidth 1 -pixels 2 -symbol $symbol
			if {[llength $wcol]} {
				set wmin [get ${uwv}(min)]
				set wmax [get ${uwv}(max)]
				set style [$object gradientstyle $basecolor $uwv 0]
				$object.g element configure e$name -weight $uwv -styles $style
			}
			if {$showregion} {
				$object.g element configure e$name -linewidth 10 -trace decreasing
			}
			set data($name,color) $color
			$object.g legend bind $name <1> [list $object elconf $name]
			incr vnum
		}
	}
	if {![isdouble $ymin]} {set ymin 0}
	if {![isdouble $ymax]} {set ymax 1}
	::$object.ends.y set [list $ymin $ymax]
	Classy::todo $object redraw
}

scrolledgraph method rearrange_roc {args} {
	private $object data colors labels graphsettings header coldata xcol wcols scols poscol negcol roccol
	if {$graphsettings(type) eq "roc"} {set roc 1} else {set roc 0}
	set showregion 0
	set linew [get graphsettings(linew) 1]
	$object rearrange_clear
	set pospos [lsearch $header [get poscol ""]]
	set negpos [lsearch $header [get negcol ""]]
	set yposs [list_remove [list_cor $header [get roccol ""]] -1]
	if {![llength $yposs]} {
		set roccol {}
		foreach field [list_sub $header -exclude [list $pospos $negpos]] {
			if {[lindex $coldata($field) 0] ne "label"} {
				lappend roccol $field
			}
		}
	} else {
		set roccol [list_sub $header $yposs]
	}
	# create axis
	::$object.ends.x set [list 0 100]
	::$object.ends.y set [list 0 100]
	if {$roc} {
		$object.g axis configure x -title "1-Specificity (%)"
		$object.g axis configure y -title "Sensitivity (%)"
		$object.g element create diagonal -label "" -xdata ::$object.ends.x -ydata ::$object.ends.y -symbol none -color gray
	} else {
		$object.g axis configure x -title "Recal (%)"
		$object.g axis configure y -title "Precision (%)"
	}
	set vnum 0
	set pos 0
	set symbols {cross circle square diamond plus splus scross triangle}
	# set dcolors {blue red gray orange yellow green violet}
	set dcolors [distinctColors [expr {[llength $header]-2}]]
	# calculate roc curve
	set negatives [lindex $coldata($negcol) 1]
	set positives [lindex $coldata($poscol) 1]
	set realpos [vector expr {sum($positives)}]
	set realneg [vector expr {sum($negatives)}]
	set todo {}
	set vnum 1
	foreach field $roccol {
		set ys [lindex $coldata($field) 1]
		$object sort $field
		set fn 0
		set tn 0
		set fp 0
		set tp 0
		set xlist {}
		set ylist {}
		set prev [$ys index 0]
		foreach p [$positives range 0 end] n [$negatives range 0 end] y [$ys range 0 end] {
			if {$y != $prev} {
				lappend ylist $yval
				lappend xlist $xval
				if {$graphsettings(roclabels)} {
					$object.g marker create text -coords [list $xval $yval] -text $prev -fill {} -anchor nw -yoffset -5
				}
				set prev $y
			}
			if {$graphsettings(rocabove)} {
				set fn [expr {$fn + $p}]
				set tn [expr {$tn + $n}]
				if {$roc} {
					# true positive rate
					set yval [expr {100.0*($realpos - $fn)/$realpos}]
					# false positive rate
					set xval [expr {100.0*($realneg - $tn)/$realneg}]
				} else {
					# precision
					if {[catch {
						set yval [expr {100.0*($realpos - $fn)/($realneg - $tn + $realpos - $fn)}]
					}]} continue
					# recall (= true positive rate)
					set xval [expr {100.0*($realpos - $fn)/$realpos}]
				}
			} else {
				set tp [expr {$tp + $p}]
				set fp [expr {$fp + $n}]
				if {$roc} {
					# true positive rate
					set yval [expr {100.0*$tp/$realpos}]
					# false positive rate
					set xval [expr {100.0*$fp/$realneg}]
				} else {
					# precision
					set yval [expr {100.0*$tp/($tp+$fp)}]
					# recall (= true positive rate)
					set xval [expr {100.0*$tp/$realpos}]
				}
			}
		}
		lappend ylist $yval
		lappend xlist $xval
		if {$graphsettings(roclabels)} {
			$object.g marker create text -coords [list $xval $yval] -text $y -fill {} -anchor nw -yoffset -5
		}
		set name $field
		set yvalv [vector create ::$object.cor.$name.yval]
		$yvalv set $ylist
		set xvalv [vector create ::$object.cor.$name.xval]
		$xvalv set $xlist
		set ys [lindex $coldata($field) 1]
		set basecolor [lindex $dcolors $vnum]
		if {$basecolor eq ""} {set basecolor blue}
		lappend data(entries) $name
		set data($name,vnum) $vnum
		set data($name,lstart) 0
		set data($name,lend) 0
		set data($name,header) $header
		$object.g element create e$name -label $name -xdata $xvalv -ydata $yvalv -symbol none
		# $object.g element configure e$name -weights ::$ws
		set color $basecolor
		$object.g element configure e$name -linewidth $linew -color $basecolor -outline $basecolor -outlinewidth 1 -pixels 2 -symbol plus
		set data($name,color) $color
		$object.g legend bind $name <1> [list $object elconf $name]
		incr vnum
	}
	set graphsettings(xmin) 0
	set graphsettings(xmax) 100
	set graphsettings(ymin) 0
	set graphsettings(ymax) 100
	Classy::todo $object redraw
}

scrolledgraph method redrawsettings {args} {
	private $object graphsettings
	set xmin [::$object.ends.x index 0]
	set xmax [::$object.ends.x index 1]
	set ymin [::$object.ends.y index 0]
	set ymax [::$object.ends.y index 1]
	if {[inlist $args resetx]} {
		set graphsettings(xmin) $xmin
		set graphsettings(xmax) $xmax
	} else {
		if {![info exists graphsettings(xmin)] || $graphsettings(xmin) eq "" || $graphsettings(xmin) < $xmin} {set graphsettings(xmin) $xmin}
		if {![info exists graphsettings(xmax)] || $graphsettings(xmax) eq "" || $graphsettings(xmax) > $xmax} {set graphsettings(xmax) $xmax}
		if {$graphsettings(xmax) <= $graphsettings(xmin)} {set graphsettings(xmax) $xmax}
		if {$graphsettings(xmax) <= $graphsettings(xmin)} {set graphsettings(xmax) [expr {$xmin+1}]}
	}
	if {[inlist $args resety]} {
		set graphsettings(ymin) $ymin
		set graphsettings(ymax) $ymax
	} else {
		if {![info exists graphsettings(ymin)] || $graphsettings(ymin) eq "" || $graphsettings(ymin) < $ymin} {set graphsettings(ymin) $ymin}
		if {![info exists graphsettings(ymax)] || $graphsettings(ymax) eq "" || $graphsettings(ymax) > $ymax} {set graphsettings(ymax) $ymax}
		if {$graphsettings(ymax) <= $graphsettings(ymin)} {set graphsettings(ymax) $ymax}
		if {$graphsettings(ymax) <= $graphsettings(ymin)} {set graphsettings(ymax) [expr {$ymin+1}]}
	}
}

scrolledgraph method redraw {args} {
puts "----------redraw $object----------"
	private $object graphsettings
	set w $object.g
	Classy::canceltodo $object redraw
	# $object _xrange
	$object redrawsettings
	$w axis configure x -min $graphsettings(xmin) -max $graphsettings(xmax) -logscale $graphsettings(xlog)
	$w axis configure y -min $graphsettings(ymin) -max $graphsettings(ymax) -logscale $graphsettings(ylog)
	if {$graphsettings(xlog)} {
		$object.scx configure -command {}
		$object.g axis configure x -scrollcommand {}
		grid forget $object.scx
	} else {
		$object.scx configure -command [list $object xview]
		$object.g axis configure x -scrollcommand [list $object xset]
		grid $object.scx -row 2 -sticky nwse
	}
	if {$graphsettings(ylog)} {
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
	private $object graphsettings
	set w $object.g
	set xmin [expr {round([$w axis cget x -min])}]
	set xmax [expr {round([$w axis cget x -max])}]
	if {[isdouble $xmin] && [isdouble $xmax]} {
		set graphsettings(xmin) $xmin
		set graphsettings(xmax) $xmax
	}
	set ymin [$w axis cget y -min]
	set ymax [$w axis cget y -max]
	catch {set ymin [format %.0f $ymin]}
	catch {set ymax [format %.0f $ymax]}
	set graphsettings(ymin) $ymin
	set graphsettings(ymax) $ymax
	Classy::todo $object _configureevent
}

scrolledgraph method _fillregion {args} {
	private $object graphsettings
	set w $object.g
	foreach type {xmin xmax ymin ymax} {
		set val [$w axis cget [string index $type 0] -[string range $type 1 end]]
		if {[isdouble $val]} {
			set val [format %.1f $val]
		}
		set graphsettings($type) $val
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
	private $object data conf graphsettings
	if {[info exists conf(entry)]} {
		set name $conf(entry)
		set data($name,color) [get conf(color) gray]
		set style [$object gradientstyle $data($name,color)]
		set color [get colors($data($name,color)-0) $data($name,color)]
		$object.g element configure e$name -linewidth $conf(linewidth) -fill $color -outline $color -color $color \
			-outlinewidth 1 -pixels 2 -symbol $conf(symbol) \
			-styles $style
	}
	# x axis settings
	if {[isint $graphsettings(xrotate)]} {
		$object.g axis configure x -rotate $graphsettings(xrotate)
	}
	if {[isint $graphsettings(xstep)]} {
		$object.g axis configure x -stepsize $graphsettings(xstep)
	}
	if {[isint $graphsettings(xsubdivisions)]} {
		$object.g axis configure x -subdivisions $graphsettings(xsubdivisions)
	}
	if {[get graphsettings(xfont) ""] eq ""} {set graphsettings(xfont) {helvetica -12}}
	$object.g axis configure x -tickfont $graphsettings(xfont)
	# y axis settings
	if {[isint $graphsettings(yrotate)]} {
		$object.g axis configure y -rotate $graphsettings(yrotate)
	}
	if {[isint $graphsettings(ystep)]} {
		$object.g axis configure y -stepsize $graphsettings(ystep)
	}
	if {[isint $graphsettings(ysubdivisions)]} {
		$object.g axis configure y -subdivisions $graphsettings(ysubdivisions)
	}
	if {[get graphsettings(yfont) ""] eq ""} {set graphsettings(yfont) {helvetica -12}}
	$object.g axis configure y -tickfont $graphsettings(yfont)
	if {[get graphsettings(legendfont) ""] eq ""} {set graphsettings(legendfont) {helvetica -12}}
	$object.g legend configure -font $graphsettings(legendfont)
	if {$graphsettings(legend)} {
		$object.g legend configure -hide yes
	} else {
		$object.g legend configure -hide no
	}
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
	private $object data graphsettings
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
	private $object data graphsettings
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
	private $object graphsettings options
	if {![llength $args]} {
		destroy $object.ps
		set file [file root [file tail [get options(-datafile) plot]]]__[$object.g axis cget y -title]_vs_[$object.g axis cget x -title]__$graphsettings(xmin)-$graphsettings(xmax).ps
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
	catch {gzclose $f}
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

scrolledgraph method changeX {} {
	private $object xcol header graphsettings
	set xcol [Classy::select "Select X axis" $header -initialvalue $xcol]
	set graphsettings(xmin) {}
	set graphsettings(xmax) {}
	$object rearrange
}

scrolledgraph method changeYs {} {
	private $object ycols header graphsettings
	set ycols [Classy::select "Select Y axis" $header -selectmode multiple -initialvalue $ycols]
	set graphsettings(ymin) {}
	set graphsettings(ymax) {}
	$object rearrange
}

scrolledgraph method changeWs {} {
	private $object wcols header
	set wcols [Classy::select "Select weight columns" $header -selectmode multiple -initialvalue $wcols]
	$object rearrange
}

scrolledgraph method changeFs {} {
	private $object fcols header
	set wcols [Classy::select "Select Element columns" $header -selectmode multiple -initialvalue $fcols]
	$object rearrange
}

scrolledgraph method changesort {} {
	private $object scols header
	set scols [get scols ""]
	set scols [Classy::select "Select sort column" $header -initialvalue $scols]
	$object rearrange
}

scrolledgraph method graphconfigure {{axis x}} {
	private $object graphsettings header poscol negcol roccol xcol ycols fcols scols wcols
	destroy $object.graphconfigure
	Classy::Dialog $object.graphconfigure -title "configure axis"
	$object.graphconfigure tab "X Axis"
	$object.graphconfigure option listbox "X column" [privatevar $object xcol] [privatevar $object header] -selectmode single \
		-browsecommand [subst {
			invoke args {
				setprivate $object graphsettings(xmin) {}
				setprivate $object graphsettings(xmax) {}
			}}]
	$object.graphconfigure option check "Scale" [privatevar $object graphsettings(xlog)] "use log scale" -- -command [list Classy::todo $object rearrange]
	$object.graphconfigure option numentry "Rotate" [privatevar $object graphsettings(xrotate)] -command [list Classy::todo $object reconf]
	$object.graphconfigure option numentry "Major stepsize" [privatevar $object graphsettings(xstep)] -command [list Classy::todo $object reconf]
	$object.graphconfigure option numentry "Subdivisions" [privatevar $object graphsettings(xsubdivisions)] -command [list Classy::todo $object reconf]
	#
	$object.graphconfigure tab "Y Axis"
	$object.graphconfigure option listbox "Y columns" [privatevar $object ycols] [privatevar $object header] -selectmode extended \
		-browsecommand [subst {
			invoke args {
				setprivate $object graphsettings(ymin) {}
				setprivate $object graphsettings(ymax) {}
			}}]
	$object.graphconfigure option check "Scale" [privatevar $object graphsettings(ylog)] "use log scale" -- -command [list Classy::todo $object rearrange]
	$object.graphconfigure option numentry "Rotate" [privatevar $object graphsettings(yrotate)] -command [list Classy::todo $object reconf]
	$object.graphconfigure option numentry "Major stepsize" [privatevar $object graphsettings(ystep)] -command [list Classy::todo $object reconf]
	$object.graphconfigure option numentry "Subdivisions" [privatevar $object graphsettings(ysubdivisions)] -command [list Classy::todo $object reconf]
	#
	$object.graphconfigure tab "Elements"
	$object.graphconfigure option listbox "Element columns" [privatevar $object fcols] [privatevar $object header] -selectmode multiple
	#
	$object.graphconfigure tab "Weights"
	$object.graphconfigure option listbox "Weight columns" [privatevar $object wcols] [privatevar $object header] -selectmode multiple
	#
	$object.graphconfigure tab "Sort"
	$object.graphconfigure option listbox "Sort column" [privatevar $object scols] [privatevar $object header] -selectmode single
	#
	$object.graphconfigure tab "Display"
	if {[get graphsettings(xfont) ""] eq ""} {set graphsettings(xfont) {helvetica -12}}
	$object.graphconfigure option font "X tick font" [privatevar $object graphsettings(xfont)] -command [list Classy::todo $object reconf]
	if {[get graphsettings(yfont) ""] eq ""} {set graphsettings(yfont) {helvetica -12}}
	$object.graphconfigure option font "Y tick font" [privatevar $object graphsettings(yfont)] -command [list Classy::todo $object reconf]
	if {[get graphsettings(legendfont) ""] eq ""} {set graphsettings(legendfont) {helvetica -12}}
	$object.graphconfigure option font "Legend font" [privatevar $object graphsettings(legendfont)] -command [list Classy::todo $object reconf]
	if {[get graphsettings(legend) ""] eq ""} {set graphsettings(legend) 1}
	$object.graphconfigure option check "Legend" [privatevar $object graphsettings(legend)] "Show legend" -- -command [list Classy::todo $object reconf]
	#
	$object.graphconfigure tab "ROC"
	if {[get graphsettings(type) ""] eq "" || ![inlist {graph roc pr} $graphsettings(type)]} {set graphsettings(type) graph}
	$object.graphconfigure option radio "Graph type" [privatevar $object graphsettings(type)] "Roc curve" roc "Precision-Recal curve" pr "Normal graph" graph -- -command [list Classy::todo $object rearrange]
	if {[get poscol ""] eq "" || ![inlist $header $poscol]} {set poscol [lindex $header end-1]}
	$object.graphconfigure option listbox "column with number of positives" [privatevar $object poscol] [privatevar $object header] -selectmode single \
		-browsecommand [subst {
			invoke args {
				setprivate $object graphsettings(xmin) {}
				setprivate $object graphsettings(xmax) {}
			}}]
	if {[get negcol ""] eq "" || ![inlist $header $negcol]} {set negcol [lindex $header end]}
	$object.graphconfigure option listbox "column with number of negatives" [privatevar $object negcol] [privatevar $object header] -selectmode single \
		-browsecommand [subst {
			invoke args {
				setprivate $object graphsettings(xmin) {}
				setprivate $object graphsettings(xmax) {}
			}}]
	if {[get roccol ""] eq "" || ![llength [list_common $header $roccol]]} {set roccol [lrange $header 1 end-2]}
	$object.graphconfigure option listbox "Cutoff columns" [privatevar $object roccol] [privatevar $object header] -selectmode multiple \
		-browsecommand [subst {
			invoke args {
				setprivate $object graphsettings(xmin) {}
				setprivate $object graphsettings(xmax) {}
			}}]
	#
	if {[get graphsettings(rocabove) ""] eq ""} {set graphsettings(rocabove) 1}
	$object.graphconfigure option check "Direction" [privatevar $object graphsettings(rocabove)] "Values above cutoff test positive" -- -command [list Classy::todo $object rearrange]
	if {[get graphsettings(roclabels) ""] eq ""} {set graphsettings(roclabels) 1}
	$object.graphconfigure option check "Labels" [privatevar $object graphsettings(roclabels)] "Show cutoffs as labels" -- -command [list Classy::todo $object rearrange]
	$object.graphconfigure add go Go [list $object rearrange]
	if {$axis eq "x"} {
		$object.graphconfigure tab "X Axis"
	}
	if {$axis eq "y"} {
		$object.graphconfigure tab "Y Axis"
	}
	update
	set xcol $xcol
	set ycols $ycols
	set fcols $fcols
	set wcols $wcols
	set scols $scols
}

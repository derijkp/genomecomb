Classy::Dialog subclass fieldsdialog

fieldsdialog method init args {
	super init
	# Create windows
	Classy::Paned $object.options.paned
	grid $object.options.paned -row 0 -column 1 -rowspan 2 -sticky nesw
	frame $object.options.buttons  \
		-borderwidth 0 \
		-relief groove \
		-height 10 \
		-width 10
	grid $object.options.buttons -row 0 -column 2 -rowspan 2 -sticky nesw
	button $object.options.buttons.addfields \
		-text -->
	grid $object.options.buttons.addfields -row 1 -column 0 -sticky nesw
	button $object.options.buttons.delfields \
		-text <--
	grid $object.options.buttons.delfields -row 2 -column 0 -sticky nesw
	button $object.options.buttons.button3 \
		-text Basic
	grid $object.options.buttons.button3 -row 4 -column 0 -sticky nesw
	button $object.options.buttons.button4 \
		-text Clear
	grid $object.options.buttons.button4 -row 5 -column 0 -sticky nesw
	button $object.options.buttons.button5 \
		-text All
	grid $object.options.buttons.button5 -row 6 -column 0 -sticky nesw
	button $object.options.buttons.button6 \
		-text Up
	grid $object.options.buttons.button6 -row 9 -column 0 -sticky nesw
	button $object.options.buttons.button7 \
		-text Down
	grid $object.options.buttons.button7 -row 10 -column 0 -sticky nesw
	frame $object.options.buttons.frame1  \
		-borderwidth 2 \
		-height 10 \
		-width 10
	grid $object.options.buttons.frame1 -row 3 -column 0 -sticky nesw
	frame $object.options.buttons.frame2  \
		-borderwidth 2 \
		-height 10 \
		-width 10
	grid $object.options.buttons.frame2 -row 8 -column 0 -sticky nesw
	button $object.options.buttons.button8 \
		-text Edit
	grid $object.options.buttons.button8 -row 7 -column 0 -sticky nesw
	grid columnconfigure $object.options.buttons 0 -uniform {}
	grid rowconfigure $object.options.buttons 0 -uniform {} -weight 1
	grid rowconfigure $object.options.buttons 1 -uniform {}
	grid rowconfigure $object.options.buttons 2 -uniform {}
	grid rowconfigure $object.options.buttons 3 -uniform {}
	grid rowconfigure $object.options.buttons 4 -uniform {}
	grid rowconfigure $object.options.buttons 5 -uniform {}
	grid rowconfigure $object.options.buttons 6 -uniform {}
	grid rowconfigure $object.options.buttons 7 -uniform {}
	grid rowconfigure $object.options.buttons 8 -uniform {}
	grid rowconfigure $object.options.buttons 9 -uniform {}
	grid rowconfigure $object.options.buttons 10 -uniform {}
	grid rowconfigure $object.options.buttons 11 -uniform {} -weight 1
	Classy::ListBox $object.options.selfields  \
		-height 4
	grid $object.options.selfields -row 0 -column 0 -rowspan 2 -sticky nesw
	Classy::ListBox $object.options.fields  \
		-height 4 \
		-width 10
	grid $object.options.fields -row 1 -column 3 -columnspan 2 -sticky nesw
	Classy::Entry $object.options.calc \
		-label {Calculated Field} \
		-width 4 \
		-combo 20
	grid $object.options.calc -row 0 -column 3 -sticky nesw
	button $object.options.button1 \
		-text Builder
	grid $object.options.button1 -row 0 -column 4 -sticky nesw
	grid columnconfigure $object.options 0 -uniform {}
	grid columnconfigure $object.options 1 -uniform {}
	grid columnconfigure $object.options 2 -uniform {}
	grid columnconfigure $object.options 3 -uniform {} -weight 1
	grid columnconfigure $object.options 4 -uniform {}
	grid columnconfigure $object.options 5 -uniform {}
	grid columnconfigure $object.options 6 -uniform {}
	grid rowconfigure $object.options 0 -uniform {}
	grid rowconfigure $object.options 1 -uniform {} -weight 1
	grid rowconfigure $object.options 2 -uniform {}
	grid rowconfigure $object.options 3 -uniform {}

	if {"$args" == "___Classy::Builder__create"} {return $object}
	# Parse this
	$object configure  \
		-title Toplevel
	$object.options.paned configure \
		-window [varsubst object {$object.options.selfields}]
	$object.options.buttons.addfields configure \
		-command [varsubst object {$object addfields}]
	$object.options.buttons.delfields configure \
		-command [varsubst object {$object delfields}]
	$object.options.buttons.button3 configure \
		-command [varsubst object {$object basicfields}]
	$object.options.buttons.button4 configure \
		-command [varsubst object {$object clear}]
	$object.options.buttons.button5 configure \
		-command [varsubst object {$object allfields}]
	$object.options.buttons.button6 configure \
		-command [varsubst object {$object fieldsup}]
	$object.options.buttons.button7 configure \
		-command [varsubst object {$object fieldsdown}]
	$object.options.buttons.button8 configure \
		-command [varsubst object {$object editfields}]
	$object.options.selfields configure \
		-filtervariable [privatevar $object filterselfields]
	$object.options.fields configure \
		-filtervariable [privatevar $object filterfields]
	$object.options.calc configure \
		-command [varsubst object {$object addcalc}]
	$object.options.button1 configure \
		-command [varsubst object {$object fieldbuilder}]
	$object add go Go [varsubst object {$object go}] default
	$object persistent set 
	# Configure initial arguments
	if {"$args" != ""} {eval $object configure $args}
# ClassyTk Finalise
$object start
	return $object
}

fieldsdialog addoption -command {command Command {}} {}

fieldsdialog method start {} {
	setprivate $object tfields {}
	$object.options.selfields configure -selectmode extended -exportselection 0
	$object.options.fields configure -selectmode extended -exportselection 0
	bind $object.options.fields <Up> "[list $object fieldsup]; break"
	bind $object.options.fields <Down> "[list $object fieldsdown]; break"
	$object.options.fields configure -browsecommand [list $object setcalc]
	$object.options.fields activate end
}

fieldsdialog method go {} {
	private $object options fields
	{*}$options(-command) $fields
}

fieldsdialog method redraw {} {
	private $object tfields fields
	set selection [$object.options.fields get]
	set active [$object.options.fields index active]
	$object.options.fields configure -content $fields
	set showtfields [list_lremove $tfields $fields]
#	set search [$object.options.search get]
#	if {$search ne "" && $search ne "*"} {
#		set poss [list_find -regexp $showtfields $search]
#		set showtfields [list_sub $showtfields $poss]
#	}
	$object.options.selfields configure -content $showtfields
	$object selectfields $selection
	$object.options.fields.list activate $active
}

fieldsdialog method tfields {fields} {
	private $object tfields
	set tfields $fields
	$object redraw
}

fieldsdialog method fields {curfields} {
	private $object fields
	set fields [list_remove $curfields {}]
	$object redraw
	$object.options.fields selection set end
	$object.options.fields activate end
}

fieldsdialog method selectfields {selfields} {
	private $object fields
	set poss [list_cor $fields $selfields]
	foreach pos $poss {
		$object.options.fields selection set $pos
	}
}

fieldsdialog method fieldadd {movefields} {
	private $object fields
	set pos [$object.options.fields index active]
	incr pos
	set fields [linsert $fields $pos {*}$movefields]
	$object redraw
	$object selectfields $movefields
	$object.options.fields activate [expr {$pos+[llength $movefields]-1}]
}

fieldsdialog method addfields {} {
	private $object fields tfields
	set movefields [$object.options.selfields get]
	if {![llength $movefields]} {
		set movefields [$object.options.selfields get 0 end]
	}
	$object fieldadd $movefields
}

fieldsdialog method delfields {} {
	private $object fields tfields
	set movefields [$object.options.fields get]
	if {![llength $movefields]} {
		set movefields [$object.options.fields get 0 end]
	}
	set fields [list_lremove $fields $movefields]
	$object redraw
}

fieldsdialog method clear {} {
	private $object fields
	set fields {}
	$object redraw
}

fieldsdialog method allfields {} {
	private $object fields tfields
	set movefields [$object.options.selfields get 0 end]
	$object fieldadd $movefields
}

fieldsdialog method basicfields {} {
	private $object tfields fields
	set poss [tsv_basicfields $tfields 6 0]
	set basicfields [list_remove [list_sub $tfields $poss] {}]
	set fields [list_union $basicfields $fields]
	$object redraw
	$object selectfields $basicfields
	$object.options.fields activate [llength $basicfields]
}

fieldsdialog method fieldsup {} {
	private $object tfields fields
	set poss [$object.options.fields curselection]
	set pos [lindex $poss 0]
	incr pos -1
	if {$pos < 0} {set pos 0}
	set movefields [list_sub $fields $poss]
	set temp [list_sub $fields -exclude $poss]
	set fields [linsert $temp $pos {*}$movefields]
	$object.options.fields configure -content $fields
	$object.options.fields activate $pos
	$object.options.fields selection set $pos [expr {$pos+[llength $movefields]-1}]
}

fieldsdialog method fieldsdown {} {
	private $object tfields fields
	set poss [$object.options.fields curselection]
	set pos [lindex $poss end]
	incr pos
	set lastfield [lindex $fields $pos]
	set movefields [list_sub $fields $poss]
	set temp [list_sub $fields -exclude $poss]
	if {$lastfield eq ""} {
		set pos end
	} else {
		set pos [lsearch $fields $lastfield]
	}
	set fields [linsert $temp $pos {*}$movefields]
	$object.options.fields configure -content $fields
	if {$pos eq "end"} {set pos [expr {[llength $fields]-[llength $movefields]}]}
	$object.options.fields activate $pos
	$object.options.fields selection set $pos [expr {$pos+[llength $movefields]-1}]
}

fieldsdialog method setcalc {args} {
	$object.options.calc nocmdset [lindex $args 0]
}

fieldsdialog method addcalc {args} {
	private $object tfields fields
putsvars args
	if {![llength $args]} {
		catch {destroy $object.calcfield}
		Classy::Dialog $object.calcfield -title "Add calculated field"
		$object.calcfield option entry "Fieldname" [privatevar $object calcfield]
		$object.calcfield option text "Code" [privatevar $object calccode] "Field Code"
		$object.calcfield add add "Add" "$object addcalc \[getprivate $object calcfield\] \[getprivate $object calccode\]" default
		$object.calcfield persistent remove add
		
	} else {
		foreach {code} $args break
		if {![regexp {^([^=]+)=} $code temp field]} {
			error "Calculated field must be of format: fieldname=code"
		}
		set pos [$object.options.fields index active]
		incr pos
		set new 1
		set curpos 0
		foreach line $fields {
			if {[regexp {^([^=]+)=} $line temp testfield] && $testfield == $field} {
				set new 0
				set pos $curpos
				break
			}
			incr curpos
		}
		if {$new} {
			set fields [linsert $fields $pos $code]
		} else {
			set fields [lreplace $fields $pos $pos $code]
		}
	}
	$object redraw
	Classy::todo $object.options.fields activate $code
}

#fieldsdialog method search {args} {
#	$object redraw
#}

fieldsdialog method editfields {} {
	private $object tfields fields
	catch {destroy $object.editfields}
	Classy::Dialog $object.editfields -title "Edit fields"
	Classy::Text $object.editfields.options.text
	pack $object.editfields.options.text -fill both -expand yes
	$object.editfields.options.text insert end [join $fields \n]
	$object.editfields add change "Change" "[list $object] fields \[split \[$object.editfields.options.text get 1.0 end\] \\n\]"
	$object.editfields persistent remove change
}

fieldsdialog method fieldbuilder {} {
	set data [split [$object.options.calc get] =]
	set p [winfo parent $object]
	private $p fieldname fieldcalc
	foreach {fieldname fieldcalc} $data break
	$p querybuilder $object
}


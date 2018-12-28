Classy::Dialog subclass sortdialog

sortdialog method init args {
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
	button $object.options.buttons.addsort \
		-text -->
	grid $object.options.buttons.addsort -row 1 -column 0 -sticky nesw
	button $object.options.buttons.delsort \
		-text <--
	grid $object.options.buttons.delsort -row 2 -column 0 -sticky nesw
	button $object.options.buttons.clear \
		-text Clear
	grid $object.options.buttons.clear -row 4 -column 0 -sticky nesw
	button $object.options.buttons.up \
		-text Up
	grid $object.options.buttons.up -row 5 -column 0 -sticky nesw
	button $object.options.buttons.down \
		-text Down
	grid $object.options.buttons.down -row 6 -column 0 -sticky nesw
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
	button $object.options.buttons.edit \
		-text Edit
	grid $object.options.buttons.edit -row 7 -column 0 -sticky nesw
	grid columnconfigure $object.options.buttons 0 -uniform {}
	grid rowconfigure $object.options.buttons 0 -uniform {} -weight 1
	grid rowconfigure $object.options.buttons 1 -uniform {}
	grid rowconfigure $object.options.buttons 2 -uniform {}
	grid rowconfigure $object.options.buttons 3 -uniform {}
	grid rowconfigure $object.options.buttons 4 -uniform {}
	grid rowconfigure $object.options.buttons 5 -uniform {}
	grid rowconfigure $object.options.buttons 8 -uniform {} -weight 1
	Classy::ListBox $object.options.selsort  \
		-height 4
	grid $object.options.selsort -row 0 -column 0 -rowspan 2 -sticky nesw
	Classy::ListBox $object.options.sort  \
		-height 4 \
		-width 10
	grid $object.options.sort -row 1 -column 3 -columnspan 2 -sticky nesw
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
		-title "Field selection"
	$object.options.paned configure \
		-window [varsubst object {$object.options.selsort}]
	$object.options.buttons.addsort configure \
		-command [varsubst object {$object addsort}]
	$object.options.buttons.delsort configure \
		-command [varsubst object {$object delsort}]
	$object.options.buttons.clear configure \
		-command [varsubst object {$object clear}]
	$object.options.buttons.up configure \
		-command [varsubst object {$object sortup}]
	$object.options.buttons.down configure \
		-command [varsubst object {$object sortdown}]
	$object.options.buttons.edit configure \
		-command [varsubst object {$object editsort}]
	$object.options.selsort configure \
		-filtervariable [privatevar $object filterselsort]
	$object.options.sort configure \
		-filtervariable [privatevar $object filtersort]
	$object add go Go [varsubst object {$object go}] default
	$object persistent set 
	# Configure initial arguments
	if {"$args" != ""} {$object configure {*}$args}
# ClassyTk Finalise
$object start
	return $object
}

sortdialog addoption -command {command Command {}} {}

sortdialog method start {} {
	setprivate $object tfields {}
	$object.options.selsort configure -selectmode extended -exportselection 0
	$object.options.sort configure -selectmode extended -exportselection 0
	bind $object.options.sort <Up> "[list $object sortup]; break"
	bind $object.options.sort <Down> "[list $object sortdown]; break"
#	$object.options.sort configure -browsecommand [list $object setcalc]
	$object.options.sort activate end
}

sortdialog method go {} {
	private $object options sort
	{*}$options(-command) $sort
}

sortdialog method redraw {} {
	private $object tfields sort
	set selection [$object.options.sort get]
	set active [$object.options.sort index active]
	$object.options.sort configure -content $sort
	set showtfields [list_lremove $tfields $sort]
	$object.options.selsort configure -content $showtfields
	$object selectsort $selection
	$object.options.sort.list activate $active
}

sortdialog method tfields {fields} {
	private $object tfields
	set tfields {}
	foreach field $fields {
		lappend tfields $field -$field
	}
	$object redraw
}

sortdialog method sort {cursort} {
	private $object sort
	set sort [list_remove $cursort {}]
	$object redraw
	$object.options.sort selection set end
	$object.options.sort activate end
}

sortdialog method selectsort {selsort} {
	private $object sort
	set poss [list_cor $sort $selsort]
	foreach pos $poss {
		$object.options.sort selection set $pos
	}
}

sortdialog method fieldadd {movesort} {
	private $object sort
	set pos [$object.options.sort index active]
	incr pos
	set sort [linsert $sort $pos {*}$movesort]
	$object redraw
	$object selectsort $movesort
	$object.options.sort activate [expr {$pos+[llength $movesort]-1}]
}

sortdialog method addsort {} {
	private $object sort
	set movesort [$object.options.selsort get]
	if {![llength $movesort]} {
		set movesort [$object.options.selsort get 0 end]
	}
	$object fieldadd $movesort
}

sortdialog method delsort {} {
	private $object sort
	set movesort [$object.options.sort get]
	if {![llength $movesort]} {
		set movesort [$object.options.sort get 0 end]
	}
	set sort [list_lremove $sort $movesort]
	$object redraw
}

sortdialog method clear {} {
	private $object sort
	set sort {}
	$object redraw
}

sortdialog method sortup {} {
	private $object sort
	set poss [$object.options.sort curselection]
	if {![llength $poss]} {set poss [$object.options.sort index active]}
	set pos [lindex $poss 0]
	incr pos -1
	if {$pos < 0} {set pos 0}
	set movesort [list_sub $sort $poss]
	set temp [list_sub $sort -exclude $poss]
	set sort [linsert $temp $pos {*}$movesort]
	$object.options.sort configure -content $sort
	$object.options.sort activate $pos
	$object.options.sort selection set $pos [expr {$pos+[llength $movesort]-1}]
	$object.options.sort see $pos
}

sortdialog method sortdown {} {
	private $object sort
	set poss [$object.options.sort curselection]
	if {![llength $poss]} {set poss [$object.options.sort index active]}
	set pos [lindex $poss end]
	incr pos
	set lastfield [lindex $sort $pos]
	set movesort [list_sub $sort $poss]
	set temp [list_sub $sort -exclude $poss]
	if {$lastfield eq ""} {
		set pos end
	} else {
		set pos [lsearch $sort $lastfield]
	}
	set sort [linsert $temp $pos {*}$movesort]
	$object.options.sort configure -content $sort
	if {$pos eq "end"} {set pos [expr {[llength $sort]-[llength $movesort]}]}
	$object.options.sort activate $pos
	$object.options.sort selection set $pos [expr {$pos+[llength $movesort]-1}]
	$object.options.sort see $pos
}

sortdialog method editsort {} {
	private $object sort
	catch {destroy $object.editsort}
	Classy::Dialog $object.editsort -title "Edit sort"
	Classy::Text $object.editsort.options.text
	pack $object.editsort.options.text -fill both -expand yes
	$object.editsort.options.text insert end [join $sort \n]
	$object.editsort add change "Change" "[list $object] sort \[split \[$object.editsort.options.text get 1.0 end\] \\n\]"
	$object.editsort persistent remove change
}

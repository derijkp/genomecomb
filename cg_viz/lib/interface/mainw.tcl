Classy::Toplevel subclass mainw
mainw method init args {
	super init
	# Create windows
	Classy::DynaTool $object.maintool  \
		-type MainTool
	grid $object.maintool -row 0 -column 0 -columnspan 3 -sticky new
	Classy::TreeWidget $object.tree \
		-width 50 \
		-height 50
	grid $object.tree -row 1 -column 0 -rowspan 4 -sticky nesw
	Classy::Paned $object.panedh \
		-orient vertical \
		-cursor sb_v_double_arrow
	grid $object.panedh -row 3 -column 2 -sticky nesw
	Classy::Paned $object.panedv
	grid $object.panedv -row 1 -column 1 -rowspan 4 -sticky nesw
	frame $object.table  \
		-borderwidth 0 \
		-relief groove \
		-height 10 \
		-width 10
	grid $object.table -row 4 -column 2 -sticky nesw
	scrollbar $object.table.sv
	grid $object.table.sv -row 1 -column 1 -sticky nesw
	scrollbar $object.table.sh \
		-orient horizontal
	grid $object.table.sh -row 2 -column 0 -sticky nesw
	frame $object.table.buttons  \
		-borderwidth 0 \
		-relief groove \
		-height 10 \
		-highlightthickness 0 \
		-width 10
	grid $object.table.buttons -row 0 -column 0 -columnspan 2 -sticky nesw
	label $object.table.buttons.label1 \
		-text Tableinfo
	grid $object.table.buttons.label1 -row 0 -column 0 -sticky nesw
	button $object.table.buttons.querybuilder \
		-text Query
	grid $object.table.buttons.querybuilder -row 0 -column 2 -sticky nesw
	Classy::Entry $object.table.buttons.query \
		-combo 20 \
		-width 4
	grid $object.table.buttons.query -row 0 -column 3 -sticky nesw
	button $object.table.buttons.button1 \
		-text Fields
	grid $object.table.buttons.button1 -row 0 -column 1 -sticky nesw
	grid columnconfigure $object.table.buttons 0 -uniform {}
	grid columnconfigure $object.table.buttons 1 -uniform {}
	grid columnconfigure $object.table.buttons 2 -uniform {}
	grid columnconfigure $object.table.buttons 3 -uniform {} -weight 1
	grid columnconfigure $object.table.buttons 4 -uniform {}
	grid columnconfigure $object.table.buttons 5 -uniform {}
	grid rowconfigure $object.table.buttons 0 -uniform {}
	grid rowconfigure $object.table.buttons 1 -uniform {}
	grid rowconfigure $object.table.buttons 2 -uniform {}
	Classy::TkTable $object.table.data \
		-rows 0 \
		-cols 1
	grid $object.table.data -row 1 -column 0 -sticky nesw
	grid columnconfigure $object.table 0 -uniform {} -weight 1
	grid columnconfigure $object.table 1 -uniform {}
	grid rowconfigure $object.table 0 -uniform {}
	grid rowconfigure $object.table 1 -uniform {} -weight 1
	grid rowconfigure $object.table 2 -uniform {}
	frame $object.canvas  \
		-borderwidth 0 \
		-relief groove \
		-height 10 \
		-highlightthickness 0 \
		-width 10
	grid $object.canvas -row 2 -column 2 -sticky nesw
	canvas $object.canvas.data \
		-height 50 \
		-highlightthickness 0 \
		-width 50
	grid $object.canvas.data -row 0 -column 0 -sticky nesw
	scrollbar $object.canvas.sv
	grid $object.canvas.sv -row 0 -column 1 -sticky nesw
	scrollbar $object.canvas.sh \
		-orient horizontal
	grid $object.canvas.sh -row 1 -column 0 -sticky nesw
	grid columnconfigure $object.canvas 0 -uniform {} -weight 1
	grid columnconfigure $object.canvas 1 -uniform {}
	grid rowconfigure $object.canvas 0 -uniform {} -weight 1
	grid rowconfigure $object.canvas 1 -uniform {}
	grid columnconfigure $object 0 -uniform {}
	grid columnconfigure $object 1 -uniform {}
	grid columnconfigure $object 2 -uniform {} -weight 1
	grid rowconfigure $object 0 -uniform {}
	grid rowconfigure $object 1 -uniform {}
	grid rowconfigure $object 2 -uniform {}
	grid rowconfigure $object 3 -uniform {}
	grid rowconfigure $object 4 -uniform {} -weight 1

	if {"$args" == "___Classy::Builder__create"} {return $object}
# ClassyTk Initialise
$object start
	# Parse this
	$object configure \
		-title [tk appname]
	$object configure \
		-destroycommand "exit"
	$object.maintool configure \
		-cmdw [varsubst object {$object}]
	$object.panedh configure \
		-window [varsubst object {$object.canvas}]
	$object.panedv configure \
		-window [varsubst object {$object.tree}]
	$object.table.sv configure \
		-command "$object.table.data yview"
	$object.table.sh configure \
		-command "$object.table.data xview"
	$object.table.buttons.query configure \
		-command [varsubst object {$object query}]
	$object.table.buttons.button1 configure \
		-command [varsubst object {$object fields}]
	$object.table.data configure \
		-xscrollcommand "$object.table.sh set" \
		-yscrollcommand "$object.table.sv set"
	$object.canvas.data configure \
		-xscrollcommand "$object.canvas.sh set" \
		-yscrollcommand "$object.canvas.sv set"
	$object.canvas.sv configure \
		-command "$object.canvas.data yview"
	$object.canvas.sh configure \
		-command "$object.canvas.data xview"
	Classy::DynaMenu attachmainmenu MainMenu $object
	# Configure initial arguments
	if {"$args" != ""} {eval $object configure $args}
	return $object
}

mainw method start {args} {
	private $object fields
#	bind $object.canvas.data <1> [list $object select %x %y]
	$object.table.buttons.querybuilder configure -command [list $object querybuilder]
	set fields {}
	Extral::event listen $object selchanged [list $object redrawselection]
	Extral::event listen $object querychanged [list $object redrawquery]
}

mainw method redrawselection {args} {
	puts redrawselection
	set var [$object.table.buttons.query cget -textvariable]
	set table [$object.tb table]
	set query [get $var ""]
	if {[string trim $query ] eq ""} {
		set query [subst {"rowid" in (select "id" from "sel_$table")}]
	} else {
		set extra [subst {"rowid" in (select "id" from "sel_$table")}]
		if {[string first $extra $query] == -1} {
			set query "($query) and $extra"
		}
	}
	$object query $query
}

mainw method redrawquery {args} {
puts "redrawquery $args"
	set len [$object.tb info len]
	set total [$object.tb info tlen]
	$object.table.buttons.label1 configure -text "[commify $len] / [commify $total] lines"
	$object.table.data configure -command [$object.table.data cget -command] -usecommand 1
	$object.table.data configure -rows [expr {$len+1}] -cols [llength [$object.tb info qfields]]
}

#mainw method select {x y} {
#putsvars x y
#	set canvas $object.canvas.data
#	set id [$canvas find closest [$canvas canvasx $x] [$canvas canvasy $y]]
#	set tagslist [list [$canvas itemcget $id -tags]]
#	set dobj [lindex $tagslist 0 0]
#	$dobj select $tagslist
#	# catch events for display later!
#}

mainw method opendb {args} {
	if {[llength $args] < 3} {
		puts stderr "format is cg viz database table refdir"
		exit 1
	}
	foreach {database table dbdir dbfarm} $args break
	catch {$object.tb destroy}
	table_monetdb new $object.tb
	$object.tb open $database $table $dbfarm
	set tables [cg_monetdb_tables $database]
	if {[lsearch $tables sel_$table] == -1} {
		cg_monetdb_sql $database [subst {create table "sel_$table" (id integer)}]
	}
	$object.tb link $object.table.data
	# default graph
	catch {$object.disp1 destroy}
	display_chr new $object.disp1 $object.canvas.data $object.tb $dbdir
	$object.disp1 redraw
	# $object.disp1 options scale 10000000
}

mainw method opentsv {args} {
	foreach {file} $args break
	catch {$object.tb destroy}
	table_tsv new $object.tb
	set file [$object.tb open $file]
	$object.tb link $object.table.data
	# default graph
#	catch {$object.disp1 destroy}
#	display_chr new $object.disp1 $object.canvas.data $object.tb $dbdir
#	$object.disp1 redraw
	# $object.disp1 options scale 10000000
	wm title $object $file
}

mainw method savetsv {file} {
	$object.tb save $file
}

mainw method open {args} {
	if {[llength $args]} {
		set file [lindex $args 0]
	} else {
		set file [Classy::selectfile -title Open -selectmode persistent]
	}
	$object opentsv $file
}

mainw method save {args} {
	if {[llength $args]} {
		set file [lindex $args 0]
	} else {
		set file [Classy::savefile -title Save]
	}
	$object savetsv $file
}

proc commify {num {sep ,}} {
    while {[regsub {^([-+]?\d+)(\d\d\d)} $num "\\1$sep\\2" num]} {}
    return $num
}

mainw method query {args} {
	set tb $object.tb
	set var [$object.table.buttons.query cget -textvariable]
	if {[llength $args]} {
		set query [lindex $args 0]
		set $var $query
	} else {
		set query [get $var ""]
	}
	$tb query $query
}

mainw method querybuilder {args} {
	set tb $object.tb
	destroy $object.querybuilder
	set var [$object.table.buttons.query cget -textvariable]
	Classy::Dialog $object.querybuilder -title Query
	$object.querybuilder option text Query $var Query
	$object.querybuilder add query "Query" [list $object query] default
}

mainw method fields {args} {
	private $object fields tfields w
	set tb $object.tb
	if {![llength $args]} {
		set tfields [$tb tfields]
		set fields [$tb info fields]
		set w $object.fieldsdialog
		catch {destroy $w}
		fieldsdialog $w
		$w fields $fields
		$w tfields $tfields
		$w configure -command [list $object fields]
	} else {
		$tb fields [lindex $args 0]
	}
}

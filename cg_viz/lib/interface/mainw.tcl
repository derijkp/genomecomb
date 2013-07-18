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
	grid $object.tree -row 1 -column 0 -rowspan 5 -sticky nesw
	Classy::Paned $object.panedh \
		-orient vertical \
		-cursor sb_v_double_arrow
	grid $object.panedh -row 3 -column 2 -sticky nesw
	Classy::Paned $object.panedv
	grid $object.panedv -row 1 -column 1 -rowspan 5 -sticky nesw
	frame $object.table  \
		-borderwidth 0 \
		-relief groove \
		-height 10 \
		-width 10
	grid $object.table -row 5 -column 2 -sticky nesw
	scrollbar $object.table.sv
	grid $object.table.sv -row 1 -column 1 -sticky nesw
	scrollbar $object.table.sh \
		-orient horizontal
	grid $object.table.sh -row 2 -column 0 -sticky nesw
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
	frame $object.buttons  \
		-borderwidth 0 \
		-relief groove \
		-height 10 \
		-highlightthickness 0 \
		-width 10
	grid $object.buttons -row 4 -column 2 -sticky nesw
	label $object.buttons.label1 \
		-text Tableinfo
	grid $object.buttons.label1 -row 0 -column 0 -sticky nesw
	button $object.buttons.querybuilder \
		-text Query
	grid $object.buttons.querybuilder -row 0 -column 2 -sticky nesw
	Classy::Entry $object.buttons.query \
		-combo 20 \
		-width 4
	grid $object.buttons.query -row 0 -column 3 -sticky nesw
	button $object.buttons.button1 \
		-text Fields
	grid $object.buttons.button1 -row 0 -column 1 -sticky nesw
	Classy::OptionMenu $object.buttons.view  \
		-list {data
summary
graph}
	grid $object.buttons.view -row 0 -column 5 -sticky nesw
	button $object.buttons.settings \
		-text Settings
	grid $object.buttons.settings -row 0 -column 4 -sticky nesw
	grid columnconfigure $object.buttons 0 -uniform {}
	grid columnconfigure $object.buttons 1 -uniform {}
	grid columnconfigure $object.buttons 2 -uniform {}
	grid columnconfigure $object.buttons 3 -uniform {} -weight 1
	grid columnconfigure $object.buttons 4 -uniform {}
	grid columnconfigure $object.buttons 5 -uniform {}
	grid columnconfigure $object.buttons 6 -uniform {}
	grid columnconfigure $object.buttons 7 -uniform {}
	grid rowconfigure $object.buttons 0 -uniform {}
	grid rowconfigure $object.buttons 1 -uniform {}
	grid rowconfigure $object.buttons 2 -uniform {}
	grid columnconfigure $object 0 -uniform {}
	grid columnconfigure $object 1 -uniform {}
	grid columnconfigure $object 2 -uniform {} -weight 1
	grid rowconfigure $object 0 -uniform {}
	grid rowconfigure $object 1 -uniform {}
	grid rowconfigure $object 2 -uniform {}
	grid rowconfigure $object 3 -uniform {}
	grid rowconfigure $object 4 -uniform {}
	grid rowconfigure $object 5 -uniform {} -weight 1

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
	$object.tree configure \
		-endnodecommand [varsubst object {$object tree_browse}] \
		-opencommand [varsubst object {$object tree_open}] \
		-closecommand [varsubst object {$object tree_close}]
	$object.panedh configure \
		-window [varsubst object {$object.canvas}]
	$object.panedv configure \
		-window [varsubst object {$object.tree}]
	$object.table.sv configure \
		-command "$object.table.data yview"
	$object.table.sh configure \
		-command "$object.table.data xview"
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
	$object.buttons.query configure \
		-command [varsubst object {$object query}]
	$object.buttons.button1 configure \
		-command [varsubst object {$object fields}]
	$object.buttons.view configure \
		-command [varsubst object {$object view}] \
		-textvariable [privatevar $object view(cur)]
	$object.buttons.settings configure \
		-command [varsubst object {$object summarybuilder}]
	Classy::DynaMenu attachmainmenu MainMenu $object
	# Configure initial arguments
	if {"$args" != ""} {eval $object configure $args}
	return $object
}

mainw method start {args} {
	private $object fields view
	set view(cur) data
#	bind $object.canvas.data <1> [list $object select %x %y]
	$object.buttons.querybuilder configure -command [list $object querybuilder]
	set fields {}
	Extral::event listen $object selchanged [list $object redrawselection]
	Extral::event listen $object querychanged [list $object redrawquery]
	$object.tree clearnode {}
	$object.table.data configure -wrap 1 -resizeborders both -bordercursor crosshair -anchor n -justify right
}

mainw method redrawselection {args} {
	puts redrawselection
	set var [$object.buttons.query cget -textvariable]
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
	$object.buttons.label1 configure -text "[commify $len] / [commify $total] lines"
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
	$object.tb link $object.table.data $object.buttons.query
	# default graph
	catch {$object.disp1 destroy}
	display_chr new $object.disp1 $object.canvas.data $object.tb $dbdir
	$object.disp1 redraw
	# $object.disp1 options scale 10000000
}

mainw method tree_close {dir} {
	catch {$object.tree clearnode $dir}
}

mainw method tree_browse {file} {
	$object open $file
}

mainw method tree_open {dir} {
putsvars dir
	catch {$object.tree clearnode $dir}
	set dirs {}
	set files {}
	foreach file [glob -nocomplain $dir/*] {
		set ext [file extension $file]
		if {$ext eq ".index"} continue
		if {[file isdirectory $file]} {
			lappend dirs $file
		} else {
			lappend files $file
		}
	}
	foreach file [lsort $dirs] {
		$object.tree addnode $dir $file -text [file tail $file]
	}
	foreach file [lsort $files] {
		$object.tree addnode $dir $file -type end -text [file tail $file]
	}
}

mainw method opentsv {args} {
	private $object view
	foreach {file} $args break
	catch {$object.tb destroy}
	table_tsv new $object.tb
	$object.tb link $object.table.data $object.buttons.query
	set file [$object.tb open $file $object]
	wm title $object $file
	# find root and start tree
	set file [file normalize $file]
	if {[file isdir $file]} {
		set root $file
	} else {
		set root [file dir $file]
	}
	lappend path $root
	while 1 {
		if {$root eq "." || $root eq "/"} {
			set root [file dir $file]
			break
		}
		if {[file exists $root/compar]} break
		if {[llength [glob -nocomplain $root/*.cgproj]]} break
		lappend path [file tail $root]
		set root [file dir $root]
	}
	$object.tree selection set node $file
	if {![llength [$object.tree selection]]} {
		$object.tree clearnode {}
		$object.tree addnode {} $root -text [file tail $root]
		$object tree_open $root
		foreach dir [list_reverse $path] {
			set root $root/$dir
			$object tree_open $root
		}
		$object.tree selection set node $file
	}
	if {[info exists tdata(sqlbackend_db)]} {
		# default graph
		catch {$object.disp1 destroy}
		display_chr new $object.disp1 $object.canvas.data $object.tb $tdata(dbdir)
		$object.disp1 redraw
		# $object.disp1 options scale 10000000
	}
	set fields [$object.tb fields]
	set sample [lindex [samples $fields] 0]
	set view(summary_rows) {chromosome}
	set view(summary_cols) [list sample {} sequenced {v}]
	set view(summary_cells) [list count]
}

mainw method savetsv {file} {
	$object.tb save $file
}

mainw method open {args} {
putsvars args
	private $object tdata
	if {[llength $args]} {
		set file [lindex $args 0]
	} else {
		set file [Classy::selectfile -title Open -selectmode persistent]
	}
	set file [file normalize $file]
	set tdata(file) $file
	$object view data
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
	set var [$object.buttons.query cget -textvariable]
	if {[llength $args]} {
		set query [lindex $args 0]
		set $var $query
	} else {
		set query [get $var ""]
	}
	$tb query $query
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

mainw method summary_redrawselect {args} {
puts summary_redrawselect
return
}

mainw method summary_recalc {args} {
	private $object summary view
	set definition [list $view(summary_rows) $view(summary_cols) $view(summary_cells)]
	if {![info exists summary] || [inlist $args querychanged] || $definition ne [get view(olddefintion) ""]} {
		set error [catch {$object.tb summary $definition} newsummary]
		if {$error} {
			error "error in summary definition (or query)\n\n$newsummary"
		} else {
			set header [list_shift newsummary]
			set summary [list $header {*}[lsort -dict $newsummary]]
			set view(olddefintion) $definition
		}
	}
	return $summary
}

mainw method savesummary {args} {
	private $object summary view
	if {[llength $args]} {
		set file [lindex $args 0]
	} else {
		set file [Classy::savefile -title "Save Summary"]
	}
	set f [open $file w]
	puts $f "\# Summary definition: [get view(olddefintion) ""]"
	csv_write $f $summary \t {}
	close $f
}

mainw method summary_redraw {args} {
puts summary_redraw
putsvars args
	private $object summary view
	set error [catch {
		set summary [$object summary_recalc {*}$args]
	} message]
	$object.summary.data configure -variable [privatevar $object summary] -variabletype list
	$object.summary.data configure -rows [llength $summary] -cols [llength [lindex $summary 0]]	
	if {$error} {error $message}
}

mainw method graph_redrawselect {args} {
	puts graph_redrawselect
	return
}

mainw method graph_redraw {args} {
	private $object view
	set error [catch {
		set summary [$object summary_recalc {*}$args]
	} message]
	$object.graph clear
	if {!$error} {
		$object.graph add $summary
	} else {
		error "error in graph definition\n\n$message"
	}
}

mainw method view {newview} {
	private $object view
	set view(cur) $newview
	catch {grid forget $object.table}
	foreach w {summary graph} {
		Extral::event remove $object.$w selchanged [list $object ${w}_redrawselect]
		Extral::event remove $object.$w querychanged [list $object ${w}_redraw]
		catch {grid forget $object.$w}
		# destroy $object.$w
	}
	switch $newview {
		data {
			grid $object.table -row 5 -column 2 -sticky nesw
		}
		summary {
			if {![winfo exists $object.summary]} {
				frame $object.summary  \
					-borderwidth 0 \
					-relief groove \
					-height 10 \
					-width 10
				Classy::TkTable $object.summary.data -rows 0 -cols 1
				scrollbar $object.summary.sv
				scrollbar $object.summary.sh -orient horizontal
				grid $object.summary.data -row 0 -column 0 -sticky nesw
				grid $object.summary.sv -row 0 -column 1 -sticky nesw
				grid $object.summary.sh -row 1 -column 0 -sticky nesw
				grid columnconfigure $object.summary 0 -weight 1
				grid rowconfigure $object.summary 0 -weight 1
				$object.summary.sv configure \
					-command "$object.summary.data yview"
				$object.summary.sh configure \
					-command "$object.summary.data xview"
				$object.summary.data configure \
					-xscrollcommand "$object.summary.sh set" \
					-yscrollcommand "$object.summary.sv set"
				$object.summary.data configure -wrap 1 -resizeborders both -bordercursor crosshair
			}
			Extral::event listen $object.summary selchanged [list $object summary_redrawselect]
			Extral::event listen $object.summary querychanged [list $object summary_redraw querychanged]
			grid $object.summary -row 5 -column 2 -sticky nesw
			$object summary_redraw
		}
		graph {
			destroy $object.graph
			scrolledgraph $object.graph
			Extral::event listen $object.graph selchanged [list $object graph_redrawselect]
			Extral::event listen $object.graph querychanged [list $object graph_redraw querychanged]
			grid $object.graph -row 5 -column 2 -sticky nesw
			$object graph_redraw
		}
	}
}

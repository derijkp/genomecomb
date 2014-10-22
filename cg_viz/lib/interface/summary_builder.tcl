# not really clean code alert:
# I do not create a separate object $object.summarybuilder, but I do use private variables to
# $object.summarybuilder to keep them separate from similar variables in querybuilder

mainw method summarybuilder_fillvalues {{samplevalues 0} args} {
	private $object.summarybuilder qvalues qvaluestext fieldsw valuesw valuew
	set w $object.summarybuilder.options.paned
	set field [lindex [$fieldsw get] 0]
	if {$field eq ""} return
	$object fieldvalues $field $samplevalues qvalues qvaluestext
	$valuew set {}
	return $qvalues
}

mainw method summarybuilder_browsefield {args} {
	private $object.summarybuilder valuesw
	Classy::todo $object summarybuilder_fillvalues
}

mainw method summarybuilder_browsevalue {args} {
	private $object.summarybuilder valuesw valuew
	$valuew set [$valuesw get]
}

mainw method summarybuilder_selvalue {value} {
	if {[isdouble $value]} {
		set insert " $value "
	} else {
		set insert " \"$value\" "
	}
	$object.summarybuilder.options.paned.query.query insert current $insert
}

mainw method summarybuilder_quotevalue {value} {
	if {[isdouble $value]} {
		return $value
	} else {
		return \"$value\"
	}
}

mainw method summarybuilder_insert {insert {join and}} {
	private $object.summarybuilder queryw
	upvar #0 [$object.buttons.query cget -textvariable] query
	set w $object.summarybuilder.options.paned
	if {[string trim $query] ne ""} {set insert "\n$join $insert"}
	$queryw insert insert $insert
}

mainw method summarybuilder_add_row {} {
	private $object.summarybuilder fieldsw valuew rowsw
	private $object view
	set fields [$fieldsw get]
	set values [list [$valuew get]]
	foreach field $fields value $values {
		$rowsw insert "end" [list $field $value]\n
	}
}

mainw method summarybuilder_add_col {} {
	private $object.summarybuilder fieldsw valuew colsw
	private $object view
	set fields [$fieldsw get]
	set values [list [$valuew get]]
	foreach field $fields value $values {
		$colsw insert "end" [list $field $value]\n
	}
}

mainw method summarybuilder_add_cell {aggr} {
	private $object.summarybuilder fieldsw valuew cellsw
	private $object view
	set fields [$fieldsw get]
	if {[string trim $view(summary_cells)] eq ""} {
		set view(summary_cells) ${aggr}([lindex $fields 0])
		return
	}
	$cellsw insert end ,${aggr}([lindex $fields 0])
	
}

mainw method summarybuilder_add {command {join and}} {
	private $object.summarybuilder fieldsw valuew
	upvar #0 [$object.buttons.query cget -textvariable] query
	set w $object.summarybuilder.options.paned
	set fields [$fieldsw get]
	set values [$valuew get]
	if {$command eq "count"} {
		set value [lindex $values 0]
		set value [$object summarybuilder_quotevalue $value]
		set insert "count\("
		foreach field $fields {
			append insert "\$$field, "
		}
		$object summarybuilder_insert $insert $join
	} elseif {$command eq "value"} {
		$object summarybuilder_insert "\$\{[join $fields "\}\$\{"]\}" $join
	} elseif {$command eq "region"} {
		destroy $object.region
		Classy::Dialog $object.region -title "Select region"
		$object.region option entry Chromosome [privatevar $object region(chr)]
		$object.region option entry Begin [privatevar $object region(begin)]
		$object.region option entry End [privatevar $object region(end)]
		$object.region add go Go "$object summarybuilder_insert region(\[getprivate $object region(chr) \]:\[getprivate $object region(begin) \]-\[getprivate $object region(end) \]) $join" default
	} elseif {$command eq "compare"} {
		private $object compare tdata
		destroy $object.compare
		Classy::Dialog $object.compare -title "Compare"
		set compare(types) {
			{compare returns sm,df,mm,un depending on type of match between two samples}
			{same same: all samples have the same genotype (does not have to be a variant) (all sequenced)}
			{sm same: variant with the same genotype in all given samples (all sequenced)}
			{df different: variant in some, reference in other (all sequenced)}
			{mm mismatch; variant in all, but different genotypes (all sequenced)}
			{un unsequenced in some samples, variant in one of the others}
		}
		$object.compare option select Comparison [privatevar $object compare(type)] [privatevar $object compare(types)]
		set compare(type) [lindex $compare(types) 0]
		set compare(samples) [cg select -n $tdata(file)]
		$object.compare option listbox Sample1 [privatevar $object compare(selsamples)] [privatevar $object compare(samples)] -selectmode multiple
		$object.compare add go Go "$object summarybuilder_insert \"\[lindex \[getprivate $object compare(type) \] 0\]\(\[join \[getprivate $object compare(selsamples) \] ,\]\)\" $join" default
	} elseif {$command in {min max lmin lmax avg}} {
		private $object compare tdata
		destroy $object.compare
		Classy::Dialog $object.compare -title "Functions"
		set compare(functypes) {
			{lmin : lmin(a1,a2,...)**: returns the minimum of a1, a2, ... a1, etc. can be a list of numbers (separated by commas, spaces or ;)}
			{lmax : lmax(a1,a2,...)**: returns the maximum of a1, a2, ... a1, etc. can be a list of numbers (separated by commas, spaces or ;)}
			{avg : returns the average of the values given. Non-number values are ignored. If no number was given, the answer will be NaN}
		}
		$object.compare option select Function [privatevar $object compare(functype)] [privatevar $object compare(functypes)]
		set compare(functype) $command
		set compare(fields) [cg select -h $tdata(file)]
		$object.compare option listbox Sample1 [privatevar $object compare(selfields)] [privatevar $object compare(fields)] -selectmode multiple
		$object.compare add go Go "$object summarybuilder_insert \"\[lindex \[getprivate $object compare(functype) \] 0\]\(\$\[join \[getprivate $object compare(selfields) \] ,\$\]\)\" $join" default
	} else {
		set insert {}
		foreach field $fields {
			set temp {}
			if {[inlist {shares} $operator]} {
				lappend insert "\$$field $operator [list $values]"
			} else { 
				foreach value $values {
					set value [$object summarybuilder_quotevalue $value]
					lappend temp "\$$field $operator $value"
				}
				if {[llength $temp] > 1} {
					lappend insert "\( [join $temp " or "] \)"
				} else {
					lappend insert "[join $temp " or "]"
				}
			}
		}
		if {$command eq "if"} {
			private $object if
			destroy $object.if
			set if(if) [join $insert " and "]
			set if(true) 1
			set if(false) 0
			Classy::Dialog $object.if -title "Conditional value (if)"
			$object.if option entry Criterion [privatevar $object if(if)]
			$object.if option entry Truevalue [privatevar $object if(true)]
			$object.if option entry FalseValue [privatevar $object if(false)]
			$object.if add go Go "$object summarybuilder_insert if(\[getprivate $object if(if) \],\[getprivate $object if(true) \],\[getprivate $object if(false) \]) $join" default
		} else {
			set insert [join $insert "\n    $command "]
			set join $command
			$object summarybuilder_insert $insert $join
		}
	}
}

mainw method summarybuilder {args} {
	private $object.summarybuilder qfields qvalues fieldsw valuesw valuew rowsw colsw cellsw qfieldsfilter valuesfilter
	private $object view
	set qfields [list_union [$object.tb qfields] [$object.tb tfields]]
	lappend qfields all sample
	foreach field [list_sub $qfields [list_find -regexp $qfields -]] {
		lappend qfields [lindex [split $field -] 0]
	}
	set qfields [list_remdup $qfields]
	set qfieldsfilter {}
	set valuesfilter {}
	destroy $object.summarybuilder
	Classy::Dialog $object.summarybuilder -title "Summary options"
	set w $object.summarybuilder.options.paned
	ttk::panedwindow $w -orient horizontal
	pack $w -fill both -expand yes
	# fields
	frame $w.fields -borderwidth 0 -highlightthickness 0
	button $w.fields.header -text "Fields" -command "[list set [privatevar $object.summarybuilder qfieldsfilter] {}]"
	pack $w.fields.header -fill x
	set fieldsw $w.fields.fields
	Classy::ListBox $w.fields.fields \
		-listvariable [privatevar $object.summarybuilder qfields] \
		-filtervariable [privatevar $object.summarybuilder qfieldsfilter] \
		-browsecommand [list $object summarybuilder_browsefield] \
		-exportselection 0 -selectmode extended
	pack $w.fields.fields -fill both -expand yes
	# values
	set w $w.pane3
	ttk::panedwindow $w -orient horizontal
	frame $w.values -borderwidth 0 -highlightthickness 0
	button $w.values.header -text "Values" -command "[list set [privatevar $object.summarybuilder valuesfilter] {}]"
	pack $w.values.header -fill x
	set valuew $w.values.value
	Classy::Entry $w.values.value
	pack $w.values.value -fill x
	frame $w.values.buttons
	pack $w.values.buttons -fill x
	button $w.values.buttons.examples -text "Get Examples" -command [list $object summarybuilder_fillvalues 1]
	pack $w.values.buttons.examples -side left
	label $w.values.label -textvariable [privatevar $object.summarybuilder qvaluestext] -justify left -anchor w
	pack $w.values.label -fill x
	set valuesw $w.values.values
	Classy::ListBox $w.values.values \
		-listvariable [privatevar $object.summarybuilder qvalues] \
		-filtervariable [privatevar $object.summarybuilder qvaluesfilter] \
		-command [list $object summarybuilder_selvalue] \
		-browsecommand [list Classy::todo $object summarybuilder_browsevalue] \
		-exportselection 0 -selectmode extended
	pack $w.values.values -fill both -expand yes
	# query
	# rows
	frame $w.query -borderwidth 0 -highlightthickness 0
	frame $w.query.summary_rows_buttons -borderwidth 0 -highlightthickness 0
	button $w.query.summary_rows_buttons.add -text "Add to row grouping (-g)" -command [list $object summarybuilder_add_row]
	pack $w.query.summary_rows_buttons.add -side left
	button $w.query.summary_rows_buttons.clear -text "Clear" -command [list set [privatevar $object view(summary_rows)] {}]
	pack $w.query.summary_rows_buttons.clear -side left
	set rowsw $w.query.summary_rows
	Classy::ScrolledText $w.query.summary_rows -textvariable [privatevar $object view(summary_rows)] -height 5
	pack $w.query.summary_rows_buttons -fill x
	pack $w.query.summary_rows -fill both -expand yes
	# cols
	frame $w.query.summary_cols_buttons -borderwidth 0 -highlightthickness 0
	button $w.query.summary_cols_buttons.add -text "Add to column grouping (-gc)" -command [list $object summarybuilder_add_col]
	pack $w.query.summary_cols_buttons.add -side left
	button $w.query.summary_cols_buttons.clear -text "Clear" -command [list set [privatevar $object view(summary_cols)] {}]
	pack $w.query.summary_cols_buttons.clear -side left
	set colsw $w.query.summary_cols
	Classy::ScrolledText $w.query.summary_cols -textvariable [privatevar $object view(summary_cols)] -height 5
	pack $w.query.summary_cols_buttons -fill x
	pack $w.query.summary_cols -fill both -expand yes
	# cells
	frame $w.query.summary_cells_buttons -borderwidth 0 -highlightthickness 0
	label $w.query.summary_cells_buttons.add -text "Data in cells"
	pack $w.query.summary_cells_buttons.add -side left
	button $w.query.summary_cells_buttons.clear -text "Clear" -command [list set [privatevar $object view(summary_cells)] {}]
	pack $w.query.summary_cells_buttons.clear -side left
	foreach aggr {count ucount min max avg sum percent gpercent distinct list} {
		button $w.query.summary_cells_buttons.$aggr -text $aggr -command [list $object summarybuilder_add_cell $aggr]
		pack $w.query.summary_cells_buttons.$aggr -side left
	}
	set cellsw $w.query.summary_cells
	Classy::Entry $w.query.summary_cells -textvariable [privatevar $object view(summary_cells)]
	pack $w.query.summary_cells_buttons -fill x
	pack $w.query.summary_cells -fill x
	# fill paned
	set w $object.summarybuilder.options.paned
	$w add $w.fields
	$w add $w.pane3
	set w $w.pane3
	$w add $w.values
	$w add $w.query
	set w $object.summarybuilder.options.paned
	$object.summarybuilder add summary "Summary" [list $object view summary]
	$object.summarybuilder add graph "Graph" [list $object view summarygraph]
#	$object.summarybuilder add scatter "Scatterplot" [list $object view summaryscatter]
	$object summarybuilder_fillvalues
	$object.summarybuilder persistent remove summary graph
	update
	set field [$w.fields.fields get 0]
	if {$field ne ""} {
		$w.fields.fields set $field
		Classy::todo $object summarybuilder_browsefield
	}
}


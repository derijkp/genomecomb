mainw method querybuilder_filteroperators {args} {
	private $object qoperators qoperatorsfilter
	set qoperators {== != < <= > >= regexp contains oneof shares }
#	if {$qoperatorsfilter ne "" && $qoperatorsfilter ne "*"} {
#		set poss [list_find -regexp $qoperators $qoperatorsfilter]
#		set qoperators [list_sub $qoperators $poss]
#	}
}

mainw method querybuilder_fillvalues {args} {
	private $object qvalues qvaluestext fieldsw valuesw valuew
	set w $object.querybuilder.options.paned
	set field [lindex [$fieldsw get] 0]
	if {$field eq ""} return
	set list [$object.tb values $field]
	if {[lindex $list end 1] ne "incomplete"} {
		set text [list "All values"]
	}
	while {![isint [lindex $list end 1]]} {
		set line [list_pop list]
		lappend text "[lindex $line 1]: [lindex $line 0]"
	}
	set qvaluestext [join $text \n]
	set qnums [list_subindex $list 1]
	set list [list_subindex $list 0]
	if {[llength $list]} {
		# found the values, do nothing else
	} elseif {[regexp ^sequenced- $field]} {
		set list {r v u}
	} elseif {[regexp _impact $field]} {
		set list {CDS RNA}
	} else {
		set list {}
	}
	set qvalues $list
	if {![inlist $qvalues [$valuew get]]} {
		$valuesw activate [lindex $qvalues 0]
		$valuesw set [lindex $qvalues 0]
	}
#	set qvalues {}
#	foreach value $list num $qnums {
#		lappend qvalues "$value ($num)"
#	}
	return $qvalues
}

mainw method querybuilder_selfield {value} {
	$object.querybuilder.options.paned.query.query insert current " \$$value "
}

mainw method querybuilder_browsefield {args} {
	Classy::todo $object querybuilder_fillvalues
}

mainw method querybuilder_browsevalue {args} {
	private $object valuesw valuew
	set w $object.querybuilder.options.paned
	$valuew set [$valuesw get]
}

mainw method querybuilder_seloperator {value} {
	$object.querybuilder.options.paned.query.query insert current " $value "
}

mainw method querybuilder_selvalue {value} {
	if {[isdouble $value]} {
		set insert " $value "
	} else {
		set insert " \"$value\" "
	}
	$object.querybuilder.options.paned.query.query insert current $insert
}

mainw method querybuilder_quotevalue {value} {
	if {[isdouble $value]} {
		return $value
	} else {
		return \"$value\"
	}
}

mainw method querybuilder_insert {insert {join and}} {
	private $object queryw
	upvar #0 [$object.buttons.query cget -textvariable] query
	set w $object.querybuilder.options.paned
	if {[string trim $query] ne ""} {set insert "\n$join $insert"}
	$queryw insert insert $insert
}

mainw method querybuilder_add {command} {
putsvars command
	private $object fieldsw operatorsw valuew
	upvar #0 [$object.buttons.query cget -textvariable] query
	set w $object.querybuilder.options.paned
	set fields [$fieldsw get]
	set operator [$operatorsw get]
	set values [$valuew get]
	if {$command eq "count"} {
		set value [lindex $values 0]
		set value [$object querybuilder_quotevalue $value]
		set insert "count\("
		foreach field $fields {
			append insert "\$$field, "
		}
		append insert "$operator $value \) > 0"
		set join and
		$object querybuilder_insert $insert $join
	} elseif {$command eq "region"} {
		destroy $object.region
		Classy::Dialog $object.region -title "Select region"
		$object.region option entry Chromosome [privatevar $object region(chr)]
		$object.region option entry Begin [privatevar $object region(begin)]
		$object.region option entry End [privatevar $object region(end)]
		$object.region add go Go "$object querybuilder_insert region(\[getprivate $object region(chr) \]:\[getprivate $object region(begin) \]-\[getprivate $object region(end) \])" default
	} elseif {$command eq "compare"} {
		private $object compare tdata
		destroy $object.compare
		Classy::Dialog $object.compare -title "Compare"
		set compare(types) {
			{same same: all samples have the same genotype (does not have to be a variant) (all sequenced)}
			{sm same: variant with the same genotype in all given samples (all sequenced)}
			{df different: variant in some, reference in other (all sequenced)}
			{mm mismatch; variant in all, but different genotypes (all sequenced)}
			{un unsequenced in some samples, variant in one of the others}
			{compare returns sm,df,mm,un depending on type of match between two samples}
		}
		$object.compare option select Comparison [privatevar $object compare(type)] [privatevar $object compare(types)]
		set compare(type) [lindex $compare(types) 0]
		set compare(samples) [cg select -n $tdata(file)]
		$object.compare option listbox Sample1 [privatevar $object compare(selsamples)] [privatevar $object compare(samples)] -selectmode multiple
		$object.compare add go Go "$object querybuilder_insert \"\[lindex \[getprivate $object compare(type) \] 0\]\(\[join \[getprivate $object compare(selsamples) \] ,\]\)\"" default
	} elseif {$command in {min max avg}} {
		private $object compare tdata
		destroy $object.compare
		Classy::Dialog $object.compare -title "Functions"
		set compare(functypes) {
			{min : min(a1,a2,...)**: returns the minimum of a1, a2, ... a1, etc. can be a list of numbers (separated by commas, spaces or ;)}
			{max : max(a1,a2,...)**: returns the maximum of a1, a2, ... a1, etc. can be a list of numbers (separated by commas, spaces or ;)}
			{avg : returns the average of the values given. Non-number values are ignored. If no number was given, the answer will be NaN}
		}
		$object.compare option select Function [privatevar $object compare(functype)] [privatevar $object compare(functypes)]
		set compare(functype) $command
		set compare(fields) [cg select -h $tdata(file)]
		$object.compare option listbox Sample1 [privatevar $object compare(selfields)] [privatevar $object compare(fields)] -selectmode multiple
		$object.compare add go Go "$object querybuilder_insert \"\[lindex \[getprivate $object compare(functype) \] 0\]\(\[join \[getprivate $object compare(selfields) \] ,\]\)\"" default
	} else {
		set insert {}
		foreach field $fields {
			set temp {}
			if {[inlist {oneof shares} $operator]} {
				lappend insert "\$$field $operator [list $values]"
			} else { 
				foreach value $values {
					set value [$object querybuilder_quotevalue $value]
					lappend temp "\$$field $operator $value"
				}
				if {[llength $temp] > 1} {
					lappend insert "\( [join $temp " or "] \)"
				} else {
					lappend insert "[join $temp " or "]"
				}
			}
		}
		set insert [join $insert "\n    $command "]
		set join $command
		$object querybuilder_insert $insert $join
	}
}

mainw method querybuilder {args} {
	private $object qfields qvalues fieldsw operatorsw valuesw valuew queryw qfieldsfilter valuesfilter
	set qfieldsfilter {}
	set valuesfilter {}
	set var [$object.buttons.query cget -textvariable]
	destroy $object.querybuilder
	Classy::Dialog $object.querybuilder -title Query
	set w $object.querybuilder.options.paned
	ttk::panedwindow $w -orient horizontal
	pack $w -fill both -expand yes
	# fields
	frame $w.fields -borderwidth 0 -highlightthickness 0
	button $w.fields.header -text "Fields" -command "[list set [privatevar $object qfieldsfilter] {}]"
	pack $w.fields.header -fill x
	set fieldsw $w.fields.fields
	Classy::ListBox $w.fields.fields \
		-listvariable [privatevar $object qfields] \
		-filtervariable [privatevar $object qfieldsfilter] \
		-command [list $object querybuilder_selfield] \
		-browsecommand [list $object querybuilder_browsefield] \
		-exportselection 0 -selectmode extended
	pack $w.fields.fields -fill both -expand yes
	# operators
	set w $w.pane2
	ttk::panedwindow $w -orient horizontal
	frame $w.operators -borderwidth 0 -highlightthickness 0
	button $w.operators.header -text "Operator" -command "[list set [privatevar $object qoperatorsfilter] {}]; [list $object querybuilder_filteroperators]"
	pack $w.operators.header -fill x
#	Classy::Entry $w.operators.filter -label Filter \
#		-textvariable [privatevar $object qoperatorsfilter] \
#		-command [list $object querybuilder_filteroperators]
#	pack $w.operators.filter -fill x
	set operatorsw $w.operators.operators
	Classy::ListBox $w.operators.operators \
		-listvariable [privatevar $object qoperators] \
		-command [list $object querybuilder_seloperator] \
		-width 6 -exportselection 0 -selectmode extended
	pack $w.operators.operators -fill both -expand yes
	# values
	set w $w.pane3
	ttk::panedwindow $w -orient horizontal
	frame $w.values -borderwidth 0 -highlightthickness 0
	button $w.values.header -text "Values" -command "[list set [privatevar $object valuesfilter] {}]"
	pack $w.values.header -fill x
	set valuew $w.values.value
	Classy::Entry $w.values.value
	pack $w.values.value -fill x
	label $w.values.label -textvariable [privatevar $object qvaluestext] -justify left -anchor w
	pack $w.values.label -fill x
	set valuesw $w.values.values
	Classy::ListBox $w.values.values \
		-listvariable [privatevar $object qvalues] \
		-filtervariable [privatevar $object qvaluesfilter] \
		-command [list $object querybuilder_selvalue] \
		-browsecommand [list $object querybuilder_browsevalue] \
		-exportselection 0 -selectmode extended
	pack $w.values.values -fill both -expand yes
	# query
	frame $w.query -borderwidth 0 -highlightthickness 0
	frame $w.query.buttons -borderwidth 0 -highlightthickness 0
	foreach command {and or count region compare min max avg} {
		button $w.query.buttons.$command -text $command -command [list $object querybuilder_add $command]
		pack $w.query.buttons.$command -side left
	}
	set queryw $w.query.query
	Classy::ScrolledText $w.query.query -textvariable $var
	pack $w.query.buttons -fill both -expand yes
	pack $w.query.query -fill both -expand yes
	# fill paned
	set w $object.querybuilder.options.paned
	$w add $w.fields
	$w add $w.pane2
	set w $w.pane2
	$w add $w.operators
	$w add $w.pane3
	set w $w.pane3
	$w add $w.values
	$w add $w.query
	set w $object.querybuilder.options.paned
	$object querybuilder_filteroperators
	$object querybuilder_fillvalues
	Classy::todo $object querybuilder_browsevalue
	$object.querybuilder add query "Query" "[list $object query]; destroy $object.querybuilder"
	$object.querybuilder add clear "Clear" "[list set $var {}]"
	update
	$w.fields.fields set [$w.fields.fields get 0]
	$operatorsw set [$operatorsw get 0]
	catch {$valuesw set [$valuesw get 0]}
	after idle $object querybuilder_browsefield
}

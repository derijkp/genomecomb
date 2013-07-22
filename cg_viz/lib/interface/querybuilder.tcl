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
		if {![llength $list]} break
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
	} else {
		$valuesw set [$valuew get]
	}
	$valuew set [$valuesw get]
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
	private $object operatorsw valuesw
	if {[$operatorsw get] eq ""} {
		catch {$operatorsw set [$operatorsw get 0]}
	}
	Classy::todo $object querybuilder_fillvalues
#	if {[$valuesw get] eq ""} {
#		catch {$valuesw set [$valuesw get 0]}
#	}
#	$object querybuilder_browsevalue
}

mainw method querybuilder_browsevalue {args} {
	private $object valuesw valuew
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

mainw method querybuilder_add {command {join and}} {
	private $object fieldsw operatorsw valuew
	upvar #0 [$object.buttons.query cget -textvariable] query
	set w $object.querybuilder.options.paned
	set fields [$fieldsw get]
	set operator [$operatorsw get]
	set values [$valuew get]
	if {$command eq "count"} {
		if {![llength $fields]} {error "Select some fields first"}
		set value [lindex $values 0]
		set value [$object querybuilder_quotevalue $value]
		set insert "count\("
		foreach field $fields {
			append insert "\$$field, "
		}
		append insert "$operator $value \) > 0"
		$object querybuilder_insert $insert $join
	} elseif {$command eq "value"} {
		if {![llength $fields]} {error "Select some fields first"}
		$object querybuilder_insert "\$\{[join $fields "\}\$\{"]\}" $join
	} elseif {$command eq "region"} {
		destroy $object.region
		Classy::Dialog $object.region -title "Select region"
		$object.region option entry Chromosome [privatevar $object region(chr)]
		$object.region option entry Begin [privatevar $object region(begin)]
		$object.region option entry End [privatevar $object region(end)]
		$object.region add go Go "$object querybuilder_insert region(\[getprivate $object region(chr) \]:\[getprivate $object region(begin) \]-\[getprivate $object region(end) \]) $join" default
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
		$object.compare add go Go "$object querybuilder_insert \"\[lindex \[getprivate $object compare(type) \] 0\]\(\[join \[getprivate $object compare(selsamples) \] ,\]\)\" $join" default
	} elseif {$command in {min max lmin lmax avg sum}} {
		if {![llength $fields]} {error "Select some fields first"}
		set insert "${command}(\$[join $fields ",\$"])"
		$object querybuilder_insert $insert $join
	} elseif {$command in {percent}} {
		if {![llength $fields]} {error "Select some fields first"}
		set insert "${command}(\$[lindex $fields 0])"
		$object querybuilder_insert $insert $join
	} else {
		if {![llength $fields]} {error "Select some fields first"}
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
			$object.if add go Go "$object querybuilder_insert if(\[getprivate $object if(if) \],\[getprivate $object if(true) \],\[getprivate $object if(false) \]) $join" default
		} else {
			set insert [join $insert "\n    $command "]
			set join $command
			$object querybuilder_insert $insert $join
		}
	}
}

mainw method querybuilder {args} {
	private $object qfields qoperators qvalues fieldsw operatorsw valuesw valuew queryw qfieldsfilter valuesfilter
	set qfields [$object.tb tfields]
	set qoperators {== != < <= > >= regexp contains shares }
	set qfieldsfilter {}
	set valuesfilter {}
	destroy $object.querybuilder
	Classy::Dialog $object.querybuilder -title Query
	if {[llength $args]} {
		$object.querybuilder configure -title "Calculated field builder"
		set fieldbuilder 1
		set fieldbuilderw [lindex $args 0]
		set var [privatevar $object fieldcalc]
		Classy::Entry $object.querybuilder.options.fieldname -label Fieldname -textvariable [privatevar $object fieldname]
		pack $object.querybuilder.options.fieldname -side top -expand yes -fill x
		set funcbuttons {value if count lmin lmax avg sum compare percent and or region isnum}
		set join ""
	} else {
		set fieldbuilder 0
		set var [$object.buttons.query cget -textvariable]
		set funcbuttons {and or count region compare lmin lmax avg sum if isnum percent}
		set join and
	}
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
	button $w.operators.header -text "Operator"
	pack $w.operators.header -fill x
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
		-browsecommand [list Classy::todo $object querybuilder_browsevalue] \
		-exportselection 0 -selectmode extended
	pack $w.values.values -fill both -expand yes
	# query
	frame $w.query -borderwidth 0 -highlightthickness 0
	frame $w.query.buttons -borderwidth 0 -highlightthickness 0
	foreach command $funcbuttons {
		button $w.query.buttons.$command -text $command -command [list $object querybuilder_add $command $join]
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
	if {$fieldbuilder} {
		$object.querybuilder add query "Insert calc field" "[list $fieldbuilderw] addcalc \"\$\{[privatevar $object fieldname]\}=\$\{$var\}\"; destroy $object.querybuilder"
	} else {
		$object.querybuilder add query "Query" "[list $object query]; destroy $object.querybuilder"
	}
	$object.querybuilder add clear "Clear" "[list set $var {}]"
	$object querybuilder_fillvalues
	update
	set field [$w.fields.fields get 0]
	if {$field ne ""} {
		$w.fields.fields set $field
		Classy::todo $object querybuilder_browsefield
	}
}

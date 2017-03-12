mainw method easyquery_refreshsel {args} {
	private $object qqdata
	update idletasks
	foreach key $args {
		catch {set qqdata(sel,$key) $qqdata(sel,$key)}
	}
}

proc easyquery_grid {w key row col} {
	if {[string index $col end] eq "+"} {
		set col [string range $col 0 end-1]
		set colweight 1
	} else {
		set colweight 0
	}
	if {[string index $row end] eq "+"} {
		set row [string range $row 0 end-1]
		set rowweight 1
	} else {
		set rowweight 0
	}
	set scol [split $col -]
	set srow [split $row -]
	foreach {col colextra} $scol break
	foreach {row rowextra} $srow break
	
	grid $w.$key -row $row -column $col -sticky nwse
	if {[isint $rowextra]} {
		grid $w.$key -rowspan [expr {$rowextra+1}]
		set rowend [expr {$row+$rowextra}]
	} else {
		set rowend $row
	}
	if {[isint $colextra]} {
		grid $w.$key -columnspan [expr {$colextra+1}]
		set colend [expr {$col+$colextra}]
	} else {
		set colend $col
	}
	grid columnconfigure $w $colend -weight $colweight
	grid rowconfigure $w $rowend -weight $rowweight
}

mainw method easyquery_list {key row col {label {}} {check {}}} {
	private $object qqdata
	if {$label eq ""} {set label $key}
	set w $object.easyquery.options.paned.dialog.frame
	destroy $w.$key
	frame $w.$key -borderwidth 0 -highlightthickness 0
	if {$check eq ""} {
		label $w.$key.label -text $label -anchor w
	} else {
		checkbutton $w.$key.label -text $label -anchor w -variable [privatevar $object qqdata(check,${key})]
	}
	if {[lindex $qqdata(list,${key}) end 1] ne "incomplete"} {
		Classy::ListBox $w.$key.list \
			-listvariable [privatevar $object qqdata(list,${key})] \
			-filtervariable [privatevar $object qqdata(filter,${key})] \
			-selvariable [privatevar $object qqdata(sel,${key})] \
			-exportselection 0 -selectmode extended
	} else {
		Classy::Text $w.$key.list \
			-textvariable [privatevar $object qqdata(sel,${key})]
	}
	pack $w.$key.label -side top -fill x
	pack $w.$key.list -side top -fill both -expand yes
	easyquery_grid $w $key $row $col
}

mainw method easyquery_num {key row col {label {}} {check {}}} {
	private $object qqdata
	if {$label eq ""} {set label $key}
	set w $object.easyquery.options.paned.dialog.frame
	Classy::NumEntry $w.$key -orient vertical -label $label -textvariable [privatevar $object qqdata(sel,${key})]
	easyquery_grid $w $key $row $col
}

mainw method easyquery_entry {key row col {label {}} {check {}}} {
	private $object qqdata
	if {$label eq ""} {set label $key}
	set w $object.easyquery.options.paned.dialog.frame
	Classy::Entry $w.$key -orient vertical -label $label -textvariable [privatevar $object qqdata(sel,${key})]
	easyquery_grid $w $key $row $col
}

proc easyquery_q {list} {
	if {[llength $list] == 1} {
		return "contains \"[lindex $list 0]\""
	} elseif {[llength $list]} {
		return "shares [list $list]"
	} else {
		return "ne \"\""
	}
}

mainw method easyquery_draw_region {args} {
	private $object qqdata
	set header [$object.tb tfields]
	set poss [tsv_basicfields $header 3]
	foreach [list qqdata(col,chromosome) qqdata(col,begin) qqdata(col,end)] [list_sub $header $poss] break
	set qqdata(list,chromosome) [$object.tb values $qqdata(col,chromosome) allif0 1000 1]
	#
	$object easyquery_entry locationstring 0 0-1+ "Location string"
	$object easyquery_list chromosome 1-2+ 0
	$object easyquery_num begin 1 1+
	$object easyquery_num end 2 1+
	$object easyquery_refreshsel chromosome
}

mainw method easyquery_do_region {args} {
	private $object qqdata
	set locationstring [get qqdata(sel,locationstring) ""]
	if {$locationstring ne ""} {
		return "region(\"$locationstring\")"
	}
	set chromosome [get qqdata(sel,chromosome) ""]
	set begin [get qqdata(sel,begin) ""]
	set end [get qqdata(sel,end) ""]
	if {$chromosome eq ""} {
		error "no chromosome selected"
	}
	if {![isint $begin]} {
		error "no begin selected"
	}
	if {![isint $end]} {
		error "no end selected"
	}
	return "region(\"$chromosome:$begin-$end\")"
}

annot_init
mainw method easyquery_draw_gene {args} {
	private $object qqdata
	set dbs {refGene intGene gencode knownGene ensGene}
	set header [$object.tb tfields]
	foreach db $dbs {
		if {![inlist $header ${db}_impact] || ![inlist $header ${db}_gene]} continue
		lappend dbs $db
		set genelist [$object.tb values ${db}_gene allif0]
		set qqdata(list,${db}) $genelist
	}
	set qqdata(list,impact) [list_remdup [list_reverse [var_impact_list]]]
	#
	$object easyquery_list impact 0+ 0+
	set num 1
	foreach db $dbs {
		if {![inlist $header ${db}_impact] || ![inlist $header ${db}_gene]} continue
		$object easyquery_list $db 0+ ${num}+ $db check
		incr num
	}
	$object easyquery_refreshsel {*}$dbs
}

mainw method easyquery_do_gene {args} {
	private $object qqdata
	set dbs {refGene intGene gencode knownGene ensGene}
	set header [$object.tb tfields]
	set query {}
	set impacts [get qqdata(sel,impact) ""]
	set impactq [easyquery_q $impacts]
	foreach db $dbs {
		if {![inlist $header ${db}_impact] || ![inlist $header ${db}_gene]} continue
		set lq {}
		set genes [get qqdata(sel,$db) ""]
		if {[llength $genes]} {
			lappend lq "\$${db}_gene [easyquery_q $genes]"
		}
		if {([llength $genes] || [get qqdata(check,$db) 0]) && [llength $impacts]} {
			lappend lq "\$${db}_impact $impactq"
		} elseif {[get qqdata(check,$db) 0]} {
			lappend lq "\$${db}_gene ne \"\""
		}
		if {[llength $lq]} {
			lappend query [join $lq " and "]
		}
	}
	if {[llength $query]} {
		return \([join $query "\)\n or \("]\)
	}
	foreach db {refGene gencode knownGene ensGene} {
		lappend query "\$${db}_impact $impactq"
	}
	return \([join $query "\) or \("]\)
}

mainw method easyquery_draw_name {args} {
	private $object qqdata
	set header [$object.tb tfields]
	set fields [list_sub $header [list_find -regexp $header snp.*_name]]
	set num 0
	foreach field $fields {
		set qqdata(list,$field) [$object.tb values $field allif0 1000 0]
		$object easyquery_list $field 0+ ${num}+
		incr num
	}
	$object easyquery_refreshsel {*}$fields
}

mainw method easyquery_do_name {args} {
	private $object qqdata
	set header [$object.tb tfields]
	set fields [list_sub $header [list_find -regexp $header snp.*_name]]
	set query {}
	foreach field $fields {
		set values [get qqdata(sel,$field) ""]
		if {[llength $values]} {
			lappend query "\$$field [easyquery_q $values]"
		}
	}
	if {![llength $query]} {error "nothing selected"}
	return \([join $query "\) or \("]\)
}

mainw method easyquery_draw_frequency {args} {
	private $object qqdata
	set header [$object.tb tfields]
	set fields [list_union {1000glow snp135_freq evs_ea_freqp evs_aa_freqp evs_ea_mfreqp evs_aa_mfreqp} [list_sub $header [list_find -regexp $header {_[mh]?(freq|var)}]]]
	set qqdata(list,freqfields) [list_common $fields $header]
	if {![info exists qqdata(sel,freqp)]} {set qqdata(sel,freqp) 5}
	if {![info exists qqdata(sel,freqnum)]} {set qqdata(sel,freqnum) 2}
	if {![info exists qqdata(sel,freqfields)]} {set qqdata(sel,freqfields) $qqdata(list,freqfields)}
	$object easyquery_list freqfields 0+ 0+ "frequency fields to check"
	$object easyquery_num freqp 1 0+ "maximum freqp (frequency of alt allele in percent)"
	$object easyquery_num freqnum 2 0+ "maximum number (only used if frequency is in numbers, as indicated by a total field)"
	$object easyquery_refreshsel freqfields
}

mainw method easyquery_do_frequency {args} {
	private $object qqdata
	set header [$object.tb tfields]
	set freqp [get qqdata(sel,freqp) ""]
	set freq [expr {$freqp/100.0}]
	set freqnum [get qqdata(sel,freqnum) ""]
	if {![isdouble $freqp]} {error "freqp is not a number"}
	set fields [get qqdata(sel,freqfields) ""]
	set query {}
	foreach field $fields {
		set total 0
		if {[regsub {_(freq|var).*$} $field _total totalfield]} {
			if {[inlist $header $totalfield]} {set total 1}
		}
		if {[regexp freqp $field]} {
			set cutoff $freqp
		} else {
			set cutoff $freq
		}
		if {!$total} {
			lappend query "lmaxd(\$$field,0) <= $cutoff"
		} elseif {![isdouble $freqnum]} {
			lappend query "(double(lmaxd(\$$field,0)))/lmax(\$$totalfield,1) <= $cutoff"
		} else {
			lappend query "(double(lmaxd(\$$field,0)))/lmax(\$$totalfield,1) <= $cutoff or lmaxd(\$$field,0) <= $freqnum"
		}
	}
	if {![llength $query]} {error "nothing selected"}
	return \([join $query "\)\n and \("]\)
}

mainw method easyquery_draw {args} {
	private $object qqdata
	update idletasks
	set w $object.easyquery.options.paned.dialog.frame
	destroy $w
	frame $w -borderwidth 0 -highlightthickness 0
	pack $w -fill both -expand yes
	set type [$object.easyquery.options.paned.qlist.queries get active]
	$object easyquery_draw_$type $w
}

mainw method easyquery_do {args} {
	private $object qqdata qqcmd
	set w $object.easyquery.options.paned.qlist.queries
	update idletasks
	set type [$w get active]
	set query [$object easyquery_do_$type]
	uplevel #0 [list {*}$qqcmd $query]
}

mainw method easyquery {{cmd {}}} {
	private $object qqueries qqueriesfilter qqueriessel qqcmd
	destroy $object.easyquery
	set qqueries {region gene frequency name}
	Classy::Dialog $object.easyquery -title Query
	set wpaned $object.easyquery.options.paned
	ttk::panedwindow $wpaned -orient horizontal
	pack $wpaned -fill both -expand yes
	# queries
	set wqlist $wpaned.qlist
	frame $wqlist -borderwidth 0 -highlightthickness 0
	label $wqlist.header -text "queries"
	pack $wqlist.header -fill x
	set queriesw $wqlist.queries
	Classy::ListBox $wqlist.queries \
		-listvariable [privatevar $object qqueries] \
		-filtervariable [privatevar $object qqueriesfilter] \
		-selvariable [privatevar $object qqueriessel] \
		-command [list $object easyquery_selfield] \
		-browsecommand [list $object easyquery_draw] \
		-exportselection 0 -selectmode extended
	pack $wqlist.queries -fill both -expand yes
	$wpaned add $wqlist
	# operators
	set wd $wpaned.dialog
	frame $wd -borderwidth 0 -highlightthickness 0 -width 400 -height 300
	$wpaned add $wd
	if {$cmd eq ""} {
		set qqcmd [list $object query]
	} else {
		set qqcmd $cmd
	}
	$object.easyquery add go "Go" [list $object easyquery_do] default
	set pos [lsearch $qqueries [get qqueriessel ""]]
	if {$pos == -1} {
		set pos 0
	}
	update idletasks
	$wqlist.queries activate $pos
	$wqlist.queries set [lindex $qqueries $pos]
}

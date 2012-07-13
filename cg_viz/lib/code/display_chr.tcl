Class subclass display_chr

proc tsv_split {c} {
	set result {}
	foreach line [split $c \n] {
		lappend result [split $line \t]
	}
	return $result
}

display_chr method init {args} {
	private $object canvas dbdir full options
	foreach {canvas src dbdir} $args break
	set full [split [cg select -sh /dev/null -f "chromosome end" [glob $dbdir/extra/reg_*_fullgenome.tsv*]] \n\t]
	array set options {
		scale 1000000
		height 100
	}
	set options(src) $src
	Extral::event listen $object querychanged [list $object redraw]
}

display_chr method destroy {} {
	Extral::event remove $object querychanged [list $object redraw]
}

display_chr method options {args} {
	private $object options
	switch [llength $args] {
		0 {return [array get options]}
		1 {return $options([lindex $args 0])}
		2 {
			set options([lindex $args 0]) [lindex $args 1]
			$object redraw
			return [lindex $args 1]
		}
		default {error "Wrong # of args"}
	}
}

display_chr method select {x y} {
	private $object canvas options
	set id [$canvas find closest [$canvas canvasx $x] [$canvas canvasy $y]]
	set tagslist [list [$canvas itemcget $id -tags]]
	set scale $options(scale)
	set src $options(src)
	set table [$src table]
	$src sql [subst {delete from "sel_$table"}]
	set ids {}
	foreach tag $tagslist {
		foreach {o type chr pos} $tag break
		set chr [string range $chr 2 end]
		set pos [string range $pos 2 end]
		if {$type ne "data"} continue
		$src sql [subst {insert into "sel_$table" select "rowid" from "query_$table" where "begin"/$scale = $pos and "chromosome" = '$chr'}]
	}
	set count [$src sql [subst {select count(*) from "sel_$table"}]]
	puts "$count selected"
	Extral::event generate selchanged $src
}

proc findticks {scale} {
	set minstep [expr {50*$scale}]
	if {$minstep < 1000000} {
		set unit M
		set labelstep [expr (($minstep+1000000)/1000000)]
		set axisstep [expr {($labelstep*1000000)/$scale}]
	} elseif {$minstep < 1000} {
		set unit K
		set labelstep [expr (($minstep+1000)/1000)]
		set axisstep [expr {($labelstep*1000)/$scale}]
	} else {
	}
	if {[string index $ticks 0] >= 5} {
		set len [string length $ticks]
		if {$len < 1000000} {}
	} else {
	}
	switch [string length $ticks] {
		
	}
	return [list $step $unit]
}

display_chr method redraw {args} {
puts "$object redraw $args"
	private $object canvas full options
	set tag $object
	set scale $options(scale)
	set height $options(height)


	regsub -all {[0-9]} $scale 0 ticks
	set ticks 50
	set src $options(src)
	#
	set qfields [$src qfields]
	set chromosome [list_common $qfields {chromosome chrom}]
	if {$chromosome eq ""} {return}
	set begin [list_common $qfields {begin start pos}]
	if {$begin eq ""} {return}
	set type [list_common $qfields {type}]
	if {$type eq ""} {return}
	set max [lmath_max [$src subsql [subst {select count(*) from "temp" group by "$chromosome","pos"}] [subst {pos="$begin"/$scale}]]]
	set list [lsort -dict [tsv_split [$src subsql [subst {select "$chromosome","pos","$type",count(*) from "temp" group by "chromosome","pos","type"}] [subst {pos="$begin"/$scale}]]]]
	# $canvas delete all
	$canvas delete $tag
	
	set start 0
	$canvas create text $start -1 -anchor e -text $max -tag [list $tag axis]
	$canvas create text $start $height -anchor e -text 0 -tag [list $tag axis]
	$canvas create line $start $height $start 0 -tag [list axis $tag]
	
	foreach {chr len} $full {
		set st($chr) $start
		set nchr [chr_clip $chr]
		set st($nchr) $start
		$canvas create line $start $height [expr {$start+($len/$scale)}] $height -tag [list $tag chr c_$chr]
		$canvas create text $start [expr {$height+1}] -anchor nw -text $chr -tag [list $tag chr text c_$chr]
		set start [expr {$start+($len/$scale)+10}]
	}
	
	array set col {snp blue del red ins green}
	set pcount 0
	set ppos -1
	set pchr {}
	list_foreach {chr pos type count} $list {
		set chr [chr_clip $chr]
		if {$chr ne $pchr || $pos != $ppos} {
			set start $st($chr)
			set pcount 0
			set ppos $pos
			set pchr $chr
		}
		incr count $pcount
		set x [expr {$start+$pos}]
		$canvas create line $x [expr {$height*($max-$pcount)/$max}] $x [expr {$height*($max-$count)/$max}] -fill [get col($type) gray] -tag [list $tag data c_$chr p_$pos]
		set pcount $count
	}
	
	$canvas raise $tag&&chr
	
	$canvas configure -scrollregion [$canvas bbox all]
	$canvas bind $tag <1> [list $object select %x %y]
}

if 0 {
	set canvas .mainw.canvas.data
	set dbdir /complgen/refseq/hg18
	display_chr new disp1 $canvas $dbdir
	disp1 redraw
}


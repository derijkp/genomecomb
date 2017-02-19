proc deffieldnames {file {deffields {chromosome begin end type ref alt}}} {
	set f [gzopen $file]
	set header [tsv_open $f]
	gzclose $f
	set poss [tsv_basicfields $header 6 0]
	set nfields {}
	set qfields {}
	set used {}
	foreach pos $poss field $deffields {
		if {$pos == -1} continue
		set ffield [lindex $header $pos]
		lappend nfields $field
		if {$ffield eq $field} {
			lappend qfields $field
		} else {
			lappend qfields "$field=\$$ffield"
		}
		lappend used $ffield
	}
	lappend qfields {*}[list_lremove $header $used]
	lappend nfields {*}[list_lremove $header $used]
	return [list $nfields $qfields]
}

proc cg_intgene {args} {
	cg_options intgene args {
	} {genefile} 1 ... {
		integrate different gene sets into one
	}
	set genesets [list $genefile {*}$args]
	set tempfile [tempfile]
	set fields {chrom start end strand cdsStart cdsEnd exonStarts exonEnds}
	set todo {}
	foreach geneset $genesets {
		foreach {header qheader} [deffieldnames $geneset {chrom start end type ref alt}] break
		set missing [list_lremove $fields $header]
		if {[llength $missing]} {
			error "file $geneset misses field(s): $missing"
		}
		if {"source" ni $qheader} {
			set source [file root [file tail [gzroot $geneset]]]
			regexp {^gene_[^_]+_([^_]+)Gene} $source temp source
			set qheader [linsert $qheader 5 "source=\"$source\""]
		}
		set tempfile [tempfile]
		cg select -f $qheader $geneset $tempfile
		lappend todo $tempfile
	}
	set tempfile [tempfile]
	cg cat -m 1 -c 0 {*}$todo | cg select -s $fields > $tempfile
	set f [open $tempfile]
	set header [tsv_open $f]
	set poss [list_cor $header $fields]
	set ploc {}
	puts [join $header \t]
	while {1} {
		if {[gets $f line] == -1} break
		set loc [list_sub [split $line \t] $poss]
		if {$loc ne $ploc} {
			puts $line
			set ploc $loc
		}
	}
	close $f
}

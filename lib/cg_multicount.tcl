#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_multicount {args} {
	set idfields {geneid genename gene exon exonid id name spliceName chromosome strand start begin end}
	cg_options multicount args {
		-idfields {
			set idfields $value
		}
	} compar_file 2
	set countfiles $args

	foreach file $countfiles {
		catch {close $a(f,$file)}
	}
	unset -nocomplain a
	set header1 {}
	unset -nocomplain idposs
	foreach file $countfiles {
		set a(f,$file) [gzopen $file]
		set header [tsv_open $a(f,$file) comment]
		set a(h,$file) $header
		set poss [list_remove [list_cor $header $idfields] -1]
		if {![info exists idposs]} {
			set idposs $poss
		} else {
			set idposs [list_common $idposs $poss]
			if {![llength $idposs]} {
				error "no id field ($idfields) shared between all files"
			}
		}
		set a(id,$file) $poss
		if {![llength $poss]} {
			error "file $file has no id field, must have at least one of: $idfields"
		}
		set a(idfields,$file) [list_sub $header $poss]
		
		set a(data,$file) [list_find -glob $header *-*]
		set a(empty,$file) [list_fill [llength $a(data,$file)] 0.0]
		set a(status,$file) [gets $a(f,$file) a(curline,$file)]
		set a(curline,$file) [split $a(curline,$file) \t]
	}
	foreach file $countfiles {
		set a(curid,$file) [list_sub $a(curline,$file) $idposs]
	}
	set header [list_sub $header $idposs]
	foreach file $countfiles {
		lappend header {*}[list_sub $a(h,$file) $a(data,$file)]
	}
	set empty [list_fill [llength $idposs] {}]
	set o [open $compar_file.temp w]
	if {$comment ne ""} {puts -nonewline $o $comment}
	puts $o [join $header \t]
	while 1 {
		set curids {}
		foreach file $countfiles {
			lappend curids $a(curid,$file)
		}
		if {[llength [list_remdup $curids]] > 1} {
			error "eror: some files have different ids: $curids"
		}
		set curid [lindex $curids 0]
		if {$curid eq $empty} break
		set line $curid
		foreach file $countfiles {
			lappend line {*}[list_sub $a(curline,$file) $a(data,$file)]
			set a(status,$file) [gets $a(f,$file) a(curline,$file)]
			set a(curline,$file) [split $a(curline,$file) \t]
			set temp [list_sub $a(curline,$file) $idposs]
			set a(curid,$file) $temp
		}
		puts $o [join $line \t]
	}
	close $o
	foreach file $countfiles {
		catch {close $a(f,$file)}
	}
	file rename -force $compar_file.temp $compar_file

}

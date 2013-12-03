#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc annotatereg {file dbfile name annotfile near {outfields {name score freq}}} {

putslog [list annotatereg $file $dbfile $name $annotfile $near $outfields]
	catch {close $f}
	set f [gzopen $file]
	set header [tsv_open $f comment]
	set poss [tsv_basicfields $header 3]
	close $f
	set fields [list_sub $header $poss]
	if {[inlist $poss -1]} {
		error "Cannot annotate $file: wrong fields"
	}
	set f [gzopen $dbfile]
	set dbposs [open_region $f dbheader]
	close $f
	set dataposs [list_cor $dbheader $outfields]
	set temp [list_find $dataposs -1]
	set nh [list_sub $outfields -exclude $temp]
	set dataposs [list_sub $dataposs -exclude $temp]
	if {[llength $nh] == 0} {
		set dataposs {}
		if {$near != -1} {
			set newh ${name}_dist
		} else {
			set newh $name
		}
	} elseif {[llength $nh] == 1} {
		set newh $name
		if {$near != -1} {
			lappend newh ${name}_dist
		}
	} else {
		set newh {}
		foreach key $nh {
			lappend newh ${name}_$key
		}
		if {$near != -1} {
			lappend newh ${name}_dist
		}
	}
	set o [open $annotfile.temp w]
	puts -nonewline $o [join [list_fill [expr {[llength [split $comment \n]]-1}] \n]]
	puts $o [join $newh \t]
	close $o
	if {[gziscompressed $file]} {
		set file "|[gzcat $file] '$file'"
	}
	# puts [list {*}[gzcat $dbfile] $dbfile | reg_annot $file {*}$poss - {*}$dbposs $near {*}$dataposs]
	if {[catch {
		exec {*}[gzcat $dbfile] $dbfile | reg_annot $file {*}$poss - {*}$dbposs $near {*}$dataposs >> $annotfile.temp 2>@ stderr
	} error]} {
		if {$error ne "child killed: write on pipe with no readers"} {error $error}
	}
	file rename -force $annotfile.temp $annotfile

}

proc annotatevar {file dbfile name annotfile {outfields {name score freq}}} {
#putsvars file dbfile name annotfile outfields
	catch {close $f}
#	set f [open $file]
#	set poss [open_region $f header]
	set f [gzopen $file]
	set header [tsv_open $f comment]
	set poss [tsv_basicfields $header 3]
	close $f
	set fields [list_sub $header $poss]
	if {[inlist $poss -1]} {
		error "Cannot annotate $file: wrong fields"
	}
	set type1pos [lsearch $header type]
	if {$type1pos == -1} {
		error "$file has no type field"
	}
	set alt1pos [lsearch $header alt]
	if {$alt1pos == -1} {
		error "$file has no alt field"
	}
	set f [open $dbfile]
	set dbposs [open_region $f dbheader]
	close $f
	set type2pos [lsearch $dbheader type]
	if {$type2pos == -1} {
		error "$dbfile has no type field"
	}
	set alt2pos [lsearch $dbheader alt]
	if {$alt2pos == -1} {
		error "$dbfile has no alt field"
	}
	set dataposs [list_cor $dbheader $outfields]
	set temp [list_find $dataposs -1]
	set nh [list_sub $outfields -exclude $temp]
	set dataposs [list_sub $dataposs -exclude $temp]
	if {[llength $nh] == 0} {
		set dataposs {}
		set newh $name
	} elseif {[llength $nh] == 1} {
		set newh $name
	} else {
		set newh {}
		foreach key $nh {
			lappend newh ${name}_$key
		}
	}
	set o [open $annotfile.temp w]
	puts -nonewline $o [join [list_fill [expr {[llength [split $comment \n]]-1}] \n]]
	puts $o [join $newh \t]
	close $o
	# puts [list var_annot $file {*}$poss $type1pos $alt1pos $dbfile {*}$dbposs $type2pos $alt2pos {*}$dataposs]
	exec var_annot $file {*}$poss $type1pos $alt1pos $dbfile {*}$dbposs $type2pos $alt2pos {*}$dataposs >> $annotfile.temp 2>@ stderr
	file rename -force $annotfile.temp $annotfile
}

proc annotatebcol {file dbfile name annotfile} {
#putslog [list annotatebcol $file $dbfile $name $annotfile]
	catch {close $f}
	set f [gzopen $file]
	set header [tsv_open $f comment]
	set poss [tsv_basicfields $header 3]
	close $f
	set fields [list_sub $header $poss]
	if {[inlist $poss -1]} {
		error "Cannot annotate $file: wrong fields"
	}
	set f [open $dbfile]
	set header [tsv_open $f]
	if {$header ne {chromosome file}} {
		error "bcol database ($dbfile) should have a header of the type: chromosome file"
	}
	set bcollist {}
	set dir [file dir [file_absolute $dbfile]]
	while {![eof $f]} {
		set line [gets $f]
		if {$line eq ""} continue
		foreach {chr dbfile} [split $line \t] break
		lappend bcollist [list $chr $dir/$dbfile]
	}
	close $f
	set bcollist [ssort -natural -index 0 $bcollist]
	set newh $name
	set o [open $annotfile.temp w]
	puts -nonewline $o [join [list_fill [expr {[llength [split $comment \n]]-1}] \n]]
	puts $o $newh
	close $o
	if {[gziscompressed $file]} {
		error "bcol_annot not supported for compressed files"
	}
	# puts "bcol_annot $file $poss [list_concat $bcollist]"
	if {[catch {
		exec bcol_annot $file {*}$poss {*}[list_concat $bcollist] >> $annotfile.temp 2>@ stderr
	} error]} {
		if {$error ne "child killed: write on pipe with no readers"} {error $error}
	}
	file rename -force $annotfile.temp $annotfile
}

proc cg_annotatedb_info {dbfile {near -1}} {
	if {[file exists [gzroot $dbfile].opt]} {
		set a [dict create {*}[file_read [gzroot $dbfile].opt]]
	} else {
		set a [dict create]
	}
	if {[dict exists $a name]} {
		set name [dict get $a name]
	} else {
		set split [split [lindex [split [file root [file tail [gzroot $dbfile]]] -] 0] _]
		set name [lindex $split end]
	}
	dict set a name $name
	set dbtype [lindex [split [file tail $dbfile] _] 0]
	dict set a dbtype $dbtype
	set f [gzopen $dbfile]
	set header [tsv_open $f comment]
	close $f
	dict set a header $header
	set newh $name
	if {$dbtype eq "gene"} {
		set outfields {}
		set dataposs {}
		set poss {}
	} elseif {$dbtype eq "var"} {
		set outfields [dict_get_default $a fields {name freq score}]
		set outfields [list_common $outfields $header]
		set dataposs [list_cor $header $outfields]
		set poss [tsv_basicfields $header 3]
	} elseif {$dbtype eq "bcol"} {
		set outfields {}
		set dataposs {}
		set poss {}
	} else {
		if {[dict exists $a fields]} {
			set outfields [dict get $a fields]
		} else {
			switch -glob $name {
				1000g* {set outfields freq}
				rmsk {set outfields name}
				evofold {set outfields score}
				rnaGene {set outfields name}
				simpleRepeat {set outfields score}
				tRNAs {set outfields name}
				default {set outfields {name freq score}}
			}
		}
		set outfields [list_common $outfields $header]
		set dataposs [list_cor $header $outfields]
		set poss [tsv_basicfields $header 3]
	}
	if {$dbtype eq "gene"} {
		set newh [list ${name}_impact ${name}_gene ${name}_descr]
	} elseif {$dbtype eq "bcol"} {
		set newh [list $name]
	} else {
		switch [llength $outfields] {
		0 {
			if {$near != -1} {
				set newh ${name}_dist
			} else {
				set newh $name
			}
		}
		1 {
			set newh $name
			if {$near != -1} {
				lappend newh ${name}_dist
			}
		}
		default {
			set newh {}
			foreach key $outfields {
				lappend newh ${name}_$key
			}
			if {$near != -1} {
				lappend newh ${name}_dist
			}
		}
		}
	}
	dict set a outfields $outfields
	dict set a newh $newh
	dict set a dataposs $dataposs
	dict set a basicfields $poss
	if {![dict exists $a fields]} {
		dict set a fields [dict get $a name]
	}
	return $a
}

proc cg_annotate {args} {
	if {([llength $args] < 3)} {
		errorformat annotate
		exit 1
	}
	set near -1
	set dbdir {}
	set pos 0
	set replace 0
	foreach {key value} $args {
		switch -- $key {
			-near {
				set near $value
			}
			-dbdir {
				set dbdir $value
			}
			-name {
				set namefield $value
			}
			-replace {
				set replace $value
			}
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	foreach {file resultfile} $args break
	set file [gztemp $file]
	set dbfiles {}
	foreach testfile [lrange $args 2 end] {
		if {[file isdir $testfile]} {
			lappend dbfiles {*}[glob -nocomplain $testfile/var_*.tsv $testfile/gene_*.tsv $testfile/reg_*.tsv $testfile/bcol_*.tsv]
		} else {
			lappend dbfiles $testfile
		}
	}
	set names {}
	set newh {}
	foreach dbfile $dbfiles {
		set dbinfo [cg_annotatedb_info $dbfile $near]
		lappend names [dict get $dbinfo name]
		lappend newh {*}[dict get $dbinfo newh]
	}
	putslog "Annotating $file"
	set f [gzopen $file]
	set poss [open_region $f header]
	catch {close $f}
	set common [list_common $header $newh]
	if {[llength $common]} {
		if {!$replace} {
			puts stderr "Error: field(s) [join $common ,] already in file"
			exit 1
		}
#		foreach name $common {
#			set skip($name) 1
#		}
	} else {
		set replace 0
	}
	set afiles {}
	foreach dbfile $dbfiles {
		putslog "Adding $dbfile"
		set dbinfo [cg_annotatedb_info $dbfile $near]
		set name [dict get $dbinfo name]
		if {[info exists namefield]} {set name $namefield}
		if {[info exists skip($name)]} {
			puts "Skipping $dbfile: $name already in file"
			continue
		}
		set dbtype [lindex [split [file tail $dbfile] _] 0]
		if {$dbtype eq "gene"} {
			if {$near != -1} {error "-near option does not work with gene dbfiles"}
			if {$dbdir eq ""} {
				set dbdir [file dir [file_absolute $dbfile]]
			}
			set genomefile [lindex [glob -nocomplain $dbdir/genome_*.ifas] 0]
			if {![file exists $genomefile]} {
				puts stderr "no genomefile (genome_*.ifas) found in $dbdir, try using the -dbdir option"
				exit 1
			}
			lappend afiles $resultfile.${name}_annot
			if {[file exists $resultfile.${name}_annot]} {
				putslog "$resultfile.${name}_annot exists: skipping scan"
				continue
			}
			set genecol [dict_get_default $dbinfo genecol name2]
			set transcriptcol [dict_get_default $dbinfo transcriptcol name]
			annotategene $file $genomefile $dbfile $name $resultfile.${name}_annot $genecol $transcriptcol
		} elseif {$dbtype eq "var"} {
			if {$near != -1} {error "-near option does not work with var dbfiles"}
			set altpos [lsearch $header alt]
			if {$altpos == -1} {
				puts "Skipping: $file has no alt field"
				continue
			}
			lappend afiles $resultfile.${name}_annot
			if {[file exists $resultfile.${name}_annot]} {
				putslog "$resultfile.${name}_annot exists: skipping scan"
				continue
			}
			set outfields [dict get $dbinfo outfields]
			annotatevar $file $dbfile $name $resultfile.${name}_annot $outfields
		} elseif {$dbtype eq "bcol"} {
			if {$near != -1} {error "-near option does not work with bcol dbfiles"}
			lappend afiles $resultfile.${name}_annot
			if {[file exists $resultfile.${name}_annot]} {
				putslog "$resultfile.${name}_annot exists: skipping scan"
				continue
			}
			set outfields [dict get $dbinfo outfields]
			annotatebcol $file $dbfile $name $resultfile.${name}_annot
		} else {
			lappend afiles $resultfile.${name}_annot
			if {[file exists $resultfile.${name}_annot]} {
				putslog "$resultfile.${name}_annot exists: skipping scan"
				continue
			}
			putslog "Adding $dbfile"
			set outfields [dict get $dbinfo outfields]
			annotatereg $file $dbfile $name $resultfile.${name}_annot.temp $near $outfields
			file rename -force $resultfile.${name}_annot.temp $resultfile.${name}_annot
		}
	}
	if {$replace} {
		cg select -f [list_lremove $header $newh] $file $resultfile.temp
		exec paste $resultfile.temp {*}$afiles > $resultfile
		file delete $resultfile.temp
	} else {
		exec paste $file {*}$afiles > $resultfile
	}
	if {[llength $afiles]} {file delete {*}$afiles}
	gzrmtemp $file
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_annotate {*}$argv
}

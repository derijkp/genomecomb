#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc annotatereg {file dbfile name annotfile near dbinfo} {
# putslog [list annotatereg $file $dbfile $name $annotfile $near $dbinfo]
	set newh [dict get $dbinfo newh]
	set dataposs [dict get $dbinfo dataposs]
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
	set o [open $annotfile.temp w]
	puts -nonewline $o [join [list_fill [expr {[llength [split $comment \n]]-1}] \n] ""]
	puts $o \t[join $newh \t]
	close $o
	if {[gziscompressed $file]} {
		set file "|[gzcat $file] '$file'"
	}
	# putslog [list {*}[gzcat $dbfile] $dbfile | reg_annot $file {*}$poss - {*}$dbposs $near {*}$dataposs]
	if {[gziscompressed $dbfile]} {
		gzcatch {
			exec {*}[gzcat $dbfile] $dbfile | reg_annot $file {*}$poss - {*}$dbposs $near {*}$dataposs >> $annotfile.temp 2>@ stderr
		}
	} else {
		exec reg_annot $file {*}$poss $dbfile {*}$dbposs $near {*}$dataposs >> $annotfile.temp 2>@ stderr
	}
	file rename -force $annotfile.temp $annotfile

}

proc annotatevar {file dbfile name annotfile dbinfo} {
#putsvars file dbfile name annotfile
	set newh [dict get $dbinfo newh]
	set dataposs [dict get $dbinfo dataposs]
	catch {close $f}
	set f [gzopen $file]
	set header [tsv_open $f comment]
	set poss [tsv_basicfields $header 3]
	gzclose $f
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
	set f [gzopen $dbfile]
	set header [tsv_open $f]
	gzclose $f
	set tempdbposs [tsv_basicfields $header 6 0]
	set dbposs [lrange $tempdbposs 0 2]
	set type2pos [lindex $tempdbposs 3]
	if {$type2pos == -1} {
		error "$dbfile has no type field"
	}
	set alt2pos [lindex $tempdbposs 5]
	if {$alt2pos == -1} {
		error "$dbfile has no alt field"
	}
	set o [open $annotfile.temp w]
	puts -nonewline $o [join [list_fill [expr {[llength [split $comment \n]]-1}] \n] ""]
	puts $o \t[join $newh \t]
	close $o
	if {[gziscompressed $file]} {
		set file "|[gzcat $file] '$file'"
	}
	# puts [list var_annot $file {*}$poss $type1pos $alt1pos $dbfile {*}$dbposs $type2pos $alt2pos {*}$dataposs]
	if {[gziscompressed $dbfile]} {
		gzcatch { 
			exec {*}[gzcat $dbfile] $dbfile | var_annot $file {*}$poss $type1pos $alt1pos - {*}$dbposs $type2pos $alt2pos {*}$dataposs >> $annotfile.temp 2>@ stderr
		}
	} else {
		# puts "var_annot $file {*}$poss $type1pos $alt1pos $dbfile {*}$dbposs $type2pos $alt2pos {*}$dataposs >> $annotfile.temp"
		exec var_annot $file {*}$poss $type1pos $alt1pos $dbfile {*}$dbposs $type2pos $alt2pos {*}$dataposs >> $annotfile.temp 2>@ stderr
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
	gzclose $f
	dict set a header $header
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
				default {set outfields {name freq score annotation}}
			}
		}
		set outfields [list_common $outfields $header]
		set dataposs [list_cor $header $outfields]
		set poss [tsv_basicfields $header 3]
	}
	if {[dict exists $a headerfields]} {
		set headerfields [dict get $a headerfields]
		set olen [llength $outfields]
		set hlen [llength $headerfields]
		if {$olen != $hlen && !($olen == 0 && $hlen == 1)} {
			error "error: number of headerfields in opt file differs from number of outfields"
		}
		set newh $headerfields
		if {$near != -1} {
			lappend newh ${name}_dist
		}
	} else {
		set newh $name
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
	set near -1
	set dbdir {}
	set pos 0
	set replace 0
	set multidb 0
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
			-multidb {
				set multidb $value
			}
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {([llength $args] < 3)} {
		errorformat annotate
		exit 1
	}
	foreach {file resultfile} $args break
	set gzfile [gztemp $file]
	set file $gzfile
	set dbfiles {}
	foreach testfile [lrange $args 2 end] {
		if {[file isdir $testfile]} {
			lappend dbfiles {*}[glob -nocomplain $testfile/var_*.tsv $testfile/gene_*.tsv $testfile/mir_*.tsv $testfile/reg_*.tsv $testfile/bcol_*.tsv]
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
	set file [tsv_convert2var $file header comment]
	set poss [tsv_basicfields $header 6 0]
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
	if {[llength $header] > 10 && [llength $dbfiles] > 4} {
		set usefile [tsv_varsfile $file]
		puts "Using varfile $usefile"
	} else {
		set usefile $file 
	}
	set tempbasefile [indexdir_file $resultfile vars.tsv ok]
	set afiles {}
	foreach dbfile $dbfiles {
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
			lappend afiles $tempbasefile.${name}_annot
			if {[file exists $tempbasefile.${name}_annot]} {
				putslog "$tempbasefile.${name}_annot exists: skipping scan"
				continue
			}
			set genecol [dict_get_default $dbinfo genecol {}]
			set transcriptcol [dict_get_default $dbinfo transcriptcol {}]
			annotategene $usefile $genomefile $dbfile $name $tempbasefile.${name}_annot $genecol $transcriptcol
		} elseif {$dbtype eq "mir"} {
			if {$near != -1} {error "-near option does not work with gene dbfiles"}
			if {$dbdir eq ""} {
				set dbdir [file dir [file_absolute $dbfile]]
			}
			set genomefile [lindex [glob -nocomplain $dbdir/genome_*.ifas] 0]
			# we do not need the genomefile for plain annotation
			set genomefile {}
			lappend afiles $tempbasefile.${name}_annot
			if {[file exists $tempbasefile.${name}_annot]} {
				putslog "$tempbasefile.${name}_annot exists: skipping scan"
				continue
			}
			set genecol [dict_get_default $dbinfo genecol name]
			set transcriptcol [dict_get_default $dbinfo transcriptcol transcript]
			set extracols [dict_get_default $dbinfo extracols status]
			annotatemir $usefile $genomefile $dbfile $name $tempbasefile.${name}_annot $genecol $transcriptcol 100 1 0 $extracols
		} elseif {$dbtype eq "var"} {
			if {$near != -1} {error "-near option does not work with var dbfiles"}
			set altpos [lsearch $header alt]
			if {$altpos == -1} {
				puts "Skipping: $file has no alt field"
				continue
			}
			lappend afiles $tempbasefile.${name}_annot
			if {[file exists $tempbasefile.${name}_annot]} {
				putslog "$tempbasefile.${name}_annot exists: skipping scan"
				continue
			}
			set outfields [dict get $dbinfo outfields]
			annotatevar $usefile $dbfile $name $tempbasefile.${name}_annot $dbinfo
		} elseif {$dbtype eq "bcol"} {
			if {$near != -1} {error "-near option does not work with bcol dbfiles"}
			lappend afiles $tempbasefile.${name}_annot
			if {[file exists $tempbasefile.${name}_annot]} {
				putslog "$tempbasefile.${name}_annot exists: skipping scan"
				continue
			}
			annotatebcol $usefile $dbfile $name $tempbasefile.${name}_annot
		} else {
			lappend afiles $tempbasefile.${name}_annot
			if {[file exists $tempbasefile.${name}_annot]} {
				putslog "$tempbasefile.${name}_annot exists: skipping scan"
				continue
			}
			putslog "Adding $dbfile"
			set outfields [dict get $dbinfo outfields]
			annotatereg $usefile $dbfile $name $tempbasefile.${name}_annot.temp $near $dbinfo
			file rename -force $tempbasefile.${name}_annot.temp $tempbasefile.${name}_annot
		}
	}
	if {$multidb} {
		cg select -f id $file $resultfile.temp2
		exec paste -d "" $resultfile.temp2 {*}$afiles > $resultfile.temp
		file delete $resultfile.temp2
		file rename -force $resultfile.temp $resultfile
	} elseif {$replace} {
		cg select -f [list_lremove $header $newh] $file $resultfile.temp2
		exec paste -d "" $resultfile.temp2 {*}$afiles > $resultfile.temp
		file delete $resultfile.temp2
		file rename -force $resultfile.temp $resultfile
	} else {
		exec paste -d "" $file {*}$afiles > $resultfile.temp
		file rename -force $resultfile.temp $resultfile
	}
	if {[llength $afiles]} {file delete {*}$afiles}
	gzrmtemp $gzfile
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


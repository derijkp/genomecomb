#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

if 0 {

	package require Tclx
	signal -restart error SIGINT
	lappend auto_path /home/peter/dev/completegenomics/lib
	append env(PATH) :/home/peter/dev/completegenomics/bin
	package require Extral

	cd /complgen/projects/dlb1
	set file dlb_compar.tsv
	set resultfile compar.tsv.temp
	set args {/complgen/refseq/hg18/reg_hg18_cytoBand.tsv /complgen/refseq/hg18/reg_hg18_simpleRepeat.tsv}
	set db 1000gCEU
	set dbfile /complgen/refseq/hg18/1000g_hg18_1000gCEU.tsv
	set db simpleRepeat
	set db phastCons28P
	set dbfile /complgen/refseq/hg18/reg_hg18_${db}.tsv
	set annotfile $file.${db}_annot
	annotate $file $dbfile $annotfile

	cg annotate dlb_compar.tsv acompar.tsv /complgen/refseq/hg18/1000g_hg18_1000gCEU.tsv
cg select -q '$type == "snp" && $more5pct != ""' -f 'chromosome begin end type more5pct 1000gCEU' acompar.tsv | less
cg select -q '$type == "snp" && $more5pct != "" && $1000gCEU < 0.05' -f 'chromosome begin end type more5pct 1000gCEU' acompar.tsv | less

	cg annotate dlb_compar.tsv acompar.tsv /complgen/refseq/hg18/reg_hg18_simpleRepeat.tsv
cg select -q '$type == "snp" && $trf != "" && $simpleRepeat == ""' -f 'chromosome begin end type trf simpleRepeat' acompar.tsv | less

	cg annotate dlb_compar.tsv acompar.tsv /complgen/refseq/hg18/reg_hg18_chainSelf.tsv
cg select -q '$type == "snp" && $selfchain != "" && $chainSelf == ""' -f 'chromosome begin end type selfchain chainSelf' acompar.tsv | less

cg select -q '$type == "snp" && $repeat != "" && $rmsk == ""' -f 'chromosome begin end type repeat rmsk' acompar.tsv | less

	cg annotate part.tsv apart.tsv /complgen/refseq/hg18/reg_*phastCons*

	cg annotate dlb_compar.tsv annotdlbcompar.tsv /complgen/refseq/hg18

cd /complgen/projects/dlb1
set file dlb_compar.tsv
set resultfile compar.tsv.temp
set args /complgen/refseq/hg18/reg_hg18_simpleRepeat.tsv
set dbfile /complgen/refseq/hg18/reg_hg18_simpleRepeat.tsv
set dbposs {0 1 2}
set poss {0 1 2}
set dataposs {4 -1}

export PATH=$PATH:/home/peter/dev/completegenomics/bin

	reg_annot dlb_compar.tsv 0 1 2 /complgen/refseq/hg18/reg_hg18_simpleRepeat.tsv 0 1 2  4 -1 > temp.tsv

cd /complgen/projects/test
cg annotate test_compar.tsv antest_compar.tsv /complgen/refseq/hg18/var_hg18_snp130.tsv

cg select -f 'chromosome begin end type ref alt snp130_name snp130_freq' antest_compar.tsv

}


proc annotate {file dbfile name annotfile near {outfields {name score freq}}} {

putslog [list annotate $file $dbfile $name $annotfile $near $outfields]
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
	set dbposs [open_region $f dbheader]
	close $f
	set dataposs [list_cor $dbheader $outfields]
	set temp [list_find $dataposs -1]
	set nh [list_sub $outfields -exclude $temp]
	set dataposs [list_sub $dataposs -exclude $temp]
	if {[llength $nh] == 0} {
		set dataposs {-1 -1}
		if {$near != -1} {
			set newh ${name}_dist
		} else {
			set newh $name
		}
	} elseif {[llength $nh] == 1} {
		lappend dataposs -1
		set newh $name
		if {$near != -1} {
			lappend newh ${name}_dist
		}
	} else {
		set dataposs [lrange $dataposs 0 1]
		set newh {}
		foreach key [lrange $nh 0 1] {
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
	# puts [list reg_annot $file {*}$poss $dbfile {*}$dbposs {*}$dataposs $near]
	exec reg_annot $file {*}$poss $dbfile {*}$dbposs {*}$dataposs $near >> $annotfile.temp 2>@ stderr
	file rename $annotfile.temp $annotfile

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
		set dataposs {-1 -1}
		set newh $name
	} elseif {[llength $nh] == 1} {
		lappend dataposs -1
		set newh $name
	} else {
		set dataposs [lrange $dataposs 0 1]
		set newh {}
		foreach key [lrange $nh 0 1] {
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

proc cg_annotate {args} {
	if {([llength $args] < 3)} {
		errorformat annotate
		exit 1
	}
	set near -1
	set dbdir {}
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-near {
				set near $value
			}
			-dbdir {
				set dbdir $value
			}
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	foreach {file resultfile} $args break
	set dbfiles [lrange $args 2 end]
	if {[file isdir [lindex $dbfiles 0]]} {
		set dbfiles [lsort -dict [list_concat [glob -nocomplain [lindex $dbfiles 0]/var_*.tsv [lindex $dbfiles 0]/gene_*.tsv [lindex $dbfiles 0]/reg_*.tsv] [lrange $dbfiles 1 end]]]
	}
	set names {}
	foreach dbfile $dbfiles {
		lappend names [lindex [split [file root [file tail $dbfile]] _] end]
	}
	puts "Annotating $file"
	if {[gzroot $file] ne $file} {
		puts stderr "annotate not supported for compressed files (yet)"
		exit 1
	}
	set f [gzopen $file]
	set poss [open_region $f header]
	catch {close $f}
	set common [list_common $header $names]
	if {[llength $common]} {
		puts "Fields [join $common ,] already in file: skipping"
		foreach name $common {
			set skip($name) 1
		}
	}
	set afiles {}
	foreach dbfile $dbfiles {
		putslog "Adding $dbfile"
		unset -nocomplain a
		if {[file exists $dbfile.opt]} {array set a [file_read $dbfile.opt]}
		if {[info exists a(name)]} {
			set name $a(name)
		} else {
			set name [lindex [split [file root [file tail [gzroot $dbfile]]] _] end]
		}
		if {[info exists skip($name)]} {
			puts "Skipping $dbfile: $name already in file"
			continue
		}
		set dbtype [lindex [split [file tail $dbfile] _] 0]
		if {$dbtype eq "gene"} {
			if {$near != -1} {error "-near option does not work with gene dbfiles"}
			if {$dbdir eq ""} {
				set dbdir [file dir [file normalize $dbfile]]
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
			set genecol [get a(genecol) name2]
			set transcriptcol [get a(transcriptcol) name]
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
			if {[info exists a(fields)]} {
				set outfields $a(fields)
			} else {
				set outfields {name freq score}
			}
			annotatevar $file $dbfile $name $resultfile.${name}_annot $outfields
		} else {
			lappend afiles $resultfile.${name}_annot
			if {[file exists $resultfile.${name}_annot]} {
				putslog "$resultfile.${name}_annot exists: skipping scan"
				continue
			}
			putslog "Adding $dbfile"
			if {[info exists a(fields)]} {
				set outfields $a(fields)
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
			annotate $file $dbfile $name $resultfile.${name}_annot.temp $near $outfields
			file rename $resultfile.${name}_annot.temp $resultfile.${name}_annot
		}
	}
	exec paste $file {*}$afiles > $resultfile
	if {[llength $afiles]} {file delete {*}$afiles}
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

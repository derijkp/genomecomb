#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

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

	cg annotate dlb_compar.tsv acompar.tsv /complgen/refseq/hg18/reg_hg18_*.tsv


}


proc annotate {file dbfile annotfile {outfields {name score freq}}} {

	set name [lindex [split [file root [file tail $dbfile]] _] end]
	catch {close $f}
	set f [open $file]
	set poss [open_region $f header]
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
	puts $o [join $newh \t]
	close $o
	exec reg_annot $file {*}$poss $dbfile {*}$dbposs {*}$dataposs >> $annotfile.temp 2>@ stderr
	file rename $annotfile.temp $annotfile

}

proc cg_annotate {args} {
	if {([llength $args] < 2)} {
		puts stderr "format is: $::base variantfile resultfile ?dbfile? ..."
		puts stderr " - adds new columns with annotation to variantfile"
		exit 1
	}
	foreach {file resultfile} $args break
	set dbfiles [lrange $args 2 end]
	set afiles {}
	foreach dbfile $dbfiles {
		set name [lindex [split [file root [file tail $dbfile]] _] end]
		if {[file exists $file.${name}_annot]} {
			puts stderr "$file.${name}_annot exists: skipping scan"
			continue
		}
		puts stderr "Adding $dbfile"
		switch -glob $name {
			1000g* {set outfields freq}
			rmsk {set outfields name}
			evofold {set outfields score}
			rnaGene {set outfields name}
			simpleRepeat {set outfields score}
			snp130 {set outfields name}
			tRNAs {set outfields name}
			default {set outfields {name freq score}}
		}
		annotate $file $dbfile $file.${name}_annot $outfields
		lappend afiles $file.${name}_annot
	}
	exec paste $file {*}$afiles > $resultfile
	file delete {*}$afiles
}

if {[info exists argv]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_annotate {*}$argv
}

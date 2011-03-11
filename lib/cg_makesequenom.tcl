#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

package require BioTcl
package require Extral

if 0 {
	cd /complgen/tests/sequenom
	set archive may2009
	set compar_file test.csv
	cg makesequenom test.csv may2009 sequenom.tsv
}


proc cg_makesequenom {args} {
	set extraseq 124
	if {([llength $args] != 3)} {
		puts stderr "format is: $::base targetsfile archive resultfile"
		puts stderr " - setup sequenom experiment for variants in sel_compar_file"
		puts stderr " - targetsfile: tab delimited file containing targets with at least following columns: chromosome begin end ref alt"
		puts stderr " - archive: ensembl version to be used (eg may2009 for ncbi36)"
		puts stderr " - resultfile: sequenom assay setup results"
		exit 1
	}
	foreach {compar_file archive resultfile} $args break
	# prepare embl downloading
	set cachedir [pwd]/cache
	file mkdir $cachedir
	catch {e destroy}
	catch {rename e {}}
	EmblFile new e
	catch {f destroy}
	catch {rename f {}}
	FastaFile new f
	#
	catch {close $f}
	set f [open $compar_file]
	set header [split [gets $f] \t]
	set poss [list_cor $header {chromosome begin end ref alt}]
	if {[lindex $poss 1] == -1} {
		lset poss 1 [lsearch $header start]
	}
	set o [open $resultfile w]
	puts $o SNP_ID\tSEQUENCE
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set sub [list_sub $line $poss]
		puts stderr $sub
		foreach {chr start end ref alt} $sub break
		regsub ^chr $chr {} chr
		set name [join [list_sub $sub {0 1 2}] -]
		set fastafile $cachedir/$name.fas
		set emblfile $cachedir/$name.embl
		set estart [expr {$start-$extraseq}]
		if {$estart < 0} {set estart 0}
		if {![file exists $emblfile]} {
			set embl [ensembl_getregion $chr $estart [expr {$end+$extraseq-1}] -archive $archive]
			file_write $emblfile $embl
		}
		if {![file exists $fastafile]} {
			set fasta [ensembl_getregion $chr $estart [expr {$end+$extraseq-1}] -archive $archive -output fasta]
			file_write $fastafile $fasta
		}
		f open $fastafile
		if {![llength [f list]]} {
			file delete $emblfile
			file delete $fastafile
			error "Could not download region for: $sub"
		}
		if {![llength [e open $emblfile]]} {
			file delete $emblfile
			error "Could not download region for: $sub"
		}
		set id [lindex [f list] 0]
		set seq [f sequence 0]
		set id [lindex [e list] 0]
		set rstart [expr {$start-$estart+1}]
		set rend [expr {$rstart+$end-$start-1}]
		set test [string range $seq $rstart $rend]
		if {[string toupper $ref] ne [string toupper $test]} {
			error "ref in vars ($ref) different from ref in genome ($test) for:\n$sub"
		}
		# mask repeats and snps
		set fts [e features $id]
		foreach ft $fts {
			foreach {type loc descr} $ft break
			set floc [lindex $loc 0]
			set start [dict get $floc start]
			if {[dict exists $floc complement]} {set complement 1} else {set complement 0}
			set end [dict get $floc end]
			if {$end < $start} {set end $start}
			if {[dict exists $floc acc]} continue
			if {[inlist {repeat_region variation} $type]} {
				set base [string range $seq $start $end]
				set seq [string_replace $seq $start $end [string tolower $base]]
			}
		}
		if {$ref eq ""} {set ref -}
		if {$alt eq ""} {set alt -}
		set list [list $ref {*}[split $alt ,]]
		set list [list_change $list {{} -}]
		set seq [string_replace $seq $rstart $rend \[[join $list /]\]]
		puts $o $name\t$seq
		flush $o
	}
	close $o
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_multireg {*}$argv
}

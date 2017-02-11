#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

if 0 {
	cd /complgen/tests/sequenom
	set dbdir /complgen/refseq/hg18
	set compar_file test.csv
	set resultfile sequenom.tsv
	cg makesequenom test.csv sequenom.tsv /complgen/refseq/hg18

	cd /complgen/tests/sequenom
	set dbdir /complgen/refseq/hg19
	set compar_file data/testvars.tsv
	cg makesequenom data/testvars.tsv sequenom.tsv /complgen/refseq/hg19
}


proc cg_makesequenom {args} {
	set extraseq 124
	set freql 0
	set freqN 0.2
	set snpdbpatterns snp
	set delsize 5
	set repeats s
	cg_options makesequenom args {
		-f - --freq {
			set freql $value
		}
		-n - --freqn {
			set freqN $value
		}
		-p - --snpdbpattern {
			set snpdbpatterns $value
		}
		-d - --delsize {
			set delsize $value
		}
		-r - --repeatmasker {
			set repeats $value
		}
	} {compar_file resultfile dbdir} 3 3
	#
	catch {close $f}; catch {close $fg}
	set fg [genome_open [lindex [glob $dbdir/genome_*.ifas] 0]]
	set f [gzopen $compar_file]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header]
	set o [open $resultfile w]
	puts $o SNP_ID\tSEQUENCE
	set dbsnpfiles [gzfiles $dbdir/var_*snp*.tsv.gz]
	set dbsnpposs {}
	foreach dbsnp $dbsnpfiles {
		set dbsnpheader [cg select -h $dbsnp]
		set temp [tsv_basicfields $dbsnpheader 4]
		lappend temp [lsearch $dbsnpheader freq]
		lappend dbsnpposs $temp
	}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set sub [list_sub $line $poss]
		putslog $sub
		foreach {chr start end type ref alt} $sub break
		set chr [chr_clip $chr]
		set name [join [list_sub $sub {0 1 2}] -]
		set estart [expr {$start-$extraseq}]
		if {$estart < 0} {set estart 0}
		set eend [expr {$end+$extraseq-1}]
		set seq [genome_get $fg $chr [expr {$estart}] [expr {$eend}]]
		set rstart [expr {$start-$estart}]
		set rend [expr {$rstart+$end-$start-1}]
		set test [string range $seq $rstart $rend]
		if {[string toupper $ref] ne [string toupper $test]} {
			error "ref in vars ($ref) different from ref in genome ($test) for:\n$sub"
		}
		set seq [genome_mask $dbdir $seq $chr [expr {$estart}] [expr {$eend}] $freql $freqN $delsize $repeats $snpdbpatterns]
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

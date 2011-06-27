#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_extracthomopolymers {args} {
	foreach file $args break
	set f [open $file]
	set o stdout
	puts $o chromosome\tbegin\tend\tbase\tsize
	while {![eof $f]} {
		set name [gets $f]
		if {$name eq ""} continue
		if {![regexp {chromosome ([0-9A-Z]+)} $name temp chr]} {
			if {![regexp {chr([0-9A-Z]+)} $name temp chr]} {
				error "no chromosome found in line $name"
			}
		}
		putslog $name\n$chr
		set seq [gets $f]
		set indices [regexp -all -inline -indices {A{6,}|G{6,}|C{6,}|T{6,}} $seq]
		putslog Writing
		list_foreach {begin end} $indices {
			puts $o $chr\t$begin\t$end\t[string index $seq $begin]\t[expr {$end-$begin+1}]
		}
	}
}

if 0 {
	set file /complgen/refseq/hg18/genome_hg18.ifas
}

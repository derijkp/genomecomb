#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_extracthomopolymers {args} {
	cg_options extracthomopolymers args {
	} genomefile 1 1 {
		extract homopolymer regions from a genome sequence (ifas)
	}
	set f [open $genomefile]
	set o stdout
	puts $o chromosome\tbegin\tend\tbase\tsize
	while {![eof $f]} {
		set name [gets $f]
		if {$name eq ""} continue
		set chr [chr_clip [string range [lindex $name 0] 1 end]]
		putslog $name\n$chr
		set seq [gets $f]
		set indices [regexp -all -inline -indices {A{6,}|G{6,}|C{6,}|T{6,}} $seq]
		putslog Writing
		list_foreach {begin end} $indices {
			puts $o $chr\t$begin\t[expr {$end+1}]\t[string index $seq $begin]\t[expr {$end-$begin+1}]
		}
	}
}

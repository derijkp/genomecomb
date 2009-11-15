cd /complgen/GS00102
set file var-GS000000078-ASM.tsv
set chrpos 3

set curchr {}

set f [open $file]

while {![eof $f]} {
	set line [slit [gets $line] \t]
	set chr [lindex $line $chrpos]
	if {$chr ne $curchr} {
		if {[info exists o]} {close $o}
		set curchr $chr
		set o [open $chr-$file a]
	}
	puts $o [join $line \t]
}
if {[info exists o]} {close $o}

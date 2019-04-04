namespace eval bio {}

proc bio_geneticcodetable {{code {}}} {
	set clen [llength $code]
	if {$clen > 1} {
		upvar #0 ::bio::transl_table_$code table
		array set table $code
		return ::bio::transl_table_$code
	} else {
		if {$clen == 0} {
			set file $::appdir/res/geneticcodes/1_Standard.txt
			set num 1
		} else {
			if [isint $code] {
				set file [lindex [glob $::appdir/res/geneticcodes/${code}_*.txt] 0]
				set num $code
			} else {
				set file [lindex [glob $::appdir/res/geneticcodes/*_$code.txt] 0]
				regexp {[0-9]+} $file num
			}
		}
		if {[string length $file] == 0} {
			error "Translation table for $code not found"
		}
		upvar #0 ::bio::transl_table_$num table
		if ![info exists table] {
			array set table [file_read $file]
		}
		return ::bio::transl_table_$num
	}
}

proc seq_translate {seq {code {}}} {
	upvar #0 [bio_geneticcodetable $code] trans
	set len [string length $seq]
	set seq [string tolower $seq]
	set result ""
	for {set i 0} {$i < $len} {incr i 3} {
		set codon [string range $seq $i [expr {$i+2}]]
		append result [get trans($codon) X]
	}
	return $result
}

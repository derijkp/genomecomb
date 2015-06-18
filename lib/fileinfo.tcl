#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc fileversion {} {
	return 0.10.0
}


# potential values in 
# key	examplevalue	meaning
# genomecomb	0.10	genomecomb version
# filetype	var	type of tsv file (e.g. var for variant list)
# split	1/0	variants with different alternative alleles are on separate lines (split = 1)
# ref	hg19	reference database/genome
# field	name 

# todo: convert vcf data
proc comment2dict {comment} {
	set result [dict create]
	set order {}
	foreach line [split [string trim $comment] \n] {
		if {[string index $line 0] ne "#"} continue
		set line [string range $line 1 end]
		set pos [string first \t $line]
		if {$pos != -1} {
			set key [string range $line 0 [expr {$pos-1}]]
			set value [string range $line [expr {$pos+1}] end]
		} else {
			set pos [string first = $line]
			if {$pos != -1} {
				set key [string range $line 0 [expr {$pos-1}]]
				set value [string range $line [expr {$pos+1}] end]
			} else {
				set key {}
				set value $line
			}
		}
		set key [string trim $key]
		lappend order $key
		dict lappend result $key $value
	}
	dict set result order [list_remdup $order]
	return $result
}

proc dict2comment {infodict} {
	if {[dict exists $infodict order]} {
		set order [dict get $infodict order]
	} else {
		set order {}
	}
	set keys [dict keys $infodict]
	set keys [list_remove $keys order]
	set keys [list_common [list_union $order $keys] $keys]
	set result {}
	foreach key $keys {
		foreach value [dict get $infodict $key] {
			append result \n\#$key\t$value
		}
	}
	return $result
}

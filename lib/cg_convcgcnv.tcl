#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc convcgcnv {srcfile dstfile} {
	set fields [cg select -h $srcfile]
	set nfields [list_concat {chr begin end {type=($calledCNVType == "-")? "del" : (($calledCNVType == "+")? "amp" : $calledCNVType)}} [list_remove $fields chr begin end calledCNVType]]
	cg select -overwrite 1 -s - -f $nfields -q {$calledCNVType != "hypervariable" && $calledCNVType != "\=" && $calledCNVType != "invariant"} $srcfile $dstfile
}

proc cg_convcgcnv {args} {
	if {[llength $args] != 2} {
		errorformat convcgcnv
	}
	foreach {srcfile dstfile} $args break
	set tempfile [filetemp $dstfile]
	convcgcnv $srcfile $tempfile
	file rename -force $tempfile $dstfile
}

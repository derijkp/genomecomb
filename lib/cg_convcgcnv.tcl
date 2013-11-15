#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc convcgcnv {srcfile dstfile} {
	set fields [cg select -h $srcfile]
	set nfields [list_concat {chr begin end {type=($calledCNVType == "-")? "del" : (($calledCNVType == "+")? "amp" : $calledCNVType)}} [list_remove $fields chr begin end calledCNVType]]
	cg select -s - -f $nfields -q {$calledCNVType != "hypervariable" && $calledCNVType != "\=" && $calledCNVType != "invariant"} $srcfile $dstfile.temp
	file rename -force $dstfile.temp $dstfile
}

proc cg_convcgcnv {args} {
	if {[llength $args] != 2} {
		errorformat convcgcnv
		exit 1
	}
	foreach {srcfile dstfile} $args break
	convcgcnv $srcfile $dstfile.temp
	file rename -force $dstfile.temp $dstfile
}

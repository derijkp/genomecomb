#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc 10x2tsv {10xdir outfile} {
	unset -nocomplain featuresa
	set featuresfile [gzfile $10xdir/features.tsv $10xdir/*features.tsv]
	set f [gzopen $featuresfile]
	set num 1
	while {[gets $f line] != -1} {
		set line [split $line \t]
		set gene [lindex $line 1]
		set featuresa($num) $line
		incr num
	}
	close $f
	unset -nocomplain barcodesa
	set barcodesfile [gzfile $10xdir/barcodes.tsv $10xdir/*barcodes.tsv]
	set f [gzopen $barcodesfile]
	set num 1
	while {[gets $f line] != -1} {
		set line [split $line \t]
		set barcode [lindex $line 0]
		regsub -- -1 $barcode {} barcode
		set barcodesa($num) $barcode
		incr num
	}
	close $f
	unset -nocomplain mtxa
	set matrixfile [gzfile $10xdir/matrix.mtx $10xdir/*matrix.mtx]
	set f [gzopen $matrixfile]
	while {[gets $f line] != -1} {
		if {[string index $line 0] ne "%"} break
	}
	if {$outfile eq "-"} {
		set o stdout
	} else {
		set o [wgzopen $outfile]
	}
	puts $o [join {geneid gene type cell count} \t]
	set num 1
	while {[gets $f line] != -1} {
		foreach {gene bc count} $line break
		puts $o [join $featuresa($gene) \t]\t$barcodesa($bc)\t$count
	}
	gzclose $o
	gzclose $f
}

proc cg_10x2tsv {args} {
	set 10xdir -
	set outfile -
	set transcripts 0
	set ignorecodon 0
	cg_options gtf2tsv args {
	} {10xdir outfile} 1 2
	proc gtf2tsv_parse_attr {attributes} {
		set a [dict create {*}[split [string_change $attributes {"; " ";"}] {;=}]]
	}
	10x2tsv $10xdir $outfile
}


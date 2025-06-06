#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc tsv210x {args} {
	set tsvfile -
	set 10xdir -
	set genefields -
	set cellbarcodefield {} 
	set countfield {}
	set round 1
	set remdups 1
	cg_options tsv210x args {
		-featurefields - -genefields {set genefields $value}
		-cellbarcodefield {set cellbarcodefield $value}
		-countfield {set countfield $value}
		-round {set round $value}
		-remdups {set remdups [true $value]}
	} {tsvfile 10xdir} 2 2
	mkdir $10xdir.temp
	set f [gzopen $tsvfile]
	set header [tsv_open $f]
	if {$genefields eq "-"} {
		set genefields [findfields $header {geneid gene}]
	}
	set geneposs {}
	foreach field $genefields {
		lappend geneposs [lsearch $header $field]
	}
	set poss [lsearch $geneposs -1]
	set geneposs [list_sub $geneposs -exclude $poss]
	set genefields [list_sub $header $geneposs]
	if {[llength $genefields] < 2} {
		error "should have at least 2 genefields/featurefields, only got: $genefields"
	}
	if {$cellbarcodefield eq ""} {
		foreach field {cell cellbarcode barcode} {
			set barcodepos [lsearch $header $field]
			if {$barcodepos != -1} break
		}
		if {$barcodepos == -1} {error "cellbarcodefield not found (tried: cell cellbarcode barcode)"}
	} else {
		set barcodepos [lsearch $header $cellbarcodefield]
		if {$barcodepos == -1} {error "cellbarcodefield $cellbarcodefield not found"}
	}
	if {$countfield eq ""} {
		foreach field {icount count} {
			set countpos [lsearch $header $field]
			if {$countpos != -1} break
		}
		if {$countpos == -1} {error "countfield not found (tried: icount count)"}
	} else {
		set countpos [lsearch $header $countfield]
		if {$countpos == -1} {error "countfield $countfield not found"}
	}
	set exclude [list $barcodepos $countpos]
	unset -nocomplain barcodea
	unset -nocomplain genea
	set barcodenum 0
	set genenum 0
	set num 0
	set om [wgzopen $10xdir.temp/matrix.mtx.temp.zst]
	set ob [wgzopen $10xdir.temp/barcodes.tsv.gz]
	set og [wgzopen $10xdir.temp/features.tsv.gz]
	unset -nocomplain dupa
	while {[gets $f line] != -1} {
		set line [split $line \t]
		set barcode [lindex $line $barcodepos]
		set count [lindex $line $countpos]
		set geneinfo [list_sub $line $geneposs]
		if {[lindex $geneinfo end] eq ""} {
			lset geneinfo end [lindex $geneinfo 0]
		}
		if {![info exists barcodea($barcode)]} {
			set barcodea($barcode) [incr barcodenum]
			puts $ob $barcode
		}
		if {![info exists genea($geneinfo)]} {
			set genea($geneinfo) [incr genenum]
			puts $og [join $geneinfo \t]\tGene\ Expression
		}
		set numgene $genea($geneinfo)
		set numbarcode $barcodea($barcode)
		if {$remdups} {
			if {[info exists dupa($numgene,$numbarcode)]} {
				continue
			} else {
				set dupa($numgene,$numbarcode) 1
			}
		}
		if {$round} {
			if {$count <= 0.01} {
				set count 0
			} elseif {$count <= 1} {
				set count 1
			} else {
				set count [expr {round($count)}]
			}
		}
		puts $om $numgene\ $numbarcode\ $count
		incr num
	}
	gzclose $om ; gzclose $ob ; gzclose $og
	close $f

	set om [wgzopen $10xdir.temp/matrix.mtx.gz]
	puts $om {%%MatrixMarket matrix coordinate integer general}
	puts $om {%metadata_json: {"format_version": 2}}
	puts $om $genenum\ $barcodenum\ $num
	set f [gzopen $10xdir.temp/matrix.mtx.temp.zst]
	fcopy $f $om
	close $f
	close $om
	file delete $10xdir.temp/matrix.mtx.temp.zst

	catch {file delete $10xdir.old}
	catch {file rename $10xdir $10xdir.old}
	file rename $10xdir.temp $10xdir
}

proc cg_tsv210x {args} {
	tsv210x {*}$args
}


#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc tsv210x {tsvfile 10xdir {genefields {geneid gene_id gene}} {cellbarcodefield {}} {countfield {}}} {
	mkdir $10xdir.temp
	set f [gzopen $tsvfile]
	set header [tsv_open $f]
	set geneposs [list_cor $header $genefields]
	set poss [lsearch $geneposs -1]
	set geneposs [list_sub $geneposs -exclude $poss]
	set genefields [list_sub $header $geneposs]
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
	while {[gets $f line] != -1} {
		set line [split $line \t]
		set barcode [lindex $line $barcodepos]
		set count [lindex $line $countpos]
		set geneinfo [list_sub $line $geneposs]
		if {![info exists barcodea($barcode)]} {
			set barcodea($barcode) [incr barcodenum]
			puts $ob $barcode
		}
		if {![info exists genea($geneinfo)]} {
			set genea($geneinfo) [incr genenum]
			puts $og [join $geneinfo \t]\tGene\ Expression
		}
		if {$count <= 0.01} {
			set count 10
		} elseif {$count <= 1} {
			set count 1
		} else {
			set count [expr {round($count)}]
		}
		puts $om $genea($geneinfo)\ $barcodea($barcode)\ $count
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
	set 10xdir -
	set tsvfile -
	set genefields {geneid gene gene_id}
	set cellbarcodefield {} 
	set countfield {}
	cg_options tsv210x args {
		-genefields {set genefields $value}
		-cellbarcodefield {set cellbarcodefield $value}
		-countfield {set countfield $value}
	} {tsvfile 10xdir} 1 2
	tsv210x $tsvfile $10xdir $genefields $cellbarcodefield $countfield
}


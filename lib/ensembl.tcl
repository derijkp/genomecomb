#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require http

proc ensembl_getregion {chr start end args} {
	set chr [chr_clip $chr]
	if {$chr eq "M"} {set chr MT}
	set url www.ensembl.org
	set species Homo_sapiens
	set fts {repeat variation gene vegagene estgene}
	set field param
	set output embl
	foreach {key value} $args {
		switch $key {
			-archive {
				if {$value ne ""} {
					set url $value.archive.ensembl.org
				} else {
					set url www.ensembl.org
				}
				set field st
			}
			-species {
				set species $value
			}
			-fts {
				set fts $value
			}
			-output {
				set output $value
			}
			default {
				error "unkown option $key"
			}
		}
	}
	set query [subst {output=$output;r=$chr:$start-$end;strand=1;genomic=soft_masked;_format=Text}]
	foreach ft $fts {
		append query \;$field=$ft
	}
	set h [http::geturl http://$url/$species/Location/Export?$query]
	set data [http::data $h]
	http::cleanup $h
	if {[regexp {^<.*Internal Server Error} $data]} {
		error "error getting ensembl region $chr-$start-$end"
	}
	return $data
}

proc cg_getembl {args} {
	global scriptname action
	if {[llength $args] < 2 && [llength $args] > 3} {
		puts stderr "format is: $scriptname $action regionsfile archive ?extraseq?"
		puts stderr " - get emblfiles"
		exit 1
	}
	foreach {regionsfile archive extraseq} $args break
	if {![isint $extraseq]} {set extraseq 50000}
	set list [file_read $regionsfile]
	set regions [split $list \n]
	list_shift regions
	list_foreach {cchr cstart cend} $regions {
		set cchr [string range $cchr 3 end]
		set filename $cchr-[expr {$cstart-$extraseq}]-[expr {$cend+$extraseq}].embl
		if {[file exists $filename]} {
			puts "$filename already done"
			continue
		}
		puts $cchr:$cstart-$cend
		set embl [ensembl_getregion $cchr [expr {$cstart-$extraseq}] [expr {$cend+$extraseq}] -archive $archive]
		file_write $filename $embl
	}
	
}

proc ncbi_getgene {id} {
	if {$id eq ""} {return ""}
	set h [http::geturl http://www.ncbi.nlm.nih.gov/gene/?term=$id]
	set data [http::data $h]
	http::cleanup $h
	regexp {<title>([^<>]+)</title>} $data temp descr
	set descr [string trim $descr]
}

proc cg_addgeneinfo {args} {
	global scriptname action
	if {[llength $args] != 1} {
		puts stderr "format is: $scriptname $action file"
		puts stderr " - add gene description"
		exit 1
	}
	foreach {filename} $args break
	set f [open $filename]
	set line [split [gets $f] \t]
	set idpos [lsearch $line geneId]
	puts [join $line \t]\tgenedescr
	while {![eof $f]} {
		set line [split [gets $f] \t]
		set id [lindex $line $idpos]
		set descr [ncbi_getgene $id]
		puts [join $line \t]\t$descr
	}
}


if 0 {
	http://www.ensembl.org/Homo_sapiens/Location/Export?output=embl;r=10:1000000-1001000;strand=1;param=similarity;param=repeat;param=genscan;param=contig;param=variation;param=marker;param=gene;param=vegagene;param=estgene;_format=Text
	wget http://may2009.archive.ensembl.org/Homo_sapiens/Location/Export?output=embl;r=1:100000-110000;strand=1;time=1258453242.93979;st=similarity;st=repeat;st=genscan;st=contig;st=variation;st=marker;st=gene;st=vegagene;st=estgene;_format=Text
	package require http
	ensembl_getregion 10 1000000 1001000 -archive may2009
	ensembl_getregion 10 1000000 1001000
}

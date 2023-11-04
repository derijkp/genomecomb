proc cg_multicount {args} {
	set potentialidfields {geneid genename gene exon exonid id name cell cellbarcode spliceName chromosome strand start begin end}
	set emptyvalue 0
	set clip_geneid_version 1
	cg_options multicount args {
		-idfields {
			set idfields $value
		}
		-empty {
			set emptyvalue $value
		}
		-clip_geneid_version {
			set clip_geneid_version $value
		}
	} compar_file 2
	set countfiles $args

	catch {gzclose $f} ; catch {gzclose $o}
	unset -nocomplain a
	set commonidfields $potentialidfields
	foreach file $countfiles {
		set f [gzopen $file]
		set header [tsv_open $f comment]
		set a(h,$file) $header
		set poss [list_remove [list_cor $header $potentialidfields] -1]
		if {![llength $poss]} {
			error "file $file has no id field, must have at least one of: $potentialidfields"
		}
		set idfields [list_sub $header $poss]
		set a(idfields,$file) $idfields
		set commonidfields [list_common $commonidfields $idfields]
		set a(chrfield,$file) [lindex [list_common $header {chromosome chr}] 0]
		set a(chrpos,$file) [lsearch $header $a(chrfield,$file)]
		set a(geneidfield,$file) [lindex [list_common $header {geneid gene_id}] 0]
		set a(geneidpos,$file) [lsearch $header $a(geneidfield,$file)]
		gzclose $f
	}
	foreach file $countfiles {catch {gzclose $a(f,$file)}}
	set newheader $commonidfields

	foreach file $countfiles {
		set rootname [file_rootname $file]
		set tempfile [tempfile].zst
		set fields $a(h,$file)
		if {$a(chrfield,$file) != ""} {
			lset fields $a(chrpos,$file) "chromosome=chr_clip(\$$a(chrfield,$file))"
		}
		if {$clip_geneid_version && $a(geneidfield,$file) != ""} {
			set temp "geneid=regsub(\$$a(geneidfield,$file),\"\\\[.\\\]\\\[0-9\\]+\\\$\",\"\")"
			lset fields $a(geneidpos,$file) $temp
		}
		cg select -stack 1 -overwrite 1 -f $fields -s $commonidfields $file $tempfile
		set a(f,$file) [gzopen $tempfile]
		set header [tsv_open $a(f,$file) comment]
		set poss [list_remove [list_cor $header $commonidfields] -1]
		set a(id,$file) $poss
		set a(dataposs,$file) [list_find -glob $header *-*]
		if {[llength $a(dataposs,$file)]} {
			set datafields [list_sub $header $a(dataposs,$file)]
		} else {
			set datafields [list_lremove $header $a(idfields,$file)]
			set a(dataposs,$file) [list_cor $header $datafields]
		}
		set a(data,$file) $datafields
		set a(dataempty,$file) [list_fill [llength $datafields] {}]
		foreach field $datafields {
			if {[regexp -- - $field]} {
				lappend newheader $field
			} else {
				lappend newheader ${field}-$rootname
			}
		}
		set a(empty,$file) [list_fill [llength $a(data,$file)] $emptyvalue]
		set a(status,$file) [gets $a(f,$file) a(curline,$file)]
		set a(curline,$file) [split $a(curline,$file) \t]
		set a(curid,$file) [list_sub $a(curline,$file) $a(id,$file)]
	}
	set empty [list_fill [llength $commonidfields] {}]

	set tempfile $compar_file.temp[gzext $compar_file]
	set o [wgzopen $tempfile]
	if {$comment ne ""} {puts -nonewline $o $comment}
	puts $o [join $newheader \t]
	while 1 {
		set curids {}
		unset -nocomplain cura
		foreach file $countfiles {
			set curid $a(curid,$file)
			lappend curids $curid
			lappend cura($curid) $file
		}
		set firstcurid [lindex [bsort $curids] 0]
		if {$firstcurid eq $empty} break
		set line $firstcurid
		foreach file $countfiles curid $curids {
			if {$curid ne $firstcurid} {
				lappend line {*}$a(empty,$file)
				continue
			}
			lappend line {*}[list_sub $a(curline,$file) $a(dataposs,$file)]
			set a(status,$file) [gets $a(f,$file) a(curline,$file)]
			set a(curline,$file) [split $a(curline,$file) \t]
			set temp [list_sub $a(curline,$file) $a(id,$file)]
			set a(curid,$file) $temp
		}
		puts $o [join $line \t]
	}
	close $o
	foreach file $countfiles {
		catch {close $a(f,$file)}
	}
	file rename -force $tempfile $compar_file
}

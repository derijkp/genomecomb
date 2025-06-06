proc annotatebcol {file dbfile name annotfile} {
# putslog [list annotatebcol $file $dbfile $name $annotfile]
# putsvars file dbfile name annotfile
	catch {close $f}
	set f [gzopen $file]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	close $f
	if {[inlist $poss -1]} {
		error "Cannot annotate $file: wrong fields"
	}
	set fields [list_sub $header $poss]
	set f [open $dbfile]
	set header [tsv_open $f]
	close $f
	if {$header ne {chromosome begin end}} {
		error "bcol database ($dbfile) should have a header of the type: chromosome begin end, old style bcols (single chr not supported)"
	}
	set binfile [gzfile $dbfile.bin]
	set newh $name
	set tempannotfile [filetemp $annotfile]
	set o [open $tempannotfile w]
	puts $o $newh
	close $o
	if {[gziscompressed $file]} {
		error "bcol_annot not supported for compressed files"
	}
	# puts stderr [list bcol_annot $file {*}$poss -1 -1 $dbfile -1]
	if {[catch {
		exec bcol_annot $file {*}$poss -1 -1 $dbfile -1 >> $tempannotfile 2>@ stderr
	} error]} {
		if {$error ne "child killed: write on pipe with no readers"} {error $error}
	}
	file rename -force -- $tempannotfile $annotfile
}

proc annotatebcolvar {file dbfile name annotfile} {
#putslog [list annotatebcol $file $dbfile $name $annotfile]
	catch {close $f}
	set f [gzopen $file]
	set header [tsv_open $f]
	close $f
	set poss [tsv_basicfields $header 6 0]
	set poss [list_sub $poss {0 1 2 3 5}]
	if {[inlist $poss -1]} {
		error "Cannot annotate $file using $dbfile: wrong fields, must contain fields chromosome,begin,end,alt"
	}
	set fields [list_sub $header $poss]
	set f [open $dbfile]
	set header [tsv_open $f]
	if {$header ne {chromosome begin end}} {
		error "bcol database ($dbfile) should have a header of the type: chromosome begin end, old style bcols (single chr not supported)"
	}
	set binfile [gzfile $dbfile.bin]
	set newh $name
	set o [open $annotfile.temp w]
	puts $o $newh
	close $o
	if {[gziscompressed $file]} {
		error "bcol_annot not supported for compressed files"
	}
	# puts stderr "bcol_annot $file $poss $dbfile -1"
	if {[catch {
		exec bcol_annot $file {*}$poss $dbfile -1 >> $annotfile.temp 2>@ stderr
	} error]} {
		if {$error ne "child killed: write on pipe with no readers"} {error $error}
	}
	file rename -force -- $annotfile.temp $annotfile
}


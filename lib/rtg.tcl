proc rtg2annotvar {file {outfile {}}} {
	if {$outfile eq ""} {
		set outfile [file dir $file]/annotvar-[file root [rzroot [file tail $file]]].tsv
	}
	catch {close $f}; catch {close $o}
	set f [rzopen $file]
	set o [open $outfile.temp w]
	puts $o [join {chromosome begin end type alleleSeq1 alleleSeq2 posterior coverage correction support_statistics reference} \t]
	while {![eof $f]} {
		set line [gets $f]
		if {[string index $line 0] ne "#"} break
		set header $line
	}
	set header [split [string range $header 1 end] \t]
	if {$header ne "name position type reference prediction posterior coverage correction support_statistics"} {
		error "$file is not a correct rtg file"
	}
	set line [split $line \t]
	set num 0
	while 1 {
		incr num
		if {![expr $num%100000]} {putslog $num}
		if {[eof $f]} break
		set line [split [gets $f] \t]
		foreach {name position het reference prediction posterior coverage correction support_statistics} $line break
		regsub ^chr $name {} chromosome
		set alleles [split $prediction :]
		if {[llength $alleles] == 1} {lappend alleles [lindex $alleles 0]}
		set begin $position
		set end [expr {$position+1}]
		puts $o "$chromosome\t$begin\t$end\tsnp\t[join $alleles \t]\t$posterior\t$coverage\t$correction\t$support_statistics\t$reference"
	}
	close $o
	close $f
	file rename $outfile.temp $outfile
}

if 0 {
	set file /complgen/rtg/GS00102-D06.snps.gz
	set dir /complgen/rtg/rtgGS102/
}

proc process_rtgsample {dir dbdir {force 0}} {
	set keepdir [pwd]
	set dir [file normalize $dir]
	cd $dir
	puts stderr "Processing sample $dir"
	set name [file tail $dir]
	set varfile [lindex [glob *snps*] 0]
	# annotated vars file
	if {$force || ![file exists annotvar-$name.tsv]} {
		puts stderr "Create annotated varfile annotvar-$name.tsv"
		rtg2annotvar $varfile temp.tsv
		cg select -s "chromosome begin" < temp.tsv > stemp.tsv
		file rename -force stemp.tsv annotvar-$name.tsv
	}
	# sample specific filters
	if {$force || ![file exists reg_cluster-$name.tsv]} {
		puts stderr "Find cluster regions for annotvar-$name.tsv"
		cg clusterregions < annotvar-$name.tsv > temp.tsv
		file rename temp.tsv reg_cluster-$name.tsv
	}
	if {$force || ![file exists fannotvar-$name.tsv]} {
		# add filterdata to annotvar
		set todo {}
		lappend todo [list cluster cl reg_cluster-$name.tsv]
		lappend todo [list trf trf $dbdir/regdb-simple_repeats.tsv]
		lappend todo [list str str $dbdir/regdb-microsatelite.tsv]
		lappend todo [list segdup sd $dbdir/regdb-segdups.tsv]
		lappend todo [list selfchain sc $dbdir/regdb-selfchain.tsv]
		lappend todo [list repeat rp $dbdir/regdb-repeatmasker.tsv]
		lappend todo [list rna rna $dbdir/regdb-rnagenes.tsv]
		foreach file [lsort -dictionary [glob $dbdir/checked*.tsv]] {
			set value [lindex [split [file root [file tail $file]] _] 1]
			if {$value eq ""} {set value checked}
			lappend todo [list checked $value $file]
		}
		annot_annotvar annotvar-$name.tsv fannotvar-$name.tsv $todo
	}
	# fake region file
	if {$force || ![file exists sreg-$name.tsv]} {
		set chrs [list_fill 22 1 1]
		lappend chrs M X Y
		set o [open sreg-$name.tsv w]
		puts $o	[join {chromosome begin end} \t]
		foreach chr $chrs {
			puts $o $chr\t1\t1000000000
		}
		close $o
	}
	cd $keepdir
}


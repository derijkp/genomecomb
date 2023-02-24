proc exons_startsends2list {starts ends {sizeVar {}}} {
	if {$sizeVar ne ""} {upvar $sizeVar size}
	set size 0
	set exons {}
	foreach begin [split $starts ,] end [split $ends ,] {
		if {$begin eq ""} continue
		set size [expr {$size + ($end - $begin)}]
		lappend exons $begin-$end
	}
	return $exons
}

proc iso_name {chromosome strand exonStarts exonEnds {sizeVar {}}} {
	if {$sizeVar ne ""} {upvar $sizeVar size}
	set size 0
	set starts [split [string trimright $exonStarts ,] ,]
	set begin [lindex $starts 0]
	set newname novelt_${chromosome}_${begin}${strand}
	set pe -1
	foreach s $starts e [split [string trimright $exonEnds ,] ,] {
		if {$s eq ""} continue
		if {$pe != -1} {
			append newname i[expr {$s-$pe}]
		}
		set pe $e
		set esize [expr {$e-$s}]
		append newname e$esize
		set size [expr {$size + $esize}]
	}
	return $newname
}

array set _gene_name_strandnamea {+ p - m . u}
proc gene_name {chromosome strand begin end} {
	return novelg_${chromosome}_[get ::_gene_name_strandnamea($strand) $strand]_${begin}_${end}
}

proc iso_combine_job {projectdir isocaller {iso_match {}}} {
	upvar job_logdir job_logdir
	# combined analysis
	cd $projectdir
	mkdir compar
	set exproot [file tail $projectdir]
	if {$isocaller eq "*"} {
		set root $exproot
	} else {
		set root ${isocaller}-$exproot
	}
	set isoformfiles {}
	foreach isoformfile [bsort [jobglob samples/*/isoform_counts-${isocaller}-*.tsv]] {
		if {[file size $isoformfile] == 0} continue
		lappend isoformfiles $isoformfile
	}
	
	if {[llength $isoformfiles]} {
		job iso_compar-isoform_counts-$root \
		-deps $isoformfiles \
		-targets {
			compar/isoform_counts-$root.tsv
		} -vars {
			isoformfiles exproot root isocaller iso_match
		} -code {
			set isoformcounts compar/isoform_counts-$root.tsv
			cg multitranscript -match $iso_match $isoformcounts {*}$isoformfiles
		}
	}
	set genefiles [bsort [jobglob samples/*/gene_counts-${isocaller}-*.tsv]]
	if {[llength $genefiles]} {
		job iso_compar-gene_counts-$root \
		-deps $genefiles \
		-targets {
			compar/gene_counts-$root.tsv
		} -vars {
			genefiles exproot root isocaller
		} -code {
			set genecounts compar/gene_counts-$root.tsv
			cg_multigene $genecounts {*}$genefiles
		}
	}
	set totalcountsfiles [bsort [jobglob samples/*/totalcounts-${isocaller}-*.tsv]]
	if {[llength $totalcountsfiles]} {
		job iso_compar-totalcounts-$root \
		-deps $totalcountsfiles \
		-targets {
			compar/totalcounts-$root.tsv
		} -vars {
			totalcountsfiles exproot root isocaller
		} -code {
			analysisinfo_combine compar/totalcounts-$root.tsv $totalcountsfiles
			cg paste {*}$totalcountsfiles > compar/totalcounts-$root.tsv
		}
	}
}

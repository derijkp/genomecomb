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
	set isoformfiles [bsort [jobglob -checkcompressed 1 samples/*/isoform_counts-${isocaller}-*.tsv]]	
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
	set genefiles [bsort [jobglob -checkcompressed 1 samples/*/gene_counts-${isocaller}-*.tsv]]
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
	set totalcountsfiles [bsort [jobglob -checkcompressed 1 samples/*/totalcounts-${isocaller}-*.tsv]]
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

proc iso_joint_job {args} {
	upvar job_logdir job_logdir
	set iso_joint {}
	set iso_match {}
	set dbdir {}
	set threads 2
	set distrreg 0
	set cleanup 1
	set reftranscripts {}
	cg_options iso_joint args {
		-reftranscripts {
			set reftranscripts [code_empty $value]
		}
		-iso_joint {
			set iso_joint $value
		}
		-iso_match {
			set iso_match $value
		}
		-threads {
			set threads $value
		}
		-distrreg {
			set distrreg [distrreg_checkvalue $value]
		}
		-dbdir - -refdir {
			set dbdir $value
		}
		-c - -cleanup {
			set cleanup $value
		}
	} {projectdir}
	set dbdir [file_absolute $dbdir]
	set refseq [refseq $dbdir]
	# combined analysis
	cd $projectdir
	mkdir compar
	set exproot [file tail $projectdir]
	foreach isocaller $iso_joint {
		set source [file_absolute compar/isoform_counts-$isocaller-$exproot.tsv]
		set ref [file_absolute compar/isoform_counts-$isocaller-$exproot-ref.tsv]
		job isojoint_ref_$isocaller -deps {
			$source $reftranscripts
		} -targets {
			$ref
		} -vars {
			source ref reftranscripts
		} -code {
			tsv_makeregular $source $ref.temp {}
			tsv_makeregular $reftranscripts $ref.temp2 {}
			cg cat -m 1 $ref.temp2 $ref.temp | cg select -s - > $ref.temp3
			set o [open $ref.temp4 w]
			set f [open $ref.temp3]
			set header [tsv_open $f comment]
			puts $o $comment[join $header \t]
			set pos [findfieldspos $header transcript]
			unset -nocomplain a
			while 1 {
				if {[gets $f line] == -1} break
				set transcript [lindex [split $line \t] $pos]
				if {[info exists a($transcript)]} continue
				puts $o $line
				set a($transcript) 1
			}
			file rename -force $ref.temp4 $ref
			file delete -force $ref.temp $ref.temp2 $ref.temp3
		}
		foreach bam [jobglob $projectdir/samples/*/*.bam $projectdir/samples/*/*.cram] {
			set preset {}
			foreach {baseisocaller preset} [split $isocaller _] break
			if {![auto_load iso_${baseisocaller}_job]} {
				error "isocaller $baseisocaller not supported"
			}
			set options {}
			if {$preset ne ""} {lappend options -preset $preset}
			set root [file_rootname $bam]
			set resultfile [file dir $bam]/isoform_counts-${isocaller}_joint-$root.tsv
			iso_${baseisocaller}_job \
				-reftranscripts $ref \
				{*}$options \
				-cleanup $cleanup \
				-distrreg $distrreg -threads $threads \
				-refseq $refseq \
				$bam $resultfile
		}
		iso_combine_job $projectdir ${isocaller}_joint $iso_match
	}
}

proc cigar2exons {cigar begin} {
	set lcigar [regsub -all {([0-9]+)([A-Z=])} $cigar {\1 \2 }]
	set exonStarts $begin ; set exonEnds {}
	set end $begin
	set exons {}
	set introns {}
	set exonsize 0
	set intronsize 0
	foreach {num op} $lcigar {
		switch $op {
			M - D - = - X {
				if {$intronsize > 0} {
					lappend introns $intronsize
					set intronsize 0
				}
				incr exonsize $num
			}
			N - D {
				if {$exonsize > 0} {
					lappend exons $exonsize
					set exonsize 0
				}
				incr intronsize $num
			}
			H - S - I - P {}
		}
	}
	if {$exonsize > 0} {
		lappend exons $exonsize
		set exonsize 0
	}
	set current $begin
	set regions [list $current]
	foreach exon $exons intron $introns {
		incr current $exon
		lappend regions $current
		if {$intron eq ""} break
		incr current $intron
		lappend regions $current
	}
	return $regions
}

proc iso_write_isoform_counts {targetisoformcountsfile regreftranscripts tcountaVar newheaderVar
	{fields {iqall iq weighedb i unique u strict s aweighed a aunique au astrict as}} {sizeaVar {}}
} {
	upvar $tcountaVar tcounta
	upvar $newheaderVar newheader
	if {$sizeaVar ne ""} {
		upvar $sizeaVar sizea
	}
	set fieldsheader [list_unmerge $fields 1 fieldsids]
	set basefields {chromosome begin end strand exonStarts exonEnds transcript gene geneid}
	# lappend remove {*}$basefields
	set newheader $basefields
	if {$regreftranscripts ne ""} {
		catch {close $f}
		set f [gzopen $regreftranscripts]
		set header [tsv_open $f comments]
		set poss [list_sub [tsv_basicfields $header 14 0] {0 1 2 6 7 8 11 12 13}]
		foreach {isopos genepos geneidpos} [lrange $poss 6 end] break
		set remove [list_sub $header $poss]
		set extrafields [list_lremove $header $remove]
		lappend poss {*}[list_cor $header $extrafields]
		set temp [list_common $extrafields $basefields]
		foreach field $temp {
			set pos [lsearch $extrafields $field]
			lset extrafields $pos ${field}_ori
		}
		lappend newheader {*}$extrafields
	}
	#
	if {$comments eq ""} {
		set comments [deindent {
			#filetype	tsv/transcriptsfile
			#fileversion	0.99
			#fields	table
			#fields	field	number	type	description
			#fields	name	1	String	Name of transcript (usually transcript_id from GTF)
			#fields	chromosome	1	String	Chromosome name
			#fields	strand	1	String	+ or - for strand
			#fields	begin	1	Integer	Transcription start position
			#fields	end	1	Integer	Transcription end position
			#fields	exonCount	1	Integer	Number of exons
			#fields	exonStarts	E	Integer	Exon start positions
			#fields	transcript	1	String	transcript id
			#fields	gene	1	String	gene name
			#fields	geneid	1	String	gene id
			#fields	cdsStart	1	Integer	Coding region start
			#fields	cdsEnd	1	Integer	Coding region end
			#fields	type	1	String	type of element
			#fields	exonEnds	E	Integer	Exon end positions
			#fields	exonFrames	E	Integer	Exon frame offsets {0,1,2}
			#fields	score	1	Float	Score
			#fields	name2	1	String	Alternate name (e.g. gene_id from GTF)
			#fields	length	1	Integer	isoform length
			#fields	exons	1	Integer	Number of exons
			#fields	category	1	String	one of the isoform categories (known, novel_in_catalog, novel_not_in_catalog, intergenic)
			#fields	associated_gene	1	String	the reference gene name
			#fields	associated_transcript	1	String	the reference transcript name
			#fields	counts	1	Integer	Number of reads mapping to isoform
		}]
	}
	catch {close $o}
	set o [wgzopen $targetisoformcountsfile]
	puts $o $comments
	puts $o [join $newheader \t]\tcategory\tsize\t[join $fieldsheader \t]
	if {$regreftranscripts ne ""} {
		set keepoutputa [array get outputa]
		while 1 {
			if {[gets $f line] == -1} break
			set line [split $line \t]
			set iso [lindex $line $isopos]
			if {$iso eq ""} continue
			if {[info exists outputa($iso)]} {
				unset outputa($iso)
			}
			# set geneid [lindex $line $geneidpos]
			set line [list_sub $line $poss]
			if {![info exists sizea($iso)]} {
				foreach {starts ends} [list_sub $line {4 5}] break
				set size 0
				foreach s [split [string trim $starts ,] ,] e [split [string trim $ends ,] ,] {
					set size [expr {$size + $e - $s}]
				}
				set sizea($iso) $size
			}
			set resultline [join $line \t]\tknown\t$sizea($iso)
			set output 0
			foreach id $fieldsids {
				if {![info exists tcounta($iso,$id)]} {
					append resultline \t0
				} elseif {![string is int $tcounta($iso,$id)]} {
					if {$tcounta($iso,$id) > 0} {set output 1}
					append resultline \t$tcounta($iso,$id)
				} else {
					if {$tcounta($iso,$id) > 0} {set output 1}
					append resultline \t$tcounta($iso,$id)
				}
			}
			if {!$output} continue
			puts $o $resultline
		}
		catch {close $f}
	}
	return $o
}

proc iso_polya {seq qstart qend} {
	if {[regexp -all T [string range $seq $qstart-10 $qstart]] > 7} {
		return -1
	} elseif {[regexp -all A [string range $seq $qend $qend+10]] > 7} {
		return 1
	} else {
		return 0
	}
}

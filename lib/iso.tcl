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
	set isoformfiles [bsort [jobglob samples/*/isoform_counts-${isocaller}-*.tsv]]	
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

proc iso_joint_job {args} {
	upvar job_logdir job_logdir
	set iso_joint {}
	set iso_match {}
	set dbdir {}
	set threads 2
	set distrreg 0
	set cleanup 1
	cg_options iso_joint args {
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
			$source
		} -targets {
			$ref
		} -vars {
			source ref
		} -code {
			set fields [cg select -h $source]
			set fields [list_sub $fields -exclude [list_find -regexp $fields ^counts]]
			cg select -overwrite 1 -f $fields $source $ref.temp
			file rename $ref.temp $ref
		}
		foreach bam [jobglob $projectdir/samples/*/*.bam $projectdir/samples/*/*.cram] {
			set preset {}
			foreach {baseisocaller preset} [split $isocaller _] break
			if {![auto_load iso_${baseisocaller}_job]} {
				error "iscaller $baseisocaller not supported"
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

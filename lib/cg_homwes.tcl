#!/bin/sh
# the next line restarts using wish \
exec cg source "$0" ${1+"$@"}

proc findfield {fields pattern args} {
	set args [list $pattern {*}$args]
	foreach pattern $args {
		set pos [lsearch -glob $fields $pattern]
		if {$pos != -1} break
	}
	lindex $fields $pos
}

proc cg_homwes {args} {
	set callers {gatk-rdsbwa- sam-rdsbwa-}
	set filterrepeats 1
	set genoqual 40
	set allowedheterozygous 1
	set homozygdensity 200
	set homozyggap 4000
	set homozygwindowsnp 20
	set snpsonly 0
	set dbdir {}
	set vcf 0
	set pos 0
	set variantsonly 0
	set samples {}
	set resultfile {}
	cg_options homwes args {
		-dbdir {
			set dbdir $value
		}
		-callers {
			set callers $value
		}
		-allowedheterozygous -- -htz {
			set allowedheterozygous $value
		}
		-filterrepeats {
			set filterrepeats $value
		}
		-genoqual {
			set genoqual $value
		}
		-density {
			set homozygdensity $value
		}
		-gap {
			set homozyggap $value
		}
		-window {
			set homozygwindowsnp $value
		}
		-variantsonly {
			set variantsonly $value
		}
		-snpsonly {
			set snpsonly $value
		}
	} {annotcomparfile samples resultfile} 1 3
	set annotcomparfile [file_absolute $annotcomparfile]
	if {$resultfile eq ""} {
		set resultfile [file root $annotcomparfile]-homwes.tsv
	}
	set resultfile [file_absolute $resultfile]
	set resultfilebase [file root $resultfile]
	set workdir $resultfilebase.work
	file mkdir $workdir
	set workbase $workdir/[file tail $resultfilebase]
	set dir [file dir $annotcomparfile]
	set tail [file tail $annotcomparfile]
	if {[file ext [gzroot $annotcomparfile]] eq ".vcf"} {
		if {$dbdir eq ""} {
			error "No -dbdir option given, this is needed for analysing"
		}
		set dbdir [dbdir $dbdir]
		set vcf 1
		set tsvfile $workdir/[file root [file tail $annotcomparfile]].tsv
		putslog "Converting vcf file $annotcomparfile to tsv $tsvfile"
		cg vcf2tsv -split 0 $annotcomparfile $tsvfile.temp 2>@ stderr
		putslog "annotating $tsvfile"
		cg annotate $tsvfile.temp $tsvfile.temp2 {*}[gzfiles $dbdir/reg_*_microsat.tsv $dbdir/reg_*_simpleRepeat.tsv]
		file rename -force $tsvfile.temp2 $tsvfile
		file delete $tsvfile.temp
		set usefile $tsvfile
	} elseif {[llength [list_common [cg select -h $annotcomparfile] {microsat simpleRepeat}]] != 2} {
		if {$dbdir eq ""} {
			error "No -dbdir option given, this is needed for annotating the file"
		}
		set tsvfile $workdir/[file root [file tail $annotcomparfile]].tsv
		putslog "annotating $tsvfile"
		cg annotate $annotcomparfile $tsvfile.temp2 {*}[gzfiles $dbdir/reg_*_microsat.tsv $dbdir/reg_*_simpleRepeat.tsv]
		file rename -force $tsvfile.temp2 $tsvfile
		file delete $tsvfile.temp
		set usefile $tsvfile
	} else {
		set usefile $annotcomparfile
	}
	cd $dir
	putslog "Gathering samples"
	set samples [list_remdup $samples]
 	set header [cg select -h $usefile]
	set hsamples {}
	foreach hsample [samples $header] {
		lappend hsamples [lindex [split $hsample -] end]
	}
	set hsamples [list_remdup $hsamples]
	if {![llength $hsamples]} {
		# putslog "Only one sample in file: not using different variant callers"
		set hassamples 0
		# set callers {}
	} else {
		set hassamples 1
	}
	if {![llength $samples]} {
		if {[llength $hsamples]} {
			set samples $hsamples
		} else {
			set samples sample1
		}
	} elseif {[llength $hsamples]} {
		set notasample [list_lremove $samples $hsamples]
		if {[llength $notasample]} {
			error "error: not a sample: $notasample"
		}
	}
	set resultfiles {}
	foreach sample $samples {
		putslog "Sample $sample"
		set sworkbase $workbase-$sample
		set query {}
		set fields {chromosome begin end type ref alt}
		set filter {}
		set caller1 [lindex $callers 0]
		if {"alleleSeq1-$caller1$sample" in $header} {
			set postfix -$caller1$sample
		} elseif {"alleleSeq1-$sample" in $header} {
			set postfix -$sample
		} elseif {"alleleSeq1" in $header} {
			set postfix ""
		} else {
			error "Could not find alleleSeq1 field for sample $sample in header (checked alleleSeq1-$caller1$sample, alleleSeq1-$sample, alleleSeq1)"
		}
		set field [findfield $header gfilter$postfix filter$postfix]
		if {$field ne ""} {
			lappend query "\$$field in {. PASS Pass pass}"
		} elseif {$postfix ne ""} {
			set field [findfield $header filter]
			if {$field ne ""} {
				lappend query "\$$field in {. PASS Pass pass}"
			} else {
				puts "warning: field \"filter$postfix\" is missing"
			}
		} else {
			puts "warning: field \"(g)filter$postfix\" is missing"
		}
		set field [findfield $header genoqual$postfix]
		if {$field ne ""} {
			lappend query "def(\$$field,100) > $genoqual"
		} else {
			puts "warning: field \"genoqual$postfix\" is missing"
		}
		if {$filterrepeats} {
			lappend query "\$microsat==\"\" && \$simpleRepeat==\"\""
		}
		if {$snpsonly} {
			lappend query {$type eq "snp"}
		}
		set dosame {}
		set callers_used {}
		foreach caller $callers {
			if {[findfield $header alleleSeq1-$caller$sample] ne "" && [findfield $header alleleSeq2-$caller$sample] ne ""} {
				lappend dosame $caller$sample
				lappend callers_used $caller
			} else {
				puts "warning: field \"alleleSeq1-$caller$sample\" or \"alleleSeq2-$caller$sample\" is missing, you can use -callers '' if the variant file does not contain data for different callers"
			}
		}
		if {[llength $dosame] > 1} {
			lappend query "same(\"[join $dosame \",\"]\")"
		}
		lappend fields "alleleSeq*-$sample=\$alleleSeq*$postfix"
		if {$variantsonly} {
			if {"sequenced$postfix" in $header} {
				lappend query "\$sequenced$postfix eq \"v\""
			} elseif {"zyg$postfix" in $header} {
				lappend query "\$zyg$postfix ni \"u r\""
			} else {
				lappend query "\$alleleSeq1$postfix != toupper(\$ref) || \$alleleSeq2$postfix != toupper(\$ref)"
			}
		} else {
			if {"sequenced$postfix" in $header} {
				lappend query "\$sequenced$postfix ne \"u\""
			} elseif {"zyg$postfix" in $header} {
				lappend query "\$zyg$postfix ni \"u\""
			}
		}
		if {$query eq ""} {set query 1}
		putslog "Quality filtering data"
		cg select -overwrite 1 \
			-q "(chr_clip(\$chromosome) ni {X Y M MT}) && [join $query " && "]" \
			-f $fields \
			$usefile ${sworkbase}-filtered.tsv
		
		# export to plink
		putslog "Export to plink"
		cg exportplink --samples $sample ${sworkbase}-filtered.tsv $sworkbase
		# 
		set c [split [string trim [file_read $sworkbase.tfam.pre]] \n]
		set o [open $sworkbase.tfam w]
		foreach line $c {
			lset line 4 1
			lset line 5 2
			puts $o [join $line \t]
		}
		close $o
		putslog "Run homozygosity mapping"
		exec plink --tfile $sworkbase --out $sworkbase --make-bed
		exec plink --bfile $sworkbase --out $sworkbase --homozyg-window-snp $homozygwindowsnp  --homozyg-window-het $allowedheterozygous --homozyg-window-threshold 0.05 --homozyg-window-missing 10 --homozyg-snp 10  --homozyg-density $homozygdensity  --homozyg-gap $homozyggap --homozyg-group
		# exec plink --tfile ${sworkbase} --missing
		
		putslog "Create result $resultfile"
		set f [open $sworkbase.hom]
		set o [open $sworkbase-hom.tsv w]
		if {!$vcf} {
			puts $o "# analysed file: $usefile"
		} else {
			puts $o "# analysed file: $annotcomparfile (converted to $usefile)"
		}
		puts $o "# samples: $samples"
		puts $o "# callers_used: $callers_used"
		puts $o "# callers_param: $callers"
		puts $o "# allowedheterozygous: $allowedheterozygous"
		puts $o "# variantsonly: $variantsonly"
		puts $o "# snpsonly: $snpsonly"
		puts $o "# dbdir: $dbdir"
		set line [gets $f]
		if {[list {*}$line] ne "FID IID PHE CHR SNP1 SNP2 POS1 POS2 KB NSNP DENSITY PHOM PHET"} {
			error "wrong header for $sworkbase.hom (plink hom output changed)"
		}
		puts $o [join {chromosome begin end DENSITY PHOM PHET FID IID PHE SNP1 SNP2 KB NSNP} \t]
		set poss [list_cor {FID IID PHE CHR SNP1 SNP2 POS1 POS2 KB NSNP DENSITY PHOM PHET} {CHR POS1 POS2 DENSITY PHOM PHET FID IID PHE SNP1 SNP2 KB NSNP}]
		while {[gets $f line] >= 0} {
			foreach {FID IID PHE CHR SNP1 SNP2 POS1 POS2 KB NSNP DENSITY PHOM PHET} $line break
			puts $o [join [list_sub $line $poss] \t]
		}
		close $o
		close $f
		file rename -force $sworkbase-hom.tsv $resultfilebase-$sample.tsv
		lappend resultfiles $resultfilebase-$sample.tsv
	}
	if {[llength $samples] == 1} {
		file rename -force $resultfilebase-$sample.tsv $resultfile
	} else {
		putslog "Making $resultfile"
		cg multireg $resultfile.temp {*}$resultfiles
		file rename -force $resultfile.temp $resultfile
	}
}

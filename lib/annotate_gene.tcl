#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc annot_init {} {
	global impact_list impact_score impact_loc annotate_type2impact
	set impact_list {
		{downstream annot downstream} {upstream annot upstream}
		{annot annot annotation_region} {A annot annotated}
		{intron intron intron}
		{reg reg regulatory} {prom reg promotor} {splice splice splice}
		{RNA transcript RNA} {RNASPLICE transcript RNA_splice}
		{UTR3 UTR3 3'UTR} {UTR3SPLICE UTR3 3'UTR_splice} {RNAEND transcript transcription_end}
		{ESPLICE splice essentialsplice}
		{CDSsilent CODING silent} {UTR5 UTR5 5'UTR} {UTR5SPLICE UTR5 5'UTR_splice} {UTR5KOZAK UTR5 5'UTR_kozak_seq} {RNASTART transcript transcription_start}
		{CDSMIS CODING missense} {CDSCOMP CODING complex} {CDSDEL CODING deletion} {CDSINS CODING insertion} {CDSSPLICE CODING splice}
		{CDSNONSENSE CODING nonsense}
		{CDSSPLICE CODING splice_deletion}
		{CDSSTOP CODING readthrough} {CDSFRAME CODING frameshift} {CDSSTART CODING startcodon} {CDSSTARTCOMP CODING startcodon} {CDSSTARTDEL CODING startcodon}
		{GENECOMP GENE genecomplex}
		{GENEDEL GENE genedeletion}
	}
	unset -nocomplain impact_score
	set num 1
	list_foreach {code location name} $impact_list {
		set impact_score($code) $num
		set impact_loc($code) $location
		incr num
	}
	array set annotate_type2impact {
		UTR UTR RNA RNA CDS CDS
		intron intron
		esplice ESPLICE
		splice splice
		promotor prom
		regulatory reg
		upstream upstream
		downstream downstream
		annotation_region annot
	}
}
annot_init

proc var_impact_list {{filter {}}} {
	if {$filter eq ""} {
		return [list_subindex $::impact_list 0]
	} else {
		set poss [list_find -regexp $filter {[*><]}]
		if {![llength $poss]} {
			return $filter
		}
		set result {}
		foreach pattern $filter {
			set temp [list_subindex $::impact_list 0]
			set first [string index $pattern 0]
			if {$first eq ">"} {
				if {[string index $pattern 1] eq "="} {
					set score [annotate_impact2score [string range $pattern 2 end]]
					incr score -1
				} else {
					set score [annotate_impact2score [string range $pattern 1 end]]
				}
				lappend result {*}[lrange $temp $score end]
			} elseif {$first eq "<"} {
				if {[string index $pattern 1] eq "="} {
					set score [annotate_impact2score [string range $pattern 2 end]]
					incr score -2
				} else {
					set score [annotate_impact2score [string range $pattern 1 end]]
					incr score -1
				}
				lappend result {*}[lrange $temp 0 $score]
			} elseif {[regexp {\*} $pattern]} {
				lappend result {*}[list_sub $temp [list_find -glob $temp $pattern]]
			} else {
				lappend result $pattern
			}
		}
		return [list_remdup $result]
	}
}

proc annotate_compl {geneobj seq} {
	if {[isint $seq]} {return $seq}
	if {![dict get $geneobj complement]} {
		return $seq
	} else {
		return [seq_complement $seq]
	}
}

proc annotate_impact2score {impact} {
	get ::impact_score($impact) 3
	
}

proc annotate_type2impact {type} {
	get ::annotate_type2impact($type) $type
}

proc annotate_getregion {refseq start end} {
	if {$start < $end} {
		return [string range $refseq $start $end]
	} else {
		return [seq_complement [string range $refseq $end $start]]
	}
}

proc switchval {var1 var2} {
	upvar $var1 v1
	upvar $var2 v2
	set temp $v1
	set v1 $v2
	set v2 $temp
}

proc annotate_complementpos {geneobj pos} {
	return [expr {[dict get $geneobj refseqlen]-$pos-1}]
}

proc calc_rpos {geneobj rsnppos} {
	set rpos [expr {$rsnppos - [dict get $geneobj rrefstart]}]
	if {$rpos >= 0} {
		return [expr {$rpos+1}]
	} else {
		return $rpos
	}
}

# set upstreamsize 2000 ; set nocds 0
proc annotatevar_gene_makegeneobj {genomef dbline dposs {upstreamsize 2000} {nocds 0}} {
	# splice 8 bases, essentialsplice 2 bases intron, upstream,downstream -> 2k bases
	foreach {dchrom dstart dend strand exonStarts exonEnds cdsStart cdsEnd transcriptname genename exonCount} [list_sub $dbline $dposs] break
	if {$nocds} {
		set cdsStart {}
		set cdsEnd {}
	}
	if {$strand eq "-"} {set complement 1} else {set complement 0}
	unset -nocomplain adata
	set adata(genomef) $genomef
	set adata(complement) $complement
	set adata(start) $dstart
	set adata(end) $dend
	set adata(transcriptname) $transcriptname
	set adata(genename) $genename
	set adata(chrom) $dchrom
	set adata(upstreamsize) $upstreamsize
	if {$cdsStart == $cdsEnd} {
		set cdsStart -1
		set cdsEnd -1
		set cds 0
	} else {
		set cds 1
	}
	set adata(cds) $cds
	set exonStarts [split [string trimright $exonStarts ,] ,]
	set exonEnds [split [string trimright $exonEnds ,] ,]
	set adata(nrexons) [llength $exonStarts]
	# use code coming from novoSNP:
	# transform to ft addressing used in novoSNP previously (ft endpos is included, so -1)

	# ftlist wil contain lines for each element of the form
	# gbegin gend type element rnabegin rnaend cdsbegin cdsend grnabegin grnaend
	# where element: exon$nr or intron$nr or pre (for variants before RNA)
	set ftlist {}
	if {!$complement} {set type upstream} else  {set type downstream}
	lappend ftlist [list -2 [expr {$dstart-1}] $type]
	set prev -1
	if {$cds} {set type UTR} else {set type RNA}
	set intronnr 0; set exonnr 0
	if {!$complement} {set nr 1; set nrstep 1} else {set nr [llength $exonStarts]; set nrstep -1}
	foreach s $exonStarts e $exonEnds {
# puts --------\n[join $ftlist \n]\n
# if {$s == 52847059} {error STOPPED}
# putsvars s e type prev
		set adata(exon$nr) [list $s [expr {$e-1}]]
		incr nr $nrstep
		if {$prev != -1} {
			lappend ftlist [list $prev [expr {$s-1}] intron]
		}
		if {($cdsStart >= $s) && ($cdsStart <= $e)} {
			set el [list $s [expr {$cdsStart-1}] $type]
			lappend ftlist $el
			set s $cdsStart
			set type CDS
		}
		if {($cdsEnd >= $s) && ($cdsEnd <= $e)} {
			set el [list $s [expr {$cdsEnd-1}] CDS]
			lappend ftlist $el
			set s $cdsEnd
			set type UTR
		}
		set prev $e
		if {$s != $e} {
			set el [list $s [expr {$e-1}] $type]
			lappend ftlist $el
		}
	}

	if {$complement} {set type upstream} else  {set type downstream}
	lappend ftlist [list [expr {$e}] 2147483647 $type]
	# add info to ftlist so data will be: begin end type exon rnabegin rnaend cdsbegin cdsend
	set adata(grefstart) 0
	set adata(rrefstart) 0
	set exon 0
	set rnabegin -1
	set rnaend -1
	set cdsbegin -1
	set cdsend -1
	set prevcheck -1
	set tempftlist {}
	set num 0
	set nums {}
	if {!$complement} {set templist $ftlist} else {set templist [list_reverse $ftlist]}
	list_foreach {s e t} $templist {
		if {[inlist {UTR RNA CDS} $t]} {
			if {!$complement} {
				if {$s != $prevcheck} {
					incr exon
				}
			} else {
				if {$e != $prevcheck} {
					incr exon
				}
			}
			if {$rnabegin == -1} {
				set rnabegin 0
			} else {
				set rnabegin [expr {$rnaend+1}]
			}
			set rnaend [expr {$rnabegin+$e-$s}]
			if {$t eq "CDS"} {
				if {$cdsbegin == -1} {
					set cdsbegin 0
					if {!$complement} {set adata(gpstart) $s} else {set adata(gpstart) $e}
					set adata(rpstart) $rnabegin
					set adata(grefstart) $adata(gpstart)
					set adata(rrefstart) $adata(rpstart)
				} else {
					set cdsbegin [expr {$cdsend+1}]
				}
				set cdsend [expr {$cdsbegin+$e-$s}]
				set adata(gpend) $e
				set adata(rpend) $rnaend
			} else {
				set cdsbegin -1
				set cdsend -1
			}
			lappend tempftlist [list $s $e $t exon$exon $rnabegin $rnaend $cdsbegin $cdsend]
			if {!$complement} {
				set prevcheck [expr {$e+1}]
			} else {
				set prevcheck [expr {$s-1}]
			}
			lappend nums $num
		} elseif {$t eq "upstream"} {
			lappend tempftlist [list $s $e $t upstream [expr {$s-$e-1}] -1 $cdsend $cdsend]
		} elseif {$t eq "downstream"} {
			set rnabegin [expr {$rnaend+1}]
			set rnaend [expr {$rnabegin+$e-$s}]
			lappend tempftlist [list $s $e $t downstream $rnabegin $rnaend $cdsend $cdsend]
		} else {
			if {$exon == 0} {set temp $t} else {set temp intron$exon}
			lappend tempftlist [list $s $e $t $temp $rnaend $rnaend $cdsend $cdsend]
		}
		incr num
	}
	set num [llength $tempftlist]
	lset tempftlist end 3 downstream
	# join $tempftlist \n
	set rnalist [list_sub $tempftlist $nums]
	if {!$complement} {
		set adata(ftlist) $tempftlist
		set adata(rnalist) $rnalist
	} else {
		set adata(ftlist) [list_reverse $tempftlist]
		set adata(rnalist) [list_reverse $rnalist]
	}
	if {[info exists adata]} {
		set adata(rrend) [lindex $adata(rnalist) end 5]
		set adata(grstart) [lindex $adata(rnalist) 0 0]
		set adata(grend) [lindex $adata(rnalist) end 1]
	} else {
		set adata(grstart) -1
		set adata(grend) -1
	}
	# join $rnalist \n
	# join $ftlist \n
	return [array get adata]
}

proc annotatevar_gene_rnaseq {geneobjVar} {
	upvar $geneobjVar geneobj
	if {[dict exists $geneobj rnaseq]} {return [dict get $geneobj rnaseq]}
	set genomef [dict get $geneobj genomef]
	set chr [dict get $geneobj chrom]
	set complement [dict get $geneobj complement]
	set result {}
	# join [dict get $geneobj rnalist] \n
	list_foreach {gs ge t el rs re ps pe} [dict get $geneobj rnalist] {
		incr ge ;# genome_get works with half open coordinates
		append result [genome_get $genomef $chr $gs $ge 0 1]
	}
	if {$complement} {
		set result [seq_complement $result]
	}
	dict set geneobj rnaseq $result
	return $result
}

proc annotatevar_gene_pseq {geneobjVar} {
	upvar $geneobjVar geneobj
	if {[dict exists $geneobj pseq]} {return [dict get $geneobj pseq]}
	set seq [string range [annotatevar_gene_rnaseq geneobj] [dict get $geneobj rpstart] [dict get $geneobj rpend]]
	set pseq [seq_translate $seq]
	dict set geneobj pseq $pseq
	return $pseq
}

proc annotatevar_gene_getrnaseq {geneobjVar qstart qend {line {}}} {
	upvar $geneobjVar geneobj
	set genomef [dict get $geneobj genomef]
	set chr [dict get $geneobj chrom]
	set complement [dict get $geneobj complement]
	if {!$complement} {
		if {[llength $line]} {
			foreach {gs ge t el rs re ps pe} $line break
			if {$qstart >= $rs && $qend <= $re} {
				set gstart [expr {$gs+$qstart-$rs}]
				set gend [expr {$gs+$qend-$rs}]
				incr gend ;# genome_get works with half open coordinates
				set result [genome_get $genomef $chr $gstart $gend 0 1]
				if {!$complement} {return $result} else {return [seq_complement $result]}
			}
		}
		set result {}
		# join [dict get $geneobj rnalist] \n
		list_foreach {gs ge t el rs re ps pe} [dict get $geneobj rnalist] {
			if {$re >= $qstart} {
				set gstart [expr {$gs+$qstart-$rs}]
				set gend [expr {$gs+$qend-$rs}]
				if {$gstart < $gs} {set gstart $gs}
				if {$gend > $ge} {set gend $ge}
				incr gend ;# genome_get works with half open coordinates
				append result [genome_get $genomef $chr $gstart $gend 0 1]
				if {$re > $qend} break
			}
		}
	} else {
		if {[llength $line]} {
			foreach {gs ge t el rs re ps pe} $line break
			if {$qstart >= $rs && $qend <= $re} {
				# for complement, take re as reference, but start from gs
				# we want to extract genomic non complement first
				set gstart [expr {$gs+$re-$qend}]
				set gend [expr {$gs+$re-$qstart}]
				incr gend ;# genome_get works with half open coordinates
				set result [genome_get $genomef $chr $gstart $gend 0 1]
				if {!$complement} {return $result} else {return [seq_complement $result]}
			}
		}
		set result {}
		# join [dict get $geneobj rnalist] \n
		list_foreach {gs ge t el rs re ps pe} [dict get $geneobj rnalist] {
			if {$rs <= $qend} {
				# for complement, take re as reference, but start from gs
				# we want to extract genomic non complement first
				set gstart [expr {$gs+$re-$qend}]
				set gend [expr {$gs+$re-$qstart}]
				if {$gstart < $gs} {set gstart $gs}
				if {$gend > $ge} {set gend $ge}
				incr gend ;# genome_get works with half open coordinates
				append result [genome_get $genomef $chr $gstart $gend 0 1]
				if {$rs < $qstart} break
			}
		}
		set result [seq_complement $result]
	}
	return $result
}

proc annotatevar_gene_getprnaseq {geneobjVar qstart qend {line {}}} {
	upvar $geneobjVar geneobj
	incr qstart [dict get $geneobj rpstart]	
	incr qend [dict get $geneobj rpstart]	
	if {[dict exists $geneobj rnaseq]} {
		return [string range [annotatevar_gene_rnaseq geneobj] $qstart $qend]
	}
	annotatevar_gene_getrnaseq geneobj $qstart $qend $line
}

# input
#   line: line for element from ftlist, contains: gbegin gend type exon rnabegin rnaend cdsbegin cdsend
# return values
#   rpos: position relative to transcript
proc annotategene_rpos {geneobj line snppos} {
	set complement [dict get $geneobj complement]
	foreach {gbegin gend eltype element rnabegin rnaend cdsbegin cdsend} $line break
	if {!$complement} {
		# position in intron/exon
		set eipos [expr {$snppos-$gbegin}]
		# position in rna
		if {$rnabegin == $rnaend} {
			# in intron
			set rpos $rnabegin
		} else {
			set rpos [expr {$rnabegin+$eipos}]
		}
	} else {
		# position in intron/exon
		set eipos [expr {$gend-$snppos}]
		# position in rna
		if {$rnabegin == $rnaend} {
			# in intron
			set rpos $rnabegin
		} else {
			set rpos [expr {$rnabegin+$eipos}]
		}
	}
	return $rpos
}

proc annotategene_element {geneobj line snppos} {
	if {$line eq "pre"} {return [list pre ]}
	set complement [dict get $geneobj complement]
	foreach {gbegin gend eltype element rnabegin rnaend cdsbegin cdsend} $line break
	if {[dict exists $geneobj $element]} {
		# for exons with UTR/CDS bounderies, we need to calculate to intron bounderies iso within small splice element
		# 
		foreach {is ie} [dict get $geneobj $element] break
		if {!$complement} {
			set elementpos [expr {$snppos-$is+1}]
		} else {
			set elementpos [expr {$ie-$snppos+1}]
		}
	} else {
		# elementpos gives position in intron/exon
		if {$eltype eq "upstream"} {
			if {!$complement} {
				set elementpos [expr {$snppos - $gend-1}]
			} else {
				set elementpos [expr {$gbegin - $snppos - 1}]
			}
		} elseif {!$complement} {
			set elementpos [expr {$snppos-$gbegin+1}]
		} else {
			set elementpos [expr {$gend-$snppos+1}]
		}
	}
	if {$elementpos > 0} {set elementpos +$elementpos}
	if {$element eq "upstream"} {set element up} elseif {$element eq "downstream"} {set element down}
	list $element $elementpos
}

proc annotategene_p_complex {snptype} {
	if {$snptype eq "del"} {
		return DEL
	} else {
		return COMP
	}
}

proc annotategene_p_loc {AZ1 prstart {AZ2 {}} {prend {}}} {
	if {$prend eq "" || $prstart == $prend} {
		incr prstart
		return p.${AZ1}${prstart}
	} else {
		incr prstart ; incr prend
		return p.${AZ1}${prstart}_${AZ2}${prend}
	}
}

proc annotategene_findreg {geneobj snppos snpend} {
	set sline [annotategene_findregc [dict get $geneobj ftlist] $snppos pre]
	set eline [annotategene_findregc [dict get $geneobj ftlist] $snpend post]
	list $sline $eline
}

proc annotategene_one_ins {geneobj loc} {
	global impact_score
	foreach {chrom snppos snpend snptype ref alt} $loc break
	set complement [dict get $geneobj complement]
	set chr [dict get $geneobj chrom]
	set dbstart [dict get $geneobj start]
	set dbend [expr {[dict get $geneobj end]-1}]
	set snp_descr {}
	# use positions before and after insert (like in del), snpend is inclusive here (not half open)
	set snpend $snppos
	set snppos [expr {$snppos - 1}]
	# reference
	if {$complement} {set strand -} else {set strand +}
	set result $strand[dict get $geneobj transcriptname]
	# get region(s)
	foreach {sline eline} [annotategene_findreg $geneobj $snppos $snpend] break
	# shift to 3'?
	if {![isint $alt] && [lindex $sline 2] eq [lindex $eline 2]} {
		set genomef [dict get $geneobj genomef]
		if {!$complement} {
			set gend [lindex $eline 1]
			while {$snpend <= $gend} {
				set from [string index $alt 0]
				if {[catch {
					set dest [genome_get $genomef $chr $snpend [expr {$snpend+1}] 0 1]
				} error]} {
					if {[string match {trying to get sequence beyond end of chromosome*} $error]} break
					error $error
				}
				if {$from ne $dest} break
				incr snppos; incr snpend
				set alt [string range $alt 1 end][string index $alt 0]
			}
			if {$snpend > $gend} {
				foreach {sline eline} [annotategene_findreg $geneobj $snppos $snpend] break
			}
		} else {
			# for complement towards 3' is going back
			set gbegin [lindex $sline 0]
			while {$snppos >= $gbegin} {
				set from [string index $alt end]
				# get base just before (current) insert
				if {[catch {
					set dest [genome_get $genomef $chr $snppos $snpend 0 1]
				} error]} {
					if {[string match {trying to get sequence beyond end of chromosome*} $error]} break
					error $error
				}
				# puts [genome_get $genomef $chr $gbegin $snpend 0 1]\ $alt
				if {$from ne $dest} break
				incr snppos -1; incr snpend -1
				set alt $from[string range $alt 0 end-1]
			}
			if {$snppos < $gbegin} {
				foreach {sline eline} [annotategene_findreg $geneobj $snppos $snpend] break
			}
		}
	}
	#
	# element
	# find start and end location line in ftlist
	if {!$complement} {
		set snp1 $snppos
		set snp2 $snpend
		set line1 $sline ; set line2 $eline
		set strand -
	} else {
		set snp1 $snpend
		set snp2 $snppos
		set line2 $sline ; set line1 $eline
		set strand +
	}
	foreach {element1 element1pos} [annotategene_element $geneobj $line1 $snp1] break
	foreach {element2 element2pos} [annotategene_element $geneobj $line2 $snp2] break
	# putsvars line1 snp1 line2 snp2 element1 element1pos element2 element2pos
	if {$element1 eq $element2} {
		if {$element1pos eq $element2pos} {
			append result	:$element1${element1pos}
		} else {
			set temp [string index $element2pos 0]
			if {$temp eq "+" && $temp eq [string index $element1pos 0]} {set element2pos [string range $element2pos 1 end]}
			append result	:$element1${element1pos}_${element2pos}
		}
	} else {
		append result	:$element1${element1pos}_$element2${element2pos}
	}
	#
	# get snp_descr
	#
	# coding DNA level
	foreach {reftype1 ref1 offset1 change} [annotategene_one_c $geneobj $chr $snp1 ins $ref $alt $line1 impact] break
	foreach {reftype2 ref2 offset2 change} [annotategene_one_c $geneobj $chr $snp2 ins $ref $alt $line2 impact2] break
	# putsvars line1 line2 reftype1 ref1 offset1 impact reftype2 ref2 offset2 impact2 change
	if {$ref1 eq $ref2 && $offset1 eq $offset2} {
		set snp_descr $reftype1$ref1${offset1}$change
	} else {
		set snp_descr $reftype1$ref1${offset1}_$ref2${offset2}$change
	}
	append result :$snp_descr
	# impact
	if {($impact eq "CDS" && [regexp UTR|RNA $impact2])
		|| ($impact eq "ESPLICE" && ([regexp UTR|CDS|RNA $impact2]))
		|| ($impact eq "ESPLICE" && $impact2 eq "splice")
		|| ($impact eq "splice" && $impact2 eq "intron")
		|| ([regexp UTR|RNA $impact] && $impact2 eq "downstream")
	} {
		set impact $impact2
		set line $line2
	} else {
		set line $line1
	}
	if {$impact eq "CDS"} {
		append result :[annotategene_one_p_ins $geneobj $snpend $ref $alt $line impact]
	}
	return [list $impact $result]
}

proc annotategene_one_del {geneobj loc} {
	global impact_score
	foreach {chrom snppos snpend snptype ref alt} $loc break
	set ref [expr {$snpend-$snppos}]
	set len $ref
	set complement [dict get $geneobj complement]
	set chr [dict get $geneobj chrom]
	set dbstart [dict get $geneobj start]
	set dbend [expr {[dict get $geneobj end]-1}]
	set snp_descr {}
	if {[isint $alt]} {
		set rlen $alt
	} else {
		set rlen [string length $alt]
	}
	#
	# reference
	set snpend [expr {$snppos+$len-1}]
	if {$complement} {set strand -} else {set strand +}
	set result $strand[dict get $geneobj transcriptname]
	# find region
	foreach {sline eline} [annotategene_findreg $geneobj $snppos $snpend] break
	# shift to 3'
	if {$snptype eq "del" && [lindex $sline 3] eq [lindex $eline 3]} {
		# shift to 3'?
		set genomef [dict get $geneobj genomef]
		if {!$complement} {
			set gend [lindex $eline 1]
			while {$snpend < $gend} {
				if {[catch {
					set from [genome_get $genomef $chr $snppos [expr {$snppos+1}] 0 1]
					set dest [genome_get $genomef $chr [expr {$snpend+1}] [expr {$snpend+2}] 0 1]
				}]} break
				if {$from ne $dest} break
				if {$from eq "N"} break
				incr snppos; incr snpend
			}
		} else {
			# for complement towards 3' is going back
			set gbegin [lindex $sline 0]
			while {$snppos > $gbegin} {
				if {[catch {
					set from [genome_get $genomef $chr [expr {$snppos-1}] $snppos 0 1]
					set dest [genome_get $genomef $chr $snpend [expr {$snpend+1}] 0 1]
				}]} break
				if {$from ne $dest} break
				if {$from eq "N"} break
				incr snppos -1; incr snpend -1
			}
		}
	}
	#
	# element
	# find start and end location line in ftlist
	if {!$complement} {
		set snp1 $snppos
		set snp2 $snpend
		set line1 $sline ; set line2 $eline
		set strand -
	} else {
		set snp1 $snpend
		set snp2 $snppos
		set line2 $sline ; set line1 $eline
		set strand +
	}
	foreach {element1 element1pos} [annotategene_element $geneobj $line1 $snp1] break
	foreach {element2 element2pos} [annotategene_element $geneobj $line2 $snp2] break
	# putsvars line1 snp1 line2 snp2 element1 element1pos element2 element2pos
	if {$element1 eq $element2} {
		if {$element1pos eq $element2pos} {
			append result	:$element1${element1pos}
		} else {
			set temp [string index $element2pos 0]
			if {$temp eq "+" && $temp eq [string index $element1pos 0]} {set element2pos [string range $element2pos 1 end]}
			append result	:$element1${element1pos}_${element2pos}
		}
	} else {
		append result	:$element1${element1pos}_$element2${element2pos}
	}
	#
	# get snp_descr
	#
	# coding DNA level
	foreach {reftype1 ref1 offset1 change} [annotategene_one_c $geneobj $chr $snp1 $snptype $ref $alt $line1 impact] break
	foreach {reftype2 ref2 offset2 change} [annotategene_one_c $geneobj $chr $snp2 $snptype $ref $alt $line2 impact2] break
	# putsvars reftype1 ref1 offset1 impact reftype2 ref2 offset2 impact2 change
	if {[annotate_impact2score $impact2] > [annotate_impact2score $impact]} {
		set impact $impact2
	}
	# putsvars reftype1 ref1 offset1 reftype2 ref2 offset2 change
	if {$ref1 eq $ref2 && $offset1 eq $offset2} {
		set snp_descr $reftype1$ref1${offset1}$change
	} else {
		set snp_descr $reftype1$ref1${offset1}_$ref2${offset2}$change
	}
	#
	set srpos [annotategene_rpos $geneobj $line1 $snp1]
	set erpos [annotategene_rpos $geneobj $line2 $snp2]
	if {$element1 == $element2} {set span 0} else {set span 1}
	set tr1 [regexp ^ex|in $element1]
	set tr2 [regexp ^ex|in $element2]
	if {($snppos <= $dbstart) && ($snpend >= $dbend)} {
		# gene completely affected
		if {$snptype eq "del"} {
			set impact GENEDEL
		} else {
			set impact GENECOMP
		}
	} elseif {!$tr1 && !$tr2} {
		# completely outside of transcript, no change to impact needed
	} elseif {!$span && [inlist {in sp es} [string range $element1 0 1]]} {
		# both start and end are in the same intron (intron splice or esplice)
		# impact already done in annotategene_one_c
	} elseif {!$span && [inlist {u d p} [string index $element1 0]]} {
		# single element, completely upstream or downstream, no further action needed
	} elseif {![dict exists $geneobj rpstart] || ($srpos > [dict get $geneobj rpend]) || ($erpos < [dict get $geneobj rpstart])} {
		# RNA deletion that does not overlap CDS
		if {!$tr1 && $tr2} {
			# del includes start or end of transcription
			set impact RNASTART
		} elseif {$tr1 && !$tr2} {
			# del includes end of transcription
			set impact RNAEND
		} elseif {$element1 != $element2} {
			# del includes at least one splicesite
			if {[dict exists $geneobj rpstart]} {
				if {$erpos < [dict get $geneobj rpstart]} {
					set impact UTR5SPLICE
					append snp_descr :p.?
				} else {
					set impact UTR3SPLICE
				}
			} else {
				set impact RNASPLICE
			}
		}
	} elseif {($srpos < [expr {[dict get $geneobj rpstart]+3}]) && ($erpos >= [dict get $geneobj rpstart])} {
		# coding deletion overlaps start codon
		set impact CDSSTART[annotategene_p_complex $snptype]
		append snp_descr :p.0?
	} elseif {$span} {
		# coding deletion overlapping splice site
		set impact CDSSPLICE
		append snp_descr :p.?
	} else {
		# coding deletion (no splice site overlap, so line1 == line2)
		append snp_descr :[annotategene_one_p_del $geneobj $snptype $srpos $ref $alt impact]
	}
	if {[regexp ^GENE $impact]} {
		set result $strand[dict get $geneobj transcriptname]
	} else {
		append result :$snp_descr
	}
	return [list $impact $result]
}

proc annotategene_one_p_snp {geneobj snppos from alt line impactVar} {
# putsvars snppos from alt line
	upvar $impactVar impact
	set genomef [dict get $geneobj genomef]
	set complement [dict get $geneobj complement]
	set start [dict get $geneobj grstart]
	set end [dict get $geneobj grend]
	set chr [dict get $geneobj chrom]
	set rpos [annotategene_rpos $geneobj $line $snppos]
	# position on protein (in nucl), starting from 0
	set pnpos [expr {$rpos - [dict get $geneobj rpstart]}]
	# position on protein (in AA), starting from 0
	set prstart [expr {$pnpos/3}]
	#
	if {$complement} {set strand -} else {set strand +}
	set codonstart [expr {3*$prstart}]
	set fromcodon [annotatevar_gene_getprnaseq geneobj $codonstart [expr {$codonstart+2}] $line]
	set cpos [expr {$pnpos-$codonstart}]
	if {$complement} {set alt [seq_complement $alt]}
	set fromAZ [seq_translate $fromcodon]
	set from [string index $fromcodon $cpos]
	set tocodon [string replace $fromcodon $cpos $cpos $alt]
	set toAZ [seq_translate $tocodon]
	if {$fromAZ eq "s"} {
		set impact CDSSTOP
		if {$fromAZ ne $toAZ} {
			set fromseq [string range [annotatevar_gene_rnaseq geneobj] [expr {[dict get $geneobj rpstart]+$codonstart+3}] end]
			set topr [seq_translate $fromseq]
			set pos [string first s $topr]
			# add 2 to pos because first changed is not included in topr, and "string first" counts from 0
			if {$pos < 0} {set tl ?} else {set tl [expr {$pos + 2}]}
			set snp_descr "p.*[expr {$prstart+1}]${toAZ}ext*$tl"
		} else {
			set snp_descr "p.*[expr {$prstart+1}]*"
		}
	} elseif {[inlist {TAA TAG TGA} $tocodon]} {
		set impact CDSNONSENSE
		set snp_descr "p.${fromAZ}[expr {$prstart+1}]*"
	} elseif {$prstart == 0} {
		set impact CDSSTART
		set snp_descr "p.${fromAZ}[expr {$prstart+1}]?"
	} elseif {$fromAZ ne $toAZ} {
		set impact CDSMIS
		set snp_descr "p.${fromAZ}[expr {$prstart+1}]${toAZ}"
	} else {
		set impact CDSsilent
		set snp_descr "p.${fromAZ}[expr {$prstart+1}]${toAZ}"
	}
	return $snp_descr
}

proc annotategene_one_p_ins {geneobj snpend from alt line impactVar} {
# putsvars snpend from alt line
	upvar $impactVar impact
	set genomef [dict get $geneobj genomef]
	set complement [dict get $geneobj complement]
	set start [dict get $geneobj grstart]
	set end [dict get $geneobj grend]
	set chr [dict get $geneobj chrom]
	set rpos [annotategene_rpos $geneobj $line $snpend]
	# position on protein (in nucl), starting from 0
	set pnpos [expr {$rpos - [dict get $geneobj rpstart]}]
	# position on protein (in AA), starting from 0
	set prstart [expr {$pnpos/3}]
	#
	set codonstart [expr {3*$prstart}]
	set cpos [expr {$pnpos-$codonstart}]
	if {$complement} {set alt [seq_complement $alt]}
	set len [string length $alt]
	set pstart [expr {$codonstart/3}]
	if {[isint $alt]} {
		set len $alt
		set fromcodon [annotatevar_gene_getprnaseq geneobj $codonstart [expr {$codonstart+2}] $line]
		set fromAZ [seq_translate $fromcodon]
		if {$fromAZ eq "s"} {
			set impact CDSSTOP
			set snp_descr p.*[expr {$prstart + 1}]Xext*?
		} elseif {[expr {$len%3}]} {
			set impact CDSFRAME
			set snp_descr p.${fromAZ}[expr {$prstart + 1}]Xfs*?
		} elseif {$cpos == 0} {
			set impact CDSINS
			set pseq [annotatevar_gene_pseq geneobj]
			set prevAZ [string index $pseq [expr {$prstart-1}]]
			set snp_descr p.${prevAZ}${prstart}_${fromAZ}[expr {$prstart + 1}]ins([expr {$alt/3}])
		} else {
			if {$prstart == 0} {
				set impact CDSSTART
				set snp_descr p.${fromAZ}[expr {$prstart + 1}]?
			} else {
				set impact CDSINS
				set snp_descr p.${fromAZ}[expr {$prstart + 1}]delins([expr {$len/3+1}])
			}
		}
		return $snp_descr
	}
	set frompr [string range [annotatevar_gene_pseq geneobj] $pstart end]
	set fromseq [string range [annotatevar_gene_rnaseq geneobj] [expr {[dict get $geneobj rpstart]+$codonstart}] end]
	set toseq [string range $fromseq 0 [expr {$cpos-1}]]$alt[string range $fromseq $cpos end]
	set topr [seq_translate $toseq]
	# shift to 3'
	set posfromdiff 0
	string_foreach f $frompr t $topr {
		if {$f eq "s"} break
		if {$t ne $f} break
		incr prstart
		incr codonstart 3
		incr numtostop -1
		incr posfromdiff
	}
	if {$prstart == 0} {
		set impact CDSSTART
		set snp_descr p.${f}1?
	} elseif {$f eq "s"} {
		set impact CDSSTOP
		if {$t eq "s"} {
			set snp_descr p.*[expr {$prstart + 1}]*
		} else {
			set toseq [string range $toseq [expr {3*$posfromdiff}] end]
			set topr [seq_translate $toseq]
			set pos [string first s $topr]
			# add 2 to pos because first changed is not included in topr, and "string first" counts from 0
			if {$pos < 0} {set tl ?} else {set tl [expr {$pos + 2}]}
			set snp_descr p.*[expr {$prstart + 1}]${t}ext*$tl
		}
	} elseif {$t eq "s"} {
		set impact CDSNONSENSE
		set snp_descr p.${f}[expr {$prstart + 1}]*
	} elseif {[expr {$len%3}]} {
		# frame shift
		set numtostop [string first s $topr]
		if {$numtostop < 0} {set tl ?} else {set tl [expr {$numtostop - $posfromdiff}]}
		set impact CDSFRAME
		set snp_descr p.${f}[expr {$prstart + 1}]${t}fs*$tl
	} else {
		set toAZ [string range $topr $posfromdiff [expr {$posfromdiff + $len/3}]]
		set pos [string first s $toAZ]
		if {$pos == -1} {
			set impact CDSINS
		} else {
			set impact CDSNONSENSE
			regsub s $toAZ * toAZ
		}
		if {$f eq [string index $toAZ end]} {
			set pseq [annotatevar_gene_pseq geneobj]
			set toAZ [string range $toAZ 0 end-1] 
			# check for dups/repeats
			set fpos $prstart
			set topos [string length $toAZ]
			while {$topos >= 0} {
				incr fpos -1
				incr topos -1
				set ftest [string index $pseq $fpos]
				set totest [string index $toAZ $topos]
				if {$ftest ne $totest} break
			}
			incr topos
			set rep [string range $toAZ $topos end]
			if {[regexp ^(${rep})+\$ $toAZ]} {
				# dup or repeat
				set repsize [string length $rep]
				# does it go back (should find at least one anyway)
				regexp (${rep})+\$ [string range $pseq 0 [expr {$prstart-1}]] temp
				set repnum [expr {([string length $temp] + [string length $toAZ])/$repsize}]
				if {$repnum == 2} {set type dup} else {set type rep\[$repnum\]}
				if {$repsize == 1} {
					set prevAZ [string index $pseq [expr {$prstart-1}]]
					set snp_descr p.${prevAZ}${prstart}$type
				} else {
					set bpos [expr {$prstart-$repsize}]
					set epos [expr {$prstart-1}]
					set prevAZ [string index $pseq $bpos]
					set snp_descr p.[string index $pseq $bpos][expr {$bpos + 1}]_[string index $pseq $epos][expr {$epos + 1}]$type
				}
			} else {
				set prevAZ [string index $pseq [expr {$prstart-1}]]
				set snp_descr p.${prevAZ}${prstart}_${f}[expr {$prstart + 1}]ins$toAZ
			}
		} else {
			set snp_descr p.${f}[expr {$prstart + 1}]delins$toAZ
			if {$prstart == 0} {append snp_descr ? ; set impact CDSSTART}
		}
	}
	return $snp_descr
}

proc annotategene_one_p_del {geneobj snptype srpos ref alt impactVar} {
# putsvars snptype srpos ref alt
	upvar $impactVar impact
	set genomef [dict get $geneobj genomef]
	set complement [dict get $geneobj complement]
#	set start [dict get $geneobj grstart]
#	set end [dict get $geneobj grend]
	set chr [dict get $geneobj chrom]
	if {[isint $ref]} {set len $ref} else {set len [string length $ref]}
	if {[isint $alt]} {set alen $alt} else {set alen [string length $alt]}
	# position on protein (in nucl), starting from 0
	set pnpos [expr {$srpos - [dict get $geneobj rpstart]}]
	# position on protein (in AA), starting from 0
	set pseq [annotatevar_gene_pseq geneobj]
	set prstart [expr {$pnpos/3}]
	set prend [expr {($pnpos+$len-1)/3}]
	#
	set codonstart [expr {3*$prstart}]
	set cpos [expr {$pnpos-$codonstart}]
	set pstart [expr {$codonstart/3}]
	if {[isint $alt]} {
		set frompr [string range $pseq $prstart $prend]
		set AZ1 [string index $frompr 0]
		set AZ2 [string index $frompr end]
		if {[regexp s $frompr]} {
			set impact CDSSTOP
			set snp_descr [annotategene_p_loc * $prstart]Xext*?
		} elseif {[expr {abs($alen-$len)%3}]} {
			set impact CDSFRAME
			set snp_descr [annotategene_p_loc $AZ1 $prstart]Xfs*?
		} elseif {$cpos == 0} {
			set impact CDS[annotategene_p_complex $snptype]
			if {$alt == 3} {
				set snp_descr [annotategene_p_loc $AZ1 $prstart $AZ2 $prend]X
			} else {
				set snp_descr [annotategene_p_loc $AZ1 $prstart $AZ2 $prend]delins([expr {$alt/3}])
			}
		} else {
			if {$prstart == 0} {
				set impact CDSSTART
				set snp_descr [annotategene_p_loc $AZ1 $prstart]?
			} else {
				set impact CDS[annotategene_p_complex $snptype]
				set snp_descr [annotategene_p_loc $AZ1 $prstart $AZ2 $prend]delins([expr {$len/3+1}])
			}
		}
		return $snp_descr
	}
	set frompr [string range [annotatevar_gene_pseq geneobj] $pstart end]
	set fromseq [string range [annotatevar_gene_rnaseq geneobj] [expr {[dict get $geneobj rpstart]+$codonstart}] end]
	set toseq [string range $fromseq 0 [expr {$cpos-1}]]$alt[string range $fromseq [expr {$cpos+$ref}] end]
	set topr [seq_translate $toseq]
	set posfromdiff 0
	string_foreach f $frompr t $topr {
		if {$f eq "s"} break
		if {$t ne $f} break
		incr prstart
		incr codonstart 3
		incr numtostop -1
		incr posfromdiff
	}
#puts [string range $fromseq 0 40]\n[string range $toseq 0 40]
#putsvars prstart cpos f t
#puts [string range $frompr 0 20]\n[string range $topr 0 20]
	if {$f eq "s"} {
		set impact CDSSTOP
		if {$t eq "s"} {
			set snp_descr [annotategene_p_loc * $prstart]*
		} else {
			set toseq [string range $toseq [expr {3*$posfromdiff}] end]
			set topr [seq_translate $toseq]
			set pos [string first s $topr]
			# add 2 to pos because first changed is not included in topr, and "string first" counts from 0
			if {$pos < 0} {set tl ?} else {set tl [expr {$pos + 2}]}
			set snp_descr [annotategene_p_loc * $prstart]${t}ext*$tl
		}
	} elseif {$t eq "s"} {
		set impact CDSNONSENSE
		set snp_descr [annotategene_p_loc $f $prstart]*
	} elseif {[expr {abs($alen-$len)%3}]} {
		# frame shift
		set numtostop [string first s $topr]
		if {$numtostop < 0} {set tl ?} else {set tl [expr {$numtostop - $posfromdiff}]}
		set impact CDSFRAME
		set snp_descr [annotategene_p_loc $f $prstart]${t}fs*$tl
	} else {
		set impact CDS[annotategene_p_complex $snptype]
		set prend [expr {$prstart + $len/3}]
		set fromp [string range $pseq $prstart $prend]
		set toend [expr {$posfromdiff + $alen/3}]
		set top [string range $topr $posfromdiff $toend]
		while {[string index $top end] eq [string index $fromp end]} {
			incr prend -1
			set top [string range $top 0 end-1]
			set fromp [string range $fromp 0 end-1]
			if {$top eq "" || $fromp eq ""} break
		}
		set snp_descr [annotategene_p_loc [string index $fromp 0] $prstart [string index $fromp end] $prend]
		if {$top ne ""} {
			if {[string length $top] == 1 && [string length $fromp] == 1} {
				append snp_descr $top
			} else {
				append snp_descr delins$top
			}
		} else {
			append snp_descr $snptype
		}
	}
	return $snp_descr
}

proc annotategene_one_c {geneobj chrom snppos snptype from alt line {impactVar {}}} {
# putsvars geneobj chrom snppos snptype from alt line
	if {$impactVar ne ""} {upvar $impactVar impact}
	foreach {gbegin gend eltype element rnabegin rnaend cdsbegin cdsend} $line break
	set complement [dict get $geneobj complement]
	set cds [dict get $geneobj cds]
	if {$cds} {set psize [expr {[dict get $geneobj rpend] - [dict get $geneobj rpstart] + 1}]}
	set rpos [annotategene_rpos $geneobj $line $snppos]
	set rrefpos [calc_rpos $geneobj $rpos]
	set impact [annotate_type2impact $eltype]
	# position on protein (in nucl)
	if {$eltype eq "UTR"} {
		if {$rrefpos > 1} {
			set impact UTR3
		} else {
			set impact UTR5
		}
	}
	if {[dict get $geneobj cds]} {set snp_descr [list c.]} else {set snp_descr [list n.]}
	if {$impact eq "UTR5"} {
		# check for kozak region
		if {$rrefpos > -7} {set impact UTR5KOZAK}
	}
	if {[regexp ^intron $element]} {
		# location for variants in introns
		if {[expr {($gend-$snppos) - ($snppos-$gbegin)}] > 0} {
			if {!$complement} {set ref ${rrefpos}} else {set ref [expr {$rrefpos + 1}]}
			set dir +
			set ipos [expr {$snppos-$gbegin+1}]
			if {$ipos <= 2} {set impact ESPLICE} elseif {$ipos <= 8} {set impact splice}
		} else {
			if {!$complement} {set ref [expr {$rrefpos + 1}]} else {set ref ${rrefpos}}
			set dir -
			set ipos [expr {$gend + 1 - $snppos}]
			if {$ipos <= 2} {set impact ESPLICE} elseif {$ipos <= 8} {set impact splice}
		}
		if {$complement} {
			if {$dir eq "+"} {set dir -} else {set dir +}
		}
		if {$cds && $ref > $psize} {set ref *[expr {$ref-$psize}]}
		lappend snp_descr ${ref}
		lappend snp_descr ${dir}$ipos
	} else {
		if {$cds && $rrefpos > $psize} {set rrefpos *[expr {$rrefpos-$psize}]}
		lappend snp_descr ${rrefpos} {}
	}
	# add change
	if {[isint $alt]} {set alt ($alt)}
	if {$snptype eq "snp"} {
		lappend snp_descr [annotate_compl $geneobj $from]>[annotate_compl $geneobj $alt]
	} else {
		if {$snptype eq "sub"} {set snptype delins}
		lappend snp_descr $snptype[annotate_compl $geneobj $alt]
	}
	return $snp_descr
}

proc annotategene_one {geneobj loc} {
# putsvars loc
	set complement [dict get $geneobj complement]
	foreach {chrom snppos snpend snptype ref alt} $loc break
	set size [expr {$snpend-$snppos}]
	if {$snptype eq "sub" && $size == 1 && [string length $alt] == 1} {set snptype snp}
	if {[inlist {del sub inv amp} $snptype] || $size > 1} {
		# treat deletions and subs separately because they need special care (can span exons, the whole annotation, etc ...)
		return [annotategene_one_del $geneobj $loc]
	} elseif {$snptype eq "ins" || $size == 0} {
		return [annotategene_one_ins $geneobj $loc]
	} else {
		# snps
		foreach line [dict get $geneobj ftlist] {
			foreach {ftstart ftend type} $line break
			if {($snppos >= $ftstart) && ($snppos <= $ftend)} break
		}
		set impact [lindex $line 2]
		if {!$complement} {set strand +} else {set strand -}
		set result $strand[dict get $geneobj transcriptname]
		foreach {element elementpos} [annotategene_element $geneobj $line $snppos] break
		append result :$element$elementpos:[join [annotategene_one_c $geneobj $chrom $snppos $snptype $ref $alt $line impact] {}]
		if {$type eq "CDS"} {
			append result :[annotategene_one_p_snp $geneobj $snppos $ref $alt $line impact]
		}
		return [list $impact $result]
	}
}

proc open_genefile {df dpossVar {genecol {}} {transcriptcol {}}} {
	upvar $dpossVar dposs
	set header [tsv_open $df comment]
	set deffields {chromosome begin end strand exonStarts exonEnds cdsStart cdsEnd transcript gene exonCount}
	set dposs [list_sub [tsv_basicfields $header 14 0] {0 1 2 6 7 8 9 10 11 12}]
	if {$transcriptcol ne ""} {
		set pos [lsearch $header $transcriptcol]
		if {$pos == -1} {
			error "error: given transcriptcol $transcriptcol not found"
		}
		lset dposs end-1 $pos
	}
	if {$genecol ne ""} {
		set pos [lsearch $header $genecol]
		if {$pos == -1} {
			error "error: given genecol $genecol not found"
		}
		lset dposs end $pos
	}
	if {-1 in $dposs} {
		set missingfields [list_sub $deffields [list_find $dposs -1]]
		set missingfields [list_remove $missingfields cdsStart cdsEnd]
		if {[llength $missingfields]} {
			error "error: gene file misses the following fields: $missingfields"
		}
	}
	lappend dposs [lsearch $header exonCount]
	return $header
}

proc annotategene {file genomefile dbfile name annotfile {genecol {}} {transcriptcol {}} {upstreamsize 2000} {limitregion {}}} {
# putsvars file genomefile dbfile name annotfile genecol transcriptcol upstreamsize limitregion
	global genomef

	annot_init
	if {[catch {eof $genomef}]} {
		set genomef [genome_open $genomefile]
	}
	catch {gzclose $f}; catch {gzclose $df}; catch {close $o};
	set f [gzopen $file]
	set header [tsv_open $f]
	if {[catch {set poss [tsv_basicfields $header]}]} {
		if {[catch {set poss [tsv_basicfields $header 4]}]} {
			set poss [tsv_basicfields $header 3]
			lappend poss -1 -1 -1
		} else {
			lappend poss -1 -1
		}
		set noref 1
	} else {
		set noref 0
	}
	set fields [list_sub $header $poss]
	set df [gzopen $dbfile]
	set header [open_genefile $df dposs $genecol $transcriptcol]
	set dbposs [lrange $dposs 0 2]
	set o [open $annotfile.temp w]
	set nh [list ${name}_impact ${name}_gene ${name}_descr]
	puts $o [join $nh \t]
	set empty [join [list_fill [llength $nh] {}] \t]
	while {![eof $df]} {
		set fdbline [split [gets $df] \t]
		set dbloc [list_sub $fdbline $dbposs]
		foreach {dbchr dbstart dbend} $dbloc break
		if {![isint $dbstart] || ![isint $dbend]} continue
		break
	}
	set dbstart [expr {$dbstart - $upstreamsize}]
	incr dbend $upstreamsize
	lset dbloc 0 $dbchr
	set dblist {}
	set counter 0
	set prevloc ""
	set prevdbloc ""
	foreach {limitchr limitstart limitend} {{} {} {}} break
	foreach {limitchr limitstart limitend} [split $limitregion :-] break
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} {
			if {[eof $f]} break
			error "error: empty line in $file"
		}
		set loc [list_sub $line $poss]
		foreach {chr start end type ref alt} $loc break
		if {$limitregion ne ""} {
			set comp [distr_chr_compare $chr $start $limitchr $limitstart $limitend]
			if {$comp > 0} break
			if {$comp != 0} continue
		}
		set ploc [lrange $loc 0 2]
		if {[lloc_compare $prevloc $ploc] > 0} {
			error "Cannot annotate because the variant file is not correctly sorted (sort correctly using \"cg select -s -\")"
		}
		set prevloc $ploc
		if {$type eq ""} {
			lset loc 3 del
			set type del
		}
		if {$start > $end} {
			error "location start > end error: $loc"
		}
		if {$noref} {
			switch $type {
				snp {
					set ref N; set alt N
					lset loc end N
				}
				default {set ref [expr {$end-$start}]}
			}
			lset loc end-1 $ref
		}
		incr counter
		if {$counter > 50000} {
			putslog $chr:$start
			set counter 0
		}
		# remove genes from dblist that are before current var in remove
		set newdblist {}
		foreach dbline $dblist {
			foreach {dc ds de} $dbline break
			set chrcompar [chr_compare $dc $chr]
			# previous chromosome, remove
			if {$chrcompar < 0} continue
			# of end of gene anotation is before current variant
			# de is already corrected for upstream
			if {$chrcompar == 0 && $de <= $start} continue
			lappend newdblist $dbline
		}
		set dblist $newdblist
		# 
		# add all overlapping to dblist
		while {![eof $df]} {
			set chrcompar [chr_compare $dbchr $chr]
			if {$chrcompar > 0} break
			if {$chrcompar == 0} {
				if {$dbstart >= $end} break
				if {$dbend >= $start} {
					lappend dblist [list $dbchr $dbstart $dbend {} $fdbline]
				}
			} else {
				set dblist {}
			}
			set ok 0
			while {![eof $df]} {
				set fdbline [split [gets $df] \t]
				set dbloc [list_sub $fdbline $dbposs]
				foreach {dbchr dbstart dbend} $dbloc break
				if {[isint $dbstart] && [isint $dbend]} {
					set ok 1
					break
				}
			}
			set pdbloc [lrange $dbloc 0 2]
			#break if eof
			if {[eof $df]} break
			if {[lloc_compare $prevdbloc $pdbloc] > 0} {
				error "Cannot annotate because the database file ($dbfile) is not correctly sorted (sort correctly using \"cg select -s -\")"
			}
			set prevdbloc $pdbloc
			lset dbloc 0 $dbchr
			if {!$ok} break
			set dbstart [expr {$dbstart - $upstreamsize}]
			incr dbend $upstreamsize
		}

		if {![llength $dblist]} {
			puts $o $empty
			continue
		}
		# join [list_subindex $dblist 4] \n\n
		# check for multiple alleles, process these separately (alist contains >1 loc)
		set malt [split $alt ,]
		if {[llength $malt] > 1} {
			set alist {}
			foreach a $malt {
				lappend alist [lreplace $loc end end $a]
			}
		} else {
			set alist [list $loc]
		}
		# go over each allele (alist), collect results per allele in ahitgenes
		set ahitgenes {}
		foreach loc $alist {
			# check for overlap, remove genes from dblist that are before current var in remove
			# annotate overlapping with annotategene_one, collect results in hitgenes
			set num 0
			set remove {}
			set hitgenes ""
			# join [list_subindex $dblist 4] \n\n
			list_foreach {dc ds de geneobj dbline} $dblist {
				set chrcompar [chr_compare $dc $chr]
				# gene cannot be before var (were removed before)
				if {$chrcompar == 0 && $ds < $end} {
					if {[catch {dict get $geneobj end}]} {
						set geneobj [annotatevar_gene_makegeneobj $genomef $dbline $dposs $upstreamsize]
						lset dblist $num 3 $geneobj
					}
					set genename [dict get $geneobj genename]
					set result [annotategene_one $geneobj $loc]
					lappend hitgenes [linsert $result 1 $genename]
				}
				incr num
			}
			# create result for this allele
			set result {}
			set impacts [list_subindex $hitgenes 0]
			set udspos [list_union [list_find $impacts upstream] [list_find $impacts downstream]]
			if {[llength $udspos] != [llength $impacts]} {
				set hitgenes [list_sub $hitgenes -exclude $udspos]
			}
			if {[llength $hitgenes] == 1} {
				set result [lindex $hitgenes 0]
			} elseif {[llength $hitgenes]} {
				if {[lsearch $impacts GENEDEL] != -1} {
					set genes [list_remdup [list_subindex $hitgenes 1]]
					set genes {}
					set descrs {}
					set impacts {}
					set prevgene {}
					list_foreach {impact gene descr} $hitgenes {
						if {$gene eq $prevgene} continue
						lappend impacts $impact
						lappend genes $gene
						if {$impact eq "GENEDEL"} {
							lappend descrs $gene:del
						} else {
							lappend descrs $descr
						}
						set prevgene $gene
					}
					set result [list [join $impacts \;] [join $genes \;] [join $descrs \;]]
				} else {
					set pos 0
					foreach h $nh {
						set list [list_subindex $hitgenes $pos]
						set nrlist [list_remdup $list]
						if {[llength $nrlist] == 1} {
							lappend result [lindex $nrlist 0]
						} else {
							lappend result [join $list \;]
						}
						incr pos
					}
				}
			}
			lappend ahitgenes $result
		}
		if {[llength $remove]} {
			set dblist [list_sub $dblist -exclude $remove]
		}
		# create final result
		if {[llength $ahitgenes] == 1} {
			set line [lindex $ahitgenes 0]
			if {[llength $line]} {
				set result [join $line \t]
			} else {
				set result $empty
			}
		} elseif {[llength $ahitgenes]} {
			set result {}
			set pos 0
			foreach h $nh {
				set list [list_subindex $ahitgenes $pos]
				set nrlist [list_remdup $list]
				if {[llength $nrlist] == 1} {
					lappend result [lindex $nrlist 0]
				} else {
					lappend result [join $list \,]
				}
				incr pos
			}
			set result [join $result \t]
		} else {
			set result $empty
		}
		puts $o $result
	}
	close $genomef
	close $o; gzclose $f;	gzclose $df
	file rename -force -- $annotfile.temp $annotfile

}

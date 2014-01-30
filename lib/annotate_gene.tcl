#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require BioTcl

proc annot_init {} {
	global snp_annot_list snp_annot_score snp_annot_loc
	set snp_annot_list {
		{downstream annot downstream} {upstream annot upstream}
		{annot annot annotation_region} {A annot annotated}
		{intron intron intron}
		{reg reg regulatory} {prom reg promotor} {splice splice splice}
		{RNA transcript RNA} {RNASPLICE transcript RNA_splice}
		{UTR3 UTR3 3'UTR} {UTR3SPLICE UTR3 3'UTR_splice} {RNAEND transcript transcription_end}
		{ESPLICE splice essentialsplice}
		{CDSsilent CODING silent} {UTR5 UTR5 5'UTR} {UTR5SPLICE UTR5 5'UTR_splice} {RNASTART transcript transcription_start}
		{CDSMIS CODING missense} {CDSDEL CODING deletion} {CDSINS CODING insertion} {CDSSPLICE CODING splice}
		{CDSNONSENSE CODING nonsense}
		{CDSDELSPLICE CODING splice deletion}
		{CDSSTOP CODING readthrough} {CDSFRAME CODING frameshift} {CDSSTART CODING startcodon} {CDSDELSTART CODING startcodon}
		{GENEDEL GENE genedeletion}
	}
	unset -nocomplain snp_annot_score
	set num 1
	list_foreach {code location name} $snp_annot_list {
		set snp_annot_score($code) $num
		set snp_annot_loc($code) $location
		incr num
	}
}

proc annotate_compl {seq} {
	global adata
	if {!$adata(complement)} {
		return $seq
	} else {
		return [seq_complement $seq]
	}
}

proc eipos eiposVar {
	upvar $eiposVar eipos
	if {![info exists eipos]} {
		return 0
	} elseif {$eipos >= 0} {
		return [expr {$eipos+1}]
	} else {
		return $eipos
	}
}

proc annotate_type2annot {type} {
	switch $type {
		esplice {set snp_annot ESPLICE}
		splice {set snp_annot splice}
		promotor {set snp_annot prom}
		regulatory {set snp_annot reg}
		intron {set snp_annot intron}
		upstream {set snp_annot upstream}
		downstream {set snp_annot downstream}
		annotation_region {set snp_annot annot}
		default {set snp_annot A}
	}
	return $snp_annot
}

proc annotate_getregion {refseq start end} {
	if {$start < $end} {
		return [string range $refseq $start $end]
	} else {
		return [seq_complement [string range $refseq $end $start]]
	}
}

proc switchval {var1 var2} {
	uplevel 1 [subst {
		set temp \$$var1
		set $var1 \$$var2
		set $var2 \$temp
	}]
}

proc annotate_complementpos {pos} {
	global adata
	return [expr {$adata(refseqlen)-$pos-1}]
}

proc calc_rpos {rsnppos} {
	global adata
	set rpos [expr {$rsnppos - $adata(rrefstart)}]
	if {$rpos >= 0} {
		return [expr {$rpos+1}]
	} else {
		return $rpos
	}
}

proc annotatevar_gene_makegeneobj {genomef dbline dposs} {
	global adata
	# splice 8 bases, essentialsplice 2 bases intron, upstream,downstream -> 2k bases
	foreach {dchrom dstart dend strand cdsStart cdsEnd exonCount exonStarts exonEnds transcriptname genename} [list_sub $dbline $dposs] break
	if {$strand eq "-"} {set complement 1} else {set complement 0}
	unset -nocomplain adata
	set adata(genomef) $genomef
	set adata(complement) $complement
	set adata(start) $dstart
	set adata(end) $dend
	set adata(transcriptname) $transcriptname
	set adata(genename) $genename
	set adata(chrom) $dchrom
	if {$cdsStart == $cdsEnd} {
		set cdsStart -1
		set cdsEnd -1
		set cds 0
	} else {
		set cds 1
	}
	set exonStarts [split [string trim $exonStarts ,] ,]
	set exonEnds [split [string trim $exonEnds ,] ,]
	set adata(nrexons) [llength $exonStarts]
	# use code coming from novoSNP:
	# transform to ft addressing used in novoSNP previously (ft endpos is included, so -1)

	set ftlist {}
	if {!$complement} {set type upstream} else  {set type downstream}
	lappend ftlist [list 0 [expr {$dstart-2001}] pre]
	lappend ftlist [list [expr {$dstart-2000}] [expr {$dstart-1}] $type]
	set prev -1
	if {$cds} {set type UTR} else {set type RNA}
	foreach s $exonStarts e $exonEnds {
# puts --------\n[join $ftlist \n]\n
# if {$s == 52847059} {error STOPPED}
# putsvars s e type prev
		if {$prev != -1} {
			# s and prev positions are not included -> -1
			set isize [expr {$s-$prev-1}]
			if {$isize >= 16} {
				lappend ftlist [list [expr {$prev}] [expr {$prev+1}] esplice]
				lappend ftlist [list [expr {$prev+2}] [expr {$prev+7}] splice]
				lappend ftlist [list [expr {$prev+8}] [expr {$s-9}] intron]
				lappend ftlist [list [expr {$s-8}] [expr {$s-3}] splice]
				lappend ftlist [list [expr {$s-2}] [expr {$s-1}] esplice]
			} elseif {$isize >= 4} {
				lappend ftlist [list [expr {$prev}] [expr {$prev+1}] esplice]
				lappend ftlist [list [expr {$prev+2}] [expr {$s-3}] splice]
				lappend ftlist [list [expr {$s-2}] [expr {$s-1}] esplice]
			} elseif {$isize >= 0} {
				lappend ftlist [list [expr {$prev}] [expr {$s-1}] esplice]
			}
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
	lappend ftlist [list [expr {$e}] [expr {$e+1999}] $type]
	lappend ftlist [list [expr {$e+2000}] 9000000000 post]
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
					set adata(gpstart) $e
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
		} else {
			if {$exon == 0} {set temp pre} else {set temp intron$exon}
			lappend tempftlist [list $s $e $t $temp $rnaend $rnaend $cdsend $cdsend]
		}
		incr num
	}
	set post [lindex $tempftlist end 3]
	set num [llength $tempftlist]
	while {$num > 0} {
		incr num -1
		set temp [lindex $tempftlist $num 3]
		if {$temp ne $post} break
		lset tempftlist $num 3 post
	}
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
		set adata(lastcodon) [expr {($adata(rrend)-1)/3}]
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

proc annotatevar_gene_getrnaseq {qstart qend {line {}}} {
	global adata
	set genomef $adata(genomef)
	set chr $adata(chrom)
	set complement $adata(complement)
	if {!$complement} {
		if {[llength $line]} {
			foreach {gs ge t el rs re ps pe} $line break
			if {$qstart >= $rs && $qend <= $re} {
				set gstart [expr {$gs+$qstart-$rs}]
				set gend [expr {$gs+$qend-$rs}]
				incr gend ;# genome_get works with half open coordinates
				set result [genome_get $genomef $chr $gstart $gend]
				if {!$complement} {return $result} else {return [seq_complement $result]}
			}
		}
		set result {}
		# join $adata(rnalist) \n
		list_foreach {gs ge t el rs re ps pe} $adata(rnalist) {
			if {$re >= $qstart} {
				set gstart [expr {$gs+$qstart-$rs}]
				set gend [expr {$gs+$qend-$rs}]
				if {$gstart < $gs} {set gstart $gs}
				if {$gend > $ge} {set gend $ge}
				incr gend ;# genome_get works with half open coordinates
				append result [genome_get $genomef $chr $gstart $gend]
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
				set result [genome_get $genomef $chr $gstart $gend]
				if {!$complement} {return $result} else {return [seq_complement $result]}
			}
		}
		set result {}
		# join $adata(rnalist) \n
		list_foreach {gs ge t el rs re ps pe} $adata(rnalist) {
			if {$rs <= $qend} {
				# for complement, take re as reference, but start from gs
				# we want to extract genomic non complement first
				set gstart [expr {$gs+$re-$qend}]
				set gend [expr {$gs+$re-$qstart}]
				if {$gstart < $gs} {set gstart $gs}
				if {$gend > $ge} {set gend $ge}
				incr gend ;# genome_get works with half open coordinates
				append result [genome_get $genomef $chr $gstart $gend]
				if {$rs < $qstart} break
			}
		}
		set result [seq_complement $result]
	}
	return $result
}

proc annotatevar_gene_getprnaseq {qstart qend {line {}}} {
	global adata
	incr qstart $adata(rpstart)	
	incr qend $adata(rpstart)	
	annotatevar_gene_getrnaseq $qstart $qend $line
}

proc annotategene_one_getsnpcoords {line snppos} {
	global adata
	set complement $adata(complement)
	foreach {gs ge t el rs re ps pe} $line break
	if {!$complement} {
		# position in intron/exon
		set eipos [expr {$snppos-$gs}]
		# position in rna
		if {$rs == $re} {
			# in intron
			set rpos $rs
		} else {
			set rpos [expr {$rs+$eipos}]
		}
	} else {
		# position in intron/exon
		set eipos [expr {$ge-$snppos}]
		# position in rna
		if {$rs == $re} {
			# in intron
			set rpos $rs
		} else {
			set rpos [expr {$rs+$eipos}]
		}
	}
	# position on protein (in nucl)
	return [list $rpos $el $eipos]
}

proc annotategene_one_del {snppos snptype ref alt} {
#putsvars snppos snptype ref alt
	global adata snp_annot_score
	set impacttype [string toupper $snptype]
	set complement $adata(complement)
	set dbstart $adata(start)
	set dbend $adata(end)
	set snp_descr {}
	if {[isint $ref]} {
		set len $ref
	} else {
		set len [string length $ref]
	}
	if {[isint $alt]} {
		set rlen $alt
	} else {
		set rlen [string length $alt]
	}
	set snpend [expr {$snppos+$len-1}]
	# find start and end location line in ftlist
	set sline pre
	set eline post
	foreach line $adata(ftlist) {
		foreach {ftstart ftend type} $line break
		if {($snppos >= $ftstart) && ($snppos <= $ftend)} {set sline $line}
		if {($snpend >= $ftstart) && ($snpend <= $ftend)} {set eline $line; break}
	}
	if {!$complement} {
		foreach {srpos sel seipos} [annotategene_one_getsnpcoords $sline $snppos] break
		foreach {erpos eel eeipos} [annotategene_one_getsnpcoords $eline $snpend] break
		set gpos [expr {$snppos - $adata(grefstart)}]
	} else {
		foreach {srpos sel seipos} [annotategene_one_getsnpcoords $eline $snpend] break
		foreach {erpos eel eeipos} [annotategene_one_getsnpcoords $sline $snppos] break
		set gpos [expr {$adata(grefstart) - $snpend}]
	}
	if {($snppos <= $dbstart) && ($snpend >= $dbend)} {
		set snp_annot GENE$impacttype
		set snp_descr $snptype
	} elseif {([string index $sel 0] ne "e") && ($sel == $eel)} {
		# both start and end are in the same intron, so just check which features are affected
		# pick out the most intersting one
		set annlist {}
		list_foreach {ts te ttype} $adata(ftlist) {
			if {[inlist {post pre} $ttype]} continue
			if {$snppos < $te && $snpend >= $ts} {
				lappend annlist $ts $te $ttype
			}
		}
		if {[llength $annlist] == 0} {
			if {(($snppos > $adata(grstart)) && ($snppos < $adata(grend))) || (($snpend > $adata(grstart)) && ($snpend < $adata(grend)))} {
				set snp_annot intron
			} else {
				set snp_annot annot
			}
		} else {
			set best_snp_annot_score -1
			foreach {tempstart tempend temptype} $annlist {
				set snp_annot [annotate_type2annot $temptype]
				if {[get snp_annot_score($snp_annot) 0] > $best_snp_annot_score} {
					set best_snp_annot $snp_annot
					set best_snp_annot_score [get snp_annot_score($snp_annot) 0]
					set dbstart $tempstart
					set dbend $tempend
				}
			}
			set snp_annot $best_snp_annot
		}
		set indel [expr {$snpend-$snppos}]
		if {!$rlen} {
			if {!$complement} {
				set snp_descr g.${gpos}_[expr {$gpos+$len-1}]$snptype$indel
			} else {
				set snp_descr g.[expr {$gpos-$len+1}]_${gpos}$snptype$indel
			}
		} else {
			if {$rlen > 70} {set alt $rlen}
			if {!$complement} {
				set snp_descr g.${gpos}_[expr {$gpos+$len-1}]del$indel-ins$alt
			} else {
				set snp_descr g.[expr {$gpos-$len+1}]_${gpos}del$indel-ins$alt
			}
		}
	} elseif {![info exists adata(rpstart)] || ($srpos > $adata(rpend)) || ($erpos < $adata(rpstart))} {
		# RNA deletion that does not overlap CDS
		if {[string index $sel 0] ne "e"} {
			# del includes start of transcription
			set snp_annot RNASTART
			set snp_descr r.[expr {- $adata(rrefstart)}]_[calc_rpos $erpos]ts
		} elseif {[string index $eel 0] ne "e"} {
			# del includes end of transcription
			set snp_annot RNAEND
			set snp_descr r.[calc_rpos $srpos]_[calc_rpos $adata(rrend)]te
		} elseif {$sel != $eel} {
			# del includes at least one splicesite
			if {[info exists adata(rpstart)]} {
				if {$erpos < $adata(rpstart)} {
					set snp_annot UTR5SPLICE
				} else {
					set snp_annot UTR3SPLICE
				}
			} else {
				set snp_annot RNASPLICE
			}
			set snp_descr r.[calc_rpos $srpos]_[calc_rpos $erpos]sd
		} else {
			if {[info exists adata(rpstart)]} {
				if {$erpos < $adata(rpstart)} {
					set snp_annot UTR5
				} else {
					set snp_annot UTR3
				}
			} else {
				set snp_annot RNA
			}
			set snp_descr r.[calc_rpos $srpos]_[calc_rpos $erpos]del[annotatevar_gene_getrnaseq $srpos $erpos]
			if {$rlen} {append snp_descr -ins$alt}
		}
	} elseif {($srpos < [expr {$adata(rpstart)+3}]) && ($erpos >= $adata(rpstart))} {
		# coding deletion overlaps start codon
		set snp_annot CDS${impacttype}START
		set snp_descr p.d?
	} elseif {$sel != $eel} {
		# coding deletion overlapping splice site
		set rps [expr {$srpos - $adata(rpstart)}]
		set rpe [expr {$erpos - $adata(rpstart)}]
		set codonpos [expr {$rps/3}]
		set codonstart [expr {3*$codonpos}]
		set fromcodon [annotatevar_gene_getprnaseq $codonstart [expr {$codonstart+2}]]
		set fromAZ [seq_translate $fromcodon]
		set snp_annot CDS${impacttype}SPLICE
		set snp_descr p.${fromAZ}[expr {$codonpos+1}]Xsd
	} else {
		# coding deletion
		set rps [expr {$srpos - $adata(rpstart)}]
		set rpe [expr {$erpos - $adata(rpstart)}]
		set codonpos [expr {$rps/3}]
		set codonstart [expr {3*$codonpos}]
		if {[expr {[expr {$len-$rlen}]%3}]} {
			# frame shift
			set fromcodon [annotatevar_gene_getprnaseq $codonstart [expr {$codonstart+2}] $sline]
			set fromAZ [seq_translate $fromcodon]
			set snp_annot CDSFRAME
			set fromAZ [seq_translate $fromcodon]
			set snp_descr p.${fromAZ}[expr {$codonpos+1}]Xfs
		} else {
			set snp_annot CDS${impacttype}
			set cpos [expr {$rps-$codonstart}]
			set endpos [expr {($rpe/3)*3 + 2}]
			set fromcodon [annotatevar_gene_getprnaseq $codonstart $endpos $sline]
			set fromAZ [seq_translate $fromcodon]
			set tocodon [string range $fromcodon 0 [expr {$cpos-1}]]$alt[string range $fromcodon [expr {$cpos+$len}] end]
			set toAZ [seq_translate $tocodon]
			if {[string index $toAZ end] eq [string index $fromAZ end]} {
				set toAZ [string range $toAZ 0 end-1]
				set fromAZ [string range $fromAZ 0 end-1]
			}
			if {[string index $toAZ 0] eq [string index $fromAZ 0]} {
				set toAZ [string range $toAZ 1 end]
				set fromAZ [string range $fromAZ 1 end]
				incr codonpos
			}
			if {[string length $fromAZ] == 1} {
				set snp_descr p.${fromAZ}[expr {$codonpos+1}]del
			} elseif {[string length $toAZ] == 0} {
				set snp_descr p.[string index $fromAZ 0][expr {$codonpos+1}]_[string index $fromAZ end][expr {$codonpos+[string length $fromAZ]}]del
			} else {
				set snp_descr p.[string index $fromAZ 0][expr {$codonpos+1}]_[string index $fromAZ end][expr {$codonpos+[string length $fromAZ]}]delins$toAZ
			}
		}
	}
	if {$adata(complement)} {set strand -} else {set strand +}
	if {$snp_annot eq "GENE${impacttype}"} {
		set p $strand$adata(transcriptname)
	} elseif {$sel eq $eel} {
		set p $strand$adata(transcriptname):$sel:$seipos
	} else {
		set p $strand$adata(transcriptname):$sel:$eeipos
		append p -$eel:$eeipos
	}
	append p :$snp_descr
	return [list $snp_annot $p]
}

proc annotategene_one_CDS {snppos snptype from to line} {
#putsvars snppos snptype from to line
	global adata snp_annot_score
	set genomef $adata(genomef)
	set complement $adata(complement)
	set start $adata(grstart)
	set end $adata(grend)
	set chr $adata(chrom)
	set snp_descr {}
	foreach {rpos el eipos} [annotategene_one_getsnpcoords $line $snppos] break
	# position on protein (in nucl)
	set ppos [expr {$rpos - $adata(rpstart)}]
	#
	set csnp_descr c.[annotate_compl $from][calc_rpos $rpos][annotate_compl $to]
	set codonpos [expr {$ppos/3}]
	set codonstart [expr {3*$codonpos}]
	set fromcodon [annotatevar_gene_getprnaseq $codonstart [expr {$codonstart+2}] $line]
	set cpos [expr {$ppos-$codonstart}]
	if {$complement} {set to [seq_complement $to]}
	set snp_annot {}
	set snp_descr {}
	switch $snptype {
		snp {
			set tocodon [string replace $fromcodon $cpos $cpos $to]
			set fromAZ [seq_translate $fromcodon]
			set toAZ [seq_translate $tocodon]
			if {[inlist {TAA TAG TGA} $tocodon]} {
				lappend snp_annot CDSNONSENSE
				lappend snp_descr "p.${fromAZ}[expr {$codonpos+1}]stop"
			} elseif {$codonpos == 0} {
				lappend snp_annot CDSSTART
				lappend snp_descr "p.${fromAZ}[expr {$codonpos+1}]${toAZ}"
			} elseif {($codonpos == $adata(lastcodon)) && ($fromAZ ne $toAZ)} {
				lappend snp_annot CDSSTOP
				lappend snp_descr "p.${fromAZ}[expr {$codonpos+1}]${toAZ}r"
			} elseif {$fromAZ ne $toAZ} {
				lappend snp_annot CDSMIS
				lappend snp_descr "p.${fromAZ}[expr {$codonpos+1}]${toAZ}"
			} else {
				lappend snp_annot CDSsilent
				lappend snp_descr "p.${fromAZ}[expr {$codonpos+1}]${toAZ}"
			}
		}
		ins {
			set len [string length $to]
			if {[expr {$len%3}]} {
				# frame shift
				lappend snp_annot CDSFRAME
				set fromAZ [seq_translate $fromcodon]
				lappend snp_descr p.${fromAZ}[expr {$codonpos+1}]Xfs
			} else {
				lappend snp_annot CDSINS
				set fromAZ [seq_translate $fromcodon]
				if {$cpos != 0} {
					set tocodon [string range $fromcodon 0 [expr {$cpos-1}]]$to[string range $fromcodon $cpos end]
					set toAZ [seq_translate $tocodon]
					if {[string index $toAZ 0] ne "$fromAZ"} {
						lappend snp_descr p.${fromAZ}[expr {$codonpos+1}]delins$toAZ
					} else {
						set nextcodon [annotatevar_gene_getprnaseq [expr {$codonstart+3}] [expr {$codonstart+5}]]
						set nextAZ [seq_translate $nextcodon]
						lappend snp_descr p.${fromAZ}[expr {$codonpos+1}]_${nextAZ}[expr {$codonpos+2}]ins[string range $toAZ 1 end]
					}
				} else {
					set prevcodon [annotatevar_gene_getprnaseq [expr {$codonstart-3}] [expr {$codonstart-1}]]
					set prevAZ [seq_translate $prevcodon]
					set toAZ [seq_translate $to]
					lappend snp_descr p.${prevAZ}[expr {$codonpos}]_${fromAZ}[expr {$codonpos+1}]ins$toAZ
				}
			}
		}
	}
	if {$adata(complement)} {set strand -} else {set strand +}
	return [list ${snp_annot} $strand$adata(transcriptname):$el:$eipos:$csnp_descr:$snp_descr]
}

proc annotategene_one_rna {snppos snptype from to line} {
	global adata snp_annot_score
	set complement $adata(complement)
	set snp_annot RNA
	set snp_descr {}
	foreach {rpos el eipos} [annotategene_one_getsnpcoords $line $snppos] break
	# position on protein (in nucl)
	set rrefpos [calc_rpos $rpos]
	if {[lindex $line 2] eq "RNA"} {
		set snp_annot RNA
	} elseif {$rrefpos > 0} {
		set snp_annot UTR3
	} else {
		set snp_annot UTR5
	}
	if {$snptype eq "snp"} {
		set snp_descr r.${rrefpos}[annotate_compl $from]>[annotate_compl $to]
	} elseif {$snptype eq "ins"} {
		if {!$complement} {
			set snp_descr r.[expr {$rrefpos-1}]_${rrefpos}ins$to
		} else {
			set snp_descr r.${rrefpos}_[expr {$rrefpos+1}]ins[seq_complement $to]
		}
	} else {
		set snp_descr ""
	}
	if {$adata(complement)} {set strand -} else {set strand +}
	set p $strand$adata(transcriptname):$el:$eipos
	return [list $snp_annot $p:$snp_descr]
}

proc annotategene_one_g {snppos snptype from to type line} {
	global adata snp_annot_score
	set complement $adata(complement)
	# non spliced annotations
	set snp_annot [annotate_type2annot $type]
	foreach {rpos el eipos} [annotategene_one_getsnpcoords $line $snppos] break
	set grefpos [expr {$snppos - $adata(grefstart)}]
	if {$snptype eq "snp"} {
		set snp_descr g.${grefpos}[annotate_compl $from]>[annotate_compl $to]
	} elseif {$snptype eq "ins"} {
		if {!$complement} {
			set snp_descr g.[expr {$grefpos-1}]_${grefpos}ins$to
		} else {
			set snp_descr g.${grefpos}_[expr {$grefpos+1}]ins[seq_complement $to]
		}
	} else {
		set snp_descr ""
	}
	if {$adata(complement)} {set strand -} else {set strand +}
	set p $strand$adata(transcriptname):$el:$eipos
	return [list $snp_annot $p:$snp_descr]
}

proc annotategene_one {loc geneobj} {
	global adata
	unset -nocomplain adata
	array set adata $geneobj
	set complement $adata(complement)
	foreach {chrom snppos snpend snptype ref alt} $loc break
	if {[inlist {del sub inv amp} $snptype]} {
		# treat deletions and subs separately because they need special care (can span exons, the whole annotation, etc ...)
		set ref [expr {$snpend-$snppos}]
		set result [annotategene_one_del $snppos $snptype $ref $alt]
	} else {
		# other annotations
		foreach line $adata(ftlist) {
			foreach {ftstart ftend type} $line break
			if {($snppos >= $ftstart) && ($snppos <= $ftend)} break
		}
		switch $type {
			CDS {
				set result [annotategene_one_CDS $snppos $snptype $ref $alt $line]
			}
			UTR - RNA {
				set result [annotategene_one_rna $snppos $snptype $ref $alt $line]
			}
			default {
				set result [annotategene_one_g $snppos $snptype $ref $alt $type $line]
			}
		}
	}
	return $result
}

proc annotategene {file genomefile dbfile name annotfile {genecol name2} {transcriptcol name}} {
	global genomef
	annot_init

	if {[catch {eof $genomef}]} {
		set genomef [genome_open $genomefile]
	}
	catch {close $f}; catch {close $df}; catch {close $o};
	set f [gzopen $file]
	set header [tsv_open $f comment]
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
	set header [tsv_open $df]
	set deffields {strand cdsStart cdsEnd exonCount exonStarts exonEnds}
	lappend deffields $transcriptcol $genecol
	set dposs [tsv_basicfields $header 3]
	lappend dposs {*}[list_cor $header $deffields]
	if {[lindex $dposs end] == -1} {
		foreach testfield {gene_name gene_id name} {
			set pos [lsearch $header $testfield]
			if {$pos != -1} {
				lset dposs end $pos
				break
			}
		}
	}
	set dbposs [lrange $dposs 0 2]
	if {[lsearch [lrange $dposs 0 end-2] -1] != -1} {
		puts stderr "error: gene file $dbfile misses the following fields: [list_sub $deffields [list_find [lrange $dposs 0 end-2] -1]]"
		exit 1
	}
	set o [open $annotfile.temp w]
	puts -nonewline $o [join [list_fill [expr {[llength [split $comment \n]]-1}] \n]]
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
	incr dbstart -2000
	incr dbend 2000
	lset dbloc 0 $dbchr
	set dblist {}
	set counter 0
	set prevloc ""
	set prevdbloc ""
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} {
			if {[eof $f]} break
			error "error: empty line in $file"
		}
		set loc [list_sub $line $poss]
		set ploc [lrange $loc 0 2]
		if {[lloc_compare $prevloc $ploc] > 0} {
			error "Cannot annotate because the variant file is not correctly sorted (sort correctly using \"cg select -s -\")"
		}
		set prevloc $ploc
		foreach {chr start end type ref alt} $loc break
		if {$type eq ""} {set type del}
		if {$start > $end} {
			puts stderr "location start > end error: $loc"
			exit 1
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
		# add all overlapping to dblist
# if {$start == 43198434} {error STOPPED}
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
			incr dbstart -2000
			incr dbend 2000
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
		# check for overlap, remove genes from dblist that are before current var
		set ahitgenes {}
		foreach loc $alist {
			set num 0
			set remove {}
			set hitgenes ""
			# join [list_subindex $dblist 4] \n\n
			list_foreach {dc ds de geneobj dbline} $dblist {
				set chrcompar [chr_compare $dc $chr]
				if {$chrcompar < 0} {
					lappend remove $num
				} elseif {$chrcompar == 0} {
					if {$de <= $start} {
						lappend remove $num
					} elseif {$ds < $end} {
						if {[catch {dict get $geneobj end}]} {
							set geneobj [annotatevar_gene_makegeneobj $genomef $dbline $dposs]
							lset dblist $num 3 $geneobj
						}
						set genename [dict get $geneobj genename]
						set result [annotategene_one $loc $geneobj]
						lappend hitgenes [linsert $result 1 $genename]
					}
				}
				incr num
			}
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
	close $o; catch {close $f};	catch {close $df}
	file rename -force $annotfile.temp $annotfile

}

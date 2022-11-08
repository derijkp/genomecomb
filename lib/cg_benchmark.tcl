proc cg_benchmarkvars {args} {
	set analyses {}
	set regionfile {}
	set refcurve 0
	set seqcode {\$sequenced-$analysis}
	set curvecode {if(lmaxd(\$genoqual-$analysis,0)>100,100,round(lmaxd(\$genoqual-$analysis,0)))}
	set refcurve_cutoffs {{} 20 30 40 50 80}
	cg_options benchmarkvars args {
		-a - -analyses {
			set analyses $value
		}
		-r - -regionfile {
			set regionfile $value
		}
		-refcurve {
			set refcurve $value
		}
		-refcurve_cutoffs {
			set refcurve_cutoffs $value
		}
		-seqcode {
			set seqcode $value
		}
		-curvecode {
			set curvecode $value
		}
	} {comparfile refanalysis resultfile} 3 3
	set resultfileall [file root [gzroot $resultfile]]-all.tsv
	if {$analyses eq ""} {
		set analyses [cg select -a $comparfile]
	}
	set analyses [list_remove $analyses $refanalysis]
	set groups {}
	set fields [split [cg select -h $comparfile] \n]
	set analysis $refanalysis
	lappend fields	"sr=[subst $seqcode]"
	if {$refcurve} {
		lappend fields	"refcurve=[subst $curvecode]"
		lappend groups refcurve *
	}
	lappend groups type *
	foreach analysis $analyses {
		lappend fields "seq-$analysis=[subst $seqcode]"
		lappend fields "curve-$analysis=[subst $curvecode]"
		lappend fields "curve-$analysis=[subst $curvecode]"
		lappend fields "compar-$analysis=\"\${sr}\$seq-$analysis\""
		lappend fields [subst {pred-$analysis=if(\$compar-$analysis in "vv","TP",\$compar-$analysis in "rv","FP",\$compar-$analysis in "vr","rFN",\$compar-$analysis in "vu","uFN","u")}]
		lappend groups pred-$analysis * curve-$analysis *
	}
	if {$regionfile ne ""} {
		exec cg regselect $comparfile $regionfile | cg select -f $fields -g $groups > $resultfileall.temp
	} else {
		cg select -f $fields -g $groups $comparfile > $resultfileall.temp
	}
	file rename -force -- $resultfileall.temp $resultfileall
	#
	set o [open $resultfile.temp w]
	puts $o [join {type analysis cutoff fscore precision recall nou_recall TP FP rFN uFN u} \t]
	set types [list_subindex [split [cg select -g type $resultfileall] \n] 0]
	set types [list_remove $types type sub]
	foreach type $types {
		cg select -overwrite 1 -q "\$type eq \"$type\"" $resultfileall $resultfileall.temp
		foreach analysis $analyses {
			foreach cutoff $refcurve_cutoffs {
				if {$cutoff eq ""} {
					set temp [cg select -g pred-$analysis -gc sum(count) $resultfileall.temp]
				} elseif {[regexp ^t $cutoff]} {
					# if we would also want to apply the cutoff to the reference
					set ucutoff [string range $cutoff 1 end]
					set temp [cg select -f [list \
						[subst {apred=if(\$refcurve < $ucutoff, "u", \$curve-$analysis >= $ucutoff,\$pred-$analysis,\$pred-$analysis in "rFN TP uFN","uFN","u")
					}]] -g apred -gc sum(count) $resultfileall.temp]
				} else {
					set temp [cg select -f [list \
						[subst {apred=if(\$curve-$analysis >= $cutoff,\$pred-$analysis,\$pred-$analysis in "rFN TP uFN","uFN","u")
					}]] -g apred -gc sum(count) $resultfileall.temp]
				}
				unset -nocomplain a
				foreach key {TP FP rFN uFN u} {
					set a($key) 0
				}
				foreach line [split $temp \n] {
					set a([lindex $line 0]) [lindex $line 1]
				}
				foreach {TP FP rFN uFN u} [list $a(TP) $a(FP) $a(rFN) $a(uFN) $a(u)] break
				if {$TP == 0 && $FP == 0} {
					set precision ?
				} else {
					set precision [format %.3f [expr {100.0*$TP/($TP+$FP)}]]
				}
				if {$TP == 0 && $rFN == 0 && $uFN == 0} {
					set recall ?
				} else {
					set recall [format %.3f [expr {100.0*$TP/($TP+$rFN+$uFN)}]]
				}
				if {$TP == 0 && $rFN == 0} {
					set nourecall ?
				} else {
					set nourecall [format %.3f [expr {100.0*$TP/($TP+$rFN)}]]
				}
				if {[catch {
					set fscore [expr {2*$precision*$recall/($precision+$recall)}]
				}]} {
					set fscore ?
				}
				puts $o [join [list $type $analysis $cutoff $fscore $precision $recall $nourecall $TP $FP $rFN $uFN $u] \t]
			}
		}
	}
	close $o
	file rename -force -- $resultfile.temp $resultfile
	file delete $resultfileall.temp
}


#!/bin/sh
# the next line restarts using wish \
exec cg source "$0" ${1+"$@"}

proc cg_compound {args} {
	set per {}
	set persample 0
	set peranalysis 0
	set geneset refGene
	set criteria {}
	set impact {>=CDSMIS}
	set limitsamples {}
	set resultfile {}
	cg_options compound args {
		-per {
			if {$value ni "sample analysis"} {error "unknown value $value for -per, must be one of: sample analysis"}
			set per $value
		}
		-criteria {
			set criteria $value
		}
		-geneset {
			set geneset $value			
		}
		-impact {
			set impact [string trim $value]
		}
		-samples {
			set persample 1
			set limitsamples [list_remdup $value]
		}
		-analyses {
			set peranalysis 1
			set limitsamples [list_remdup $value]
		}
		-dbdir {
			set dbdir $value
		}
	} {varfile resultfile} 1 2
	if {$per eq ""} {
		if {$persample} {
			if {$peranalysis} {
				error "cannot give both -samples and analyses option"
			}
			set per sample
		} else {
			set per analysis
		}
	}
	set varfile [file_absolute $varfile]
	if {$resultfile eq ""} {
		set resultfile [file root $varfile]-compound.tsv
	}
	set resultfile [file_absolute $resultfile]
	set resultgenelist [file root [gzroot $resultfile]]-genelist.tsv
	set tempfile [tempfile]
	if {$criteria eq ""} {
		if {$per eq "analysis"} {
			set criteria {$sequenced == "v"}
		} else {
			error "for analysis per sample, -criteria must be specified"
		}
	}
	# putsvars tempfile per criteria geneset impact limitsamples varfile resultfile resultgenelist

	#
	# make tempfile with hqv, and query for compound to $resultgenelist
	#
	# make tempfile with hqv
	set f [gzopen $varfile]
	set header [tsv_open $f]
	gzclose $f
	set fields [list *]
	set neededfields {}
	set tokens [tsv_select_tokenize $header $criteria neededfields]
	if {$per eq "sample"} {
		set samples [listsamples $header]
		if {[llength $limitsamples]} {
			set samples [list_common $limitsamples $samples]
			if {[llength $samples] != [llength $limitsamples]} {
				error "some samples given in -samples are not in the file $varfile: [list_remove $limitsamples $samples]"
			}
		}
		foreach sample $samples {
			set temp [tsv_select_replacevars $tokens $header $sample]
			set code [tsv_select_saggr_detokenize $temp $header neededfields missing]
			lappend fields "hqv-$sample=if($code,1,0)"
		}
		lappend fields "transcripts=transcripts(\"$geneset\",\"$impact\")"
		cg select -overwrite 1 -f $fields -q {
			scount($hqv == 1) > 0
		} $varfile $tempfile
	} else {
		set samples [listanalyses $header]
		if {[llength $limitsamples]} {
			set samples [list_common $limitsamples $samples]
			if {[llength $samples] != [llength $limitsamples]} {
				error "some analyses given in -analyses are not in the file $varfile: [list_remove $limitsamples $samples]"
			}
		}
		foreach sample $samples {
			set temp [tsv_select_replacevars $tokens $header $sample]
			set code [tsv_select_saggr_detokenize $temp $header neededfields missing]
			lappend fields "hqv-$sample=if($code,1,0)"
		}
		lappend fields "transcripts=transcripts(\"$geneset\",\"$impact\")"
		# cg select -overwrite 1 -f $fields $varfile $tempfile
		cg select -overwrite 1 -f $fields -q {
			acount($hqv == 1) > 0
		} $varfile $tempfile
	}
#	set fields {}
#	lappend fields "hqv-*=if($criteria,1,0)"
#	lappend fields "transcripts=transcripts(\"$geneset\",\"$impact\")"
#	# make tempfile with hqv
#	if {$per eq "sample"} {
#		cg select -overwrite 1 -f [list * {*}$fields] -q {
#			scount($hqv == 1) > 0
#		} $varfile $tempfile
#	} else {
#		cg select -overwrite 1 -f [list * {*}$fields] -q {
#			acount($hqv == 1) > 0
#		} $varfile $tempfile
#	}
	# look for compound vars
	exec cg select -g {
		analysis * -transcripts
	} -gc {
		hqv 1 count
	} $tempfile | cg select -sh /dev/null -q {$count-1 > 1} {*}[compresspipe $resultgenelist] > $resultgenelist.temp
	file rename -force -- $resultgenelist.temp $resultgenelist

	#
	# parse $resultgenelist for info on coumpounds
	#
	set temp [exec {*}[gzcat $resultgenelist] $resultgenelist]
	unset -nocomplain a
	unset -nocomplain transcripta
	unset -nocomplain samplea
	foreach line [split [string trim $temp] \n] {
		foreach {sample transcript num} [split $line \t] break
		set samplea($sample) 1
		set transcripta($transcript) 1
		set a($transcript,$sample) $num
	}
	# set samples [bsort [array names samplea]]
	set transcripts [bsort [array names transcripta]]

	#
	# write $resultfile
	# go over tempfile with hqv, and add compound info cols
	#
	catch {close $f} ; catch {close $o}
	set f [gzopen $tempfile]
	set header [tsv_open $f]
	if {$per eq "sample"} {
		set samples [listsamples $header]
	} else {
		set samples [listanalyses $header]
	}
	if {[llength $limitsamples]} {
		set samples [list_common $limitsamples $samples]
	}
	set poss {}
	foreach sample $samples {
		lappend poss [lsearch $header hqv-$sample]
	}
	# set tpos [lsearch $header transcripts]
	set tempresultfile $resultfile.temp[gzext $resultfile]
	set o [wgzopen $tempresultfile]
	set newheader $header
	foreach sample $samples {
		lappend newheader compound-$sample
	}
	puts $o [join $newheader \t]
	while 1 {
		if {[gets $f line] == -1} break
		set line [split $line \t]
		set ltranscripts {}
		set temp [list_pop line]
		foreach t [split $temp {,; }] {
			if {[info exists transcripta($t)]} {
				lappend ltranscripts $t
			}
		}
		lappend line $ltranscripts
		set found 0
		foreach hqv [list_sub $line $poss] sample $samples {
			set compound {}
			if {$hqv} {
				foreach t $ltranscripts {
					if {[info exists a($t,$sample)]} {
						set compound $a($t,$sample)
						set found 1
						break
					}
				}
			}
			lappend line $compound
		}
		if {$found} {
			puts $o [join $line \t]
		}

	}
	close $o
	close $f
	file rename -force -- $tempresultfile $resultfile

}



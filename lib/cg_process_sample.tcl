#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg {args} {
	# puts "cg $args"
	eval exec cg $args 2>@ stderr
}

proc process_sample {dir destdir dbdir {force 0}} {
	set keepdir [pwd]
	set dir [file normalize $dir]
	set destdir [file normalize $destdir]
	set dbdir [file normalize $dbdir]
	file mkdir $destdir
	cd $destdir
	set name [file tail $destdir]
	putslog "Processing sample $dir -> $destdir"
	# sort files
	if {[file exists $destdir/variants-$name.tsv]} {
		# prepared dir, not coming from CG
		if {![file exists $destdir/annotvar-$name.tsv]} {
			file link -symbolic $destdir/annotvar-$name.tsv $destdir/variants-$name.tsv
		}
	} else {
		if {$force || ![file exists [gzfile svar-$name.tsv]] || ![file exists $destdir/info.txt]} {
			set varfile [glob $dir/ASM/var-*-ASM*.tsv*]
			if {[llength $varfile] != 1} {error "could not identify varfile"}
			if {[llength $varfile] > 1} {error "could not identify varfile"}
			set f [gzopen $varfile]
			set info {}
			while {![eof $f]} {
				set line [gets $f]
				if {[string index $line 0] ne "#"} break
				lappend info $line
			}
			catch {close $f}
			if {[file exists $destdir/info.txt]} {
				set test [split [file_read $destdir/info.txt] \n]
				if {$info ne $test} {
					error "$destdir already has info.txt that differs from data in the source $dir"
				}
			}
			file_write $destdir/info.txt [join $info \n]
			putslog "Sort var file ($varfile)"
			cg select -s "chromosome begin end varType" $varfile svar-$name.tsv.temp
			file rename -force svar-$name.tsv.temp svar-$name.tsv
		}
		if {$force || ![file exists [gzfile sgene-$name.tsv]]} {
			set genefile [list_lremove [glob $dir/ASM/gene-*-ASM*.tsv*] [glob -nocomplain $dir/ASM/gene-var-summary-*-ASM*.tsv*]]
			if {[llength $genefile] != 1} {error "could not identify genefile"}
			putslog "Sort gene file ($genefile)"
			cg select -s "chromosome begin end" $genefile sgene-$name.tsv.temp
			file rename -force sgene-$name.tsv.temp sgene-$name.tsv
		}
		# annotated vars file
		if {$force || ![file exists [gzfile annotvar-$name.tsv]]} {
			putslog "Create annotated varfile annotvar-$name.tsv"
			# set file svar-$name.tsv; set genefile sgene-$name.tsv; set outfile temp.tsv
			cg var2annot svar-$name.tsv sgene-$name.tsv annotvar-$name.tsv.temp
			file rename -force annotvar-$name.tsv.temp annotvar-$name.tsv
		}
		if {$force || ![file exists [gzfile sreg-$name.tsv]]} {
			set regfile [glob -nocomplain $dir/ASM/reg-*-ASM*.tsv*]
			if {[file exists $regfile]} {
				putslog "Sort region file ($regfile)"
				cg select -s "chromosome begin end" $regfile sreg-$name.tsv.temp
				file rename -force sreg-$name.tsv.temp sreg-$name.tsv
			} else {
				putslog "Extract sreg-$name.tsv from svar-$name.tsv"
				cg select -q {$varType != "no-call" && $varType != "no-ref"} -f "chromosome begin end" svar-$name.tsv sreg-$name.tsv.temp
				cg regjoin sreg-$name.tsv.temp > temp2.tsv
				file rename -force temp2.tsv sreg-$name.tsv
			}
		}
		# sample specific filters
		if {$force || ![file exists [gzfile reg_refcons-$name.tsv]]} {
			putslog "Find refcons regions for var-$name.tsv"
			cg refconsregions svar-$name.tsv > reg_refcons-$name.tsv.temp
			file rename -force reg_refcons-$name.tsv.temp reg_refcons-$name.tsv
		}
		if {$force || ![file exists [gzfile reg_nocall-$name.tsv]]} {
			putslog "Find partial no-call regions for var-$name.tsv"
			if {[catch {
				nocallregions svar-$name.tsv reg_nocall-$name.tsv.temp
			}]} {
				puts stderr "Could not make reg_nocall-$name.tsv (old version files ?)"
			} else {
				file rename -force reg_nocall-$name.tsv.temp reg_nocall-$name.tsv
			}
		}
	}
	# make coverage files
	file mkdir coverage
	set files [lsort -dict [glob -nocomplain $dir/ASM/REF/coverage*]]
	foreach file $files {
		set chr [lindex [split $file -] 1]
		regsub ^chr $chr {} chr
		set header [cg select -h $file]
		foreach posfield {offset pos} {
			if {[lsearch $header $posfield] != -1}  break
		}
		if {$posfield == -1} {
			exiterror "No position/offset field found in $file"
		}
		foreach covfield {uniqueSequenceCoverage coverage} {
			if {[lsearch $header $covfield] != -1}  break
		}
		if {$covfield == -1} {
			exiterror "No coverage/uniqueSequenceCoverage field found in $file"
		}
		set other [list_remove $header $posfield $covfield]
		foreach field $other {
			set base coverage/$field-$name-$chr
			if {![file exists $base.bcol]} {
				putslog "Making $base.bcol"
				if {[catch {
					exec [catprog $file] $file | cg bcol make -p $posfield -t s -n -1 $base $field <  $file
				} e]} {
					exec [catprog $file] $file | cg bcol make -p $posfield -t i -n -1 $base $field <  $file
				}
				cg razip $base.bcol.bin
			}
		}
		set base coverage/coverage-$name-$chr
		if {![file exists $base.bcol]} {
			putslog "Making $base.bcol"
			exec [catprog $file] $file | cg bcol make -p $posfield -t su $base $covfield
			cg razip $base.bcol.bin
		}
	}
	# copy extra info
	if {$force || ![file exists CNV]} {
		putslog "Copying CNV"
		file delete -force CNV.temp
		catch {
			file copy $dir/ASM/CNV CNV.temp
			file rename CNV.temp CNV
		} result
		putslog $result
	}
	if {$force || ![file exists SV]} {
		putslog "Copying SV"
		file delete -force SV.temp
		catch {
			file copy $dir/ASM/SV SV.temp
			file rename SV.temp SV
		} result
		putslog $result
	}
	set dstfile cgsv-$name.tsv
	if {$force || ![file exists $dstfile]} {
		putslog "Convert SV: making $dstfile"
		set list [glob -nocomplain SV/allJunctionsBeta-*.tsv* SV/annotatedJunctionsAlpha-*.tsv*]
		if {[llength $list]} {
			set srcfile [lindex $list 0]
			cg convcgsv $srcfile $dstfile
		}
	}
	set dstfile cgcnv-$name.tsv
	if {$force || ![file exists $dstfile]} {
		putslog "Convert SV: making $dstfile"
		set list [glob -nocomplain CNV/cnvSegmentsBeta-*.tsv* CNV/cnvSegmentsAlpha-*.tsv*]
		if {[llength $list]} {
			set srcfile [lindex $list 0]
			cg convcgcnv $srcfile $dstfile
		}
	}
#	if {$force || ![file exists EVIDENCE]} {
#		putslog "Copying EVIDENCE"
#		file delete -force EVIDENCE.temp
#		catch {
#			file copy $dir/ASM/EVIDENCE .
#			file rename EVIDENCE.temp EVIDENCE
#		} result
#		putslog $result
#	}
	if {$force || ![file exists [gzfile reg_cluster-$name.tsv]]} {
		putslog "Find cluster regions for annotvar-$name.tsv"
		cg clusterregions < annotvar-$name.tsv > reg_cluster-$name.tsv.temp
		file rename -force reg_cluster-$name.tsv.temp reg_cluster-$name.tsv
	}
	if {$force || ![file exists [gzfile reg_ns-$name.tsv]]} {
		putslog "Find regions with N's for annotvar-$name.tsv"
		cg select -f {chromosome begin end} -q {$alleleSeq1 ~ /[N?]/ || $alleleSeq2 ~ /[N?]/} < annotvar-$name.tsv > reg_ns-$name.tsv.temp
		file rename -force reg_ns-$name.tsv.temp reg_ns-$name.tsv
	}
	if {$force || ![file exists [gzfile reg_lowscore-$name.tsv]]} {
		set header [cg select -h annotvar-$name.tsv]
		if {[llength [list_common $header {totalScore1 totalScore2}]] == 2} {
			putslog "Find regions with lowscores for annotvar-$name.tsv"
			cg select -f {chromosome begin end} -q {$totalScore1 < 60 || $totalScore2 < 60} < annotvar-$name.tsv > reg_lowscore-$name.tsv.temp
			file rename -force reg_lowscore-$name.tsv.temp reg_lowscore-$name.tsv
		}
	}
	if {$force || ![file exists [gzfile fannotvar-$name.tsv]]} {
		# add filterdata to annotvar
		set todo {}
		if {[file exists reg_refcons-$name.tsv]} {
			lappend todo [list refcons rc reg_refcons-$name.tsv]
		}
		# lappend todo [list nocall nc reg_nocall-$name.tsv]
		lappend todo [list cluster cl reg_cluster-$name.tsv]
		# lappend todo [list trf trf $dbdir/regdb-simple_repeats.tsv]
		# lappend todo [list str str $dbdir/regdb-microsatelite.tsv]
		# lappend todo [list segdup sd $dbdir/regdb-segdups.tsv]
		# lappend todo [list selfchain sc $dbdir/regdb-selfchain.tsv]
		# lappend todo [list repeat rp $dbdir/regdb-repeatmasker.tsv]
		# lappend todo [list rna rna $dbdir/regdb-rnagenes.tsv]
		# lappend todo [list more5pct m5 $dbdir/regdb-1000genomesmore5pct.tsv]
		# lappend todo [list more1pct m1 $dbdir/regdb-1000genomesmore1pct.tsv]
		# foreach file [lsort -dictionary [glob -nocomplain $dbdir/checked*.tsv]] {
		# 	set value [lindex [split [file root [file tail $file]] _] 1]
		# 	if {$value eq ""} {set value checked}
		# 	lappend todo [list checked $value $file]
		# }
		annot_annotvar annotvar-$name.tsv fannotvar-$name.tsv $todo $destdir
	}
	# coverage
	if {$force || ![file exists reg-$name.covered]} {
		putslog "Coverage of sequenced regions"
		cg covered sreg-$name.tsv > reg-$name.covered.temp
		file rename -force reg-$name.covered.temp reg-$name.covered
	}
	if {$force || ![file exists filteredrefcons-$name.covered] && [file exists reg_refcons-$name.tsv]} {
		putslog "Coverage of refcons region"
		cg regsubtract sreg-$name.tsv reg_refcons-$name.tsv > filteredrefcons-$name.tsv.temp
		file rename -force filteredrefcons-$name.tsv.temp filteredrefcons-$name.tsv
		cg covered filteredrefcons-$name.tsv > filteredrefcons-$name.covered.temp
		file rename -force filteredrefcons-$name.covered.temp filteredrefcons-$name.covered
	}
	if {$force || ![file exists [gzfile filteredns-$name.tsv]]} {
		putslog "Coverage of ns region"
		cg regsubtract sreg-$name.tsv reg_ns-$name.tsv > filteredns-$name.tsv.temp
		file rename -force filteredns-$name.tsv.temp filteredns-$name.tsv
	}
	if {$force || ![file exists filteredns-$name.covered]} {
		putslog "Making filteredns-$name.covered"
		cg covered filteredns-$name.tsv > filteredns-$name.covered.temp
		file rename -force filteredns-$name.covered.temp filteredns-$name.covered
	}
	if {$force || ![file exists [gzfile filteredlowscore-$name.tsv]] && [file exists reg_lowscore-$name.tsv]} {
		putslog "Coverage of lowscore region"
		cg regsubtract sreg-$name.tsv reg_lowscore-$name.tsv > filteredlowscore-$name.tsv.temp
		file rename -force filteredlowscore-$name.tsv.temp filteredlowscore-$name.tsv
	}
	if {$force || ![file exists filteredlowscore-$name.covered] && [file exists filteredlowscore-$name.tsv]} {
		putslog "Making filteredlowscore-$name.covered"
		cg covered filteredlowscore-$name.tsv > filteredlowscore-$name.covered.temp
		file rename -force filteredlowscore-$name.covered.temp filteredlowscore-$name.covered
	}
	if {$force || ![file exists [gzfile histo-refcons-$name.tsv]] && [file exists reg_refcons-$name.tsv]} {
		putslog "Making histo-refcons-$name.covered"
		cg reghisto reg_refcons-$name.tsv > histo-refcons-$name.tsv.temp
		file rename -force histo-refcons-$name.tsv.temp histo-refcons-$name.tsv
	}
	if {$force || ![file exists [gzfile filteredcluster-$name.tsv]]} {
		putslog "Coverage of clusters region"
		cg regsubtract sreg-$name.tsv reg_cluster-$name.tsv > filteredcluster-$name.tsv.temp
		file rename -force filteredcluster-$name.tsv.temp filteredcluster-$name.tsv
	}
	if {$force || ![file exists filteredcluster-$name.covered]} {
		putslog "Making filteredcluster-$name.covered"
		cg covered filteredcluster-$name.tsv > filteredcluster-$name.covered.temp
		file rename -force filteredcluster-$name.covered.temp filteredcluster-$name.covered
	}
#	if {$force || ![file exists [gzfile reg_below20-$name.tsv]]} {
#		putslog "Make region file reg_below20-$name.tsv"
#		set files [gzfiles coverage/coverageRefScore-*.tsv coverage/coverage-*.bcol]
#		if {[llength $files] < 24} {
#			putslog "WARNING: only [llength $files] coverage files found"
#		}
#		file delete reg_below20-$name.tsv.temp
#		cg regextract -above 0 20 {*}$files > reg_below20-$name.tsv.temp
#		file rename reg_below20-$name.tsv.temp reg_below20-$name.tsv
#	}
#	if {$force || ![file exists reg_below20-$name.covered]} {
#		putslog "Make reg_below20-$name.covered"
#		cg covered [gzfile reg_below20-$name.tsv] > reg_below20-$name.covered.temp
#		file rename -force reg_below20-$name.covered.temp reg_below20-$name.covered
#	}
#	if {$force || ![file exists reg_below20-$name.tsv.gz] || [file exists reg_below20-$name.tsv]} {
#		catch {file delete reg_below20-$name.tsv.gz}
#		exec bgzip reg_below20-$name.tsv
#	}
#	if {$force || ![file exists [gzfile reg_above100-$name.tsv]]} {
#		putslog "Make region file reg_above100-$name.tsv"
#		set files [gzfiles coverage/coverageRefScore-*.tsv]
#		if {[llength $files] < 24} {
#			putslog "WARNING: only [llength $files] coverage files found"
#		}
#		file delete reg_above100-$name.tsv.temp
#		cg regextract -above 1 100 {*}$files > reg_above100-$name.tsv.temp
#		file rename -force reg_above100-$name.tsv.temp reg_above100-$name.tsv
#	}
#	if {$force || ![file exists reg_above100-$name.covered]} {
#		putslog "Make reg_above100-$name.covered"
#		cg covered [gzfile reg_above100-$name.tsv] > reg_above100-$name.covered.temp
#		file rename reg_above100-$name.covered.temp reg_above100-$name.covered
#	}
#	if {$force || ![file exists reg_above100-$name.tsv.gz] || [file exists reg_above100-$name.tsv]} {
#		catch {file delete reg_above100-$name.tsv.gz}
#		exec bgzip reg_above100-$name.tsv
#	}
	file_write FINISHED ""
	putslog "Finished $destdir"
	cd $keepdir
}

proc process_compare_checkfield {file field} {
	set f [open $file]
	set header [tsv_open $f]
	close $f
	inlist $header $field
}

proc process_indexcompress {file} {
	set ext [file extension $file]
	if {$ext eq ".gz"} {
		gunzip $file
		set file [file root $file]
	}
	set f [gzopen $file]
	set header [tsv_open $f]
	foreach field {offset end1} {
		set fpos [lsearch $header $field]
		if {$fpos != -1} break
	}
	if {$fpos == -1} {error "no column offset, end1 in file $file"}
	if {([gets $f] eq "") && [eof $f]} return
	close $f
	if {![file exists $file.${field}_index]} {
		tsv_index $field $file
	}
	putslog "Compressing $file"
	if {![inlist {.rz .bgz} $ext]} {
		exec bgzip -c $file > $file.gz.temp
		file rename -force $file.gz.temp $file.gz
		file delete $file
	}
}

proc cg_process_indexcompress {args} {
	global scriptname action
	if {[llength $args] != 1} {
		puts stderr "format is: $scriptname $action file"
		puts stderr " - makes index, and compresses to bgzip"
		exit 1
	}
	foreach {file} $args break
	process_indexcompress $file
}

proc cg_process_sample {args} {
	if {([llength $args] < 3) || ([llength $args] > 4)} {
		errorformat process_sample
		exit 1
	}
	process_sample {*}$args
}

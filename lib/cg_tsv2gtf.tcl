proc cg_tsv2gtf {args} {
	set file {}
	set outfile {}
	set upstream 0
	set nocds 0
	set genecol {}
	set addgene 0
	cg_options tsv2gtf args {
		-upstream {set upstream $value}
		-nocds {set nocds $value}
		-genecol {set genecol $value}
		-addgene {set addgene $value}
	}  {file outfile} 0 2
	if {[llength $args] > 2} {
		errorformat tsv2gtf
	}
	catch {sclose $f} ; catch {sclose $o}
	if {$file eq ""} {
		set f stdin
	} else {
		set f [gzopen $file]
	}
	if {$outfile eq ""} {
		set o stdout
	} else {
		set o [open $outfile.temp w]
	}
	set header [open_genefile $f dposs $genecol]
	set regposs [tsv_basicfields $header 3]
	set sourcepos [lsearch $header source]
	set scorepos [lsearch $header score]
	set rest [list_sub $header -exclude $dposs]
	set restposs [list_cor $header $rest]
	set gchr {}
	unset -nocomplain genea
	while {![eof $f]} {
		set attr {}
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		if {$sourcepos == -1} {
			set source genepred2tsv
		} else {
			set source [lindex $line $sourcepos]
		}
		if {$scorepos == -1} {
			set score .
		} else {
			set score [lindex $line $scorepos]
		}
		set geneobj [annotatevar_gene_makegeneobj {} $line $dposs 2000 $nocds]
		set ftlist [dict get $geneobj ftlist]
		foreach {chr begin end} [list_sub $line $regposs] break
		set gene [dict get $geneobj genename]
		set transcript [dict get $geneobj transcriptname]
		set compl [dict get $geneobj complement]
		if {!$compl} {set strand +} else {set strand -}
		set rest [list_sub $line $restposs]
		set attr "gene_id \"$gene\"\; transcript_id \"$transcript\"\;"
		#
		if {$addgene} {
			if {![info exists genea($gene)]} {
				set genea($gene) [list $chr $source gene [expr {$begin+1}] $end . $strand . "gene_id \"$gene\"\;"]
			} else {
				set gend [lindex $genea($gene) 4]
				if {$end > $gend} {lset genea($gene) 4 $end}
			}
		}
		#
		puts $o [join [list $chr $source transcript [expr {$begin+1}] $end $score $strand . $attr] \t]
		if {$upstream} {
			set dline [lindex $ftlist 0]
			set begin [expr {[lindex $dline 1]-$upstream+1}]
			set end [expr {[lindex $dline 1]+1}]
			puts $o [join [list $chr $source upstream $begin $end $score $strand . $attr] \t]
		}
		set list [lrange $ftlist 1 end]
		set types [list_subindex $list 2]
		set poss [list_find -regexp $types {^UTR|CDS|RNA$}]
		set list [list_sub $list $poss]
		set types [list_subindex $list 2]
		set noexon 0
		set cdsposs [list_find $types CDS]
		set remove {}
		if {$compl} {
			set adjstartcds [lindex $cdsposs 0]
			set startfix 3
			if {$adjstartcds ne ""} {
				set temp [lindex $list $adjstartcds]
				foreach {begin end} $temp break
				set size [expr {$end - $begin}]
				if {$size <= 3} {
					set remove $adjstartcds
					set adjstartcds [lindex $cdsposs 1]
					set startfix [expr {3-$size}]
				}
			}
			#
			set adjendcds {}
			set endfix 0
		} else {
			set adjstartcds {}
			set startfix 0
			#
			set adjendcds [lindex $cdsposs end]
			set endfix -3
			if {$adjendcds ne ""} {
				set temp [lindex $list $adjendcds]
				foreach {begin end} $temp break
				set size [expr {$end - $begin + 1}]
				if {$size <= 3} {
					set remove $adjendcds
					set adjendcds [lindex $cdsposs end-1]
					set endfix [expr {$size - 3}]
				}
			}
		}
		set num 0
		set framespos 0
		foreach dline $list ndline [lrange $list 1 end] {
			foreach {begin end type element rna_start rna_end protein_start protein_end} $dline break
			incr begin
			incr end
			set nbegin [lindex $ndline 0]
			set ntype [lindex $ndline 0]
			if {$type eq $ntype || $nbegin != $end} {
				if {$noexon} {
					set noexon 0
				} else {
					puts $o [join [list $chr $source exon $begin $end $score $strand . $attr] \t]
				}
			} else {
				set nend [lindex $ndline 1]
				incr nend
				puts $o [join [list $chr $source exon $begin $nend $score $strand . $attr] \t]
				set noexon 1
			}
			if {$type eq "CDS"} {
				if {$num == $adjstartcds} {
					incr begin $startfix
				}
				if {$num == $adjendcds} {
					incr end $endfix
				}
				if {$num != $remove} {
					set frame [expr {3-$protein_start%3}]
					if {$frame == 3} {set frame 0}
					puts $o [join [list $chr $source $type $begin $end $score $strand $frame $attr] \t]
				}
			}
			incr num
		}
		if {$upstream} {
			set dline [lindex $ftlist end]
			set end [lindex $dline 0]
			set begin [expr {$end+$upstream+1}]
			puts $o [join [list $chr $source downstream $begin $end $score $strand . $attr] \t]
		}
	}
	if {$addgene} {
		foreach gene [array names genea] {
			puts $o [join $genea($gene) \t]
		}
	}
	if {$o ne "stdout"} {
		close $o
		file rename -force -- $outfile.temp $outfile
	}
	if {$f ne "stdout"} {gzclose $f}
}

if 0 {
	set file tmp/genetmp.tsv
	set outfile tmp/out.tsv
	set args [list $file $outfile]
}


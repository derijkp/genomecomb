proc cg_bamreorder {args} {
	set method genomecomb
	set threads 1
	cg_options bamreorder args {
		-method {
			if {$value ni {picard genomecomb}} {error "bamreorder: unsupported -method $value"}
			set method $value
		}
		-threads {
			set threads $value
		}
	} {src dest ref}
	if {![file exists $ref.fai]} {
		exec samtools faidx $ref
	}
	set index [csv_parse [file_read $ref.fai] \t]
	set refchrs [list_subindex $index 0]
	set reflens [list_subindex $index 3]
	# check bam file to see if we need to rename some contigs
	set bamheader [exec samtools view -H $src]
	if {$method eq "picard"} {
		set bamchrs {}
		set bamlens {}
		foreach line [split $bamheader \n] {
			if {[regexp {SN:([^ ]+)\tLN:([0-9]+)} $line temp chr len]} {
				lappend bamchrs $chr
				lappend bamlens $len
			}
		}
		# putsvars bamchrs bamlens
		set torename {}
		foreach chr $bamchrs len $bamlens {
			set pos [lsearch $refchrs $chr]
			set newchr $chr
			if {$pos == -1} {
				if {[regexp ^chr $chr]} {
					set newchr [string range $chr 3 end]
					set pos [lsearch $refchrs $newchr]
				} else {
					set newchr chr$chr
					set pos [lsearch $refchrs $newchr]
				}
			}
			if {$pos == -1} {
				puts stderr "WARNING: contig $chr was not found in ref"
			} else {
				set reflen [lindex $reflens $pos]
				if {$len != $reflen} {
					puts stderr "WARNING: dropping contig $chr: has different size from $newchr in ref"
					set newchr _rem_$chr
				}
				if {$newchr ne $chr} {
					lappend torename SN:$chr SN:$newchr
				}
			}
		}
		if {[llength $torename]} {
			set bamheader [string_change $bamheader $torename]
			set tempfile [tempfile]
			file_write $tempfile $bamheader
			set tempsrc $dest.temp
			exec samtools reheader $tempfile $src > $tempsrc
			exec samtools index $tempsrc
		} else {
			set tempsrc $src
		}
		if {[file extension $ref] ne ".fa"} {
			set nref [tempfile].fa
			mklink $ref $nref
			mklink $ref.fai $nref.fai
			mklink $ref.dict $nref.dict
		} else {
			set nref $ref
		}
		if {![file exists $nref.dict]} {
			picard CreateSequenceDictionary R=$nref O=$nref.dict
		}
		picard ReorderSam R=$nref I=$tempsrc O=$dest.temp2 ALLOW_INCOMPLETE_DICT_CONCORDANCE=true VALIDATION_STRINGENCY=LENIENT 2>@ stderr >@ stdout
		if {$nref ne $ref} {file delete $nref}
		exec samtools index $dest.temp2
		file delete $dest.temp  $dest.temp.[indexext $dest]
		file rename $dest.temp2 $dest
		file rename $dest.temp2.[indexext $dest] $dest.[indexext $dest]
	} else {
		# new bamheader
		unset -nocomplain chrsa
		set newheader {}
		unset -nocomplain post
		foreach line [split $bamheader \n] {
			if {[regexp ^@SQ $line]} {
				if {![regexp {SN:([^ \t]+)} $line temp chr]} {
					error "@SQ line without SN: $line"
				}
				set chrsa($chr) $line
				set post {}
			} elseif {[info exists post]} {
				lappend post $line
			} else {
				lappend newheader $line
			}
		}
		set todo {}
		foreach chr $refchrs reflen $reflens {
			set cchr [chr_clip $chr]
			if {[info exists chrsa($cchr)]} {
				set bamchr $cchr
			} elseif {[info exists chrsa(chr$cchr)]} {
				set bamchr chr$cchr
			} else {
				puts stderr "WARNING: reference contig $chr was not present in bam"
				continue
			}
			set line $chrsa($bamchr)
			if {[regexp {LN:([0-9]+)} $line temp bamlen]} {
				if {$bamlen != $reflen} {
					puts stderr "WARNING: dropping contig $bamchr: has different size from $chr in ref"
					continue
				}
			}
			set line [string_change $line [list SN:$bamchr SN:$chr]]
			lappend newheader $line
			lappend todo $bamchr
			unset chrsa($bamchr)
		}
		set left [array names chrsa]
		if {[llength $left]} {
			puts stderr "WARNING: dropping bam contigs [join $left ,] : not found in ref"
			
		}
		lappend newheader {*}$post
		set o [open [list | samtools view -b - > $dest.temp] w]
		puts $o [join $newheader \n]
		foreach chr $todo {
			set f [open [list | samtools view -@ $threads $src $chr]]
			fcopy $f $o
			close $f
		}
		close $o
		exec samtools index $dest.temp
		file rename -force $dest.temp $dest
		file rename -force $dest.temp.[indexext $dest] $dest.[indexext $dest]
	}
}

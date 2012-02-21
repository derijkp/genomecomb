proc multiexists {pattern {chrs {0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y}}} {
	foreach chr $chrs {
		if {![file exists [subst $pattern]]} {
			return 0
		}
	}
	return 1
}

proc cg_map2bam {readfile mapfile reffile outfile} {
	set name [file tail $outfile]
	set chrs {0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y}
	if {[multiexists $outfile-chr\$chr.bam $chrs]} {
		putslog "skipping files $outfile: all bams exist already"
		return
	}
	set scratchdir [scratchdir]
	set scratchbase $scratchdir/temp_[file tail $outfile]
	set samfile $scratchbase.sam
	putslog "Making sam from $mapfile"
	if {[catch {
		exec cgatools map2sam -r $readfile -m $mapfile -s $reffile --mate-sv-candidates -o $samfile.temp
	} error]} {
		regsub -all {IND_[^=]+=[0-9]+} $error {} error
		set error [string trim $error]
		if {$error ne {}} {error $error}
	}
	file rename $samfile.temp $samfile
	if {[catch {
		set header [exec samtools view -S -H $samfile 2> /dev/null]
	} error]} {
		if {$error ne {[samopen] SAM header is present: 25 sequences.}} {error $error}
	}
	set header [split $header \n]
	set len [llength $header]
	set poss [list_find -glob $header @SQ*]
	set list [list_sub $header $poss]
	set list [lsort -dict -index 1 $list]
	set pos [lindex $poss 0]
	set header [lreplace [list_sub $header -exclude $poss] $pos -1 {*}$list]
	putslog "Distributing sam to chromosomes for $name"
	exec tail -n +[expr {($len+1)}] $samfile | distrsam $scratchbase
	file delete $samfile
	putslog "Making and sorting sam per chromosome for $name"
	foreach chr $chrs {
		file_write $scratchbase-chr$chr.sam [join $header \n]\n
		exec gnusort8 -T $scratchdir -t \t -N -s -k3 -k4 $scratchbase-$chr >> $scratchbase-chr$chr.sam
		file delete $scratchbase-$chr
	}
	putslog "Making bams for $name"
	foreach chr $chrs {
		if {[catch {
			exec samtools view -S -h -b -o $outfile-chr$chr.bam.temp $scratchbase-chr$chr.sam
		} error]} {
			if {$error ne {[samopen] SAM header is present: 25 sequences.}} {error $error}
		}
		file delete $scratchbase-chr$chr.sam
		catch {exec samtools index $scratchbase-chr$chr.sam}
		file rename -force $outfile-chr$chr.bam.temp $outfile-chr$chr.bam
	}
}

proc cg_mergebam {outfile args} {
	if {![llength $args] > 0} {
		puts stderr "format is: cg mergebam outfile bamfile1 ?bamfile2? ..."
		exit 1
	}
	putslog "Merging $outfile"
	exec samtools merge $outfile.temp {*}$args
	file rename -force $outfile.temp $outfile
	putslog "Delete old files for building $outfile ($args)"
	file delete {*}$args
}

proc cg_cg2bam {args} {
	set ::submit_direct 0
	while 1 {
		set key [lindex $args 0]
		if {$key eq "-direct"} {
			set ::submit_direct [lindex $args 1]
			set args [lrange $args 2 end]
		} else {
			break
		}
	}
	set len [llength $args]
	if {($len < 3) || ($len > 3)} {
		puts stderr "format is: cg cg2bam cgdir destdir dbdir"
		exit 1
	}	
	foreach {cgdir destprefix dbdir} $args break
	set destdir [file dir $destprefix]
	set chrs {0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y}
	set chrstodo {}
	foreach chr $chrs {
		if {[file exists $destprefix-chr$chr.bam]} continue
		lappend chrstodo $chr
	}
	set reffile [glob $dbdir/*.crr]
	set deps {}
	set files [glob $cgdir/MAP/*/mapping_*]
	set todo {}
	file mkdir -force $destdir/temp
	foreach mapfile $files {
		set tail [file tail $mapfile]
		regexp {mapping_([^.]*)} $tail temp name
		set readfile [lindex [glob [file dir $mapfile]/reads_$name*] 0]
		set outfile $destdir/temp/$name
		if {![multiexists $outfile-chr\$chr.bam $chrstodo]} {
			lappend deps [submit -io 1 cg map2bam $readfile $mapfile $reffile $outfile]
		} else {
			putslog "Skipping $outfile: already done"
		}
		lappend todo $outfile-chr\$chr.bam
	}
	set fdeps {}
	foreach chr $chrstodo {
		set outfile $destprefix-chr$chr.bam
		if {![file exists $outfile]} {
			lappend fdeps [submit -deps $deps -io 1 cg mergebam $outfile {*}[subst $todo]]
		} else {
			putslog "Skipping $outfile: already done"
		}
	}
	putslog "finished $destprefix"
	return $fdeps
}

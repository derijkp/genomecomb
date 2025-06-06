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
#	if {[multiexists $outfile-chr\$chr.bam $chrs]} {
#		putslog "skipping files $outfile: all bams exist already"
#		return
#	}
	set scratchdir [scratchdir]
	set scratchbase $scratchdir/temp_[file tail $outfile]
	set samfile $scratchbase.sam
	putslog "Making sam from $mapfile"
	set temptarget [filetemp $samfile]
	if {[catch {
		exec cgatools map2sam -r $readfile -m $mapfile -s $reffile --mate-sv-candidates -o $temptarget
	} error]} {
		regsub -all {IND_[^=]+=[0-9]+} $error {} error
		set error [string trim $error]
		if {$error ne {}} {error $error}
	}
	file rename -force -- $temptarget $samfile
	# read header
	if {[catch {
		set header [catch_exec samtools view --no-PG -S -H $samfile 2> /dev/null]
	} error]} {
		if {$error ne {[samopen] SAM header is present: 25 sequences.}} {error $error}
	}
	set header [split $header \n]
	set len [llength $header]
	set poss [list_find -glob $header @SQ*]
	set list [list_sub $header $poss]
	set list [bsort -index 1 $list]
	set pos [lindex $poss 0]
	set header [lreplace [list_sub $header -exclude $poss] $pos -1 {*}$list]
	putslog "Distributing sam to chromosomes for $name"
	exec tail -n +[expr {($len+1)}] $samfile | distr2chr $scratchbase 2
	file delete $samfile
	putslog "Making and sorting sam per chromosome for $name"
	foreach chr $chrs {
		file_write $scratchbase-chr$chr.sam [join $header \n]\n
		exec gnusort8 -T $scratchdir -t \t -N -s -k3,3 -k4,4 -k1,1 -k2,2 $scratchbase-$chr >> $scratchbase-chr$chr.sam
		file delete $scratchbase-$chr
	}
	putslog "Making bams for $name"
	foreach chr $chrs {
		set temptarget [filetemp $outfile-chr$chr.bam]
		if {[catch {
			catch_exec samtools view --no-PG -S -h -b -o $temptarget $scratchbase-chr$chr.sam
		} error]} {
			if {$error ne {[samopen] SAM header is present: 25 sequences.}} {error $error}
		}
		file delete $scratchbase-chr$chr.sam
		catch {bam_index $scratchbase-chr$chr.sam}
		file rename -force -- $temptarget $outfile-chr$chr.bam
	}
}

proc cg_mergebam {outfile args} {
	if {![llength $args] > 0} {
		error "format is: cg mergebam outfile bamfile1 ?bamfile2? ..."
	}
	putslog "Merging $outfile"
	set temptarget [filetemp $outfile]
	catch_exec samtools merge $temptarget {*}$args
	file rename -force -- $temptarget $outfile
	putslog "Delete old files for building $outfile ($args)"
	file delete {*}$args
}

proc cg_cg2bam {args} {
	set len [llength $args]
	if {($len < 3) || ($len > 3)} {
		error "format is: cg cg2bam cgdir destdir dbdir"
	}	
	foreach {cgdir destprefix dbdir} $args break
	set dbdir [dbdir $dbdir]
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
	file mkdir $destdir/temp
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

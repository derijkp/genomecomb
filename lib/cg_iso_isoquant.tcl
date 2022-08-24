proc sqanti3 {src dest refseq gtfannotation} {
	mkdir $dest
	set lines [lindex [cg select -g all $src] end]
	if {$lines == 0} {
		file_write $dest/[file tail $dest]_classification.txt "isoform	chrom	strand	length	exons	structural_category	associated_gene	associated_transcript	ref_length	ref_exons	diff_to_TSS	diff_to_TTS	diff_to_gene_TSS	diff_to_gene_TTS	subcategory	RTS_stage	all_canonical	min_sample_cov	min_cov	min_cov_pos	sd_cov	FL	n_indels	n_indels_junc	bite	iso_exp	gene_exp	ratio_exp	FSM_class	coding	ORF_length	CDS_length	CDS_start	CDS_end	CDS_genomic_start	CDS_genomic_end	predicted_NMD	perc_A_downstream_TTS	seq_A_downstream_TTS	dist_to_cage_peak	within_cage_peak	dist_to_polya_site	within_polya_site	polyA_motif	polyA_dist	ORF_seq	ratio_TSS\n"
		file_write $dest/[file tail $dest]_corrected.genePred {}
		file_write $dest/[file tail $dest]_corrected.faa {}
		file_write $dest/[file tail $dest]_corrected.fasta {}
		file_write $dest/[file tail $dest]_corrected.gtf {}
		file_write $dest/[file tail $dest]_corrected.gtf.cds.gff {}
		file_write $dest/[file tail $dest]_junctions.txt {}
	} else {
		catch_exec sqanti3_qc.py \
			$src \
			$gtfannotation \
			$refseq \
			-d $dest \
			-o [file tail $dest] \
			--report skip --skipORF
	}
	cg genepred2tsv \
		$dest/[file tail $dest]_corrected.genePred \
		$dest/[file tail $dest]_corrected.genepred.tsv
	set tempfileg [tempfile]
	set tempfilec [tempfile]
	cg select -overwrite 1 -s name $dest/[file tail $dest]_corrected.genepred.tsv $tempfileg
	cg select -overwrite 1 -s isoform $dest/[file tail $dest]_classification.txt $tempfilec
	set fg [open $tempfileg]
	set gheader [tsv_open $fg gcomment]
	set fc [open $tempfilec]
	set cheader [tsv_open $fc ccomment]
	set crem {isoform chrom strand exons 
		dist_to_cage_peak within_cage_peak dist_to_polya_site within_polya_site polyA_motif polyA_dist ratio_TSS ref_length ref_exons
		diff_to_TSS diff_to_TTS diff_to_gene_TSS diff_to_gene_TTS
		min_sample_cov min_cov min_cov_pos sd_cov FL n_indels n_indels_junc bite iso_exp gene_exp ratio_exp
	}
	set rmposs [list_cor $cheader $crem]
	set o [open $dest/isoforms-[file tail $dest].isoforms.tsv.temp w]
	puts $o [deindent {
		#filetype	tsv/transcriptsfile
		#fileversion	0.99
		#fields	table
		#fields	field	number	type	description
		#fields	name	1	String	Name of gene (usually transcript_id from GTF)
		#fields	chromosome	1	String	Chromosome name
		#fields	strand	1	String	+ or - for strand
		#fields	begin	1	Integer	Transcription start position
		#fields	end	1	Integer	Transcription end position
		#fields	cdsStart	1	Integer	Coding region start
		#fields	cdsEnd	1	Integer	Coding region end
		#fields	exonCount	1	Integer	Number of exons
		#fields	exonStarts	E	Integer	Exon start positions
		#fields	exonEnds	E	Integer	Exon end positions
		#fields	score	1	Float	Score
		#fields	name2	1	String	Alternate name (e.g. gene_id from GTF)
		#fields	cdsStartStat	1	String	Status of CDS start annotation (none, unknown, incomplete, or complete)
		#fields	cdsEndStat	1	String	Status of CDS end annotation (none, unknown, incomplete, or complete)
		#fields	exonFrames	E	Integer	Exon frame offsets {0,1,2}
		#fields	length	1	Integer	isoform length
		#fields	exons	1	Integer	Number of exons
		#fields	structural_category	1	String	one of the isoform categories absed on the best matching reference transcript (https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-output-explanation)
		#fields	associated_gene	1	String	the reference gene name
		#fields	associated_transcript	1	String	the reference transcript name
		#fields	ref_length	1	Integer	reference transcript length
		#fields	ref_exons	1	Integer	reference transcript number of exons
		#fields	diff_to_TSS	1	Integer	distance of query isoform 5' start to reference transcript start end. Negative value means query starts downstream of reference
		#fields	diff_to_TTS	1	Integer	distance of query isoform 3' end to reference annotated end site. Negative value means query ends upstream of reference
		#fields	diff_to_gene_TSS	1	Integer	distance of query isoform 5' start to the closest start end of any transcripts of the matching gene
		#fields	diff_to_gene_TTS	1	Integer	distance of query isoform 3' end to the closest end of any transcripts of the matching gene
		#fields	subcategory	1	String	additional splicing categorization, separated by semi-colons
		#fields	RTS_stage	1	String	TRUE if one of the junctions could be a RT switching artifact
		#fields	all_canonical	1	String	TRUE if all junctions have canonical splice sites
		#fields	min_sample_cov	1	String	sample with minimum coverage
		#fields	min_cov	1	Integer	minimum junction coverage based on short read STAR junction output file. NA if no short read given
		#fields	min_cov_pos	1	Integer	the junction that had the fewest coverage. NA if no short read data given
		#fields	sd_cov	1	Integer	standard deviation of junction coverage counts from short read data. NA if no short read data given
		#fields	FL	1	Integer	FL count associated with this isoform per sample if --fl_count is provided, otherwise NA
		#fields	n_indels	1	Integer	total number of indels based on alignment
		#fields	n_indels_junc	1	Integer	number of junctions in this isoform that have alignment indels near the junction site (indicating potentially unreliable junctions)
		#fields	bite	1	String	TRUE if contains at least one "bite" positive SJ
		#fields	iso_exp	1	Integer	short read expression for this isoform if --expression is provided, otherwise NA
		#fields	gene_exp	1	Integer	short read expression for the gene associated with this isoform (summing over all isoforms) if --expression is provided, otherwise NA
		#fields	ratio_exp	1	Integer	ratio of iso_exp to gene_exp if --expression is provided, otherwise NA
		#fields	FSM_class	1	String	classifies the transcript according to the expression of other isoforms in the gene to which the transcript belongs
		#fields	coding	1	String	coding or non_coding transcript
		#fields	ORF_length	1	Integer	predicted ORF length
		#fields	CDS_length	1	Integer	predicted CDS length
		#fields	CDS_start	1	Integer	CDS start
		#fields	CDS_end	1	Integer	CDS end
		#fields	CDS_genomic_start	1	Integer	genomic coordinate of the CDS start. If on - strand, this coord will be greater than the end
		#fields	CDS_genomic_end	1	Integer	genomic coordinate of the CDS end. If on - strand, this coord will be smaller than the start
		#fields	predicted_NMD	1	String	TRUE if there's a predicted ORF and CDS ends at least 50bp before the last junction; FALSE if otherwise. NA if non-coding
		#fields	perc_A_downstream_TTS	1	Float	percent of genomic "A"s in the downstream 20 bp window. If this number if high (say > 0.8), the 3' end site of this isoform is probably not reliable
		#fields	seq_A_downstream_TTS	1	String	sequence of the downstream 20 bp window
		#fields	dist_to_cage_peak	1	Integer	distance to closest TSS based on CAGE Peak data
		#fields	within_cage_peak	1	String	TRUE if the PacBio transcript start site is within a CAGE Peak
		#fields	dist_to_polya_site	1	Integer	if --polyA_motif_list is given, shows the location of the last base of the hexamer. Position 0 is the putative poly(A) site. This distance is hence always negative because it is upstream
		#fields	within_polya_site	1	String	
		#fields	polyA_motif	1	String	if --polyA_motif_list is given, shows the top ranking polyA motif found within 50 bp upstream of end
		#fields	polyA_dist	1	Integer	: if --polyA_motif_list is given, shows the location of the last base of the hexamer. Position 0 is the putative poly(A) site. This distance is hence always negative because it is upstream
		#fields	ORF_seq	1	String	ORF sequence
		#fields	ratio_TSS	1	Float	Using Short-Read data, we measure the mean coverage of the 100bp upstream and downstream a reported TSS.
		#fields	counts	1	Integer	Number of reads mapping to isoform
		#fields	tpm	1	Float	Transcripts per million (number of reads mapping nomralized to 1m reads total)
	}]
	# some checks need to be build in later
	puts $o [join [list {*}$gheader {*}[list_sub $cheader -exclude $rmposs]] \t]
	while 1 {
		if {[gets $fg gline] == -1} break
		if {[gets $fc cline] == -1} break
		set cline [split $cline \t]
		if {[lindex [split $gline \t] 0] ne [lindex $cline 0]} {
			error "error: mismatch at:\n$gline\n$cline"
		}
		puts $o $gline\t[join [list_sub $cline -exclude $rmposs] \t]
	}
	close $o ; close $fg ; close $fc
	file rename -force $dest/isoforms-[file tail $dest].isoforms.tsv.temp $dest/isoforms-[file tail $dest].isoforms.tsv
}

proc cg_iso_isoquant_genecounts {isoformcounts genecounts {genefield associated_gene}} {
	set tempfile [tempfile]
	cg select -overwrite 1 \
		-g $genefield \
		-gc {distinct(structural_category),distinct(subcategory),sum(counts*)} \
		$isoformcounts $tempfile
	catch {close $f} ; catch {close $o}
	set f [open $tempfile]
	set header [tsv_open $f]
	set oheader {gene gene_type}
	foreach field [lrange $header 3 end] {
		regsub sum_ $field {} field
		lappend oheader $field
	}
	set o [open $genecounts w]
	puts $o [join $oheader \t]
	while {[gets $f line] != -1} {
		set line [split $line \t]
		foreach {gene category sub} $line break
		set novel [regexp {^novel} $gene]
		if {!$novel && ![regexp _ $gene]} {
			set genetype known
		} elseif {[regexp fusion $category]} {
			set genetype fusion
		} elseif {!$novel} {
			set genetype overlap
		} elseif {[regexp antisense $category]} {
			set genetype novel_antisense
		} elseif {[regexp intergenic $category]} {
			set genetype novel_intergenic
		} elseif {[regexp genic_intron $category]} {
			set genetype novel_genic_intron
		} else {
			set genetype novel
		}
		puts $o $gene\t$genetype\t[join [lrange $line 3 end] \t]
	}
	close $o
	close $f
}

proc convert_isoquant_add_ambig {varVar ambig} {
	upvar $varVar var
	if {![info exists var]} {
		if {$ambig} {
			set var [expr {1.0/$ambig}]
		} else {
			set var 1
		}
	} else {
		if {$ambig} {
			set var [expr {$var + 1.0/$ambig}]
		} else {
			set var [expr {$var + 1}]
		}
	}
}

proc convert_isoquant_reggenedb {genedbtsv samregions refseq reggenedbtsv reggenedb} {
	set regfile [tempfile].tsv
	distrreg_reg2tsv $regfile $samregions $refseq
	cg regselect $genedbtsv $regfile > $reggenedbtsv
	# if empty
	set f [open $reggenedbtsv]
	tsv_open $f
	set read [gets $f line]
	close $f
	if {$read == -1} {
		set o [open $reggenedbtsv w]
		puts $o [join {chromosome begin end name gene strand cdsStart cdsEnd exonCount exonStarts exonEnds source transcript_name gene_id gene_name} \t]
		puts $o [join {chr1 1 2 dummy1 dummy2 + 1 1 1 1, 2, dummy3 dummy4 dummy5 dummy6} \t]
		close $o
		file_write [file dir $reggenedbtsv]/info_empty_known ""
	}
	cg tsv2gtf -genecol gene_id -addgene 1 $reggenedbtsv $reggenedb
	set tempgenedb [tempfile].db
}

proc convert_isoquant_add {varVar} {
	upvar $varVar var
	if {![info exists var]} {
		set var 1
	} else {
		incr var
	}
}

proc convert_isoquant {isodir destdir sample refseq genedb genedbtsv} {
	set read_assignmentsfile [gzfile $isodir/*.read_assignments.tsv]
	set isoformcountsfile $destdir/isoform_counts-isoquant-$sample.tsv
	set strictpct 90

	# model isos
	file delete -force $destdir/sqanti3-isoquant_models-$sample
	set src [gzfile $isodir/*.transcript_models.gtf]
	set dest $destdir/sqanti3-isoquant_models-$sample
	sqanti3 $src $dest $refseq $genedb

	# known isos
	unset -nocomplain sizea
	unset -nocomplain transcriptidsa
	array set transcriptidsa [split [cg select -hc 1 -g isoform_id $read_assignmentsfile] \n\t]
	set knownfile $destdir/transcripts_known-$sample.genepred.tsv
	set f [gzopen $genedbtsv]
	set header [tsv_open $f]
	set namepos [lsearch $header name]
	set startspos [lsearch $header exonStarts]
	set endspos [lsearch $header exonEnds]
	set o [open $knownfile w]
	puts $o [join $header \t]
	while {[gets $f line] != -1} {
		set split [split $line \t]
		set name [lindex $line $namepos]
		if {[info exists transcriptidsa($name)]} {
			set starts [lindex $split $startspos]
			set ends [lindex $split $endspos]
			set size 0
			foreach begin [split $starts ,] end [split $ends ,] {
				if {$begin eq ""} continue
				set size [expr {$size + ($end - $begin)}]
			}
			set sizea($name) $size
			puts $o $line
		}
	}
	close $o
	close $f
	cg tsv2gtf $knownfile [file root $knownfile].gtf
	sqanti3 [file root $knownfile].gtf $destdir/sqanti3-isoquant_known-$sample $refseq $genedb

	# check which reads have ambiguous mappings
	unset -nocomplain ambiga
	array set ambiga [split [cg select -hc 1 -g read_id $read_assignmentsfile | cg select -q {$count > 1}] \n\t]
	#
	# get convert data from sqanti3 models file
	unset -nocomplain converta
	unset -nocomplain outputa
	catch {close $f}
	set f [open [gzfile $destdir/sqanti3-isoquant_models-$sample/isoforms-sqanti3-isoquant_models-$sample.isoforms.tsv]]
	set header [tsv_open $f]
	set poss [list_cor $header {name subcategory associated_transcript exonStarts exonEnds}]
	set namepos [lindex $poss 0]
	while {[gets $f line] != -1} {
		set line [split $line \t]
		foreach {name subcategory associated_transcript starts ends} [list_sub $line $poss] break
		if {$subcategory in "reference_match mono-exon full_splice_match"} {
			set converta($name) $associated_transcript
			lset line $namepos $associated_transcript
			set name $associated_transcript
		}
		set outputa($name) $line
		set size 0
		foreach begin [split $starts ,] end [split $ends ,] {
			if {$begin eq ""} continue
			set size [expr {$size + ($end - $begin)}]
		}
		set sizea($name) $size
	}
	close $f
	# reads	
	# load model transcripts read info
	unset -nocomplain read2isoa
	catch {close $f} ; set f [gzopen [gzfile $isodir/*.transcript_model_reads.tsv]]
	gets $f
	while {[gets $f line] != -1} {
		foreach {read transcript} $line break
		lappend read2isoa($read) $transcript
	}
	gzclose $f
	catch {close $f} ; catch {close $o} ; 
	unset -nocomplain tcounta
	unset -nocomplain gcounta
	for {set p 0} {$p <= 100} {incr p} {set pcta($p) 0}
	#
	# make read_assignments
	# ---------------------
	# set temptarget $destdir/read_assignments-isoquant-$sample.tsv.temp
	# set target $destdir/read_assignments-isoquant-$sample.tsv
	set temptarget $destdir/read_assignments-isoquant-$sample.tsv.temp
	set target $destdir/read_assignments-isoquant-$sample.tsv
	set f [gzopen $read_assignmentsfile]
	set prev {}
	while {[gets $f line] != -1} {
		if {[string index $line 0] ne {#}} break
		set prev $line
	}
	set header [split [string range $prev 1 end] \t]
	set newheader {read_id chromosome begin end strand exonStarts exonEnds aligned_size}
	lappend newheader {*}[list_remove $header read_id chromosome chr strand exons]
	set poss [list_cor $header $newheader]
	lset poss 1 [lsearch $header chr]
	set exonspos [lsearch $header exons]
	set chrpos [lsearch $newheader chromosome]
	set readpos [lsearch $newheader read_id]
	set isopos [lsearch $newheader isoform_id]
	set genepos [lsearch $newheader gene_id]
	set infopos [lsearch $newheader additional_info]
	set beginpos [lsearch $newheader begin]
	set endpos [lsearch $newheader end]
	set strandpos [lsearch $newheader strand]
	set exonStartspos [lsearch $newheader exonStarts]
	set exonEndspos [lsearch $newheader exonEnds]
	set assignment_typepos [lsearch $newheader assignment_type]
	set aligned_sizepos [lsearch $newheader aligned_size]
#	set nposs [list_cor $newheader {
#		chromosome begin end strand exonStarts exonEnds 
#		read_id isoform_id gene_id additional_info assignment_type aligned_size
#	}]
	#
	# set o [wgzopen $temptarget w {} 4]
	set o [open $temptarget w]
	puts $o [join $newheader \t]\tambiguity\tcovered_pct\tpolya\tclassification\tclosest_known	
	while {[gets $f line] != -1} {
		set line [split $line \t]
		set exons [lindex $line $exonspos]
		set line [list_sub $line $poss]
		set read [lindex $line $readpos]
		set exonStarts {}
		set exonEnds {}
		set size 0
		foreach exon [split $exons ,] {
			foreach {s e} [split $exon -] break
			lappend exonStarts $s
			lappend exonEnds $e
			set size [expr {$size + ($e - $s)}]
		}
		set chr [lindex $line $chrpos]
		set strand [lindex $line $strandpos]
		set begin [lindex $exonStarts 0]
		set end [lindex $exonEnds end]
		set gene [lindex $line $genepos]
		if {[info exists ambiga($read)]} {set ambig $ambiga($read)} else {set ambig 0}
		set assignment_type [lindex $line $assignment_typepos]
		lset line $beginpos $begin
		lset line $endpos $end
		lset line $exonStartspos [join $exonStarts ,]
		lset line $exonEndspos [join $exonEnds ,]
		lset line $aligned_sizepos $size
		set closest_known [lindex $line $isopos]
# if {$read eq "509e19b0-14bb-4cec-95de-b491324817c7"} {error stop}
		if {$assignment_type eq "inconsistent" && [info exists read2isoa($read)]} {
			set expnum [get ambiga($read) 1] ; if {$expnum == 0} {set expnum 1}
			set iso [lindex $read2isoa($read) 0]
			if {[llength $read2isoa($read)] > 1} {set read2isoa($read) [lrange $read2isoa($read) 1 end]}
			if {$iso eq "*"} {
				set iso $closest_known
			} elseif {[info exists converta($iso)]} {
				set iso $converta($iso)
			}
			lset line $isopos $iso
		} else {
			set iso $closest_known
		}
		set info [lindex $line $infopos]
		set newinfo {}
		set polya {}
		set classification {}
		foreach temp [split $info \;] {
			set temp [string trim $temp]
			set pos [string first = $temp]
			if {$pos != -1} {
				set key [string range $temp 0 [expr {$pos-1}]]
				set value [string range $temp [expr {$pos+1}] end]
				if {$key eq "PolyA"} {
					set polya $value
				} elseif {$key eq "Classification"} {
					set classification $value
				} else {
					lappend newinfo $temp
				}
			} elseif {$temp ne ""} {
				lappend newinfo $temp
			}
			lset line $infopos $newinfo
		}
		if {[info exists sizea($iso)]} {
			set pct [expr {100*$size/$sizea($iso)}]
			if {$pct > 100} {set pct 100}
			incr pcta($pct)
		} else {
			set pct 0
		}
		puts $o [join $line \t]\t$ambig\t$pct\t$polya\t$classification\t$closest_known
		# counts
		convert_isoquant_add_ambig tcounta($iso,t) $ambig
		if {!$ambig} {
			# unique
			convert_isoquant_add tcounta($iso,u)
			if {$pct >= $strictpct} {
				# strict
				convert_isoquant_add tcounta($iso,s)
			}
		}
		if {$polya eq "True"} {
			convert_isoquant_add_ambig tcounta($iso,a) $ambig
			if {!$ambig} {
				# unique
				convert_isoquant_add tcounta($iso,au)
				if {$pct >= $strictpct} {
					# strict
					convert_isoquant_add tcounta($iso,as)
				}
			}
		}
		if {![info exists gcounta($gene)]} {
			set gcounta($gene) 0
		} else {
			incr gcounta($gene)
		}
	}
	gzclose $o
	gzclose $f
	file rename -force $temptarget $target

	# make isoform_counts file
	# ------------------------
	# read isoquant counts
	unset -nocomplain countsa
	set f [gzopen [gzfile $isodir/*.transcript_counts.tsv]]
	gets $f
	array set countsa [read $f]
	gzclose $f
	set f [gzopen [gzfile $isodir/*.transcript_model_counts.tsv]]
	gets $f
	array set countsa [read $f]
	gzclose $f

	catch {close $f} ; catch {close $o}
	set f [open [gzfile $destdir/sqanti3-isoquant_known-$sample/isoforms-sqanti3-isoquant_known-$sample.isoforms.tsv]]
	set header [tsv_open $f comments]
	set o [open $isoformcountsfile.temp w]
	set isopos [lsearch $header name]
	puts $o $comments
	puts $o [join $header \t]\tcounts-$sample\tcounts_weighed-isoquant-$sample\tcounts_unique-isoquant-$sample\tcounts_strict-isoquant-$sample\tcounts_aweighed-isoquant-$sample\tcounts_aunique-isoquant-$sample\tcounts_astrict-isoquant-$sample
	while {[gets $f line] != -1} {
		set split [split $line \t]
		set iso [lindex $split $isopos]
		if {[info exists outputa($iso)]} {
			unset outputa($iso)
		}
		set t [get tcounta($iso,t) 0]
		if {[get countsa($iso) 0] < 1 && $t < 1} continue
		puts $o $line\t[get countsa($iso)]\t[format %.2f $t]\t[get tcounta($iso,u) 0]\t[get tcounta($iso,s) 0]\t[format %.2f [get tcounta($iso,a) 0]]\t[get tcounta($iso,au) 0]\t[get tcounta($iso,as) 0]
	}
	foreach iso [array names outputa] {
		set line [join $outputa($iso) \t]
		puts $o $line\t[get countsa($iso)]\t[format %.2f [get tcounta($iso,t) 0]]\t[get tcounta($iso,u) 0]\t[get tcounta($iso,s) 0]\t[format %.2f [get tcounta($iso,a) 0]]\t[get tcounta($iso,au) 0]\t[get tcounta($iso,as) 0]
	}
	close $o
	close $f
	cg select -overwrite 1 -s - $isoformcountsfile.temp $isoformcountsfile.temp2
	file rename -force $isoformcountsfile.temp2 $destdir/isoform_counts-isoquant-$sample.tsv
	#
	# make gene_counts file
	# ---------------------
	# file delete $isoformcountsfile.temp
	set o [open $destdir/gene_counts-isoquant-$sample.tsv w]
	puts $o [join [list gene count-$sample] \t]
	foreach gene [array names gcounta] {
		puts $o $gene\t$gcounta($gene)
	}
	close $o

}


proc iso_isoquant_job {args} {
	# putslog [list iso_isoquant_job {*}$args]
	set cmdline "[list cd [pwd]] \; [list cg iso_isoquant {*}$args]"
	global appdir
	set distrreg chr
	set refseq {}
	set skips {}
	set genedbtsv {}
	set sqanti 1
	set compar joint
	set threads 8
	set quantification all
	upvar job_logdir job_logdir
	cg_options iso_isoquant args {
		-refseq {
			set refseq $value
		}
		-genedb {
			set genedbtsv $value
		}
		-distrreg {
			if {$value eq "regionfile"} {
				set distrreg regionfile
			} else {
				set distrreg [distrreg_checkvalue $value]
			}
		}
		-quantification {
			set quantification $value
		}
		-compar {
			set compar $value
		}
		-threads {
			set threads $value
		}
		-skip {
			lappend skips -skip $value
		}
	} {projectdir}
	# hanging problems with threads, run single
	set threads 1
	set projectdir [file_absolute $projectdir]
	set refseq [refseq $refseq]
	set genedbtsv [lindex [bsort [gzfiles [file dir $refseq]/extra/*gencode*.tsv]] end]
	if {![file exists $genedbtsv]} {
		set genedb [lindex [bsort [glob [file dir $refseq]/extra/*gencode*.gtf]] end]
		if {$genedb eq ""} {
			set genedb [lindex [bsort [glob [file dir $refseq]/*.gtf]] end]
		}
		set genedbtsv [tempfile].tsv
		cg_gtf2tsv $genedb $genedbtsv
	}
	# set iso_isoquantdir [findiso_isoquant]
	cd $projectdir
	set sampledirs [glob samples/*]
	job_logfile $projectdir/iso_isoquant_[file tail $projectdir] $projectdir $cmdline \
		{*}[versions iso_isoquant dbdir zstd os]

	# analysis per sample
	set regions [list_remove [distrreg_regs $distrreg $refseq] unaligned]
	foreach sampledir $sampledirs {
		putsvars sample
		cd $projectdir/$sampledir
		set sample [file tail $sampledir]
		set bam [lindex [jobglob map-sminimap*.bam map-*.bam] 0]
		set rootname [file_rootname $bam]
		set isofiles {}
		set genefiles {}
		set readfiles {}
		foreach region $regions {
			set regdir isoquant-$rootname/isoquant-$rootname-$region
			job isoquant-$sample-$region -mem 40G -cores $threads -deps {
				$bam $refseq $genedbtsv
			} -targets {
				$regdir
			} -vars {
				sampledir sample bam refseq regdir region genedbtsv quantification threads
			} -code {
				mkdir $regdir.temp
				# region bamfile
				# set tempbam [tempfile].bam
				set tempbam $regdir.temp/regali.bam
				set samregions [samregions $region $refseq]
				exec samtools view -h -b -1 $bam {*}$samregions > $tempbam
				exec samtools index $tempbam
				# region gene file
				set reggenedbtsv [tempfile].tsv
				set reggenedb [tempfile].gtf
				convert_isoquant_reggenedb $genedbtsv $samregions $refseq $reggenedbtsv $reggenedb
				set tempgenedb [tempfile].db
				exec isoquant3_gtf2db --complete_genedb --input $reggenedb --output $tempgenedb
				exec isoquant3 \
					--threads $threads \
					--transcript_quantification $quantification \
					--gene_quantification $quantification \
					--reference $refseq \
					--bam $tempbam \
					--genedb $tempgenedb \
					--data_type nanopore \
					--keep_tmp \
					-o $regdir.temp 2>@ stderr >@ stdout
				file rename $regdir.temp $regdir
			}
			lappend isofiles $regdir/isoform_counts-isoquant-$sample.tsv
			lappend genefiles $regdir/gene_counts-isoquant-$sample.tsv
			lappend readfiles $regdir/read_assignments-isoquant-$sample.tsv
			job isoquant-convert-$sample-$region -cores 1 -mem 10g -deps {
				$regdir $refseq $genedbtsv
			} -targets {
				$regdir/isoform_counts-isoquant-$sample.tsv
				$regdir/gene_counts-isoquant-$sample.tsv
				$regdir/read_assignments-isoquant-$sample.tsv
			} -vars {
				sampledir sample refseq regdir region genedbtsv
			} -code {
				set isodir $regdir/00_regali
				set destdir $regdir
				set reggenedbtsv [tempfile].tsv
				set reggenedb [tempfile].gtf
				set samregions [samregions $region $refseq]
				convert_isoquant_reggenedb $genedbtsv $samregions $refseq $reggenedbtsv $reggenedb
				convert_isoquant $isodir $destdir $sample $refseq $reggenedb $reggenedbtsv
			}
		}
		cd $projectdir/$sampledir
		set sample [file tail $sampledir]
		set bam [lindex [jobglob map-sminimap*.bam map-*.bam] 0]
		set rootname [file_rootname $bam]
		#
		# combine
		set isofiles {}
		set genefiles {}
		set readfiles {}
		set missing {}
		foreach region $regions {
			set regdir isoquant-$rootname/isoquant-$rootname-$region
			set isofile $regdir/isoform_counts-isoquant-$sample.tsv
			set genefile $regdir/gene_counts-isoquant-$sample.tsv
			set readfile $regdir/read_assignments-isoquant-$sample.tsv
			lappend isofiles $isofile
			lappend genefiles $genefile
			lappend readfiles $readfile
			if {![file exists $isofile]} {
				puts "$isofile missing"
				lappend missing $isofile
			}
			if {![file exists $genefile]} {
				puts "$genefile missing"
				lappend missing $genefile
			}
			if {![file exists $readfile]} {
				puts "$readfile missing"
				lappend missing $readfile
			}
		}
		job isoquant-join-$sample -cores 1 -deps [list {*}$isofiles {*}$genefiles {*}$readfiles] -targets {
			isoform_counts-isoquant-$sample.tsv
			gene_counts-isoquant-$sample.tsv
			read_assignments-isoquant-$sample.tsv
			totalcounts-isoquant-$sample.tsv
		} -vars {
			isofiles genefiles readfiles sampledir sample refseq regdir region genedbtsv genedb
		} -code {
			# cg cat {*}$isofiles > isoform_counts-isoquant-$sample.tsv.temp
			set o [open isoform_counts-isoquant-$sample.tsv.temp w]
			set f [open [lindex $isofiles 0]]
			set header [tsv_open $f comments]
			close $f
			puts $o $comments[join $header \t]
			unset -nocomplain donea
			foreach isofile $isofiles {
				set f [open $isofile]
				set cheader [tsv_open $f comments]
				if {$cheader ne $header} {
					error "header of $isofile differs from header of [lindex $isofiles 0]"
				}
				set pos [lsearch $cheader name]
				while {[gets $f line] != -1} {
					set split [split $line \t]
					set name [lindex $split $pos]
					if {[info exists donea($name)]} {
						if {[regexp {^(transcript_[0-9]+)(.*)$} $name temp pre post]} {
							set num 2
							while 1 {
								set test $pre-$num$post
								if {![info exists donea($test)]} break
								incr num
							}
							set name $test
						} else {
							set num 2
							while 1 {
								set test $name-$num
								if {![info exists donea($test)]} break
								incr num
							}
							set name $test
						}
						lset split $pos $name
						puts $o [join $split \t]
					} else {
						puts $o $line
					}
					set donea($name) 1
				}
				close $f
			}
			close $o
			file rename -force isoform_counts-isoquant-$sample.tsv.temp isoform_counts-isoquant-$sample.tsv
			#
			cg cat {*}$genefiles > gene_counts-isoquant-$sample.tsv.temp
			file rename gene_counts-isoquant-$sample.tsv.temp gene_counts-isoquant-$sample.tsv
			#
			cg cat {*}$readfiles | cg zst > read_assignments-isoquant-$sample.tsv.temp.zst
			file rename read_assignments-isoquant-$sample.tsv.temp.zst read_assignments-isoquant-$sample.tsv.zst
			# 
#			set temp [cg select -g all -gc {sum(count*-*)} isoform_counts-isoquant-$sample.tsv]
#			set temp [split $temp \n]
#			set header [split [list_shift temp] \t]
#			set header [list_regsub ^sum_ $header {}]
#			set o [open totalcounts-isoquant-$sample.tsv.temp w]
#			puts $o [join $header \t]
#			puts $o [join $temp \n]
#			close $o
			cg select -overwrite 1 -g all -gc {sum(count*-*)} isoform_counts-isoquant-$sample.tsv | cg select -rf all > totalcounts-isoquant-$sample.tsv.temp
			file rename -force totalcounts-isoquant-$sample.tsv.temp totalcounts-isoquant-$sample.tsv
		}
	}

	if {$compar eq "0"} return
	# combined analysis
	cd $projectdir
	mkdir compar
	set exproot [file tail $projectdir]
	set isoformfiles [bsort [jobglob samples/*/isoform_counts-isoquant-*.tsv]]
	set totalcountsfiles [bsort [jobglob samples/*/totalcounts-isoquant-*.tsv]]
	job iso_isoquant_compar-$exproot {*}$skips \
	-cores $threads \
	-deps [list_concat $isoformfiles $totalcountsfiles] \
	-targets {
		compar/isoform_counts-isoquant-$exproot.genepred.tsv
		compar/gene_counts-isoquant-$exproot.tsv
		compar/totalcounts-isoquant-$exproot.tsv
	} -vars {
		isoformfiles totalcountsfiles exproot
	} -code {
		analysisinfo_write [lindex $isoformfiles 0] $target
		analysisinfo_write [lindex $isoformfiles 0] compar/gene_counts-isoquant-$exproot.genepred.tsv
		analysisinfo_write [lindex $totalcountsfiles 0] compar/totalcounts-isoquant-$exproot.tsv
		set isoformcounts compar/isoform_counts-isoquant-$exproot.genepred.tsv
		cg multitranscript $isoformcounts {*}$isoformfiles
		set genecounts compar/gene_counts-isoquant-$exproot.tsv
		cg_iso_isoquant_genecounts $isoformcounts $genecounts
		cg paste {*}[bsort $totalcountsfiles] > compar/totalcounts-isoquant-$exproot.tsv
	}
}

proc cg_iso_isoquant {args} {
	set args [job_init {*}$args]
	iso_isoquant_job {*}$args
	job_wait
}


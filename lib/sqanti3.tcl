# (more) generic sqanti3 annotation
# originally made for isoquant analysis, but sqanti3 missed several annotations
# so in the end used other way to annotate there
# 
proc sqanti3 {src dest refseq gtfannotation {rename 0} {geneconvaVar {}}} {
	if {$geneconvaVar ne ""} {upvar $geneconvaVar geneconva}
	mkdir $dest
	set f [gzopen $src]
	while 1 {
		set read [gets $f line]
		if {$read == -1} break
		if {[string index $line 0] ne {#}} break
	}
	if {$read == -1} {
		close $f
		file_write $dest/[file tail $dest]_classification.txt "isoform	chrom	strand	length	exons	structural_category	associated_gene	associated_transcript	ref_length	ref_exons	diff_to_TSS	diff_to_TTS	diff_to_gene_TSS	diff_to_gene_TTS	subcategory	RTS_stage	all_canonical	min_sample_cov	min_cov	min_cov_pos	sd_cov	FL	n_indels	n_indels_junc	bite	iso_exp	gene_exp	ratio_exp	FSM_class	coding	ORF_length	CDS_length	CDS_start	CDS_end	CDS_genomic_start	CDS_genomic_end	predicted_NMD	perc_A_downstream_TTS	seq_A_downstream_TTS	dist_to_cage_peak	within_cage_peak	dist_to_polya_site	within_polya_site	polyA_motif	polyA_dist	ORF_seq	ratio_TSS\n"
		file_write $dest/[file tail $dest]_corrected.genePred {}
		file_write $dest/[file tail $dest]_corrected.faa {}
		file_write $dest/[file tail $dest]_corrected.fasta {}
		file_write $dest/[file tail $dest]_corrected.gtf {}
		file_write $dest/[file tail $dest]_corrected.gtf.cds.gff {}
		file_write $dest/[file tail $dest]_junctions.txt {}
	} else {
		set ok 1
		while 1 {
			if {[lindex [split $line \t] 6] ni "+ -"} {
				set ok 0
				break
			}
			set read [gets $f line]
			if {$read == -1} break
		}
		close $f
		if {!$ok} {
			# sqanti3 gives an error on transcripts with undetermined strand (strand . is not actually correct in gtf)
			# replace with + to make it work (although not neccesarily correct) anyway
			set tempfile [tempfile].gtf
			catch {close $f} ; catch {close $o}
			set f [open $src]
			set o [open $tempfile w]
			while 1 {
				set read [gets $f line]
				if {$read == -1} break
				if {[string index $line 0] ne "#"} break
				puts $o $line
			}
			while 1 {
				set split [split $line \t]
				if {[lindex $split 6] ni "+ -"} {
					lset split 6 +
					set line [join $split \t]
				}
				puts $o $line
				set read [gets $f line]
				if {$read == -1} break
			}
			close $o ; close $f
			set src $tempfile
		}
		file delete -force $dest
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
		#fields	polyA_dist	1	Integer	if --polyA_motif_list is given, shows the location of the last base of the hexamer. Position 0 is the putative poly(A) site. This distance is hence always negative because it is upstream
		#fields	ORF_seq	1	String	ORF sequence
		#fields	ratio_TSS	1	Float	Using Short-Read data, we measure the mean coverage of the 100bp upstream and downstream a reported TSS.
		#fields	counts	1	Integer	Number of reads mapping to isoform
		#fields	tpm	1	Float	Transcripts per million (number of reads mapping nomralized to 1m reads total)
	}]
	if {$rename} {
		set renameposs [list_cor $gheader {chromosome strand begin exonStarts exonEnds gene}]
	}
	# some checks need to be build in later
	if {!$rename} {
		puts $o [join [list {*}$gheader {*}[list_sub $cheader -exclude $rmposs]] \t]
	} else {
		puts $o [join [list name size oriname {*}[lrange $gheader 1 end] {*}[list_sub $cheader -exclude $rmposs]] \t]
	}
	while 1 {
		if {[gets $fg gline] == -1} break
		if {[gets $fc cline] == -1} break
		set cline [split $cline \t]
		set splitgline [split $gline \t]
		if {[lindex $splitgline 0] ne [lindex $cline 0]} {
			error "error: mismatch at:\n$gline\n$cline"
		}
		if {!$rename} {
			puts $o $gline\t[join [list_sub $cline -exclude $rmposs] \t]
		} else {
			foreach {chromosome strand begin exonStarts exonEnds gene} [list_sub $splitgline $renameposs] break
			set newname [iso_name $chromosome $strand $exonStarts $exonEnds size]
			if {[info exists geneconva($gene)]} {
				lset splitgline [lindex $renameposs end] $geneconva($gene)
				set gline [join $splitgline \t]
			}
			puts $o $newname\t$size\t$gline\t[join [list_sub $cline -exclude $rmposs] \t]
		}
	}
	close $o ; close $fg ; close $fc
	file rename -force $dest/isoforms-[file tail $dest].isoforms.tsv.temp $dest/isoforms-[file tail $dest].isoforms.tsv
}


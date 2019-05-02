if 0 {
set build hg19
set dest /complgen/refseq/
set path $dest/tmp
set url ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz
}

proc cg_download_exac {args} {
	set url ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz
	set keep 0
	cg_options download_exac args {
		-k {set keep $value}
	} {resultfile url} 1 2 {
		download data from exac
	}
	set tempdir $resultfile.temp
	file mkdir $tempdir
	set tail [file tail $url]
	if {![regexp {[0-9][0-9.]+[0-9]} $tail version]} {
		if {![regexp {[0-9]+} $tail version]} {
			error "could not get version from $tail"
		}
	}
	wgetfile $url $tempdir/$tail
	set tsvfile $tempdir/[file root $tail].tsv
	set tempresult $tempdir/[file tail resultfile]
	if {![file exists $tsvfile]} {
		set populations {AFR AMR EAS FIN NFE OTH SAS}
		set nheader {chromosome begin end type ref alt}
		foreach pop $populations {
			lappend nheader "freqp_${pop}=if(\$AN_$pop < 20,\"-\",format(\"%.3f\",100.0*\$AC_$pop/\$AN_$pop))"
			lappend nheader "mfreqp_${pop}=if(\$AN_$pop < 20,\"-\",format(\"%.3f\",100.0*(def(\$Hom_$pop,0) + def(\$Hemi_$pop,0))/\$AN_$pop))"
			lappend nheader "nchrom_${pop}=\$AN_$pop"
		}
		lappend nheader filter
		cg vcf2tsv -split 1 -sort 0 $tempdir/$tail | cg select -f $nheader > $tsvfile.temp
		file rename -force $tsvfile.temp $tsvfile
	}
	cg select -s - $tsvfile | cg collapsealleles > $tempresult
	# opt
	file_write $tempresult.opt "fields\t{freqp_AFR mfreqp_AFR freqp_AMR mfreqp_AMR freqp_EAS mfreqp_EAS freqp_FIN mfreqp_FIN freqp_NFE mfreqp_NFE freqp_OTH mfreqp_OTH freqp_SAS mfreqp_SAS}\n"
	# info
	file_write $tempresult.info [deindent [subst {
		ExAC variants (version $version)
		=============
		
		Download info
		-------------
		dbname	evs (ExAC variants)
		version	$version
		citation	Exome Aggregation Consortium (ExAC), Cambridge, MA (URL: http://exac.broadinstitute.org) \[[clock format [clock seconds] -format "%b %Y"]\]
		website	http://exac.broadinstitute.org
		source	$url
		time	[timestamp]
		
		About ExAC
		----------
		The Exome Aggregation Consortium (ExAC) is a coalition of investigators
		seeking to aggregate and harmonize exome sequencing data from a variety of
		large-scale sequencing projects, and to make summary data available for
		the wider scientific community.
		
		The data set provided on this website spans 60,706 unrelated individuals
		sequenced as part of various disease-specific and population genetic
		studies. We have removed individuals affected by severe pediatric disease,
		so this data set should serve as a useful reference set of allele
		frequencies for severe disease studies. All of the raw data from these
		projects have been reprocessed through the same pipeline, and jointly
		variant-called to increase consistency across projects.
		
		Data Usage
		----------
		All data here are released under a Fort Lauderdale Agreement for the
		benefit of the wider biomedical community. You can freely download and
		search the data, and use it for publications focused on specific sets of
		variants (for instance, assessing the frequency of a set of candidate
		causal variants observed in a collection of rare disease patients).
		However, we ask that you not publish global (genome-wide) analyses of
		these data until after the ExAC flagship paper has been published,
		estimated to be in early 2015. If you're uncertain which category
		your analyses fall into, please email us.
		
		The aggregation and release of summary data from the exomes collected by
		the Exome Aggregation Consortium has been approved by the Partners IRB
		(protocol 2013P001477, \u201cGenomic approaches to gene discovery in rare
		neuromuscular diseases\u201d).
		
		Citation
		--------
		We anticipate the publication of a flagship paper describing this data set
		in early 2015. Until that is published, please cite: Exome Aggregation
		Consortium (ExAC), Cambridge, MA (URL: http://exac.broadinstitute.org)
		\[[clock format [clock seconds] -format "%b %Y"]\].
		
		Annotation fields added
		-----------------------
		For each population following values are given
		  freqp_POP: frequency (as percentage!) of alt aleles in population POP
		  mfreqp_POP: frequency (as percentage!) of homozygous (or hemizygous) alt aleles in population POP
		  nchrom_POP: nr of sequenced chromosomes in population POP
		
		for each of the following populations
		AFR African/African American
		AMR American
		EAS East Asian
		FIN Finnish
		NFE Non-Finnish European
		OTH Other
		SAS South Asian
		
		If less than 20 chromosomes are sequenced, the frequency will be put as -
		
		more information on http://exac.broadinstitute.org/about
	}]]
	file rename -force $tempresult.opt [gzroot $resultfile].opt
	file rename -force $tempresult.info [gzroot $resultfile].info
	compress $tempresult $resultfile
	if {!$keep} {file delete -force $tempdir}
}


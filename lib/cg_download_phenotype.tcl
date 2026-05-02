proc cg_download_phenotype {args} {
	cg_options download_phenotype args {
	} {resultfile build}
	if {$build eq "hg19"} {
		set nbuild GRCh37
	} elseif {$build eq "hg38"} {
		set nbuild GRCh38
	} else {
		error "build $build is not supported (only hg19 and hg38)"
	}
	set tempdir $resultfile.temp
	file mkdir $tempdir
	wgetfile https://www.ebi.ac.uk/gene2phenotype/api/panel/all/download $tempdir/gene2phenotype.csv
	cg csv2tsv $tempdir/gene2phenotype.csv | cg select -f {
		{gene symbol} {disease name} {disease mim} confidence
	} -nh {
		gene phenotype_description disease_mim confidence
	} > $tempdir/gene2phenotype.tsv

	# clinvar
	set url ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_$nbuild/clinvar.vcf.gz
	wgetfile $url $tempdir/clinvar_$nbuild.vcf.gz
	cg vcf2tsv $tempdir/clinvar_$nbuild.vcf.gz $tempdir/clinvar_$nbuild.tsv

	# make final
	unset -nocomplain a
	# process hsapiens_gene_ensembl
	set f [open $tempdir/gene2phenotype.tsv]
	set header [tsv_open $f]
	while 1 {
		if {[gets $f line] == -1} break
		set line [split $line \t]
		foreach {gene pheno} $line break
		if {$gene eq "" || $pheno eq ""} continue
		set pheno [string tolower $pheno]
		set a($gene) [list_union [get a($gene) ""] [list $pheno]]
	}
	close $f
	# process clinvar

	catch {gzclose $f}
	set f [gzopen $tempdir/clinvar_$nbuild.tsv]
	set header [tsv_open $f]
	set poss [list_cor $header {GENEINFO CLNDN}]
	while 1 {
		if {[gets $f line] == -1} break
		set line [split $line \t]
		foreach {name pheno} [list_sub $line $poss] break
		set gene [lindex [split $name :] 0]
		if {$gene eq ""} continue
		set pheno [string_change $pheno {_ { } {\x2c} {-}}]
		set pheno [string tolower $pheno]
		set pheno [split $pheno |]
		set pheno [list_lremove $pheno {{not specified} {not provided}}]
		if {![llength $pheno]} continue
		set a($gene) [list_union [get a($gene) ""] $pheno]
	}
	gzclose $f
	set o [open $tempdir/phenotype.tsv w]
	puts $o "gene\tphenotype_description"
	foreach gene [bsort [array names a]] {
		foreach pheno $a($gene) {
			puts $o $gene\t$pheno
		}
	}
	close $o
	# info
	file_write [gzroot $resultfile].info [subst [deindent {
		phenotype
		=========
		
		Download info
		-------------
		dbname	phenotype
		version	[timestamp]
		website	https://www.ncbi.nlm.nih.gov/clinvar/ , http://www.ensembl.org
		source	$url, https://www.ebi.ac.uk/gene2phenotype/api/panel/all/download
		time	[timestamp]
		
		Description
		-----------
		Gene-phenotype data file
		
		These gene-phenotype correlations are extracted from the ensembl gene database 
		using biomart combined with those found in the clinvar database.
	}]]\n
	compress $tempdir/phenotype.tsv $resultfile
}

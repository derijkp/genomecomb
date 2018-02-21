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
	# hsapiens_gene_ensembl
	cg_download_mart $tempdir/ens_phenotype.tsv hsapiens_gene_ensembl gene_ensembl_config {hgnc_symbol phenotype_description}
	# clinvar
	set url ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_$nbuild/clinvar.vcf.gz
	wgetfile $url $tempdir/clinvar_hg19.vcf.gz
	cg vcf2tsv $tempdir/clinvar_hg19.vcf.gz $tempdir/clinvar_hg19.tsv
	unset -nocomplain a
	# process hsapiens_gene_ensembl
	set f [open $tempdir/ens_phenotype.tsv]
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
	set f [gzopen $tempdir/clinvar_hg19.tsv]
	set header [tsv_open $f]
	set poss [list_cor $header {GENEINFO CLNDBN}]
	while 1 {
		if {[gets $f line] == -1} break
		set line [split $line \t]
		foreach {name pheno} [list_sub $line $poss] break
		set gene [lindex [split $name :] 0]
		if {$gene eq ""} continue
		set pheno [string_change $pheno {_ { } {\x2c} {-}}]
		set pheno [string tolower $pheno]
		set pheno [split $pheno |,]
		set pheno [list_lremove $pheno {{not specified} {not provided}}]
		if {![llength $pheno]} continue
		set a($gene) [list_union [get a($gene) ""] $pheno]
	}
	gzclose $f
	set o [open $tempdir/phenotype.tsv w]
	puts $o "gene\tphenotype_description"
	foreach gene [lsort -dict [array names a]] {
		foreach pheno $a($gene) {
			puts $o $gene\t$pheno
		}
	}
	close $o
	# info
	regsub -all \n\t\t {
		Gene-phenotype data file
		
		These gene-phenotype correlations are extracted from the ensembl gene database 
		using biomart combined with those found in the clinvar database.
	} \n temp
	file_write [gzroot $resultfile].info [string trim $temp]
file_write [gzroot $resultfile].info [subst [string trim {
phenotype
=========

Download info
-------------
dbname	phenotype
version	[timestamp]
website	https://www.ncbi.nlm.nih.gov/clinvar/ , http://www.ensembl.org
source	$url, http://www.ensembl.org/biomart
time	[timestamp]

Description
-----------
Gene-phenotype data file

These gene-phenotype correlations are extracted from the ensembl gene database 
using biomart combined with those found in the clinvar database.
}]]
	if {[file extension $resultfile] eq ".lz4"} {
		cg lz4 -i 1 $tempdir/phenotype.tsv
		file rename -force $tempdir/phenotype.tsv.lz4 $resultfile
		file rename -force $tempdir/phenotype.tsv.lz4.lz4i $resultfile.lz4i
	} else {
		file rename -force $tempdir/phenotype.tsv $resultfile
	}
}

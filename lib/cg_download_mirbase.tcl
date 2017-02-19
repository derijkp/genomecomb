proc cg_download_mirbase {resultfile species {release 21}} {
	set resultfile [file_absolute $resultfile]
	file delete $resultfile
	file mkdir $resultfile.index
	set gff3file $resultfile.index/$species.gff3
	set structfile $resultfile.index/miRNA.str
	set genomefile [glob [file dir $resultfile]/genome_*.ifas]
	wgetfile ftp://mirbase.org/pub/mirbase/$release/genomes/$species.gff3 $gff3file
	wgetfile ftp://mirbase.org/pub/mirbase/$release/miRNA.str.gz $structfile.gz
	wgetfile ftp://mirbase.org/pub/mirbase/$release/README $resultfile.info.temp
file_write $resultfile.info [subst [string trimleft {
= miRBase =

== Download info ==
dbname	mirbase
version	$release
citation	Kozomara, Ana, and Sam Griffiths-Jones. 2014. miRBase: Annotating High Confidence microRNAs Using Deep Sequencing Data. Nucleic Acids Research 42 (D1): D68\u201373. doi:10.1093/nar/gkt1181.
website	http://mirbase.org/
source	ftp://mirbase.org/pub/mirbase/$release/genomes/$species.gff3
time	[timestamp]

}]]
	exec cat $resultfile.info.temp >> $resultfile.info
	file delete $resultfile.info.temp
	cg unzip $structfile.gz
	convertmirbase $gff3file $resultfile $genomefile $structfile
	file delete -force $resultfile.index
}

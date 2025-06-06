proc cg_download_mirbase {resultfile species {release 21}} {
	set resultfile [file_absolute $resultfile]
	file delete $resultfile
	set workdir [workdir $resultfile]
	set gff3file $workdir/$species.gff3
	set structfile $workdir/miRNA.str
	set genomefile [glob [file dir $resultfile]/genome_*.ifas]
	# set base ftp://mirbase.org/pub/mirbase
	# set base https://www.mirbase.org/ftp
	# wgetfile $base/$release/genomes/$species.gff3 $gff3file
	# wgetfile $base/$release/miRNA.str.gz $structfile.gz
	# wgetfile $base/$release/README $resultfile.info.temp
	set base https://www.mirbase.org/download
	wgetfile $base/$species.gff3 $gff3file
	wgetfile https://github.com/derijkp/genomecomb/releases/download/0.109.0/miRNA.str.gz $structfile.gz
	wgetfile $base/README/ $resultfile.info.temp
file_write [gzroot $resultfile].info [subst [string trimleft {
= miRBase =

== Download info ==
dbname	mirbase
version	$release
citation	Kozomara, Ana, and Sam Griffiths-Jones. 2014. miRBase: Annotating High Confidence microRNAs Using Deep Sequencing Data. Nucleic Acids Research 42 (D1): D68\u201373. doi:10.1093/nar/gkt1181.
website	http://mirbase.org/
source	ftp://mirbase.org/pub/mirbase/$release/genomes/$species.gff3
time	[timestamp]

}]]
	exec cat $resultfile.info.temp >> [gzroot $resultfile].info
	file delete -force $resultfile.info.temp
	file delete -force $structfile
	cg unzip $structfile.gz
	convertmirbase $gff3file $resultfile $genomefile $structfile
	file delete -force $workdir
}

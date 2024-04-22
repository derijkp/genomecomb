proc cg_install {args} {
	set bindir {}
	set refdir {}
	cg_options annotate args {
		-refdir {set refdir $value}
		-bindir {set bindir $value}
	} {} 0 ...
	unset -nocomplain typea
	unset -nocomplain urla
	unset -nocomplain knowna
	foreach {item type url} {
		hg38	ref	https://genomecomb.bioinf.be/download/refdb_hg38-0.109.0.tar.gz
		hg38-cadd	ref	https://genomecomb.bioinf.be/download/refdb_hg38-cadd-0.109.0.tar.gz
		hg38-minimap2	ref	https://genomecomb.bioinf.be/download/refdb_hg38-minimap2-0.109.0.tar.gz
		hg38-star	ref	https://genomecomb.bioinf.be/download/refdb_hg38-star-0.109.0.tar.gz
		liftover	ref	https://genomecomb.bioinf.be/download/refdb_liftover-0.106.0.tar.gz
		mm10	ref	https://genomecomb.bioinf.be/download/refdb_mm10-0.106.0.tar.gz
		mm10-minimap2	ref	https://genomecomb.bioinf.be/download/refdb_mm10-minimap2-0.106.0.tar.gz
		mm10-star	ref	https://genomecomb.bioinf.be/download/refdb_mm10-star-0.106.0.tar.gz
		sacCer3	ref	https://genomecomb.bioinf.be/download/refdb_sacCer3-0.106.0.tar.gz
		sacCer3-minimap2	ref	https://genomecomb.bioinf.be/download/refdb_sacCer3-minimap2-0.106.0.tar.gz
		sacCer3-star	ref	https://genomecomb.bioinf.be/download/refdb_sacCer3-star-0.106.0.tar.gz
		ce11	ref	https://genomecomb.bioinf.be/download/refdb_ce11-0.106.0.tar.gz
		ce11-minimap2	ref	https://genomecomb.bioinf.be/download/refdb_ce11-minimap2-0.106.0.tar.gz
		dm6	ref	https://genomecomb.bioinf.be/download/refdb_dm6-0.106.0.tar.gz
		dm6-minimap2	ref	https://genomecomb.bioinf.be/download/refdb_dm6-minimap2-0.106.0.tar.gz

		clair3	bin	https://genomecomb.bioinf.be/download/extra/clair3-1.0.4-linux-x86_64.tar.gz
		cutesv	bin	https://genomecomb.bioinf.be/download/extra/cutesv-1.0.11-linux-x86_64.tar.gz
		dirR	bin	https://genomecomb.bioinf.be/download/extra/dirR-4.2.1-linux-x86_64.tar.gz
		flair	bin	https://genomecomb.bioinf.be/download/extra/flair-2.0-linux-x86_64.tar.gz
		flames	bin	https://genomecomb.bioinf.be/download/extra/flames-c1_413e09c-linux-x86_64.tar.gz
		flye	bin	https://genomecomb.bioinf.be/download/extra/flye-2.9.2-linux-x86_64.tar.gz
		gatk3	bin	https://genomecomb.bioinf.be/download/extra/GATK-3.8.1.0-gf15c1c3ef-java.tar.gz
		gatk	bin	https://genomecomb.bioinf.be/download/extra/gatk-4.1.8.1-java.tar.gz
		hisat2	bin	https://genomecomb.bioinf.be/download/extra/hisat2-2.2.1-linux-x86_64.tar.gz
		isoquant	bin	https://genomecomb.bioinf.be/download/extra/isoquant-3.3.0-linux-x86_64.tar.gz
		java-1.8	bin	https://genomecomb.bioinf.be/download/extra/java-1.8.0-openjdk-1.8.0.275.b01-0-linux-x86_64.tar.gz
		longshot	bin	https://genomecomb.bioinf.be/download/extra/longshot-0.4.1-linux-x86_64.tar.gz
		lumpy	bin	https://genomecomb.bioinf.be/download/extra/lumpy-0.3.1-linux-x86_64.tar.gz
		manta	bin	https://genomecomb.bioinf.be/download/extra/manta-1.6.0-linux-x86_64.tar.gz
		medaka	bin	https://genomecomb.bioinf.be/download/extra/medaka-1.4.4-linux-x86_64.tar.gz
		minimap2	bin	https://genomecomb.bioinf.be/download/extra/minimap2-2.24-linux-x86_64.tar.gz
		nanopolish	bin	https://genomecomb.bioinf.be/download/extra/nanopolish-0.13.2-linux-x86_64.tar.gz
		nextflow	bin	https://genomecomb.bioinf.be/download/extra/nextflow-22.10.0-linux-x86_64.tar.gz
		java	bin	https://genomecomb.bioinf.be/download/extra/openjdk-22-linux-x86_64.tar.gz
		picard	bin	https://genomecomb.bioinf.be/download/extra/picard-2.21.3-java.tar.gz
		python3	bin	https://genomecomb.bioinf.be/download/extra/python3-3.9-linux-x86_64.tar.gz
		scywalker	bin	https://genomecomb.bioinf.be/download/extra/scywalker-0.108.0-Linux-x86_64.tar.gz
		sniffles	bin	https://genomecomb.bioinf.be/download/extra/sniffles-2.2-linux-x86_64.tar.gz
		sqanti3	bin	https://genomecomb.bioinf.be/download/extra/sqanti3-4.2-linux-x86_64.tar.gz
		star	bin	https://genomecomb.bioinf.be/download/extra/STAR-2.7.9a_2021-06-25-linux-x86_64.tar.gz
		strelka	bin	https://genomecomb.bioinf.be/download/extra/strelka-2.9.10-linux-x86_64.tar.gz
		qorts	bin	https://genomecomb.bioinf.be/download/extra/QoRTs-1.3.6.tar.gz
		modkit	bin	https://genomecomb.bioinf.be/download/extra/modkit-0.2.4-linux-x86_64.tar.gz

		ont	preset	{minimap2 sniffles cutesv clair3 modkit}
		ontr	preset	{minimap2 isoquant clair3}
		srs	preset	{gatk gatk3 picard java-1.8 strelka manta lumpy java}
		rseq	preset	{star qorts gatk gatk3 java}
		scywalker	preset	{minimap2 isoquant dirR}
	} {
		set typea($item) $type
		set urla($item) $url
		lappend knowna($type) $item
	}
	set errorlist "presets: $knowna(preset)\nbinaries: $knowna(bin)\nreferences: $knowna(ref)"
	if {![llength $args]} {
		error "no item given, must be one or more of:\n$errorlist"
	}

	if {$bindir eq ""} {
		set bindir $::genomecombdir/extra
	}
	if {$refdir eq ""} {
		set refdir [file dir $::genomecombdir]/refseq
	}
	set items {}
	foreach item $args {
		set type $typea($item)
		if {$type eq "preset"} {
			foreach item $urla($item) {
				if {![info exists urla($item)]} {
					error "unknown item $item, must be one of:\n$errorlist"
				}
				lappend items $item
			}
		} else {
			if {![info exists urla($item)]} {
				error "unknown item $item, must be one of:\n$errorlist"
			}
			lappend items $item
		}
	}
	set items [list_remdup $items]
	foreach item $items {
		set type $typea($item)
		set url $urla($item)
		if {$type eq "bin"} {
			puts "Downloading $item ($type) from $url"
			mkdir $bindir
			cd $bindir
			set tempfile [tempdir]/[file tail $url]
			wgetfile $urla($item) $tempfile
			puts "Unpacking $item ($type) in $bindir from $url"
			exec tar xvzf $tempfile
			file delete $tempfile
		} else {
			puts "Downloading $item ($type) from $url"
			mkdir $refdir
			cd $refdir
			set tempfile [tempdir]/[file tail $url]
			wgetfile $urla($item) $tempfile
			puts "Installing $item ($type) in $refdir from $url"
			exec tar xvzf $tempfile
			file delete $tempfile
		}
	}
}

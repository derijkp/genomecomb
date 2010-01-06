package require http

proc ensembl_getregion {chr start end args} {
	set url www.ensembl.org
	set species Homo_sapiens
	set fts {repeat variation gene vegagene estgene}
	set field param
	foreach {key value} $args {
		switch $key {
			-archive {
				if {$value ne ""} {
					set url $value.archive.ensembl.org
				} else {
					set url www.ensembl.org
				}
				set field st
			}
			-species {
				set species $value
			}
			-fts {
				set fts $value
			}
			default {
				error "unkown option $key"
			}
		}
	}
	set query [subst {output=embl;r=$chr:$start-$end;strand=1;_format=Text}]
	foreach ft $fts {
		append query \;$field=$ft
	}
	set h [http::geturl http://$url/$species/Location/Export?$query]
	set data [http::data $h]
	http::cleanup $h
	if {[regexp ^> $data]} {
		error "error getting ensembl region $chr-$start-$end"
	}
	return $data
}

if 0 {
	http://www.ensembl.org/Homo_sapiens/Location/Export?output=embl;r=10:1000000-1001000;strand=1;param=similarity;param=repeat;param=genscan;param=contig;param=variation;param=marker;param=gene;param=vegagene;param=estgene;_format=Text
	wget http://may2009.archive.ensembl.org/Homo_sapiens/Location/Export?output=embl;r=1:100000-110000;strand=1;time=1258453242.93979;st=similarity;st=repeat;st=genscan;st=contig;st=variation;st=marker;st=gene;st=vegagene;st=estgene;_format=Text
	package require http
	ensembl_getregion 10 1000000 1001000 -archive may2009
	ensembl_getregion 10 1000000 1001000
}

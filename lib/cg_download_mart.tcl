proc cg_download_mart {file dataset config attributes} {
	set attributesxml "<Attribute name=\"[join $attributes "\"/><Attribute name=\""]\"/>"
	set query [string trim [regsub -all \n [subst {
		query=<!DOCTYPE Query>
			<Query client=\"true\" processor=\"TSV\" limit=\"-1\" header=\"1\">
			<Dataset name=\"$dataset\" interface=\"default\">
			$attributesxml
			</Dataset></Query>
	}] {}]]
	puts "Downloading $file"
	if {[catch {
		exec wget --quiet -O $file.temp http://central.biomart.org/biomart/martview/results?$query
	}]} {
		exec wget --quiet -O $file.temp http://www.ensembl.org/biomart/martservice?$query
	}
	set header [list_change $attributes {external_gene_name gene name_1006 go_name namespace_1003 go_domain go_linkage_type go_evidence}]
	cg select -overwrite 1 -nh $header $file.temp $file.temp2
	cg select -overwrite 1 -q {$gene ne ""} $file.temp2 $file.temp3[file extension $file]
	file rename -force -- $file.temp3[file extension $file] $file
	file delete $file.temp $file.temp2
}

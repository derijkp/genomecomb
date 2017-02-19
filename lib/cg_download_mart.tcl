proc cg_download_mart {file dataset config attributes} {
	set attributesxml "<Attribute name=\"[join $attributes "\"/><Attribute name=\""]\"/>"
	set query "query=<!DOCTYPE Query><Query client=\"true\" processor=\"TSV\" limit=\"-1\" header=\"1\"><Dataset name=\"$dataset\" config=\"$config\">$attributesxml</Dataset></Query>"
	if {[catch {
		exec wget --quiet -O $file.temp http://central.biomart.org/biomart/martview/results?$query
	}]} {
		exec wget --quiet -O $file.temp http://www.ensembl.org/biomart/martservice?$query
	}
	set header [list_change $attributes {hgnc_symbol gene}]
	cg select -nh $header $file.temp $file.temp2
	cg select -q {$gene ne ""} $file.temp2 $file.temp3
	file rename -force $file.temp3 $file
	file delete $file.temp $file.temp2
}

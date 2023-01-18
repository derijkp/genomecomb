#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

# todo:
# genotype in haploid calls (Y chromosome)

proc cg_gff2sft {args} {
	cg_gff2tsv {*}$args
}

proc cg_gff2tsv {args} {
	set filename -
	set outfile -
	set transcripts 0
	set ignorecodon 0
	cg_options gtf2tsv args {
		-transcripts {
			set transcripts $value
		}
		-ignorecodon {
			set ignorecodon $value
		}
	} {filename outfile} 0 2
	proc gtf2tsv_parse_attr {attributes} {
		set a [dict create {*}[split [string_change $attributes {"; " ";"}] {;=}]]
	}
	gtf2tsv $filename $outfile $transcripts $ignorecodon
}

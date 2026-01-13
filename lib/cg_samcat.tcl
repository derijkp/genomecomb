proc cg_samcat {args} {
	set checkheaders 1
	cg_options samcat args {
		-checkheaders {set checkheaders $value}
	} {} 1
	foreach sam $args {
		if {[gziscompressed $sam]} {
			set header [exec cg zcat $sam | samtools view --no-PG -H]
		} else {
			set header [exec samtools view --no-PG -H $sam]
		}
		if {![info exists refheader]} {
			set refheader $header
			puts $header
		} else {
			if {$checkheaders && $header ne "$refheader"} {
				puts stderr "cannot samcat: differences in header"
				exit 1
			}
		}
		if {[gziscompressed $sam]} {
			exec cg zcat $sam | samtools view --no-PG >@ stdout
		} else {
			exec samtools view --no-PG $sam >@ stdout
		}
	}
}

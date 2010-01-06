proc ucsc2region {ucsc_file} {
	catch {close $f}
	set f [open $ucsc_file]
	set temp [string range [gets $f] 1 end]
	set header [split $temp \t]
	set poss [list_cor $header {chrom chromStart chromEnd name}]
	puts [join [list chromosome start end name] \t]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		foreach {chrom chromStart chromEnd name} [list_sub $line $poss] break
		puts [join [list $chrom $chromStart $chromEnd $name] \t]
	}
	close $f
}

if 0 {
# ------------------------------------------------------------------------------
lappend auto_path ~/dev/completegenomics/lib
package require Extral
package require Tclx
signal -restart error SIGINT

cd /data/db
set ucsc_file _data_db_ucsc-simple_repeats.tsv

}

proc map2besthits {file outfile} {
	set o [open $outfile w]
	set f [opencgifile $file header]
	set num 1
	while {![eof $f]} {
		incr num
		if {![expr {$num % 100000}]} {putslog $num}
		set list [readcgimap $f]
		foreach {flags chr start end strand side} [lindex $list 0] break
		puts $o $chr\t[expr {($start+$end)/2}]
		if {$side eq "l"} {
			set pos [lsearch [list_subindex $list 5] r]
			if {$pos == -1} continue
			set line [lindex $list $pos]
			foreach {flags chr start end side} $line break
			puts $o $chr\t[expr {($start+$end)/2}]
		}
	}
	close $o
	close $f
}

proc cg_cg2cnv {args} {
	if {[llength $argv] != 2} {
		puts stderr "format is: $scriptname $action dir outdir"
		exit 1
	}
	foreach {dir outdir} $args break
	cd $dir
	file mkdir $outdir
	set files [glob */mapping.tsv.gz]
	set base [file tail [file root $dir]]_hits
	foreach file $files {
		puts $file
		exec zcat $file | map2besthit | distr2chr $outdir/$base
	}
	exec touch $outdir/FINISHED
}

proc cg_cnv-seq {args} {
	if {[llength $argv] != 3} {
		puts stderr "format is: $scriptname $action refdir testdir outdir"
		exit 1
	}
	foreach {ref test outdir} $args break
	set ref [file normalize $ref]
	set test [file normalize $test]
	file mkdir $outdir
	cd $outdir
	foreach {chr size} {
	        1 249250621 
	        2 243199373 
	        3 199501827 
	        4 191273063
	        5 180857866
	        6 170899992
	        7 158821424
	        8 146274826
	        9 140273252
	        10 135374737
	        11 134452384
	        12 132349534
	        13 114142980
	        14 106368585
	        15 100338915
	        16 88827254 
	        17 78774742
	        18 76117153
	        19 63811651
	        20 62435964
	        21 46944323
	        X 154913754
	        Y 57772954 
	} {
		puts $chr
		set error [catch {
			exec cnv-seq.pl --genome-size=$size --ref=$ref/MAP_hits-$chr --test=$test/MAP_hits-$chr >@ stdout
		} result]
		puts $result
		set cnvout [glob MAP_hits-$chr-vs-MAP_hits-$chr.log2-*.cnv]
		set cnvout [file root $cnvout]
		set f [open "| R --vanilla --slave >@ stdout" w]
		puts $f "suppressMessages(library(cnv))"
		puts $f "data<-read.delim(\"$cnvout.cnv\")"
		puts $f "sink(file=\"$cnvout.result\")"
		puts $f "cnv.print(data)"
		puts $f "sink()" 
		puts $f "sink(file=\"$cnvout.summary\")"
		puts $f "cnv.summary(data)"
		puts $f "sink()" 
		catch {close $f} result
		puts $result
	}
}

proc cnvmedian {} {

# raw sata
	set outfile [gzroot $filename].cov
	set f [open "| zless $filename"]
	set o [open $outfile w]
	puts $o "pos\tcoverage"
	set header [tsv_open $f]
	set num 0
	while {![eof $f]} {
		incr num
		if {$num > 1000000} break
		if {![expr {$num%100000}]} {putslog $num}
		set line [getnotempty $f]
		puts $o [lindex $line 0]\t[lindex $line 2]
	}
	close $o
	close $f	

# code
package require Tclx
signal -restart error SIGINT
lappend auto_path ~/dev/completegenomics/lib
package require Extral
cd /complgen
cd /media/passport/complgen
set filename GS00102/coverageRefScore-20-GS000000078-ASM.tsv.gz
set filename GS00103/coverageRefScore-20-GS000000079-ASM.tsv.gz
set outfile [gzroot $filename].med
set windowsize 5000
catch {close $o}
catch {close $f}
	set f [open "| zless $filename"]
	set o [open $outfile w]
	puts $o "offset\tmedcoverage"
	set header [tsv_open $f]
	set list {}
	unset -nocomplain plist
	set line [getnotempty $f]
	set start [lindex $line 0]
	set start [expr {$start-$start%$windowsize}]
	set end [expr {$start + $windowsize -1}]
	set list [list [lindex $line 2]]
	set num 0
	while {![eof $f]} {
		incr num
		if {![expr {$num%100000}]} {putslog $num}
		while {![eof $f]} {
			if {![llength $line]} continue
			set refpos [lindex $line 0]
			if {$refpos >= $end} break
			lappend list [list [lindex $line 2]]
			set line [split [gets $f] \t]
		}
		if {[info exists plist]} {
			set temp [list_concat $plist $list]
			if {[llength $temp] > $windowsize} {
				set temp [lsort -integer $temp]
				set v [lindex $temp [expr {[llength $temp]/2}]]
				if {![isint $v]} {set v 0}
				puts $o $start\t$v
			}
		}
		set plist $list
		set list {}
		incr start $windowsize
		incr end $windowsize
	}

	close $o
	close $f

}

if 0 {

eval cg2cnv $argv

set dir /mnt/extra/CompleteGenomics/GS00102-DNA-D06/MAP/
set outdir /home/MOLGEN/peterdr/complgen/cnv-GS00102-DNA-D06

cgicnv /mnt/extra/CompleteGenomics/GS00102-DNA-D06/MAP/ /home/MOLGEN/peterdr/complgen/cnv-GS00102-DNA-D06

map2besthit < mapping.tsv | distr2chr distr

lappend auto_path ~/dev/completegenomics/lib
package require Extral
set file /complgen/GS00103/MAP/GS000005015-MAP/mapping.tsv.gz
set outfile /complgen/cnv/79.hits

map2besthits $file $outfile


set file /data/peter/complgendata/CGI-DATA-sample/GS00028-DNA-C01/MAP/GS000005323-MAP-sample/mapping-100000.tsv
time {map2besthits $file /data/peter/complgendata/tbesthits.tsv}
1743500 microseconds per iteration

time {exec ~/dev/completegenomics/bin/map2besthit < $file > /data/peter/complgendata/besthits.tsv}




}

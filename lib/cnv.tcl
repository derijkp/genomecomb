proc map2besthits {file outfile} {
	set o [open $outfile w]
	set f [opencgifile $file header]
	set num 1
	while {![eof $f]} {
		incr num
		if {![expr {$num % 100000}]} {puts $num}
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

proc cg2cnv {dir outdir} {
	cd  $dir
	file mkdir $outdir
	set files [glob */mapping.tsv.gz]
	set base [file tail [file root $dir]]_hits
	foreach file $files {
		puts $file
		exec zcat $file | map2besthit | distr2chr $outdir/$base
	}
	exec touch $outdir/FINISHED
}

proc cnvseq {ref test outdir} {
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

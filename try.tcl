lappend auto_path ~/dev/completegenomics/lib
package require Extral




scp ~/dev/completegenomics/bin/* 143.169.30.14:/home/peter/dev/completegenomics/bin
scp ~/dev/completegenomics/bin/* peterdr@ristretto:bin
ssh peterdr@ristretto




cd /data/peter/complgendata/
~/dev/completegenomics/bin/map2besthit < mapping-100000.tsv > mapping.hits
 ~/dev/completegenomics/bin/distr2chr distr < mapping.hits

~/dev/completegenomics/bin/map2besthit < mapping-100000.tsv | ~/dev/completegenomics/bin/distr2chr distr

cat mapping.tsv.gz | ~/dev/completegenomics/bin/map2besthit | distr2chr distr



zcat mapping.tsv.gz | ~/dev/completegenomics/bin/map2besthit | less



set file /complgen/GS00102/var-GS000000078-ASM.tsv
set genefile /complgen/GS00102/sgene-GS000000078-ASM.tsv
set outfile /complgen/GS00102/annotvar-GS000000078-ASM.tsv


var2annot $file $genefile $outfile





	set file /complgen/var-test.tsv
	set genefile /complgen/GS00102/sgene-GS000000078-ASM.tsv
	set outfile /complgen/annotvar-test.tsv
catch {close $f1}
catch {close $g}
	set f1 [opencgifile $file header]
		set list [readlocus $f1]
	set cur [readonevar $f1]
set flocus 12552239
	set g [opencgifile $genefile header]
	if {[split $header \t] ne "index locus haplotype chromosome begin end varType reference call xRef geneId mrnaAcc proteinAcc orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef"} {
		error "header error in $genefile"
	}
	set curgene [readonegene $g]
		set glocus [lindex $curgene 0]
		while {[isint $glocus] && ($glocus < $flocus)} {
			set curgene [readonegene $g]
			set glocus [lindex $curgene 0]
if {$glocus == 4311} break
putsvars glocus
		}


set cur1 {962269 2 2090 2090 snp A G 126 126 {}}
set cur2 {1110504 2 2153 2153 snp A T 120 62 {}}
set r1 file8
set r2 file9

4303    1       1101320 1101320 snp     A       C       53      40              100132716       XM_001724630.1  XP_001724682.1  -       EXON    0       Y       COMPATIBLE      335     111     R       R       R       

putsvars cur1 cur2

seek $f2 2013648400

wc diff.tsv
906085
grep -v ref-inconsistent diff.tsv | wc
835622
grep -v ref-inconsistent diff.tsv | grep -P '\tsnp\t' | wc
561545 (other = 274077)
grep -v 'ref-inconsistent' diff.tsv | grep -v '\?' | wc
525874
grep -v 'ref-inconsistent' diff.tsv | grep -v '\?' |  grep -P '\tsnp\t' | wc
399994 (other = 125880)

same
wc same.tsv
2998478
grep -v ref-inconsistent same.tsv | wc
2963002
grep -v ref-inconsistent same.tsv | grep -P '\tsnp\t' | wc
2642820 (other 320182)
grep -v 'ref-inconsistent' same.tsv | grep -v '\?' | wc
2826136
grep -v 'ref-inconsistent' same.tsv | grep -v '\?' |  grep -P '\tsnp\t' | wc
2584269 (other = 241867)
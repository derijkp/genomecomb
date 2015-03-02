set build hg19
set dest /complgen/refseq/
set path $dest/tmp
set url ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz
proc downloaddb_exac {path build {url {}}} {
	if {$url eq ""} {
		set url ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz
	}
	file mkdir ${path}/tmp/$build
	cd ${path}/tmp/$build
	# exec -ignorestderr wget -c --tries=45 --directory-prefix=${path}/tmp/$build/exac http://exac.broadinstitute.org/about
	exec -ignorestderr wget -c --tries=45 --directory-prefix=${path}/tmp/$build/exac $url
	set file [file tail $url]
	set tsvfile [file root $file].tsv
	cg vcf2tsv $file $tsvfile
	#
	catch {close $f} ; catch {close $o}
	set f [open $tsvfile]
	set header [tsv_open $f]
	set bposs [tsv_basicfields $header]
	set populations {AFR AMR EAS FIN NFE OTH SAS}
	set o [open var_${build}_exac.tsv.temp w]
	set nheader [list {*}[list_sub $header $bposs] ]
	foreach pop $populations {
		lappend nheader freqp_${pop} mfreqp_${pop} nchrom_${pop}
	}
	lappend nheader filter
	puts $o [join $nheader \t]
	set altpos [lsearch $bposs 5]
	set list {}
	foreach pop $populations {
		lappend list AC_$pop Hom_$pop Hemi_$pop AN_$pop
	}
	set poss [list_cor $header $list]
	set restposs [list_cor $header filter]
	while {[gets $f line] >= 0} {
		set line [split $line \t]
		set alt [lindex $line $altpos]
		set result [list_sub $line $bposs]
		foreach {acs homs hemis total} [list_sub $line $poss] pop $populations {
			set freqs {} ; set mfreqs {};
			if {$total >= 20} {
				foreach ac [split $acs ,] hom [split $homs ,] hemi [split $hemis ,] {
					if {$hemi eq ""} {set hemi 0}
					set freq [expr {100.0*$ac/$total}]
					lappend freqs [trimformat %.3f $freq]
					set mfreq [expr {100.0*(2*$hom+$hemi)/$total}]
					lappend mfreqs [trimformat %.3f $mfreq]
				}
				lappend result [join $freqs ,] [join $mfreqs ,] $total
			} else {
				lappend result - - 0
			}
		}
		lappend result {*}[list_sub $line $restposs]
		puts $o [join $result \t]
	}
	close $o
	close $f
# opt
file_write $path/$build/extra/var_${build}_exac.tsv.opt "fields\t{freqp_AFR mfreqp_AFR freqp_AMR mfreqp_AMR freqp_EAS mfreqp_EAS freqp_FIN mfreqp_FIN freqp_NFE mfreqp_NFE freqp_OTH mfreqp_OTH freqp_SAS mfreqp_SAS}\n"
# info
file_write $path/$build/extra/var_${build}_exac.tsv.info [string trim [subst {
ExAC variants
=============
About ExAC
----------
The Exome Aggregation Consortium (ExAC) is a coalition of investigators
seeking to aggregate and harmonize exome sequencing data from a variety of
large-scale sequencing projects, and to make summary data available for
the wider scientific community.

The data set provided on this website spans 60,706 unrelated individuals
sequenced as part of various disease-specific and population genetic
studies. We have removed individuals affected by severe pediatric disease,
so this data set should serve as a useful reference set of allele
frequencies for severe disease studies. All of the raw data from these
projects have been reprocessed through the same pipeline, and jointly
variant-called to increase consistency across projects.

Data Usage
----------
All data here are released under a Fort Lauderdale Agreement for the
benefit of the wider biomedical community. You can freely download and
search the data, and use it for publications focused on specific sets of
variants (for instance, assessing the frequency of a set of candidate
causal variants observed in a collection of rare disease patients).
However, we ask that you not publish global (genome-wide) analyses of
these data until after the ExAC flagship paper has been published,
estimated to be in early 2015. If you're uncertain which category
your analyses fall into, please email us.

The aggregation and release of summary data from the exomes collected by
the Exome Aggregation Consortium has been approved by the Partners IRB
(protocol 2013P001477, \u201cGenomic approaches to gene discovery in rare
neuromuscular diseases\u201d).

Citation
--------
We anticipate the publication of a flagship paper describing this data set
in early 2015. Until that is published, please cite: Exome Aggregation
Consortium (ExAC), Cambridge, MA (URL: http://exac.broadinstitute.org)
\[[clock format [clock seconds] -format "%b %Y"]\].

Annotation fields added
-----------------------
For each population following values are given
  freqp_POP: frequency (as percentage!) of alt aleles in population POP
  mfreqp_POP: frequency (as percentage!) of homozygous (or hemizygous) alt aleles in population POP
  nchrom_POP: nr of sequenced chromosomes in population POP

for each of the following populations
AFR African/African American
AMR American
EAS East Asian
FIN Finnish
NFE Non-Finnish European
OTH Other
SAS South Asian

If less than 20 chromosomes are sequenced, the frequency will be put as -

more information on http://exac.broadinstitute.org/about
}]]
	file rename -force var_${build}_exac.tsv.temp $path/$build/extra/var_${build}_exac.tsv
}


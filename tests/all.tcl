#!/bin/sh
# the next line restarts using wish \
exec cg source "$0" "$@"

source analysis.tcl
source annot.tcl
source bcol.tcl
source cg2tsv.tcl
source clip.tcl
source compar.tcl
source convert.tcl
source genome_seq.tcl
source libext.tcl
source lift.tcl
source mirannot.tcl
source misc.tcl
source pmulticompar.tcl
source process_multicompar.tcl
source reg.tcl
source remap.tcl
source select.tcl
source select_group.tcl
source select_saggregates.tcl
source tsv.tcl
source val.tcl
source varia.tcl
source vcf2tsv.tcl
source map.tcl
source job.tcl

# next ones take longer, use larger data in genomecomb.testdata
source queries.tcl
source reports.tcl
source bam.tcl
source homwes.tcl

# long time, run separately
# source process.tcl

# not really used
# source mselect.tcl

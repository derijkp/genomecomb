#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source analysis.tcl
source annot.tcl
source bcol.tcl
source cg2tsv.tcl
source clip.tcl
source compar.tcl
source convert.tcl
source genome_seq.tcl
source job.tcl
source libext.tcl
source lift.tcl
source mirannot.tcl
source misc.tcl
source pmulticompar.tcl
source process_multicompar.tcl
source queries.tcl
source reg.tcl
source remap.tcl
source select.tcl
source select_group.tcl
source select_saggregates.tcl
source tsv.tcl
source val.tcl
source varia.tcl
source vcf2tsv.tcl

# long time, run separately
# source process.tcl

# not really used
# source mselect.tcl

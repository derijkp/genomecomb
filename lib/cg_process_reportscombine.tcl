proc pget {var args} {
	upvar $var data
	foreach field $args {
		set name [lindex [bsort [array names data $field]] 0]
		if {$name ne ""} {
			return $data($name)
		}
	}
	return {}
}

proc report_getpattern {pattern data {default {}}} {
	if {[regexp $pattern $data temp result]} {
		return $result
	} else {
		return $default
	}
}

proc ppercent {value ref} {
	if {![isint $value] || ![isint $ref] || $ref == 0} {
		return ""
	} else {
		return [format %.2f [expr {100.0*$value/$ref}]]
	}
}

proc pgetsum {var p1 p2} {
	upvar $var data
	set v1 [pget data $p1]
	set v2 [pget data $p2]
	if {$v1 eq ""} {return $v2}
	if {$v2 eq ""} {return $v1}
	if {[catch {expr {$v1+$v2}} result]} {
		return $v1,$v2
	} else {
		return $result
	}
}

proc report_htmltable {tablename table} {
	set html "<table id=\"$tablename\" class=\"sort\">\n"
	set header [lindex $table 0]
	append html "<theader>\n<tr><th class=\"sort-header\">#</th><th class=\"sort-header\">[join $header {</th><th class="sort-header">}]</th></tr>\n</theader>\n<tbody>\n"
	set num 1
	foreach line [lrange $table 1 end] {
		append html <tr><td>$num</td><td>[join $line {</td><td>}]</td></tr>\n
		incr num
	}
	append html </tbody></table>\n
	append html [deindent [subst {
		<script>
		  new Tablesort(document.getElementById('$tablename'));
		</script>
	}]]
}

proc fastqc_readtable {file pattern {headerVar {}}} {
	if {$headerVar ne ""} {upvar $headerVar header}
	set f [open $file]
	while {[gets $f line] != -1} {
		if {[regexp $pattern $line]} break
	}
	set header [split [string range [gets $f] 1 end] \t]
	set result {}
	while {[gets $f line] != -1} {
		if {[regexp {END_MODULE} $line]} break
		set line [split $line \t]
		lappend result $line
	}
	close $f
	return $result
}

proc report_css {} {
return {
<style>
/* ------------------------------------------
  Reset ?
  http://meyerweb.com/eric/tools/css/reset/
  v2.0 | 20110126
  License: none (public domain)
---------------------------------------------*/
html, body, div, span, applet, object, iframe,
h1, h2, h3, h4, h5, h6, p, blockquote, pre,
a, abbr, acronym, address, big, cite, code,
del, dfn, em, img, ins, kbd, q, s, samp,
small, strike, strong, sub, sup, tt, var,
b, u, i, center,
dl, dt, dd, ol, ul, li,
fieldset, form, label, legend,
table, caption, tbody, tfoot, thead, tr, th, td,
article, aside, canvas, details, embed,
figure, figcaption, footer, header, hgroup,
menu, nav, output, ruby, section, summary,
time, mark, audio, video {
	margin: 0;
	padding: 0;
	border: 0;
	font-size: 100%;
	font: inherit;
	vertical-align: baseline;
}
/* HTML5 display-role reset for older browsers */
article, aside, details, figcaption, figure,
footer, header, hgroup, menu, nav, section {display: block;}
body {line-height: 1;}
ol, ul {list-style: none;}
blockquote, q {quotes: none;}
blockquote:before, blockquote:after,
q:before, q:after {content: ''; content: none;}
/* tables still need 'cellspacing="0"' in the markup */
table {border-collapse: collapse; border-spacing: 0;}
/**/
/* remember to define focus styles. Hee Haw */
:focus { outline: 0; }
*, *:before, *:after {
	-moz-box-sizing:border-box;
	     box-sizing:border-box;
}

body {
	margin:40px;
	font-family:'Helvetica Neue', Helvetica, Arial, sans-serif;
	font-size:12px;
	line-height:18px;
	color:#303030;
	background-color:#fafafa;
	-webkit-font-smoothing:antialiased;
}
h1,h2,h3,h4,h5 {
	font-weight:bold;
	display:block;
	margin:0 0 10px;
	break-after: avoid-page;
	page-break-after: avoid;
}
h1 {
	font-size:32px;
	margin:0 0 20px;
	display:block;
	font-weight:normal;
	text-shadow:0 1px 0 #fff;
}
h1 span.description {
	color:#6d6d6d;
}
h2 { font-size:18px; line-height:24px; margin:20px 0 10px;}
h3 { font-size:15px; }
ul { margin:0 0 20px; }
li {
	margin-left:30px;
	margin-bottom:3px;
}
ul li { list-style:disc; }
ol li { list-style:decimal; }
p {
	widows: 2;
	orphans: 2;
}
strong {
	font-weight:bold;
}
.center {
	text-align:center;
}
/* ------------------------------------------
mostly based on tablesort css
---------------------------------------------*/
table {
	background:#fff;
	max-width:100%;
	border-spacing:0;
	margin:10px 0;
	border:1px solid #ddd;
	border-collapse:separate;
	*border-collapse:collapsed;
	-webkit-box-shadow:0 0 4px rgba(0,0,0,0.10);
	   -moz-box-shadow:0 0 4px rgba(0,0,0,0.10);
	        box-shadow:0 0 4px rgba(0,0,0,0.10);
}
table th,
	table td {
		padding:8px;
		line-height:18px;
		text-align:left;
		border-top:1px solid #ddd;
	}
	table th {
		background:#eee;
		background:-webkit-gradient(linear, left top, left bottom, from(#f6f6f6), to(#eee));
		background:-moz-linear-gradient(top, #f6f6f6, #eee);
		text-shadow:0 1px 0 #fff;
		font-weight:bold;
		vertical-align:bottom;
	}
	table td {
		vertical-align:top;
	}
	table thead:first-child tr th,
	table thead:first-child tr td {
		border-top:0;
	}
	table tbody + tbody {
		border-top:2px solid #ddd;
	}
	table th + th,
	table td + td,
	table th + td,
	table td + th {
		border-left:1px solid #ddd;
	}
	table thead:first-child tr:first-child th,
	table tbody:first-child tr:first-child th,
	table tbody:first-child tr:first-child td {
		border-top:0;
	}

/*tablesort specific styling*/
th.sort-header::-moz-selection { background:transparent; }
th.sort-header::selection      { background:transparent; }
table th.sort-header:before {
	content:'';
	float:right;
	margin-top:7px;
	border-width:0 4px 4px;
	border-style:solid;
	border-color:#404040 transparent;
	visibility:hidden;
}
table th.sort-header:hover:before {
	visibility:visible;
}
table th.sort-up:before,
table th.sort-down:before,
table th.sort-down:hover:before {
	visibility:visible;
	opacity:0.4;
}
table th.sort-up:before {
	border-bottom:none;
	border-width:4px 4px 0;
}
.heading {
	margin-top:90px;
}
.options {
	margin:10px 0 30px 15px;
}
.options h3 {
	display:block;
	padding-top:10px;
	margin-top:20px;
}
.options h3:first-child {
	border:none;
	margin-top:0;
}

@media print {
	@page { 
		size: auto;   /* auto is the initial value */ 
		/* this affects the margin in the printer settings */ 
		margin: 1.5cm;  
	} 
	body {
		margin:0;
		font-family:'Helvetica Neue', Helvetica, Arial, sans-serif;
		font-size:10pt;
		line-height:15pt;
		color:#000000;
		background-color:#ffffff;
		-webkit-font-smoothing:antialiased;
	}
	a { page-break-inside:avoid }
	blockquote { page-break-inside: avoid; }
	h1, h2, h3, h4, h5, h6 {
		break-after: avoid-page;
		page-break-after:avoid; 
		page-break-inside:avoid
	}
	img {
		page-break-inside:avoid; 
		page-break-after:avoid;
	}
	table, pre { page-break-inside:avoid }
	ul, ol, dl  { page-break-before:avoid }
	a:link, a:visited, a {
		background: transparent;
		color: #520;
		font-weight: bold;
		text-decoration: underline;
		text-align: left;
	}
	a { page-break-inside:avoid }
	a[href^=http]:after { content:" < " attr(href) "> "; }
	$a:after > img { content: "";}
	article a[href^="#"]:after { content: ""; }
	a:not(:local-link):after { content:" < " attr(href) "> "; }
	/* Avoid chart divs being cutoff by having a pagebreak in them. */
	.nobreak {
		page-break-inside: avoid;
	}
}

/*-----------------------------------
  Markup free clearing
  Details: http: //perishablepress.com/press/2009/12/06/new-clearfix-hack
-------------------------------------*/
.clearfix:after {
	content: '.';
	display: block;
	height: 0;
	clear: both;
	visibility: hidden;
}
* html .clearfix { height: 1%; } /* IE6 */
*:first-child + html .clearfix { min-height: 1%; } /* IE7 */
</style>
}
}

proc reportscombine_singlecell {reportdirs dataVar} {
#putsvars reportdirs
#error stop
	upvar $dataVar data
	
	set html "\n<h2>Single cell overview</h2>\n"
	# singlecell summary
	set fields {
		nrcells
		mean_readcounts
		median_readcounts
		mean_genes_percell
		median_genes_percell
		mean_isoforms_percell
		median_isoforms_percell
		total_reads
		barcoded_reads
		pct_barcoded_reads
		validbarcoded_reads
		pct_validbarcoded_reads
		total_umis
		barcoded_umis
		pct_barcoded_umis
		validbarcoded_umis
		pct_validbarcoded_umis
		rawgenecount
		pct_rawgenecount
		filteredgenecount
		pct_filteredgenecount
		filteredisoformcount
		pct_filteredisoformcount
	}
	set table [list [list sample {*}$fields]]
	foreach dir $reportdirs {
		set sample [file tail [file dir $dir]]
		set file [glob -nocomplain $dir/report_singlecell-$sample.tsv]
		if {![file exists $file]} continue
		unset -nocomplain a
		array set a [split [string trim [cg select -f {parameter value} $file]] \n\t]
		set line {}
		lappend line $sample
		foreach field $fields {
			set value [get a(sc_$field) ?]
			if {[regexp \\. $value]} {
				set value [format %.2f $value]
			}
			lappend line $value
		}
		lappend table $line
	}
	append html [report_htmltable singlecellinfo $table]
	append html "\n<h2>Single cell celltyping overview (nrcells)</h2>\n"
	#
	# celltyping
	unset -nocomplain samplea
	unset -nocomplain groupa
	unset -nocomplain sourcea
	unset -nocomplain a
	foreach dir $reportdirs {
		set sample [file tail [file dir $dir]]
		set file [glob -nocomplain $dir/report_singlecell-cellgrouping-$sample.tsv]
		if {![file exists $file]} continue
		foreach {sample source group value} [split [string trim [cg select -sh /dev/null -f {sample source group value} $file]] \n\t] {
			set a($sample,$source,$group) $value
			set samplea($sample) 1
			set groupa($group) 1
			set sourcea($source) 1
		}
	}
	set groups [bsort [array names groupa]]
	set samples [bsort [array names samplea]]
	set sources [bsort [array names sourcea]]
	set table [list [list source sample {*}$groups]]
	foreach source $sources {
		foreach sample $samples {
			set line [list $source $sample]
			foreach group $groups {
				lappend line [get a($sample,$source,$group) 0]
			}
			lappend table $line
		}
	}
	append html [report_htmltable celltyping $table]
	#
	# isoform pct_covered
	set tabledata {}
	foreach dir $reportdirs {
		set sampledir [file dir $dir]
		set sample [file tail $sampledir]
		set file [glob -nocomplain $sampledir/reports/singlecell-isoform_pctcovered-$sample.tsv]
		if {![file exists $file]} continue
		set xs {}
		set ys {}
		set f [gzopen $file]
		set header [tsv_open $f]
		while {[gets $f line] != -1} {
			foreach {x y} [split $line \t] break
			if {$x == 0} continue
			lappend xs $x
			lappend ys $y
		}
		gzclose $f
		lappend tabledata $sample $xs $ys
	}
	append html [plotly isoform_pct_covered $tabledata "Read coverage of isoforms" "pct of isoform covered" "nr of umi corrected reads"]\n
	#
	# isoforms
	unset -nocomplain xsa
	unset -nocomplain ysa
	foreach dir $reportdirs {
		set sample [file tail [file dir $dir]]
		set file [glob -nocomplain $dir/singlecell-pseudobulk_isoforms-$sample.tsv]
		if {![file exists $file]} continue
		foreach {sample source group nrisoforms count} [split [string trim [cg select -sh /dev/null -f {sample source group nrisoforms count} $file]] \n\t] {
			set element $source-$group-$sample
			lappend xsa($element) $nrisoforms
			lappend ysa($element) $count
		}
	}
	set tabledata {}
	foreach element [bsort [array names xsa]] {
		lappend tabledata $element $xsa($element) $ysa($element)
	}
	append html [plotly nrisoforms $tabledata "Number of genes with given nr of isoforms" "nr of isoforms" "nr of genes"]\n
	return $html
}

proc reportscombine_table_yield {samples dataVar} {
	upvar $dataVar data
	set html "\n<h2>Yield and quality overview</h2>\n"
	set table [list {
		sample numreads numbases-mb pct-unique-reads
		fw-qual-mean fw-qual-stdev rev-qual-mean rev-qual-stdev
	}]
	foreach sample $samples {
		set line {}
		lappend line $sample
		set numreads [pgetsum data $sample,fastq-stats,fw_numreads $sample,fastq-stats,rev_numreads]
		lappend line $numreads
		set data($sample,numbases) [pgetsum data $sample,fastq-stats,fw_total_bases $sample,fastq-stats,rev_total_bases]
		lappend line [catchdef {expr {$data($sample,numbases)/1000000.0}} ""]
		lappend line [pget data $sample,pct_pf_unique_reads]
		foreach field {
			fw_qual_mean fw_qual_stdev rev_qual_mean rev_qual_stdev
		} {
			lappend line [pget data $sample,fastq-stats,$field]
		}
		lappend table $line
	}
	append html [report_htmltable yieldandquality $table]
	return $html
}

proc reportscombine_table_composition {samples dataVar} {
	upvar $dataVar data
	set html "\n<h2>Sequence composition</h2>\n"
	set table [list {
		sample 
		fw-pct-A rev-pct-A fw-pct-C rev-pct-C 
		fw-pct-G rev-pct-G fw-pct-T rev-pct-T fw-pct-N rev-pct-N
	}]
	foreach sample $samples {
		set line {}
		lappend line $sample
		foreach field {
			fw_pct_A rev_pct_A fw_pct_C rev_pct_C 
			fw_pct_G rev_pct_G fw_pct_T rev_pct_T fw_pct_N rev_pct_N
		} {
			lappend line [pget data $sample,fastq-stats,$field]
		}
		lappend table $line
	}
	append html [report_htmltable composition $table]
	return $html
}

proc reportscombine_table_alignment {alignments dataVar dbdir} {
	upvar $dataVar data
	set html "\n<h2>Alignment overview</h2>\n"
	set table [list {
		sample alignment numalignments numreads mapped-reads pct_mapped-reads duplicate-reads properly-paired-reads
		num-mapped-bases-mb pct-mapped-bases pct-mapped-ontarget
	}]
	foreach alignment $alignments {
		set sample [lindex $alignment end]
		set line [list $sample [join [lrange $alignment 0 end-1] -]]
		set alignment [join $alignment -]
		lappend line [pget data "$alignment,flagstat_reads,in total"]
		set numreads [pget data "$alignment,flagstat_reads,primary" "$alignment,samstats_summary,raw_total_sequences"]
		lappend line $numreads
		set mapped_reads [pget data "$alignment,flagstat_reads,primary mapped" "$alignment,samstats_summary,reads_mapped"]
		lappend line $mapped_reads
		lappend line [ppercent $mapped_reads $numreads]
		lappend line [pget data $alignment,flagstat_reads,duplicates]
		lappend line [pget data "$alignment,flagstat_reads,properly paired"]
		# mapping
		set numbases [pget data $sample,numbases]
		set mappedbases_cigar [pget data $alignment,samstats_summary,bases_mapped_cigar]
		set mappedbases [pget data $alignment,samstats_summary,bases_mapped]
		if {$mappedbases eq ""} {
			set mappedbases [pgetsum data $alignment,numbases_offtarget $alignment,numbases_ontarget]
		}
		set numbases_ontarget [pget data $alignment,numbases_ontarget]
		lappend line [catchexpr {$mappedbases/1000000.0}]
		lappend line [ppercent $mappedbases_cigar $numbases]
		lappend line [ppercent $numbases_ontarget $mappedbases]
		lappend table $line
	}
	append html [report_htmltable alignmentoverview $table]\n
	append html "\n<h2>Alignment target overview</h2>\n"
	set table [list {
		sample alignment
		target-region-mb avg-depth-ontarget pct-target-1X pct-target-2X pct-target-10X pct-target-20X pct-target-30X
	}]
	set notargetregs {}
	foreach alignment $alignments {
		set sample [lindex $alignment end]
		set line [list_reverse $alignment]
		set alignment [join $alignment -]
		set targetbases [pget data $alignment,histodepth,targetbases]
		if {$targetbases eq ""} {
			# there was no target defined, we'll take sequencedgenome as target
			if {![info exists tfile_targetbases]} {
				set tfile [gzfile $dbdir/extra/reg_*_sequencedgenome.tsv]
				if {[file exists $tfile]} {
					set tfile_targetbases [lindex [cg covered $tfile] end]
				} else {
					set tfile_targetbases ""
				}
			}
			set targetbases $tfile_targetbases
			if {$targetbases eq ""} {
				set targetbases 0
				lappend line 0
				lappend line {}
				foreach num {1 2 10 20 30} {lappend line {}}
			} else {
				lappend notargetregs $alignment
				# with no target defined, all bases are indicated offtarget
				set numbases_offtarget [pget data $alignment,numbases_offtarget]
				lappend line [expr {$targetbases/1000000.0}]
				set avg_target_depth [catchexpr {double($numbases_offtarget)/$targetbases}]
				lappend line [catchformat %.2f $avg_target_depth]
				foreach num {1 2 10 20 30} {
					set value [pget data $alignment,histodepth,offtarget_bases_${num}X]
					if {$value ne ""} {
						set value [catchformat %.4f [catchexpr {$value*100.0/$targetbases}]]
					}
					lappend line $value
				}
			}
		} else {
			lappend line [catchexpr {$targetbases/1000000.0}]
			set avg_target_depth [catchexpr {double($numbases_ontarget)/$targetbases}]
			lappend line [catchformat %.2f $avg_target_depth]
			foreach num {1 2 10 20 30} {
				set value [pget data $alignment,histodepth,pct_target_bases_${num}X]
				if {$value eq "" && $targetbases != 0} {
					set value [pget data $alignment,histodepth,ontarget_bases_${num}X]
					if {$value ne ""} {
						set value [catchformat %.4f [catchexpr {$value*100.0/$targetbases}]]
					}
				}
				lappend line $value
			}
		}
		lappend table $line
	}
	if {[llength $notargetregs]} {
		append html "no target region was supplied for [join $notargetregs ,].\n\n"
		append html "The full \"sequenced genome\" region in $tfile was used as target region for calulations in the table\n"
	}
	append html [report_htmltable alignmenttarget $table]
}

proc reportscombine_table_variant {alignments dataVar vcallers} {
	upvar $dataVar data
	set table [list {
		sample alignment varcaller vars qvars qvars-refcoding qvars-zyg-m qvars-zyg-t qvars-TiTv-ratio
	}]
	foreach vcaller $vcallers {
		set sample [lindex $vcaller end]
		set line [list_reverse $vcaller]
		set vcaller [join $vcaller -]
		foreach field {vars qvars qvars_refcoding qvars_zyg_m qvars_zyg_t} {
			lappend line [pget data $vcaller,genomecomb,$field]
		}
		set titv_i [pget data $vcaller,genomecomb,qvars_titv_i]
		set snps [pget data $vcaller,genomecomb,qvars_type_snp]
		if {![isint $titv_i] || ![isint $snps]} {
			lappend line ""
		} else {
			set titv [catchexpr {double($titv_i)/($snps-$titv_i)}]
			lappend line [catchformat %.2f $titv]
		}
		lappend table $line
	}
	set html "\n<h2>Variants overview</h2>\n"
	append html [report_htmltable alignmentoverview $table]\n
}

proc reportscombine_depth {histofiles depthaVar} {
	upvar $depthaVar deptha
	# depth histo chart
	# -----------------
	set html {}
	set depthdata {}
	foreach name [bsort [array names deptha]] {
		set xs [lrange [list_subindex $deptha($name) 0] 1 end]
		set ys [lrange [list_subindex $deptha($name) 1] 1 end]
		lappend depthdata $name $xs $ys
	}
	# make chart
	append html [plotly depthontarget $depthdata "Number of on-target bases with given depth" "depth" "number of positions" 100]\n
	# offtarget depth histo chart
	# ---------------------------
	set depthdata {}
	foreach name [bsort [array names deptha]] {
		set xs [lrange [list_subindex $deptha($name) 0] 1 end]
		set ys [lrange [list_subindex $deptha($name) 2] 1 end]
		lappend depthdata $name $xs $ys
	}
	# make chart
	append html [plotly depthofftarget $depthdata "Number of off-target bases with given depth" "depth" "number of positions" 100]\n
}


proc report_fastc_perposqual {dir files} {
	set qdata {}
	set maxx 0
	set maxy 40
	set colors [plotly_colors [llength $files]]
	foreach file $files color $colors {
		set name [file_analysis [file dir $file]]
		set table [fastqc_readtable $file {Per base sequence quality} header]
		if {$header ne {Base Mean Median {Lower Quartile} {Upper Quartile} {10th Percentile} {90th Percentile}}} {
			continue
		}
		set base {}
		foreach value [list_subindex $table [list_find $header Base]] {
			set value [lindex [split $value -] 0]
			if {$value > $maxx} {set maxx $value}
			lappend base $value
		}
		set medians [list_subindex $table [list_find $header Median]]
		if {![catch {
			set max [lmath_max $medians]
		}]} {
			if {$max > $maxy} {set maxy $max}
		}
		set p10 [list_subindex $table [list_find $header {10th Percentile}]]
		if {![catch {
			set max [lmath_max $p10]
		}]} {
			if {$max > $maxy} {set maxy $max}
		}
		lappend qdata [plotly_element $name line \
			line \{[subst {color: $color}]\} \
			x \[[join $base ,]\] y \[[join $medians ,]\]
		]
		lappend qdata [plotly_element ${name}_p10 line \
			line \{[subst {dash: 'dot', color: $color}]\} \
			x \[[join $base ,]\] y \[[join $p10 ,]\]
		]
	}
	if {![llength $qdata]} {return {}}
	set result [string_change [deindent {
		<div class = "nobreak">
		<h2>Fastqc <DIR> read quality per position (median and 10th percentile)</h2>
		<div id="<DIVID>" class="nobreak" style="height: 400px; width: 100%"></div>
		<script>
		var chart_<DIVID> = Plotly.plot(
		document.getElementById('<DIVID>'),
		[
		<DATA>
		],{
			xaxis: {
				title: 'read position'
			},
			yaxis: {
				title: 'quality'
			},
			margin: { t: 0 },
			paper_bgcolor: 'rgba(0,0,0,0)',
			plot_bgcolor: 'rgba(0,0,0,0)',
			shapes: [{
				type: 'rect',
				xref: 'x', yref: 'y',
				x0: 0, y0: 0, x1: <MAXX>, y1: 20,
				fillcolor: 'red', opacity: 0.2,
				line: {width: 0}				
			},{
				type: 'rect',
				xref: 'x', yref: 'y',
				x0: 0, y0: 20, x1: <MAXX>, y1: 28,
				fillcolor: 'orange', opacity: 0.2,
				line: {width: 0}				
			},{
				type: 'rect',
				xref: 'x', yref: 'y',
				x0: 0, y0: 28, x1: <MAXX>, y1: 40,
				fillcolor: 'green', opacity: 0.2,
				line: {width: 0}				
			}]
		},
		{responsive: true,
		  modeBarButtonsToAdd: [{
		    name: 'tosvg',
		    icon: Plotly.Icons.camera,
		    click: function(gd) {
		      Plotly.downloadImage(gd, {format: 'svg'})
		    }
		  }]
		}
		);
		</script>
		</div>
	}] [list \
		<DIR> $dir \
		<DIVID> fastqc_${dir}_perposqual \
		<DATA> [join $qdata ,\n] \
		<MAXX> $maxx \
	]]
	return $result
}

proc report_fastc_perseqqual {dir files} {
	set qdata {}
	set maxx 0
	set maxy 40
	set colors [plotly_colors [llength $files]]
	foreach file $files color $colors {
		set name [file_analysis [file dir $file]]
		set table [fastqc_readtable $file {Per sequence quality scores} header]
		set xpos [lsearch $header Quality]
		set ypos [lsearch $header Count]
		set xs [list_subindex $table $xpos]
		set ys [list_subindex $table $ypos]
		lappend qdata [plotly_element $name line \
			line \{[subst {color: $color}]\} \
			x \[[join $xs ,]\] y \[[join $ys ,]\]
		]
	}
	if {![llength $qdata]} {return {}}
	set result [string_change [deindent {
		<div class = "nobreak">
		<h2>Fastqc <DIR> reads per quality</h2>
		<div id="<DIVID>" class="nobreak" style="height: 400px; width: 100%"></div>
		<script>
		var chart_<DIVID> = Plotly.plot(
		document.getElementById('<DIVID>'),
		[
		<DATA>
		],{
			xaxis: {
				title: 'quality'
			},
			yaxis: {
				title: 'number of reads'
			},
			margin: { t: 0 },
			paper_bgcolor: 'rgba(0,0,0,0)',
			plot_bgcolor: 'rgba(0,0,0,0)',
		},
		{responsive: true,
		  modeBarButtonsToAdd: [{
		    name: 'tosvg',
		    icon: Plotly.Icons.camera,
		    click: function(gd) {
		      Plotly.downloadImage(gd, {format: 'svg'})
		    }
		  }]
		}
		);
		</script>
		</div>
	}] [list \
		<DIVID> fastqc_${dir}_perseqqual \
		<DATA> [join $qdata ,\n] \
		<DIR> $dir \
	]]
	return $result
}

proc report_fastc_chart {dir files id title pattern xfield yfield} {
	set qdata {}
	set maxx 0
	set maxy 40
	set colors [plotly_colors [llength $files]]
	foreach file $files color $colors {
		set name [file_analysis [file dir $file]]
		set table [fastqc_readtable $file $pattern header]
		set xpos [lsearch $header $xfield]
		set ypos [lsearch $header $yfield]
		if {$xpos == -1 || $ypos == -1} continue
		set xs [list_subindex $table $xpos]
		if {![isdouble [lindex $xs]]} {
			set xs [list_regsub -- {-.*} $xs {}]
		}
		set ys [list_subindex $table $ypos]
		lappend qdata [plotly_element $name line \
			line \{[subst {color: $color}]\} \
			x \[[join $xs ,]\] y \[[join $ys ,]\]
		]
	}
	if {![llength $qdata]} {return {}}
	set divid fastqc_${dir}_$id
	set data [join $qdata ,\n]
	subst -nocommands [deindent {
		<div class = "nobreak">
		<h2>$title</h2>
		<div id="$divid" class="nobreak" style="height: 400px; width: 100%"></div>
		<script>
		var chart_$divid = Plotly.plot(
		document.getElementById('$divid'),
		[
		$data
		],{
			xaxis: {
				title: '$xfield'
			},
			yaxis: {
				title: '$yfield'
			},
			margin: { t: 0 },
			paper_bgcolor: 'rgba(0,0,0,0)',
			plot_bgcolor: 'rgba(0,0,0,0)',
		},
		{responsive: true,
		  modeBarButtonsToAdd: [{
		    name: 'tosvg',
		    icon: Plotly.Icons.camera,
		    click: function(gd) {
		      Plotly.downloadImage(gd, {format: 'svg'})
		    }
		  }]
		}
		);
		</script>
		</div>
	}]
}

proc reportscombine_fastqc {fastqcfiles} {
	#
	# gather fastqc data
	# ------------------
	set fastqctable {}
	lappend fastqctable [list \
		sample dir gc per-base-quality per-sequence-quality \
		per-base-content per-sequence-gc per-base-N \
		length-distribution sequence-duplication overrepresented-seq \
	]
	unset -nocomplain fastqcfilesa
	foreach file $fastqcfiles {
		set name [file tail [file dir $file]]
		if {[regexp ^fastqc_fw $name]} {
			set dir forward
		} elseif {[regexp ^fastqc_rev $name]} {
			set dir reverse
		} else {
			set dir {}
		}
		lappend fastqcfilesa($dir) $file
		set c [file_read $file]
		set sample [file_sample [file dir $file]]
		set line [list $sample $dir]
		set pct-gc {}
		lappend line [report_getpattern {%GC	([0-9.]+)} $c]
		lappend line [report_getpattern {>>Per base sequence quality	([a-z]+)} $c]
		lappend line [report_getpattern {>>Per sequence quality scores	([a-z]+)} $c]
		lappend line [report_getpattern {>>Per base sequence content	([a-z]+)} $c]
		lappend line [report_getpattern {>>Per sequence GC content	([a-z]+)} $c]
		lappend line [report_getpattern {>>Per base N content	([a-z]+)} $c]
		lappend line [report_getpattern {>>Sequence Length Distribution	([a-z]+)} $c]
		lappend line [report_getpattern {>>Sequence Duplication Levels	([a-z]+)} $c]
		lappend line [report_getpattern {>>Overrepresented sequences	([a-z]+)} $c]
		lappend fastqctable $line
	}
	set html {}
	#
	# fastqc overview
	# ---------------
	append html "\n<h2>fastqc overview</h2>\n"
	append html [report_htmltable fastqcoverview $fastqctable]\n
	append html </script>\n</div>\n
	#
	# fastqc charts
	# -------------
	#
	# fastqc quality per pos
	# ----------------------
	foreach dir [lsort [array names fastqcfilesa]] {
		append html [report_fastc_perposqual $dir $fastqcfilesa($dir)]\n
	}
	#
	# fastqc per seq quality
	# ----------------------
	foreach dir [lsort [array names fastqcfilesa]] {
		# set html [report_fastc_perseqqual $dir $fastqcfilesa($dir)]
		append html [report_fastc_chart $dir $fastqcfilesa($dir) \
			perseqqual \
			"Fastqc $dir reads per quality" \
			{Per sequence quality scores} \
			Quality Count \
		]\n
	}
	#
	# fastqc per seq gc content
	# -------------------------
	foreach dir [lsort [array names fastqcfilesa]] {
		append html [report_fastc_chart $dir $fastqcfilesa($dir) \
			perseqgc \
			"Fastqc $dir reads GC content" \
			{Per sequence GC content} \
			{GC Content} Count \
		]\n
	}
	#
	# sequence length
	# ---------------
	foreach dir [lsort [array names fastqcfilesa]] {
		append html [report_fastc_chart $dir $fastqcfilesa($dir) \
			seqlen \
			"Fastqc $dir reads length distribution" \
			{Sequence Length Distribution} \
			{Length} Count \
		]\n
	}
	return $html
}

proc reportscombine_samstats_quality {qualityfiles {name qualityffq} {title "Quality distribution"}} {
	if {![llength $qualityfiles]} {return {}}
	set chartdata {}
	set xmax 30
	foreach file $qualityfiles {
		if {![file exists $file]} continue
		set tail [file tail $file]
		set prefix {}
		if {[regsub ^unaligned $tail {} tail]} {
			append prefix unaligned
		} elseif {[regsub ^aligned $tail {} tail]} {
			append prefix aligned
		}
		if {[regexp ^samstats_FFQ $tail]} {
			append prefix fw-
		} elseif {[regexp ^samstats_LFQ $tail]} {
			append prefix rv-
		}
		set analysis $prefix[file_analysis $file]
		set f [gzopen $file]
		set header [tsv_open $f]
		if {[eof $f]} continue
		set xs {}
		set ys {}
		while {[gets $f line] != -1} {
			foreach {x y} [split $line \t] break
			if {$y > 10 && $x > $xmax && $x <= 100} {set xmax $x}
			lappend xs $x
			lappend ys $y
		}
		close $f
		lappend chartdata $analysis $xs $ys
		
	}
	#
	# make chart
	plotly $name $chartdata $title "quality" "number of bases" $xmax
}

proc reportscombine_samstats_gc {rfiles {name samstatsgc} {title "GC distribution"}} {
	if {![llength $rfiles]} {return {}}
	set chartdata {}
	foreach file $rfiles {
		if {![file exists $file]} continue
		set tail [file tail $file]
		set prefix {}
		if {[regsub ^unaligned $tail {} tail]} {
			append prefix unaligned
		} elseif {[regsub ^aligned $tail {} tail]} {
			append prefix aligned
		}
		if {[regexp ^samstats_GCF $tail]} {
			append prefix fw-
		} elseif {[regexp ^samstats_GCL $tail]} {
			append prefix rv-
		}
		set analysis $prefix[file_analysis $file]
		set f [gzopen $file]
		set header [tsv_open $f]
		if {[eof $f]} continue
		set xs {}
		set ys {}
		while {[gets $f line] != -1} {
			foreach {x y} [split $line \t] break
			lappend xs $x
			lappend ys $y
		}
		close $f
		lappend chartdata $analysis $xs $ys
		
	}
	#
	# make chart
	plotly $name $chartdata $title "% GC" "number of reads" 100
}

proc reportscombine_samstats_readlength {statsrlfiles} {
	if {![llength $statsrlfiles]} {return ""}
	set chartdata {}
	set xmax 50
	foreach file $statsrlfiles {
		if {![file exists $file]} continue
		set analysis [file_analysis $file]
		set f [gzopen $file]
		set header [tsv_open $f]
		set c [split [string trim [read $f]] \n]
		close $f
		set poss [list_cor $header {read_length count}]
		set xs {}
		set ys {}
		foreach line $c {
			set line [list_sub [split $line \t] $poss]
			foreach {x y} $line break
			if {$x > $xmax && $y > 10} {set xmax $x}
			lappend xs $x
			lappend ys $y
		}
		lappend chartdata $analysis $xs $ys
		
	}
	#
	# make chart
	plotly readlength $chartdata "Readlength distribution" "readlength" "number of sequences" $xmax
}

proc reportscombine_samstats_uareadlength {uastatsrlfiles} {
	if {![llength $uastatsrlfiles]} {return ""}
	set chartdata {}
	set xmax 50
	foreach file $uastatsrlfiles {
		if {![file exists $file]} continue
		set prefix {}
		regsub samstats [lindex [split [file tail $file] _] 0] {} prefix
		set analysis $prefix-[file_analysis $file]
		set f [gzopen $file]
		set header [tsv_open $f]
		set c [split [string trim [read $f]] \n]
		close $f
		set poss [list_cor $header {read_length count}]
		set xs {}
		set ys {}
		foreach line $c {
			set line [list_sub [split $line \t] $poss]
			foreach {x y} $line break
			if {$x > $xmax && $y > 10} {set xmax $x}
			lappend xs $x
			lappend ys $y
		}
		lappend chartdata $analysis $xs $ys
		
	}
	#
	# make chart
	plotly uareadlength $chartdata "Readlength distribution (aligned and unaligned)" "readlength" "number of sequences" $xmax
}

proc reportscombine_remove_hlongshot {files} {
	unset -nocomplain a
	foreach file $files {
		set a([file_rootname $file]) $file
	}
	foreach name [array names a hlongshot-*] {
		if {[info exists a([string range $name 10 end])]} {
			unset a($name)
		}
	}
	set result {}
	foreach name [bsort [array names a]] {
		lappend result $a($name)
	}
	return $result
}

proc process_reportscombine_job {args} {
	upvar job_logdir job_logdir
	set overwrite 0
	set dbdir {}
	cg_options process_reportscombine args {
		-experimentname {
			set experimentname $value
		}
		-overwrite {
			set overwrite 1
		}
		-dbdir {
			set dbdir $value
		}
	} {destdir reportsdir} 2 ... {
		combines reports directories (per sample) into one overview reports directory
	}
	catch {set dbdir [dbdir $dbdir]}
	set destdir [file_absolute $destdir]
	set reportstodo [list $reportsdir {*}$args]
	if {![info exists job_logdir]} {
		set_job_logdir $destdir/log_jobs
	}
	if {![info exists experimentname]} {
		set experimentname [file tail $destdir]
		if {$experimentname eq "reports"} {
			set experimentname [file tail [file dir $destdir]]
		}
	}
	#
	# combine hsmetrics (if found)
	set deps {}
	foreach dir $reportstodo {
		lappend deps {*}[jobglob -checkcompressed 1 $dir/hsmetrics-*.hsmetrics $dir/reports/hsmetrics-*.hsmetrics]
		if {[file extension $dir] eq ".hsmetrics"} {
			lappend deps {*}[jobglob $dir]
		}
	}
	if {[llength $deps]} {
		set target $destdir/report_hsmetrics-${experimentname}.tsv
		set target2 $destdir/${experimentname}_hsmetrics_report.tsv
		job reportscombine_hsmetrics-$experimentname -deps $deps -targets {$target $target2} -code {
			file mkdir [file dir $target]
			cg cat -m 1 -c 0 {*}$deps > $target.temp
			cg select -rc 1 $target.temp $target.temp2
			file rename -force -- $target.temp2 $target
			file delete $target.temp
			mklink $target $target2
		}
	}

	# combine other reports
	set reportdirs {}
	set reports {}
	set histofiles {}
	set fastqcfiles {}
	set statsrlfiles {}
	set uastatsrlfiles {}
	set statsqualityfiles {}
	set uastatsqualityfiles {}
	set statsgcfiles {}
	set uastatsgcfiles {}
	set singlecellfiles {}
	foreach dir $reportstodo {
		if {[file extension [gzroot $dir]] eq ".txt"} {
			if {[regexp ^fastqc $dir]} {
				lappend fastqcfiles $dir
			}
		} elseif {[file extension [gzroot $dir]] eq ".tsv"} {
			if {[regexp ^histodepth- $dir]} {
				lappend histofiles $dir
			} else {
				lappend reports $dir
			}
			lappend reportdirs [file dir $dir]
		} else {
			set temp [jobglob -checkcompressed 1 $dir/reports/report_*.tsv]
			if {[llength $temp]} {
				lappend reportdirs $dir/reports
				lappend reports {*}$temp
			} else {
				lappend reportdirs $dir
				lappend reports {*}[jobglob -checkcompressed 1 $dir/report_*.tsv]
			}
			lappend histofiles {*}[jobglob -checkcompressed 1 $dir/histodepth-*.tsv $dir/reports/histodepth-*.tsv]
			lappend fastqcfiles {*}[jobglob -checkcompressed 1 $dir/fastqc_*.fastqc/fastqc_data.txt $dir/reports/fastqc_*.fastqc/fastqc_data.txt]
			lappend statsrlfiles {*}[jobglob -checkcompressed 1 $dir/samstats_RL-*.tsv.zst]
			lappend uastatsrlfiles {*}[jobglob -checkcompressed 1 $dir/alignedsamstats_RL-*.tsv.zst $dir/unalignedsamstats_RL-*.tsv.zst]
			lappend statsqualityfiles {*}[jobglob -checkcompressed 1 $dir/samstats_FFQs-*.tsv.zst $dir/samstats_LFQs-*.tsv.zst]
			lappend uastatsqualityfiles {*}[jobglob -checkcompressed 1 \
				$dir/alignedsamstats_FFQs-*.tsv.zst $dir/unalignedsamstats_FFQs-*.tsv.zst \
				$dir/alignedsamstats_LFQs-*.tsv.zst $dir/unalignedsamstats_LFQs-*.tsv.zst \
			]
			lappend statsqualityfiles {*}[jobglob -checkcompressed 1 $dir/samstats_FFQs-*.tsv.zst $dir/samstats_LFQs-*.tsv.zst]
			lappend uastatsqualityfiles {*}[jobglob -checkcompressed 1 \
				$dir/alignedsamstats_FFQs-*.tsv.zst $dir/unalignedsamstats_FFQs-*.tsv.zst \
				$dir/alignedsamstats_LFQs-*.tsv.zst $dir/unalignedsamstats_LFQs-*.tsv.zst \
			]
			lappend statsgcfiles {*}[jobglob -checkcompressed 1 $dir/samstats_GCF-*.tsv.zst $dir/samstats_GCL-*.tsv.zst]
			lappend uastatsgcfiles {*}[jobglob -checkcompressed 1 \
				$dir/alignedsamstats_GCF-*.tsv.zst $dir/unalignedsamstats_GCF-*.tsv.zst \
				$dir/alignedsamstats_GCL-*.tsv.zst $dir/unalignedsamstats_GCL-*.tsv.zst \
			]
			lappend singlecellfiles {*}[jobglob -checkcompressed 1 \
				$dir/*singlecell*.tsv \
			]
		}
	}
	set reports [bsort [list_remdup $reports]]
	set histofiles [reportscombine_remove_hlongshot $histofiles]
	set fastqcfiles [bsort [list_remdup $fastqcfiles]]
	set statsrlfiles [bsort [list_remdup $statsrlfiles]]
	set uastatsrlfiles [bsort [list_remdup $uastatsrlfiles]]
	set statsqualityfiles [bsort [list_remdup $statsqualityfiles]]
	set uastatsqualityfiles [bsort [list_remdup $uastatsqualityfiles]]
	set statsgcfiles [bsort [list_remdup $statsgcfiles]]
	set uastatsgcfiles [bsort [list_remdup $uastatsgcfiles]]
	set singlecellfiles [bsort [list_remdup $singlecellfiles]]
	set deps [list_concat $reports $histofiles $fastqcfiles $statsrlfiles $uastatsrlfiles $statsqualityfiles $statsgcfiles $uastatsgcfiles $singlecellfiles]

	if {[llength $deps]} {
		set target $destdir/report_stats-${experimentname}.tsv
		set target2 $destdir/report_summarytable-${experimentname}.tsv
		set target3 $destdir/report-${experimentname}.html
		if {$overwrite} {
			file delete $target $target2 $target3
		}
		job reportscombine_stats-$experimentname -deps $deps -vars {
			reportdirs reports histofiles fastqcfiles experimentname dbdir 
			statsrlfiles statsqualityfiles uastatsqualityfiles uastatsrlfiles statsgcfiles uastatsgcfiles singlecellfiles
		} -targets {
			$target $target2 $target3
		} -code {
			# make report_stats
			# =================
			file mkdir [file dir $target]
			cg cat -m -c 0 {*}[bsort $reports] > $target.temp
			cg select -rc 1 $target.temp $target.temp2
			file rename -force -- $target.temp2 $target
			file delete $target.temp
			# make report_summarytable
			unset -nocomplain data
			unset -nocomplain analysesa
			set f [open $target]
			set header [tsv_open $f]
			while {[gets $f line] != -1} {
				foreach {analysis source parameter value} [split $line \t] break
				set analysesa($analysis) 1
				set parameter [string trim $parameter]
				set data($analysis,$source,$parameter) $value
			}
			close $f
			#
			# make report_summarytable
			# ========================
			set o [open $target2.temp w]
			puts $o [join {
				analysis sample numreads 
				pf_reads pct_pf_reads pf_unique_reads pct_pf_unique_reads pf_mapped pct_pf_aligned_reads
				targetbases	pct_target_bases_2X pct_target_bases_10X pct_target_bases_20X pct_target_bases_30X
				covered_total qvars qvars_target qvars_refcoding
			} \t]
			set analyses [bsort [array names analysesa]]
			foreach analysis $analyses {
				# only report for "endpoints" analysis (data of partial will be included in endpoints)
				if {[regexp {^cov5-} $analysis]} continue
				if {[regexp {^hlongshot-} $analysis]} {
					if {[inlist $analyses [string range $analysis 10 end]]} continue
				}
				if {[llength [array names a *-$analysis]]} continue
				set sample [lindex [split $analysis -] end]
				set bamname [join [lrange [split $analysis -] end-1 end] -]
				set resultline [list $analysis $sample]
				# numreads
				set numreads {}
				set fw_numreads [get data($sample,fastq-stats,fw_numreads) {}]
				set rev_numreads [get data($sample,fastq-stats,rev_numreads) {}]
				if {[isint $fw_numreads] && [isint $rev_numreads]} {
					set numreads [expr {$fw_numreads+$rev_numreads}]
				}
				if {$numreads eq ""} {
					set numreads [get data($bamname,samstats_summary,raw_total_sequences) {}]
				}
				lappend resultline $numreads
				set data($sample,numreads) $numreads
				# bam based stats
				if {[info exists "data($bamname,flagstat_reads,primary)"]} {
					set pf_reads [get "data($bamname,flagstat_reads,primary)"]
				}  else {
					set pf_reads [expr {
						[get "data($bamname,flagstat_reads,in total)" -1] - [get "data($bamname,flagstat_reads,secondary)" 0] - [get "data($bamname,flagstat_reads,supplementary)" 0]
					}]
					if {$pf_reads < 0} {set pf_reads {}}
				}
				set pf_duplicates [get "data($bamname,flagstat_reads,duplicates)" {}]
				if {[info exists "data($bamname,flagstat_reads,primary mapped)"]} {
					set pf_mapped [get "data($bamname,flagstat_reads,primary mapped)"]
				}  else {
					set pf_mapped [expr {
						[get "data($bamname,flagstat_reads,mapped)" -1] - [get "data($bamname,flagstat_reads,secondary)" 0] - [get "data($bamname,flagstat_reads,supplementary)" 0]
					}]
				}
				if {$pf_mapped < 0} {set pf_mapped {}}
				set pf_properlypaired [get "data($bamname,flagstat_reads,properly paired)" {}]
				if {[isint $numreads] && [isint $pf_reads]} {
					set pct_pf_reads [ppercent $pf_reads $numreads]
					set pf_unique_reads [expr {($pf_reads - $pf_duplicates)}]
					set pct_pf_unique_reads [ppercent $pf_unique_reads $pf_reads]
					set pct_pf_aligned_reads [ppercent $pf_mapped $pf_reads]
				} else {
					set pct_pf_reads {}
					set pf_unique_reads {}
					set pct_pf_unique_reads {}
					set pct_pf_aligned_reads {}
				}
				lappend resultline $pf_reads $pct_pf_reads $pf_unique_reads $pct_pf_unique_reads $pf_mapped $pct_pf_aligned_reads
				# keep for further use
				if {$pf_reads ne ""} {set data($sample,pf_reads) $pf_reads}
				if {$pct_pf_reads ne ""} {set data($sample,pct_pf_reads) $pct_pf_reads}
				if {$pf_unique_reads ne ""} {set data($sample,pf_unique_reads) $pf_unique_reads}
				if {$pct_pf_unique_reads ne ""} {set data($sample,pct_pf_unique_reads) $pct_pf_unique_reads}
				if {$pf_mapped ne ""} {set data($sample,pf_mapped) $pf_mapped}
				if {$pct_pf_aligned_reads ne ""} {set data($sample,pct_pf_aligned_reads) $pct_pf_aligned_reads}
				# rest
				lappend resultline [get data($bamname,histodepth,targetbases) {}]
				lappend resultline [get data($bamname,histodepth,pct_target_bases_2X) {}]
				lappend resultline [get data($bamname,histodepth,pct_target_bases_10X) {}]
				lappend resultline [get data($bamname,histodepth,pct_target_bases_20X) {}]
				lappend resultline [get data($bamname,histodepth,pct_target_bases_30X) {}]
				lappend resultline [get data($analysis,genomecomb,covered_total) {}]
				lappend resultline [get data($analysis,genomecomb,qvars) {}]
				lappend resultline [get data($analysis,genomecomb,qvars_target) {}]
				lappend resultline [get data($analysis,genomecomb,qvars_refcoding) {}]
				puts $o [join $resultline \t]
			}
			close $o
			file rename -force -- $target2.temp $target2
			#
			# make html report
			# ================

			catch {close $o} ; set o [open $target3 w]
			puts $o [string_change [string trim [deindent {
				<!DOCTYPE html>
				<html>
				<head>
				<meta charset="UTF-8">
				<title><PAGETITLE></title>
				<script src="https://cdnjs.cloudflare.com/ajax/libs/tablesort/4.1.0/tablesort.min.js"></script>
				<script src="https://cdnjs.cloudflare.com/ajax/libs/tablesort/4.1.0/src/sorts/tablesort.number.js"></script>
				<script src="https://cdnjs.cloudflare.com/ajax/libs/tablesort/4.1.0/src/sorts/tablesort.dotsep.js"></script>
				<script src="https://cdnjs.cloudflare.com/ajax/libs/tablesort/4.1.0/src/sorts/tablesort.date.js"></script>
				<script src="https://cdnjs.cloudflare.com/ajax/libs/tablesort/4.1.0/src/sorts/tablesort.monthname.js"></script>
				<script src="https://cdnjs.cloudflare.com/ajax/libs/tablesort/4.1.0/src/sorts/tablesort.filesize.js"></script>
				<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
				<REPORT_CSS>
				<style>
					text {
						font: 12px sans-serif;
					}
					.js-plotly-plot .plotly .modebar {
						left: 50%; 
						transform: translate(-50%,-100%);
					}
				</style>
				</head>
				<body>
				<h1 class="title"><PAGETITLE></h1>
			}]] [list \
				<PAGETITLE> "NGS reports exp $experimentname" \
				<REPORT_CSS> [report_css]\
			]]
			#
			# gather analysis data
			# --------------------
			foreach file $deps {
				set expname ""
				set dirs [file split [file dir $file]]
				if {[lindex $dirs end] eq "reports"} {
					if {[lindex $dirs end-2] eq "samples"} {
						set expname [lindex $dirs end-3]
					} else {
						set expname [lindex $dirs end-2]
					}
				} elseif {[lindex $dirs end-1] eq "samples"} {
					set expname [lindex $dirs end-2]
				}
				set expa($expname) 1
			}
			set experiments [list_remove [array names expa] ""]
			set samples [list_regsub -all {.*-} $analyses {}]
			set samples [bsort [list_remdup $samples]]
			set alignments {}
			foreach temp [array names data *,histodepth,offtarget_bases_20X] {
				regsub ,histodepth,offtarget_bases_20X\$ $temp {} temp
				lappend alignments [split $temp -]
			}
			foreach temp [array names data {*,flagstat_reads,in total}] {
				regsub {,flagstat_reads,in total$} $temp {} temp
				lappend alignments [split $temp -]
			}
			foreach temp [array names data {*,samstats_summary,bases_mapped_cigar}] {
				regsub {,samstats_summary,bases_mapped_cigar$} $temp {} temp
				lappend alignments [split $temp -]
			}
			# remove hlongshot alignments if without exists
			set temp {}
			foreach alignment [bsort -index end [list_remdup $alignments]] {
				if {[lindex $alignment 0] eq "hlongshot"} {
					if {[inlist $analyses [join [lrange $alignment 1 end] -]]} continue
				}
				lappend temp $alignment
			}
			set alignments [list_remdup $temp]
			set vcallers {}
			foreach temp [array names data *,qvars] {
				regsub ,genomecomb,qvars\$ $temp {} temp2
				if {$temp2 eq {}} continue
				lappend vcallers [split $temp2 -]
			}
			set vcallers [bsort -index end [bsort [list_remdup $vcallers]]]
			# gather depth data
			# -----------------
			unset -nocomplain deptha
			foreach file [bsort $histofiles] {
				set analysis [file_analysis $file]
				set f [open $file]
				set header [tsv_open $f]
				set c [split [string trim [read $f]] \n]
				close $f
				set poss [list_cor $header {depth ontarget offtarget}]
				set depthtable [list {depth ontarget offtarget}]
				set numbases_offtarget 0
				set numbases_ontarget 0
				foreach line $c {
					set line [list_sub $line $poss]
					lappend depthtable $line
					foreach {depth ontarget offtarget} $line break
					incr numbases_offtarget [expr {$depth * $offtarget}]
					incr numbases_ontarget [expr {$depth * $ontarget}]
				}
				set deptha($analysis) $depthtable
				set data($analysis,numbases_offtarget) $numbases_offtarget
				set data($analysis,numbases_ontarget) $numbases_ontarget
				
			}
			#
			# intro
			# -----
			puts $o "\n<h2>Experiment info</h2>"
			if {[llength $experiments] > 1} {
				set exp [subst { of experiments [join $experiments ", "]}]
			} elseif {[llength $experiments] == 1} {
				set exp [subst { of experiment $experiments}]
			} else {
				set exp {}
			}
			puts $o [subst {
				This report was generated by genomecomb version [version genomecomb].
				It combines overview results on [llength $samples] samples$exp
			}]
			#
			# Yield and quality
			# -----------------
			puts $o [reportscombine_table_yield $samples data]
			# Sequence composition
			# --------------------
			# puts $o [reportscombine_table_composition $samples data]
			#
			# Alignment overview
			# ------------------
			#
			# alignment target table
			# ----------------------
			puts $o [reportscombine_table_alignment $alignments data $dbdir]
			#
			# singlecell info table
			# ---------------------
			if {[llength $singlecellfiles]} {
				puts $o [reportscombine_singlecell $reportdirs data]
			}
			#
			# Variant info table
			# ------------------
			puts $o [reportscombine_table_variant $alignments data $vcallers]
			#
			# gender prediction
			# ------------------
			set table [list {
				sample predicted-gender Y-over-X-ratio normalised-Y-over-X-ratio pct-heterozygous-Y
			}]
			foreach sample $samples {
				set line [list $sample]
				foreach field {predgender pg_yxratio pg_yxnratio pg_pctheterozygous} {
					lappend line [pget data $sample,genomecomb,$field]
				}
				if {![llength [list_remove [lrange $line 1 end] {}]]} continue
				lappend table $line
			}
			if {[llength $table] > 1} {
				puts $o "\n<h2>Gender prediction overview</h2>"
				puts $o [report_htmltable genderprediction $table]\n
			}
			# depth histo
			# -----------
			if {[llength $histofiles]} {
				puts $o [reportscombine_depth $histofiles deptha]\n
			}
			#
			# fastqc info
			# -----------
			if {[llength $fastqcfiles]} {
				puts $o [reportscombine_fastqc $fastqcfiles]
			}
			# samstats quality chart
			# ----------------------
			# gather seq size data
			if {[llength $statsqualityfiles]} {
				puts $o [reportscombine_samstats_quality $statsqualityfiles qualityffq "Quality distribution"]
			}
			# samstats ua quality chart
			# ----------------------
			# gather seq size data
			if {[llength $uastatsqualityfiles]} {
				puts $o [reportscombine_samstats_quality $uastatsqualityfiles uaqualityffq "Quality distribution (aligned and unaligned)"]
			}
			# readlength chart
			# --------------
			# gather seq size data
			if {[llength $statsrlfiles]} {
				puts $o [reportscombine_samstats_readlength $statsrlfiles]
			}
			# ua readlength chart
			# --------------
			# gather seq size data
			if {[llength $uastatsrlfiles]} {
				puts $o [reportscombine_samstats_uareadlength $uastatsrlfiles]
			}
			if {![llength $fastqcfiles]} {
				# samstats GC chart
				# ----------------------
				# gather seq size data
				if {[llength $statsgcfiles]} {
					puts $o [reportscombine_samstats_gc $statsgcfiles samstatsGC "GC distribution"]
				}
				# samstats ua GC chart
				# ----------------------
				# gather seq size data
				if {[llength $uastatsgcfiles]} {
					puts $o [reportscombine_samstats_gc $uastatsgcfiles uasamstatsGC "GC distribution (aligned and unaligned)"]
				}
			}
			#
			# footer
			# ------
			puts $o [string_change [string trim [deindent {
				<footer class="footer">
				<p>Made by genomecomb <GENOMECOMBVERSION></p>
				</footer>
				</body>
				</html>
			}]] [list <GENOMECOMBVERSION> [version genomecomb]]]
			close $o

			puts "Created $target3"
		}
	}
}

proc cg_process_reportscombine {args} {
	set args [job_init {*}$args]
	process_reportscombine_job {*}$args
	job_wait
}



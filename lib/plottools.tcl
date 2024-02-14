proc plotly_element {name type args} {
	set result "\{name: \"${name}\",\ntype: '$type'"
	foreach {key value} $args {
		append result ",\n$key: $value"
	}
	append result \n\}
	return $result
}

proc plotly_colors {num} {
	set list {'#1f77b4' '#ff7f0e' '#2ca02c' '#d62728' '#9467bd' '#8c564b' '#e377c2' '#7f7f7f' '#bcbd22' '#17becf'}
	set result $list
	while {$num > [llength $result]} {
		lappend result {*}$list
	}
	lrange $result 0 [expr {$num-1}]
}

# chartdata: elementname xs ys elementname2 xs2 ys2 ...
proc plotly {name chartdata title xaxis yaxis {xmax {}} {ymax {}}} {
	set resulthtml "\n<div class = \"nobreak\">\n<h2>$title</h2>\n"
	append resulthtml "<div id=\"$name\" class=\"nobreak\" style=\"height: 400px; width: 100%\"></div>\n"
	set fmtdata {}
	if {$xmax ne ""} {
		set xmax "range: \[0,$xmax\]"
	}
	if {$ymax ne ""} {
		set ymax "range: \[0,$ymax\]"
	}
	foreach {elementname xs ys} $chartdata {
		lappend fmtdata [plotly_element $elementname line \
			x \[[join $xs ,]\] y \[[join $ys ,]\]
		]
	}
	append resulthtml <script>\n
	append resulthtml [string_change [deindent {
		var chart_@NAME@ = Plotly.plot(
		document.getElementById('@NAME@'),
		[
		<DATA>
		],{
			xaxis: {
				title: '@XAXIS@',
				@XMAX@
			},
			yaxis: {
				title: '@YAXIS@'
				@YMAX@
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
	}] [list \
		@NAME@ $name <DATA> [join $fmtdata ,\n] @XAXIS@ $xaxis @YAXIS@ $yaxis @XMAX@ $xmax @YMAX@ $ymax \
	]]
	append resulthtml \n</script>\n</div>
	return $resulthtml
}

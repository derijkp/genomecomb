
proc cg_viz_transcripts {args} {
	set referencetranscript {}
	set countfields {}
	set countfieldpatterns ^count
	set rescale 1
	set textsize 10
	set panelwidths {3 1}
	set width 397
	set height 210
	set proportions 0
	set round_digits {}
	cg_options viz_transcripts args {
		-rescale {
			set rescale $value
		}
		-countfields {
			set countfields $value
		}
		-countfieldpatterns {
			set countfieldpatterns $value
		}
		-textsize {
			set textsize $value
		}
		-panelwidths {
			if {[llength $panelwidths] != 2 || ![isint [lindex $panelwidths 0]] || ![isint [lindex $panelwidths 1]]} {
				error "option -panelwidths should be a list of 2 numbers"
			}
			set panelwidths $value
		}
		-sortcol {
			set sortcol $value
		}
		-proportions {
			set proportions $value
		}
		-round_digits {
			set round_digits $value
		}
		-width {
			set width $value
		}
		-height {
			set height $value
		}
	} {isoform_counts_file gene output_file}
	set isoform_counts_file [file_absolute $isoform_counts_file]
	set output_file [file_absolute $output_file]

	set tempiso [tempfile].tsv
	set tempgtf [file root $tempiso].gtf
	file delete $tempiso
	cg select -q [subst {\$gene eq "$gene"}] $isoform_counts_file $tempiso
	cg tsv2gtf $tempiso $tempgtf
	if {$countfields eq ""} {
		set countfields {}
		set header [cg select -h $tempiso]
		foreach pattern $countfieldpatterns {
			foreach field $header {
				set out [regexp -inline $pattern $field]
				if {$out eq ""} continue
				regsub -all -- - $field . field
				if {[llength $out] > 1} {
					set name [join [lrange $out 1 end] .]
				} else {
					set name $field
					regsub counts_ $field {} name
				}
				lappend countfields $name $field
			}
		}
		# join $countfields \n
	}
	set cmd [subst {
		gtffile="$tempgtf"
		tsvfile="$tempiso"
		output_file="$output_file"
		rescale=$rescale
		textsize=$textsize
		panelwidths=c([join $panelwidths ,])
		height=$height
		width=$width
	}]
	append cmd {
		library(dplyr)
		library(tidyr)
		library(ggplot2)
		library(ggtranscript)
		library(ggpubr)
		
		tsv=read.table(tsvfile,header = T, sep = '\t')
		# sort transcripts on counts_weight
	}
	set temp {}
	foreach countfield $countfields {
		lappend temp [string_change $countfield {+ .}]
	}
	set countfields $temp
	if {![info exists sortcol]} {
		set sortcol [lindex $countfields 1]
	}
	if {$sortcol ne ""} {
		append cmd [subst -nocommands {
			transcripts = tsv\$transcript[order(tsv[["$sortcol"]])]
		}]
	}
	if {$proportions} {
		foreach {shortname countfield} $countfields {
			append cmd [subst -nocommands {
				tsv[["$countfield"]] = 100*tsv[["$countfield"]]/sum(tsv[["$countfield"]])
			}]
		}
	}
	if {$round_digits != ""} {
		foreach {shortname countfield} $countfields {
			append cmd [subst -nocommands {
				tsv[["$countfield"]] = round(tsv[["$countfield"]],$round_digits)
			}]
		}
	}
	append cmd {
		gtf=rtracklayer::import(gtffile)
		gtf = gtf %>% dplyr::as_tibble()
		exons = gtf %>% dplyr::filter(type == "exon")
		exons$novel = grepl("^novel",exons$transcript_id)
		if (!rescale) {
			# want to add locations of all potential exons everywhere
			diffs = to_diff(
			    exons = exons,
			    ref_exons = exons,
			    group_var = "transcript_id"
			)
			t = exons %>% 
			    ggplot(aes(
			        xstart = start,
			        xend = end,
			        y = transcript_id
			    )) +
			    geom_range(
			    ) +
			    geom_intron(
			        data = to_intron(exons, "transcript_id"),
			        aes(strand = strand), 
			        arrow.min.intron.length = 300
			    ) +
			    geom_range(
			        data = diffs,
			        aes(fill = diff_type,color=diff_type),
			        alpha = 0.2,
				show.legend = FALSE
			    ) +
			    scale_y_discrete(limits = transcripts,position = "right",label=function(x) abbreviate(x, minlength=18))
		} else {
			rescaled <- shorten_gaps(
			  exons = exons, 
			  introns = to_intron(exons, "transcript_id"), 
			  group_var = "transcript_id"
			)
			# shorten_gaps() returns exons and introns all in one data.frame()
			# let's split these for plotting 
			rescaled_exons <- rescaled %>% dplyr::filter(type == "exon") 
			rescaled_introns <- rescaled %>% dplyr::filter(type == "intron") 
			# want to add locations of all potential exons everywhere
			rescaled_diffs = to_diff(
			    exons = rescaled_exons,
			    ref_exons = rescaled_exons,
			    group_var = "transcript_id"
			)
			
			rescaled_exons$transcript_id
			
			t = rescaled_exons %>% 
			    ggplot(aes(
			        xstart = start,
			        xend = end,
			        y = transcript_id
			    )) +
			    geom_range(
			    ) +
			    geom_intron(
			        data = rescaled_introns,
			        aes(strand = strand), 
			        arrow.min.intron.length = 300
			    ) +
			    geom_range(
			        data = rescaled_diffs,
			        aes(fill = diff_type,color=diff_type),
			        alpha = 0.2,
				show.legend = FALSE
			    ) +
			    scale_y_discrete(limits = transcripts,position = "right",label=function(x) abbreviate(x, minlength=18))
		}
		# counts
		temp = data.frame(
			transcript_id=tsv$transcript,
	}
	set count_types {}
	set rcounts {}
	foreach {shortname countfield} $countfields {
		lappend count_types $shortname
		lappend rcounts [subst -nocommands {		$shortname=tsv[["$countfield"]]}]
	}
	append cmd [join $rcounts ,\n]
	append cmd [subst {
		)
		count_types = c("[join $count_types \",\"]")
	}]
	if {$proportions} {
		set output percent
	} else {
		set output count
	}
	append cmd [string_change {
		counts = temp %>%
		  as_tibble() %>%
		  pivot_longer(-transcript_id, names_to = "count_type", values_to = "@OUTPUT@")
		
		c = ggplot(counts, aes(y=transcript_id, x=count_type)) +
		  geom_tile(aes(fill = @OUTPUT@)) +
		  geom_text(aes(label = round(@OUTPUT@, 2)),size=textsize/.pt) +
		  scale_fill_gradient(low = "white", high = "steelblue") +
		  scale_y_discrete(limits = transcripts,label=function(x) return("")) +
		  scale_x_discrete(limits = count_types) + 
		  theme(axis.title.y = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
		
		p = ggarrange(t+rremove("ylab"),c+rremove("xlab"),align="h",widths=panelwidths)
		ggsave(output_file,p,width=width,height=height,units="mm")
	} [list @OUTPUT@ $output]]

	R $cmd
}

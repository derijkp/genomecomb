proc igvopen {} {
	# set port 60151
	# find free port
	for {set port 49152} {$port <= 65535} {incr port} {
		if {![catch {socket -server Connect $port} sid]} break
	}
	close $sid
	set igvpid [catch_exec igv -p $port &]
	set count 30
	while {[incr count -1]} {
		after 2000
		if {![catch {
			set igvport [socket 127.0.0.1 $port]
		}]} break
	}
	if {![info exists igvport]} {error "Could not connect to IGV (port $port)"}
	fconfigure $igvport -blocking 0
	return [list igv $igvpid $igvport]
}

proc igvclose {igv} {
	foreach {temp igvpid igvport} $igv break
	close $igvport
	exec kill $igvpid
}

proc igv {igv args} {
	foreach {temp igvpid igvport} $igv break
	puts $igvport $args
	flush $igvport
	gets $igvport
}

#	IGV control Commands: https://software.broadinstitute.org/software/igv/PortCommands
#	
#	new
#	Create a new session.  Unloads all tracks except the default genome annotations.
#	
#	load file 	
#	Loads data or session files.  Specify a comma-delimited list of full paths or URLs.
#	
#	Note: for Google Genomics readgroup sets the id is the "file",  specify format=ga4gh (version 2.3.80 and greater only).  For example
#	load CMvnhpKTFhCjz9_25e_lCw format=ga4gh
#	
#	collapse trackName
#	Collapses a given trackName. trackName is optional, however, and if it is not supplied all tracks are collapsed.
#	
#	echo
#	Writes "echo" back to the response.  (Primarily for testing)
#	
#	exit
#	Exit (close) the IGV application.
#	
#	expand trackName
#	Expands a given trackName. trackName is optional, however, and if it is not supplied all tracks are expanded.
#	
#	genome genomeIdOrPath
#	Selects a genome by id, or loads a genome (or indexed fasta) from the supplied path.
#	
#	goto locus or listOfLoci
#	Scrolls to a single locus or a space-delimited list of loci. If a list is provided, these loci will be displayed in a split screen view.  Use any syntax that is valid in the IGV search box.
#	
#	goto all
#	Scrolls to a whole genome view.
#	
#	region chr start end
#	Defines a region of interest bounded by the two loci (e.g., region chr1 100 200).
#	
#	maxPanelHeight height
#	Sets the number of vertical pixels (height) of each panel to include in image. Images created from a port command or batch script are not limited to the data visible on the screen. Stated another way, images can include the entire panel not just the portion visible in the scrollable screen area. The default value for this setting is 1000, increase it to see more data, decrease it to create smaller images.
#	
#	setSleepInterval ms
#	Sets a delay (sleep) time in milliseconds.  The sleep interval is invoked between successive commands.
#	
#	snapshotDirectory path
#	Sets the directory in which to write images.
#	
#	snapshot filename
#	Saves a snapshot of the IGV window to an image file.  If filename is omitted, writes a PNG file with a filename generated based on the locus.  If filename is specified, the filename extension determines the image file format, which must be .png, .jpg, or .svg.
#	
#	sort option locus
#	Sorts an alignment track by the specified option.  Recognized values for the option parameter are: base, position, strand, quality, sample,  readGroup, AMPLIFICATION, DELETION, EXPRESSION, SCORE, and MUTATION_COUNT.  The locus option can define a single position, or a range.  If absent sorting will be perfomed based on the region in view, or the center position of the region in view, depending on the option.
#	
#	squish trackName
#	Squish a given trackName. trackName is optional, and if it is not supplied all annotation tracks are squished.
#	
#	viewaspairs trackName
#	Set the display mode for an alignment track to "View as pairs".  trackName is optional.
#	
#	preference key value
#	Temporarily set the preference named key to the specified value. This preference only lasts until IGV is shut down.

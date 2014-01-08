invoke {} {
	global select_functions select_functions_insert
	unset -nocomplain select_functions
	unset -nocomplain select_functions_insert
	list_foreach {type name insert syntax descr} [split [string trim {
		calc	exp	field	{exp(arg)} {exponential of arg.}
		calc	fmod	field_value	{fmod(x,y)} {floating-point remainder of the division of x by y.}
		calc	isqrt	field	{isqrt(arg)} {Computes the integer part of the square root of arg.}
		calc	log	field	{log(arg)} {natural logarithm of arg. Arg must be a positive value.}
		calc	log10	field	{log10(arg)} {base 10 logarithm of arg. Arg must be a positive value.}
		calc	pow	field_value	{pow(x,y)} {Computes the value of x raised to the power y.}
		calc	sqrt	field	{sqrt(arg)} {The argument may be any non-negative numeric value.}
		calc	min	fields	{min(a1,a2,...)} {returns the minimum of a1, a2, ... min will return an error if one of the values is not a number. Use lmin if some values are list of numbers, or not numbers.}
		calc	max	fields	{max(a1,a2,...)} {returns the maximum of a1, a2, ... max will return an error if one of the values is not a number. Use lmax if some values are list of numbers, or not numbers.}
		calc	avg	fields	{avg(value,...)} {returns the average of the values given. Non-number values are ignored.}
		calc	rand	empty	{rand()} {Returns a pseudo-random floating-point value in the range (0,1). The generator algorithm is a simple linear congruential generator that is not cryptographically secure. Each result from rand completely determines all future results from subsequent calls to rand, so rand should not be used to generate a sequence of secrets, such as one-time passwords. The seed of the generator is initialized from the internal clock of the machine or may be set with the srand function.}
		calc	srand	field	{srand(arg)} {The arg, which must be an integer, is used to reset the seed for the random number generator of rand. Returns the first random number (see rand) from that seed. Each interpreter has its own seed.}
		calc	acos	field	{acos(arg)} {arc cosine of arg, in the range [0,pi] radians. Arg should be in the range [-1,1].}
		calc	asin	field	{asin(arg)} {arc sine of arg, in the range [-pi/2,pi/2] radians. Arg should be in the range [-1,1].}
		calc	atan	field	{atan(arg)} {arc tangent of arg, in the range [-pi/2,pi/2] radians.}
		calc	atan2	field_value	{atan2(y,x)} {arc tangent of y/x, in the range [-pi,pi] radians. x and y cannot both be 0. If x is greater than 0, this is equivalent to âatan [expr {y/x}]â.}
		calc	cos	field	{cos(arg)} {cosine of arg, measured in radians.}
		calc	cosh	field	{cosh(arg)} {hyperbolic cosine of arg. If the result would cause an overflow, an error is returned.}
		calc	hypot	field_value	{hypot(x,y)} {Computes the length of the hypotenuse of a right-angled triangle âsqrt [expr {x*x+y*y}]â.}
		calc	sin	field	{sin(arg)} {sine of arg, measured in radians.}
		calc	sinh	field	{sinh(arg)} {hyperbolic sine of arg. If the result would cause an overflow, an error is returned.}
		calc	tan	field	{tan(arg)} {tangent of arg, measured in radians.}
		calc	tanh	field	{tanh(arg)} {hyperbolic tangent of arg.}
		number	between	field_values	{between(value,{min max}) or between(value,min,max)} {true of value is >= min and <= max (e.g. "between($begin,1000,2000)")}
		number	isnum	field	{isnum(value)} {true if value is a valid number}
		number	def	field_value	{def(value,default)} {if value is not a number, it returns the given default, otherwise value}
		number	percent	field	{percent(value)} {returns a fraction as a percent }
		convert	ceil	field	{ceil(arg)} {smallest integral floating-point value (i.e. with a zero fractional part) not less than arg. The argument may be any numeric value.}
		convert	floor	field	{floor(arg)} {largest integral floating-point value (i.e. with a zero fractional part) not greater than arg. The argument may be any numeric value.}
		convert	round	field	{round(arg)} {If arg is an integer value, returns arg, otherwise converts arg to integer by rounding and returns the converted value.}
		convert	abs	field	{abs(arg)} {absolute value of arg.}
		convert	double	field	{double(arg)} {The argument may be any numeric value, If arg is a floating-point value, returns arg, otherwise converts arg to floating-point and returns the converted value. May return Inf or -Inf when the argument is a numeric value that exceeds the floating-point range.}
		convert	entier	field	{entier(arg)} {The argument may be any numeric value. The integer part of arg is determined and returned. The integer range returned by this function is unlimited, unlike int and wide which truncate their range to fit in particular storage widths.}
		convert	int	field	{int(arg)} {The argument may be any numeric value. The integer part of arg is determined, and then the low order bits of that integer value up to the machine word size are returned as an integer value. For reference, the number of bytes in the machine word are stored in tcl_platform(wordSize).}
		convert	bool	field	{bool(arg)} {Accepts any numeric value, or any string acceptable to string is boolean, and returns the corresponding boolean value 0 or 1. Non-zero numbers are true. Other numbers are false. Non-numeric strings produce boolean value in agreement with string is true and string is false.}
		convert	wide	field	{wide(arg)} {The argument may be any numeric value. The integer part of arg is determined, and then the low order 64 bits of that integer value are returned as an integer value.}
		logical	if	special {if(condition,true,?condition2,true2, ...?false}	{if condition is true, the value for "true" will be returned, otherwise the last parameter (false) is returned}
		bio	region	special {region("chromosome:begin-end",...)}	{is true for any region in dataset that overlaps the given regions. Can also be given as region(chromosome,begin,end,...).}
		bio	hovar	special {hovar(samplename)} {true if the given sample is a homozygous variant. This is equivalent to ($sequenced-samplename == "v" && $alleleSeq1-samplename == $alleleSeq2-samplename)}
		bio	compare	special	{compare(samplename1,samplename2, ...)} {compares the variant in the given samples, and returns one of: sm (variant with the same genotype in all given samples, with all sequenced), df (different: variant in some, reference in other, with all sequenced), mm (mismatch; variant in all, but different genotypes, with all sequenced), un (unsequenced in some samples, variant in one of the others)}
		bio	same	special	{same(sample1,sample2, ...)} {same: all samples have the same genotype (does not have to be a variant) (all sequenced)}
		bio	sm	special	{sm(sample1,sample2, ...)} {same: variant with the same genotype in all given samples (all sequenced)}
		bio	df	special	{df(sample1,sample2, ...)} {different: variant in some, reference in other (all sequenced)}
		bio	mm	special	{mm(sample1,sample2, ...)} {mismatch; variant in all, but different genotypes (all sequenced)}
		bio	un	special	{un(sample1,sample2, ...)} {unsequenced in some samples, variant in one of the others}
		string	oneof	field_values	{oneof($field,value1,value2,...)} {returns true if the given field is equal to one of the values}
		string	regexp	field_value	{regexp(value,pattern)} {true if value matches the regular expression given by **pattern**}
		string	matches	field_value	{matches(value,pattern)} { true if value matches the glob pattern given by **pattern** (using wildcards * for anything, ? for any single character, [chars] for any of the characters in chars)}
		string	concat	fields	{concat(value,...)} {makes one long string by appending all values.}
		list	vector	fields	{vector(value1,value2, ...)} {creates a vector from a number of values. If some elements are vectors themselves, they will be concatenated}
		list	lindex	field_value	{lindex(vector, position)} {the value of the element at the given **position** in the list. The first element is at position 0!}
		list	llen	field	{llen(vector)} {number of elements in the vector}
		list	lmin	fields	{lmin(vector, ...)} {the minimum of the list of numbers in vector(s). A default value (NaN or not a number) is given for non-numeric characters (-); any comparison with NaN is false.}
		list	lmax	fields	{lmax(vector, ...)} {the maximum of the vector. A default value (NaN or not a number) is given for non-numeric characters (-); any comparison with NaN is false.}
		list	lmind	fields	{lmind(vector, ..., def)} {same as lmin, but you can set the default value for non-numeric characters is given as the last parameter}
		list	lmaxd	fields	{lmaxd(vector, ..., def)} {same as lmax, but you can set the default value for non-numeric characters is given as the last parameter}
		list	lminpos	fields	{lminpos(vector, ...)} {position (within the index) of the minimum value. If more than one vector is given, the position of the minimum of all vectors is given}
		list	lmaxpos	fields	{lmaxpos(vector, ...)} {position (within the index) of the maximum value. If more than one vector is given, the position of the maximum of all vectors is given}
		list	lsum	fields	{lsum(vector, ...)} {the sum of the list of numbers in vector(s). Non numeric values are ignored. If no numeric value is present in the vectors, NaN (not a number) will be returned; any comparison with NaN is false.}
		list	lavg	fields	{lavg(vector, ...)} {the average of the vector. Non numeric values are ignored. If no numeric value is present in the vectors, NaN (not a number) will be returned; any comparison with NaN is false.}
		list	lstdev	fields	{lstdev(vector, ...)} {the standard deviation of the vector. Non numeric values are ignored. If no numeric value is present in the vectors, NaN (not a number) will be returned; any comparison with NaN is false.}
		list	lmedian	fields	{lmedian(vector, ...)} {the median of the vector. Non numeric values are ignored. If no numeric value is present in the vectors, NaN (not a number) will be returned; any comparison with NaN is false.}
		list	lmode	fields	{lmode(vector, ...)} {the mode (element that is most abundant) of the vector. The result can be a new vector (if multiple values occur at the same count)}
		list	contains	field_value	{contains(vector, value)} {true if **vector** contains **value**. This can also be used as an operator: vector contains value}
		list	shares	field_value	{shares(vector, valuelist)} {true if **vector** and the list in **valuelist** (a SPACE separated list!) share a value. This can also be used as an operator: vector shares valuelist}
		list	lone	field	{lone(vector)} {true if one of elements of the vector is true}
		list	lall	field	{lall(vector)} {true if all elements of the vector are true}
		list	lcount	field	{lcount(vector)} {number of elements in vector that are true}
		list	vdistinct	fields	{vdistinct(vector, ...)} {returns a vector in which each element in one of the vectors occurs only once}
		list	vabs	field	{vabs(vector)} {returns vector of absolute values of given vector}
		list	vavg	fields	{vavg(vector1,vector2,...)} {returns vector with average value for each position in the vector}
		list	vmax	fields	{vmax(vector1,vector2,...)} {returns vector with maximum value for each position in the vector}
		list	vmin	fields	{vmin(vector1,vector2,...)} {returns vector with minimum value for each position in the vector}
		list	vdef	field_value	{vdef(vector,default)} {returns the given vector, but with all non numbers replaced by default}
		list	vif	special	{vif(condition,true,?condition2,true2, ...?false)} {like if, but conditions, true1, ... and false may be vectors, and a vector is returned}
		sample_aggregates	scount	saggr	{scount(condition)} {number of samples for which **condition** is true}
		sample_aggregates	slist	saggr	{slist(?condition?,value)} {returns a (comma separated) list with results of value for each sample for which (if given) **condition** is true}
		sample_aggregates	sdistinct	saggr	{sdistinct(?condition?,value)} {returns a non-redundant (comma separated) list of the results of value for each sample for which (if given) **condition** is true}
		sample_aggregates	smin	saggr	{smin(?condition?,value)} {returns the minimum of results of value for each sample for which (if given) **condition** is true}
		sample_aggregates	smax	saggr	{smax(?condition?,value)} {returns the maximum of results of value for each sample for which (if given) **condition** is true}
		sample_aggregates	ssum	saggr	{ssum(?condition?,value)} {returns the sum of results of value for each sample for which (if given) **condition** is true}
		sample_aggregates	savg	saggr	{savg(?condition?,value)} {returns the average of results of value for each sample for which (if given) **condition** is true}
		sample_aggregates	sstdev	saggr	{sstdev(?condition?,value)} {returns the standard deviation of results of value for each sample for which (if given) **condition** is true}
		sample_aggregates	smedian	saggr	{smedian(?condition?,value)} {returns the median of results of value for each sample for which (if given) **condition** is true}
		sample_aggregates	smode	saggr	{smode(?condition?,value)} {returns the mode of results of value for each sample for which (if given) **condition** is true}
		sample_aggregates	spercent	saggr	{spercent(condition1,condition2)} {returns 100.0*(number of samples for which condition1 and condition2 are true)/(number of samples for which condition1 is true)}
	}] \n] {
		lappend select_functions($type) "$name  -  $syntax  -  $descr"
		set select_functions_insert($name) $insert
	}
}

mainw method querybuilder_fillvalues {args} {
	private $object qvalues qvaluestext fieldsw valuesw valuew
	set w $object.querybuilder.options.paned
	set field [lindex [$fieldsw get] 0]
	if {$field eq ""} return
	set list [$object.tb values $field]
	if {[lindex $list end 1] ne "incomplete"} {
		set text [list "All values"]
	}
	while {![isint [lindex $list end 1]]} {
		if {![llength $list]} break
		set line [list_pop list]
		lappend text "[lindex $line 1]: [lindex $line 0]"
	}
	set qvaluestext [join $text \n]
	set qnums [list_subindex $list 1]
	set list [list_subindex $list 0]
	if {[llength $list]} {
		# found the values, do nothing else
	} elseif {[regexp ^sequenced- $field]} {
		set list {r v u}
	} elseif {[regexp _impact $field]} {
		set list {CDS RNA}
	} else {
		set list {}
	}
	set qvalues $list
	if {![inlist $qvalues [$valuew get]]} {
		$valuesw activate [lindex $qvalues 0]
		$valuesw set [lindex $qvalues 0]
	} else {
		$valuesw set [$valuew get]
	}
	$valuew set [$valuesw get]
#	set qvalues {}
#	foreach value $list num $qnums {
#		lappend qvalues "$value ($num)"
#	}
	return $qvalues
}

mainw method querybuilder_selfield {value} {
	$object.querybuilder.options.paned.query.query insert current " \$$value "
}

mainw method querybuilder_browsefield {args} {
	private $object operatorsw valuesw
	if {[$operatorsw get] eq ""} {
		catch {$operatorsw set [$operatorsw get 0]}
	}
	Classy::todo $object querybuilder_fillvalues
#	if {[$valuesw get] eq ""} {
#		catch {$valuesw set [$valuesw get 0]}
#	}
#	$object querybuilder_browsevalue
}

mainw method querybuilder_browsevalue {args} {
	private $object valuesw valuew
	$valuew set [$valuesw get]
}

mainw method querybuilder_seloperator {value} {
	$object.querybuilder.options.paned.query.query insert current " $value "
}

mainw method querybuilder_selvalue {value} {
	if {[isdouble $value]} {
		set insert " $value "
	} else {
		set insert " \"$value\" "
	}
	$object.querybuilder.options.paned.query.query insert current $insert
}

mainw method querybuilder_quotevalue {value} {
	if {[isdouble $value]} {
		return $value
	} else {
		return \"$value\"
	}
}

mainw method querybuilder_quotevalues {list} {
	set result {}
	foreach value $list {
		if {[isdouble $value]} {
			lappend result $value
		} else {
			lappend result \"$value\"
		}
	}
	return $result
}

mainw method querybuilder_insert {insert {join and}} {
	private $object queryw
	upvar #0 [$object.buttons.query cget -textvariable] query
	set w $object.querybuilder.options.paned
	if {[string trim $query] ne ""} {set insert "\n$join $insert"}
	$queryw insert insert $insert
}

mainw method querybuilder_makecondition {fields operator values {join " or "}} {
	set insert {}
	foreach field $fields {
		set temp {}
		if {[inlist {shares in ni} $operator]} {
			lappend insert "\$$field $operator [list $values]"
		} else { 
			foreach value $values {
				set value [$object querybuilder_quotevalue $value]
				lappend temp "\$$field $operator $value"
			}
			if {[llength $temp] > 1} {
				lappend insert "\( [join $temp $join] \)"
			} else {
				lappend insert "[join $temp $join]"
			}
		}
	}
	return $insert
}

mainw method querybuilder_add {command {join and}} {
	private $object fieldsw operatorsw valuew queryw
	upvar #0 [$object.buttons.query cget -textvariable] query
	set w $object.querybuilder.options.paned
	set fields [$fieldsw get]
	set operator [$operatorsw get]
	set values [$valuew get]
	if {$command eq "count"} {
		if {![llength $fields]} {error "Select some fields first"}
		set value [lindex $values 0]
		set value [$object querybuilder_quotevalue $value]
		set insert "count\("
		foreach field $fields {
			append insert "\$$field, "
		}
		append insert "$operator $value \) > 0"
		$object querybuilder_insert $insert $join
	} elseif {$command eq "field"} {
		if {![llength $fields]} {error "Select some fields first"}
		$queryw insert insert "\$\{[join $fields "\},\$\{"]\}"
	} elseif {$command eq "value"} {
		if {![llength $values]} {error "Select some values first"}
		set temp {}
		foreach v $values {
			if {[isdouble $v]} {lappend temp $v} else {lappend temp \"$v\"}
		}
		$queryw insert insert "[join $temp ,]"
	} elseif {$command eq "cond"} {
		if {![llength $values]} {error "Select some values first"}
		set value [$object querybuilder_quotevalue [lindex $values 0]]
		$queryw insert insert " $operator $value"
	} elseif {$command eq "region"} {
		destroy $object.region
		Classy::Dialog $object.region -title "Select region"
		$object.region option entry Chromosome [privatevar $object region(chr)]
		$object.region option entry Begin [privatevar $object region(begin)]
		$object.region option entry End [privatevar $object region(end)]
		$object.region add go Go "$object querybuilder_insert region(\[getprivate $object region(chr) \]:\[getprivate $object region(begin) \]-\[getprivate $object region(end) \]) $join" default
	} elseif {$command in {compare same sm df mm un}} {
		private $object compare tdata
		destroy $object.compare
		Classy::Dialog $object.compare -title "Compare"
		set compare(types) {
			{compare returns sm,df,mm,un depending on type of match between two samples}
			{same same: all samples have the same genotype (does not have to be a variant) (all sequenced)}
			{sm same: variant with the same genotype in all given samples (all sequenced)}
			{df different: variant in some, reference in other (all sequenced)}
			{mm mismatch; variant in all, but different genotypes (all sequenced)}
			{un unsequenced in some samples, variant in one of the others}
		}
		$object.compare option select Comparison [privatevar $object compare(type)] [privatevar $object compare(types)]
		set compare(type) [lindex $compare(types) [lsearch [list_subindex $compare(types) 0] $command]]
		set compare(samples) [cg select -n $tdata(file)]
		$object.compare option listbox Sample1 [privatevar $object compare(selsamples)] [privatevar $object compare(samples)] -selectmode multiple -filtervariable compare_filter
		update idletasks
		set compare(selsamples) [samples $fields]
#		$object.compare add go Go "$object querybuilder_insert \"\[lindex \[getprivate $object compare(type) \] 0\]\(\"\[join \[getprivate $object compare(selsamples) \] \",\"\]\)\" $join" default
		$object.compare add go Go [list invoke {object} {
			set type [lindex [getprivate $object compare(type)] 0]
			set samples [getprivate $object compare(selsamples)]
			$object querybuilder_insert "${type}(\"[join $samples \",\"]\")"
			$object.compare destroy
		} $object]
	} elseif {$command in {if vif}} {
		if {![llength $fields]} {error "Select some fields first"}
		if {$command eq "vif"} {
			if {![regexp ^@ $operator]} {set operator @$operator}
			set condjoin " vor "
		} else {
			set condjoin " or "
		}
		set insert [join [$object querybuilder_makecondition $fields $operator $values $condjoin] " and "]
		private $object if
		destroy $object.if
		set if(if) $insert
		set if(true) 1
		set if(false) 0
		Classy::Dialog $object.if -title "Conditional value (if)"
		$object.if option entry Criterion [privatevar $object if(if)]
		$object.if option entry Truevalue [privatevar $object if(true)]
		$object.if option entry FalseValue [privatevar $object if(false)]
		$object.if add go Go "$object querybuilder_insert if(\[getprivate $object if(if) \],\[getprivate $object if(true) \],\[getprivate $object if(false) \]) $join" default
	} elseif {$command in {and or}} {
		if {![llength $fields]} {error "Select some fields first"}
		set insert [join [$object querybuilder_makecondition $fields $operator $values] "\n    $command "]
		set join $command
		$object querybuilder_insert $insert $join
	} else {
		global select_functions_insert
		if {![llength $fields]} {error "Select some fields first"}
		set type [get select_functions_insert($command) fields]
		if {$type eq "field"} {
			# single field functions
			set insert {}
			foreach field $fields {
				lappend insert "${command}(\$$field)"
			}
			set insert [join $insert " $join "]
		} elseif {$type eq "field_values"} {
			# single field functions
			set values [$object querybuilder_quotevalues $values]
			set insert {}
			foreach field $fields {
				lappend insert "${command}(\$$field,[join $values ,])"
			}
			set insert [join $insert " $join "]
		} elseif {$type eq "field_value"} {
			# single field functions
			set values [$object querybuilder_quotevalues $values]
			set insert {}
			foreach field $fields {
				set temp {}
				foreach value $values {
					lappend temp "${command}(\$$field,$value)"
				}
				if {[llength $temp] > 1} {
					lappend insert \([join $temp " or "]\)
				} else {
					lappend insert [join $temp " or "]
				}
			}
			set insert [join $insert " $join "]
		} elseif {$type eq "fields_cond"} {
			# multi field functions with condition
			set values [$object querybuilder_quotevalues $values]
			set insert {}
			foreach value $values {
				lappend insert "${command}(\$[join $fields ,\$]) $operator $value"
			}
			set insert [join $insert " $join "]
		} elseif {$type eq "field_cond"} {
			# multi field functions with condition
			set values [$object querybuilder_quotevalues $values]
			set insert {}
			foreach field $fields {
				set temp {}
				foreach value $values {
					lappend temp "${command}(\$field) $operator $value"
				}
				if {[llength $temp] > 1} {
					lappend insert \([join $temp " or "]\)
				} else {
					lappend insert [join $temp " or "]
				}
			}
			set insert [join $insert " $join "]
		} elseif {$type eq "cond"} {
			set condition [$object querybuilder_makecondition $fields $operator $values]
			set insert ${command}($condition)
		} elseif {$type eq "empty"} {
			set condition [$object querybuilder_makecondition $fields $operator $values]
			set insert ${command}()
		} elseif {$type eq "sample_aggregates" || $command in {scount slist sdistinct smin smax ssum savg sstdev smedian smode spercent}} {
			# sample aggregates
			set sfields {}
			set ssamples {}
			foreach field $fields {
				set pos [string first - $field]
				if {$pos == -1} {
					lappend sfields $field
				} else {
					lappend sfields [string range $field 0 [expr {$pos-1}]]
					lappend ssamples [string range $field [expr {$pos+1}] end]
				}
			}
			set condition [$object querybuilder_makecondition $sfields $operator $values]
			if {[llength $ssamples] == 1} {
				set pre "\$sample == \"[list $ssamples]\" and "
			} elseif {[llength $ssamples]} {
				set pre "\$sample in [list $ssamples] and "
			} else {
				set pre ""
			}
			if {![llength $condition] > 1} {
				set condition "${pre}$condition"
			} else {
				set condition "${pre}([join $condition " and "])"
			}
			if {$command eq "scount"} {
				set insert "${command}($condition)"
			} else {
				set insert "${command}($condition, )"
			}
			$object querybuilder_insert $insert $join
			tk::TextSetCursor $queryw insert-1c
			return
		} else {
			# multi field functions (fields)
			set insert ${command}(\$[join $fields ,\$])
		}
		$object querybuilder_insert $insert $join
	}
}

mainw method querybuilder_function_add {function {join and}} {
	private $object fieldsw operatorsw valuew functions
putsvars function
	$object querybuilder_add $function $join
#	set fields [$fieldsw get]
#	set operator [$operatorsw get]
#	set values [$valuew get]
#	$object querybuilder_insert ${function}(\$[join $fields ,\$]) $join
}

mainw method querybuilder_functions_type {type} {
	private $object functions
	global select_functions
	set select_functions(filter) {}
	if {$type eq "all"} {
		set temp {}
		foreach type [array names select_functions] {
			lappend temp {*}$select_functions($type)
		}
		set temp [lsort -dict $temp]
	} else {
		set temp $select_functions($type)
	}
	set functions(available) $temp
}

mainw method querybuilder_functions {{join and} {starttype all}} {
	private $object functions
	global select_functions
	destroy $object.functions
	Classy::Dialog $object.functions -title "Functions"
	set temp {}
	set types [list_union [list_common {calc number string list convert bio logical sample_aggregates} [array names select_functions]] [array names select_functions]]
	lappend temp all [list $object querybuilder_functions_type all]
	foreach type $types {
		lappend temp $type [list $object querybuilder_functions_type $type]
	}
	$object.functions option buttons Type {*}$temp
	$object querybuilder_functions_type $starttype
	set functions(sel) [lindex $functions(available) 0]
	$object.functions option listbox Functions [privatevar $object functions(sel)] \
		[privatevar $object functions(available)] -selectmode single \
		-command [list invoke {object join sel} {
			$object querybuilder_function_add [lindex $sel 0 0] $join
			$object.functions destroy
		} $object $join] -filtervariable select_functions(filter)
}

mainw method querybuilder_fieldsreset {type} {
	private $object qfieldsfilter qfields
	if {$type eq "samplefields"} {
		set list {}
		foreach field [$object.tb tfields] {
			if {[regexp {^([^-]+)-} $field temp field]} {
				lappend list $field
			}
		}
		set qfields [list_remdup $list]
	} else {
		set qfields [$object.tb tfields]
	}
	set qfieldsfilter {}
}

mainw method querybuilder {args} {
	private $object qfields qoperators qvalues fieldsw operatorsw valuesw valuew queryw qfieldsfilter valuesfilter
	set qfields [$object.tb tfields]
	set qoperators {== != < <= > >= matches regexp contains in ni shares}
	set qfieldsfilter {}
	set valuesfilter {}
	destroy $object.querybuilder
	Classy::Dialog $object.querybuilder -title Query
	if {[llength $args]} {
		$object.querybuilder configure -title "Calculated field builder"
		set fieldbuilder 1
		set fieldbuilderw [lindex $args 0]
		set var [privatevar $object fieldcalc]
		Classy::Entry $object.querybuilder.options.fieldname -label Fieldname -textvariable [privatevar $object fieldname]
		pack $object.querybuilder.options.fieldname -side top -expand yes -fill x
		set funcbuttons {_blank field value _blank if count lmin lmax avg sum compare percent and or region isnum}
		set join ""
	} else {
		set fieldbuilder 0
		set var [$object.buttons.query cget -textvariable]
		set funcbuttons {_blank and or field value cond _blank count region compare lmin lmax avg sum if isnum percent}
		set join and
	}
	set w $object.querybuilder.options.paned
	ttk::panedwindow $w -orient horizontal
	pack $w -fill both -expand yes
	# fields
	frame $w.fields -borderwidth 0 -highlightthickness 0
	button $w.fields.header -text "Fields" -command "[list $object querybuilder_fieldsreset fields]"
	pack $w.fields.header -fill x
	button $w.fields.samplefields -text "Sample fields" -command "[list $object querybuilder_fieldsreset samplefields]"
	pack $w.fields.samplefields -fill x
	set fieldsw $w.fields.fields
	Classy::ListBox $w.fields.fields \
		-listvariable [privatevar $object qfields] \
		-filtervariable [privatevar $object qfieldsfilter] \
		-command [list $object querybuilder_selfield] \
		-browsecommand [list $object querybuilder_browsefield] \
		-exportselection 0 -selectmode extended
	pack $w.fields.fields -fill both -expand yes
	# operators
	set w $w.pane2
	ttk::panedwindow $w -orient horizontal
	frame $w.operators -borderwidth 0 -highlightthickness 0
	button $w.operators.header -text "Operator"
	pack $w.operators.header -fill x
	set operatorsw $w.operators.operators
	Classy::ListBox $w.operators.operators \
		-listvariable [privatevar $object qoperators] \
		-command [list $object querybuilder_seloperator] \
		-width 6 -exportselection 0 -selectmode extended
	pack $w.operators.operators -fill both -expand yes
	# values
	set w $w.pane3
	ttk::panedwindow $w -orient horizontal
	frame $w.values -borderwidth 0 -highlightthickness 0
	button $w.values.header -text "Values" -command "[list set [privatevar $object valuesfilter] {}]"
	pack $w.values.header -fill x
	set valuew $w.values.value
	Classy::Entry $w.values.value
	pack $w.values.value -fill x
	label $w.values.label -textvariable [privatevar $object qvaluestext] -justify left -anchor w
	pack $w.values.label -fill x
	set valuesw $w.values.values
	Classy::ListBox $w.values.values \
		-listvariable [privatevar $object qvalues] \
		-filtervariable [privatevar $object qvaluesfilter] \
		-command [list $object querybuilder_selvalue] \
		-browsecommand [list Classy::todo $object querybuilder_browsevalue] \
		-exportselection 0 -selectmode extended
	pack $w.values.values -fill both -expand yes
	# query
	frame $w.query -borderwidth 0 -highlightthickness 0
	frame $w.query.buttons -borderwidth 0 -highlightthickness 0
	button $w.query.buttons.clear -text clear -command "[list set $var {}]"
	pack $w.query.buttons.clear -side left
	button $w.query.buttons.functions -text functions -command [list $object querybuilder_functions $join]
	pack $w.query.buttons.functions -side left
	set sep 1
	foreach command $funcbuttons {
		if {$command eq "_blank"} {
			frame $w.query.buttons.sep$sep -width 8
			pack $w.query.buttons.sep$sep -side left
			incr sep
		} else {
			button $w.query.buttons.$command -text $command -command [list $object querybuilder_add $command $join]
			pack $w.query.buttons.$command -side left
		}
	}
	set command string
	button $w.query.buttons.$command -text Functions -command [list $object querybuilder_functions $join]
	
	set queryw $w.query.query
	Classy::ScrolledText $w.query.query -textvariable $var
	pack $w.query.buttons -fill x -expand no
	pack $w.query.query -fill both -expand yes
	# fill paned
	set w $object.querybuilder.options.paned
	$w add $w.fields
	$w add $w.pane2
	set w $w.pane2
	$w add $w.operators
	$w add $w.pane3
	set w $w.pane3
	$w add $w.values
	$w add $w.query
	set w $object.querybuilder.options.paned
	if {$fieldbuilder} {
		$object.querybuilder add query "Insert calc field" "[list $fieldbuilderw] addcalc \"\$\{[privatevar $object fieldname]\}=\$\{$var\}\"; destroy $object.querybuilder"
	} else {
		$object.querybuilder add query "Query" "[list $object query]; destroy $object.querybuilder"
	}
	$object.querybuilder add clear "Clear" "[list set $var {}]"
	$object querybuilder_fillvalues
	update
	set field [$w.fields.fields get 0]
	if {$field ne ""} {
		$w.fields.fields set $field
		Classy::todo $object querybuilder_browsefield
	}
}

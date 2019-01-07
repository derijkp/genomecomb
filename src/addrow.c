/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE64_SOURCE 1

#include "cg.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include "tools.h"
#include "tools_bcol.h"
#include "debug.h"
#include "gztools.h"

int main(int argc, char *argv[]) {
	DString *line = NULL;
	unsigned int cnt, row = 0;
	if ((argc != 1)) {
		fprintf(stderr,"Format is: addrow\n");
		exit(EXIT_FAILURE);
	}
	line = DStringNew();
/*
	skip_header(stdin,line,&numfields,&pos);
	DStringputs(line,stdout);
	fprintf(stdout,"\tROW\n");
 */
	while(1) {
		cnt = DStringGetLine(line,stdin);
		if (cnt == -1) break;
		DStringputs(line,stdout);
		fprintf(stdout,"\t%d\n",row);
		row++;
	}
	exit(EXIT_SUCCESS);
}

/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include "tools.h"
#include "debug.h"

int main(int argc, char *argv[]) {
	DString *line = DStringNew();
	char *endptr;
	int read;
	long int value;
	int pos=0;
	if ((argc != 1)) {
		fprintf(stderr,"Format is: test\n");
		exit(EXIT_FAILURE);
	}
	while ((read = DStringGetLine(line, stdin)) != -1) {
/*
		value = strtol(line->string,&endptr,10);
		if ((*endptr == '\0' && *(line->string) != '\0') && value > 20) {
			fprintf(stdout,"%d\n",pos);
		}
*/
		if (line->size == 1 && line->string[0] == 'u') {
			fprintf(stdout,"%d\n",pos);
		}
		pos++;
	}
	if (line) {DStringDestroy(line);}
	exit(EXIT_SUCCESS);
}

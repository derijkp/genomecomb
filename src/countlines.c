/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "debug.h"

int main(int argc, char *argv[]) {
	uint64_t count;
	char c;
	if (argc != 1) {
		fprintf(stderr,"Format is: countlines\n");
		exit(EXIT_FAILURE);
	}
	count = 0;
	while (1) {
		c = getc_unlocked(stdin);
		if (c == '\n') {
			count++;
		} else if (c == EOF) {
			break;
		}
	}
	fprintf(stdout,"%" PRIu64 "\n",count);
	exit(EXIT_SUCCESS);
}

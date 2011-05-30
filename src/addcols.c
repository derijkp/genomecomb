/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
	char *line = NULL,*scan;
	size_t len = 0;
	ssize_t read;
	int cols=-1,num,pos=0;
	if (argc != 2) {
		fprintf(stderr,"Format is: map2besthit numcols");
		exit(EXIT_FAILURE);
	}
	cols = atoi(argv[1]);
	while ((read = getline(&line, &len, stdin)) != -1) {
		pos++;
		num = cols - 1;
		scan = line;
		while(*scan) {
			if (*scan++ == '\t') {num--;}
		}
		*(scan-1) = '\0';
		if (num < 0) {
			fprintf(stderr,"too many columns at line %d: %s\n", pos, line);
			exit(1);
		}
		if (num) {
			fprintf(stdout,"%s", line);
			while (num--) {
				fprintf(stdout,"\t");
			}
			fprintf(stdout,"\n");
		} else {
			fprintf(stdout,"%s\n", line);
		}
	}
	if (line) {
		free(line);
	}
	exit(EXIT_SUCCESS);
}

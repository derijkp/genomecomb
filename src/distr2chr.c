
#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

int main(int argc, char *argv[]) {
	FILE *o[26];
	char *line = NULL,*linepos = NULL;
	size_t len = 0;
	ssize_t read;
	char chromosome[10];
	int i,nchr,col = 0,count;
	if ((argc != 2)&&(argc != 3)) {
		fprintf(stderr,"Format is: distr2chr output_pre ?col?");
		exit(EXIT_FAILURE);
	}
	if (argc == 3) {
		col = atoi(argv[2]);
	}
	line = malloc(strlen(argv[1])+5);
	for (i=0 ; i < 23 ; i++) {
		sprintf(line,"%s-%d",argv[1],i);
		o[i] = fopen64(line,"a");
	}
	sprintf(line,"%s-X",argv[1]);
	o[i++] = fopen64(line,"a");
	sprintf(line,"%s-Y",argv[1]);
	o[i++] = fopen64(line,"a");
	sprintf(line,"%s-M",argv[1]);
	o[i++] = fopen64(line,"a");
	while ((read = getline(&line, &len, stdin)) != -1) {
		if (col == 0) {
			sscanf(line,"%9s\t",chromosome);
		} else {
			linepos = line;
			count = col;
			while (*linepos && count) {
				if (*linepos == '\t') {
					count--;
				}
				linepos++;
			}
			if (*linepos) {
				sscanf(linepos,"%9s\t",chromosome);
			} else {
				chromosome[0] = '0';
				chromosome[1] = '\0';
			}
		}
		nchr = chromosomenum(chromosome);
		fprintf(o[nchr],"%s", line);
	}
	for (i=1 ; i < 26 ; i++) {
		fclose(o[i]);
	}
	if (line) {
		free(line);
	}
	exit(EXIT_SUCCESS);
}

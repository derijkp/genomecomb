
#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *line = NULL;
static size_t len;

int main(int argc, char *argv[]) {
	char *chr, *linepos = NULL, *scanpos = NULL;
	ssize_t read;
	int poscol,valuecol,above,shift,maxcol,pos,value,cutoff;
	int status,begin=0,count;
	if ((argc != 7)) {
		fprintf(stderr,"Format is: getregions chromosome poscol valuecol cutoff above shift");
		exit(EXIT_FAILURE);
	}
	chr = argv[1];
	poscol = atoi(argv[2]);
	valuecol = atoi(argv[3]);
	cutoff = atoi(argv[4]);
	above = atoi(argv[5]);
	shift = atoi(argv[6]);
	pos = 0 - shift;
	maxcol = poscol ; if (valuecol > maxcol) {maxcol = valuecol;}
	getline(&line, &len, stdin);
	while ((read = getline(&line, &len, stdin)) != -1) {
		if (line[0] == '\0') continue;
		if (line[0] != '#' && line[0] != '>') break;
	}
	begin = -1;
	status = 0;
	while (1) {
		linepos = line;
		count = 0;
		while (*linepos && (count <= maxcol)) {
			if (*linepos == '\t' || (linepos == line)) {
				if (*linepos == '\t') {
					scanpos = linepos+1;
				} else {
					scanpos = linepos;
				}
				if (count == poscol) {
					sscanf(scanpos,"%d",&pos);
				} else if (count == valuecol) {
					sscanf(scanpos,"%d",&value);
				}
				count++;
			}
			linepos++;
		}
		if (count > maxcol) {
			if (status) {
				if ((above && (value <= cutoff)) || (!above && (value >= cutoff))) {
					fprintf(stdout,"%s\t%d\t%d\n", chr, begin+shift, pos+shift);
					status = 0;
				}
			} else {
				if ((above && (value > cutoff)) || (!above && (value < cutoff))) {
					begin = pos;
					status = 1;
				}
			}
		}
		if ((read = getline(&line, &len, stdin)) == -1) break;
	}
	if (status) {
		fprintf(stdout,"%s\t%d\t%d\n", chr, begin+shift, pos+1+shift);
	}
	if (line) {
		free(line);
	}
	exit(EXIT_SUCCESS);
}

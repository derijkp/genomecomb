
#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int chromosomenum(char *chromosome) {
	int i;
	if (strncmp(chromosome,"chr",3) == 0) {
		chromosome += 3;
	} else if (strncmp(chromosome,"CHR",3) == 0) {
		chromosome += 3;
	}
	if (*chromosome == 'M') {
		return 95;
	} else if (*chromosome == 'X') {
		return 96;
	} else if (*chromosome == 'Y') {
		return 97;
	} else {
		i = atoi(chromosome);
		if (i < 0 || i > 22) {
			fprintf(stderr,"wrong chromosome %s",chromosome);
			exit(EXIT_FAILURE);
		}
		return i;
	}
}

static char *line;
static size_t len;

int get_region(FILE *f1, int chr1pos, int start1pos, int end1pos, int max1, char *chromosome1, int *start1, int *end1) {
	char *linepos = NULL,*scanpos = NULL;
	ssize_t read;
	int count;
	while ((read = getline(&line, &len, f1)) != -1) {
		linepos = line;
		count = 0;
		while (*linepos && (count <= max1)) {
			if (*linepos == '\t' || (linepos == line)) {
				if (*linepos == '\t') {
					scanpos = linepos+1;
				} else {
					scanpos = linepos;
				}
				if (count == chr1pos) {
					sscanf(scanpos,"%s\t",chromosome1);
				} else if (count == start1pos) {
					sscanf(scanpos,"%d",start1);
				} else if (count == end1pos) {
					sscanf(scanpos,"%d",end1);
				}
				count++;
			}
			linepos++;
		}
		if (count > max1) break;
	}
	if (read == -1) {return 1;}
	return 0;
}

int main(int argc, char *argv[]) {
	FILE *f1;
	char *line = NULL;
	char chromosome1[10],keepchromosome[10];
	int64_t size=0,totsize=0;
	int chr1pos,start1pos,end1pos,max1;
	int nchr1=0,start1,end1;
	int error,curchr=0,nextpos=0;
	if ((argc != 5)) {
		fprintf(stderr,"Format is: covered file1 chrpos1 startpos1 endpos1");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64(argv[1],"r");
	chr1pos = atoi(argv[2]);
	start1pos = atoi(argv[3]);
	end1pos = atoi(argv[4]);
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	fprintf(stdout,"chromosome\tbases\n");
	getline(&line, &len, f1);
	while (1) {
		error = get_region(f1,chr1pos,start1pos,end1pos,max1,chromosome1,&start1,&end1);
		if (error) break;
		nchr1 = chromosomenum(chromosome1);
		if (curchr == 0) {
			strcpy(keepchromosome,chromosome1);
			curchr = 1;
		}
		if (nchr1 > curchr) {
			totsize += size;
			fprintf(stdout,"%s\t%lld\n", keepchromosome, size);
			size = 0;
			curchr = nchr1;
			nextpos = 0;
			strcpy(keepchromosome,chromosome1);
		}
		if (start1 >= nextpos) {
			fprintf(stderr, "%s-%d\n",chromosome1,start1);
			fflush(stderr);
			nextpos += 10000000;
		}
		size += end1 - start1;
	}
	totsize += size;
	fprintf(stdout,"%s\t%lld\n", keepchromosome, size);
	fprintf(stdout,"\ntotal\t%lld\n", totsize);
	fclose(f1);
	if (line) {
		free(line);
	}
	exit(EXIT_SUCCESS);
}

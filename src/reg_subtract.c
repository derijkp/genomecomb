
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
	FILE *f1,*f2;
	char *line = NULL,*linepos = NULL;
	char chromosome1[10],chromosome2[10];
	int chr1pos,start1pos,end1pos,chr2pos,start2pos,end2pos,max1,max2;
	int nchr1=0,start1,end1,nchr2=0,start2,end2;
	int error,error2,curchr=0,nextpos=0;
	if ((argc != 9)) {
		fprintf(stderr,"Format is: reg_subtract file1 chrpos1 startpos1 endpos1 file2 chrpos2 startpos2 endpos2");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64(argv[1],"r");
	chr1pos = atoi(argv[2]);
	start1pos = atoi(argv[3]);
	end1pos = atoi(argv[4]);
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	f2 = fopen64(argv[5],"r");
	chr2pos = atoi(argv[6]);
	start2pos = atoi(argv[7]);
	end2pos = atoi(argv[8]);
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	fprintf(stdout,"chromosome\tbegin\tend\n");
	getline(&line, &len, f2);
	getline(&line, &len, f1);
	error = get_region(f1,chr1pos,start1pos,end1pos,max1,chromosome1,&start1,&end1);
	if (error) {exit(0);}
	nchr1 = chromosomenum(chromosome1);
	error2 = get_region(f2,chr2pos,start2pos,end2pos,max2,chromosome2,&start2,&end2);
	nchr2 = chromosomenum(chromosome2);
	while (1) {
/*
fprintf(stdout,"----- %d\t%s\t%d\t%d\n",1,chromosome1,start1,end1);
fprintf(stdout,"--------- %d\t%s\t%d\t%d\n",2,chromosome2,start2,end2);
*/
		if (nchr1 > curchr) {
			curchr = nchr1;
			nextpos = 0;
		}
		if (start1 >= nextpos) {
			fprintf(stderr, "%s-%d\n",chromosome1,start1);
			fflush(stderr);
			nextpos += 10000000;
		}
		linepos = line;
		if ((nchr2 < nchr1) || ((nchr2 == nchr1) && (end2 <= start1))) {
			error2 = get_region(f2,chr2pos,start2pos,end2pos,max2,chromosome2,&start2,&end2);
			if (error2)  {
				fprintf(stdout,"%s\t%d\t%d\n", chromosome1,start1,end1);
				break;
			}
			nchr2 = chromosomenum(chromosome2);
		} else if ((nchr1 < nchr2) || ((nchr1 == nchr2) && (end1 <= start2))) {
			fprintf(stdout,"%s\t%d\t%d\n", chromosome1,start1,end1);
			error = get_region(f1,chr1pos,start1pos,end1pos,max1,chromosome1,&start1,&end1);
			if (error) break;
			nchr1 = chromosomenum(chromosome1);
		} else {
			if (start2 > start1) {
				fprintf(stdout,"%s\t%d\t%d\n", chromosome1,start1,start2);
			}
			if (end2 >= end1) {
				error = get_region(f1,chr1pos,start1pos,end1pos,max1,chromosome1,&start1,&end1);
				if (error) break;
				nchr1 = chromosomenum(chromosome1);
			} else {
				start1 = end2;
				error2 = get_region(f2,chr2pos,start2pos,end2pos,max2,chromosome2,&start2,&end2);
				if (error2)  {
					fprintf(stdout,"%s\t%d\t%d\n", chromosome1,start1,end1);
					break;
				}
				nchr2 = chromosomenum(chromosome2);
			}
		}
	}
	while (!get_region(f1,chr1pos,start1pos,end1pos,max1,chromosome1,&start1,&end1)) {
		fprintf(stdout,"%s\t%d\t%d\n", chromosome1,start1,end1);
	}
	fclose(f1);
	fclose(f2);
	if (line) {
		free(line);
	}
	exit(EXIT_SUCCESS);
}

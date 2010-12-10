
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

/*static char *line;*/
/*static size_t len;*/

int get_region(
	char **line,size_t *len,
	FILE *f1, int chr1pos, int start1pos, int end1pos, int max1, char *chromosome1, int *start1, int *end1,
	int data1pos, int data2pos, char**data1, char **data2
) {
	char *linepos = NULL,*scanpos = NULL;
	ssize_t read;
	int count;
	while ((read = getline(line, len, f1)) != -1) {
		linepos = *line;
		count = 0;
		while (*linepos && (count <= max1)) {
			if (*linepos == '\t' || (linepos == *line)) {
				if (*linepos == '\t') {
					scanpos = linepos+1;
					*linepos = '\0';
				} else {
					scanpos = linepos;
				}
				if (count == chr1pos) {
					sscanf(scanpos,"%s\t",chromosome1);
				} else if (count == start1pos) {
					sscanf(scanpos,"%d",start1);
				} else if (count == end1pos) {
					sscanf(scanpos,"%d",end1);
				} else if (count == data1pos) {
					*data1 = scanpos;
				} else if (count == data2pos) {
					*data2 = scanpos;
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
	size_t len1=0,len2=0;
	char *line1 = NULL,*line2 = NULL,*data1 = NULL,*data2 = NULL;
	char chromosome1[10],chromosome2[10];
	int chr1pos,start1pos,end1pos,chr2pos,start2pos,end2pos,data1pos,data2pos,max1,max2;
	int nchr1=0,start1,end1,nchr2=0,start2,end2;
	int error2,curchr=0,nextpos=0;
	if ((argc != 11)) {
		fprintf(stderr,"Format is: reg_annot file1 chrpos1 startpos1 endpos1 file2 chrpos2 startpos2 endpos2 data1pos data2pos");
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
	data1pos = atoi(argv[9]);
	data2pos = atoi(argv[10]);
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	if (data1pos > max2) {max2 = data1pos;} ; if (data2pos > max2) {max2 = data2pos;} ;
	max2+=1;
	getline(&line1, &len1, f2);
	getline(&line1, &len1, f1);
	error2 = get_region(&line2,&len2,f2,chr2pos,start2pos,end2pos,max2,chromosome2,&start2,&end2,data1pos,data2pos,&data1,&data2);
	nchr2 = chromosomenum(chromosome2);
	while (!get_region(&line1,&len1,f1,chr1pos,start1pos,end1pos,max1,chromosome1,&start1,&end1,-1,-1,NULL,NULL)) {
/*
fprintf(stdout,"----- %d\t%s\t%d\t%d\n",1,chromosome1,start1,end1);
fprintf(stdout,"--------- %d\t%s\t%d\t%d\n",2,chromosome2,start2,end2);
*/
		nchr1 = chromosomenum(chromosome1);
		if (nchr1 > curchr) {
			curchr = nchr1;
			nextpos = 0;
		}
		if (start1 >= nextpos) {
			fprintf(stderr, "%s-%d\n",chromosome1,start1);
			fflush(stderr);
			nextpos += 50000000;
		}
		while (!error2 && ((nchr2 < nchr1) || ((nchr2 == nchr1) && (end2 <= start1)))) {
			error2 = get_region(&line2,&len2,f2,chr2pos,start2pos,end2pos,max2,chromosome2,&start2,&end2,data1pos,data2pos,&data1,&data2);
			if (error2)  {break;}
			nchr2 = chromosomenum(chromosome2);
		}
		if (error2 || (nchr1 < nchr2) || ((nchr1 == nchr2) && (end1 <= start2))) {
			if (data2pos == -1) {
				fprintf(stdout,"\n");
			} else {
				fprintf(stdout,"\t\n");
			}
		} else {
			if (data1pos == -1) {
				fprintf(stdout,"1\n");
			} else if (data2pos == -1) {
				fprintf(stdout,"%s\n", data1);
			} else {
				fprintf(stdout,"%s\t%s\n", data1, data2);
			}
		}
	}
	fclose(f1);
	fclose(f2);
	if (line1) {free(line1);}
	if (line1) {free(line2);}
	exit(EXIT_SUCCESS);
}

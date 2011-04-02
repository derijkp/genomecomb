
#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

int get_tab(
	char **line,size_t *len,
	FILE *f1, int max, char **result
) {
	char *linepos = NULL;
	ssize_t read;
	int count;
	while ((read = getline(line, len, f1)) != -1) {
		linepos = *line;
		count = 0;
		result[count] = linepos;
		while (*linepos) {
			if (*linepos == '\0') break;
			if (*linepos == '\n') {
				*linepos = '\0';
				break;
			} else if (*linepos == '\t') {
				*linepos = '\0';
				count++;
				if (count > max) break;
				result[count] = linepos+1;
			}
			linepos++;
		}
		if (count >= max) break;
	}
	if (read == -1) {return 1;}
	return 0;
}

int main(int argc, char *argv[]) {
	FILE *f1,*f2;
	char **result1=NULL,**result2=NULL,**resultkeep=NULL,**resulttemp=NULL;
	size_t len1=0,len2=0,lenkeep=0,lentemp=0;
	char *line1 = NULL,*line2 = NULL,*linekeep = NULL,*linetemp = NULL,*data1 = NULL,*data2 = NULL;
	char chromosome1[10],chromosome2[10];
	int chr1pos,start1pos,end1pos,chr2pos,start2pos,end2pos,data1pos,data2pos,max1,max2;
	int nchrkeep=0,startkeep,endkeep,near,near2;
	int nchr1=0,start1,end1,nchr2=0,start2,end2;
	int error2,curchr=0,nextpos=0,datanear=-1;
	if (argc == 12) {
		datanear = atoi(argv[11]);
	} else if ((argc != 11)) {
		fprintf(stderr,"Format is: reg_annot file1 chrpos1 startpos1 endpos1 file2 chrpos2 startpos2 endpos2 data1pos data2pos ?datanear?");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64(argv[1],"r");
	chr1pos = atoi(argv[2]);
	start1pos = atoi(argv[3]);
	end1pos = atoi(argv[4]);
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	result1 = (char **)malloc((max1+1)*sizeof(char *));
	f2 = fopen64(argv[5],"r");
	chr2pos = atoi(argv[6]);
	start2pos = atoi(argv[7]);
	end2pos = atoi(argv[8]);
	data1pos = atoi(argv[9]);
	data2pos = atoi(argv[10]);
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	if (data1pos > max2) {max2 = data1pos;} ; if (data2pos > max2) {max2 = data2pos;} ;
	result2 = (char **)malloc((max2+1)*sizeof(char *));
	resultkeep = (char **)malloc((max2+1)*sizeof(char *));
	getline(&line1, &len1, f1);
	getline(&line2, &len2, f2);
	error2 = get_tab(&line2,&len2,f2,max2,result2);
	sscanf(result2[chr2pos],"%9s\t",chromosome2);
	sscanf(result2[start2pos],"%d",&start2);
	sscanf(result2[end2pos],"%d",&end2);
	if (data1pos != -1) {data1 = result2[data1pos];}
	if (data2pos != -1) {data2 = result2[data2pos];}
	nchr2 = chromosomenum(chromosome2);
	while (!get_tab(&line1,&len1,f1,max1,result1)) {
		sscanf(result1[chr1pos],"%9s\t",chromosome1);
		sscanf(result1[start1pos],"%d",&start1);
		sscanf(result1[end1pos],"%d",&end1);
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
			/* keep data of previous */
			/* to avoid allocating new memory everytime, reuse linekeep and associated data */
			nchrkeep = nchr2; startkeep = start2; endkeep = end2;
			linetemp = linekeep;
			lentemp = lenkeep;
			resulttemp = resultkeep;
			linekeep = line2;
			lenkeep = len2;
			resultkeep = result2;
			line2 = linetemp;
			len2 = lentemp;
			result2 = resulttemp;
			/* get new line */
			error2 = get_tab(&line2,&len2,f2,max2,result2);
			if (error2)  {
				nchr2 = -1;
				break;
			}
			sscanf(result2[chr2pos],"%9s\t",chromosome2);
			sscanf(result2[start2pos],"%d",&start2);
			sscanf(result2[end2pos],"%d",&end2);
			if (data1pos != -1) {data1 = result2[data1pos];}
			if (data2pos != -1) {data2 = result2[data2pos];}
			nchr2 = chromosomenum(chromosome2);
		}
		if (error2 || (nchr1 < nchr2) || ((nchr1 == nchr2) && (end1 <= start2))) {
			if (datanear != -1) {
				near = datanear;
				if (nchr1 == nchr2) {
					near = start2-end1;
					if (data1pos != -1) {data1 = result2[data1pos];}
					if (data2pos != -1) {data2 = result2[data2pos];}
				}
				if (nchr1 == nchrkeep) {
					near2 = start1-endkeep;
					if (near2 < near) {
						near = near2;
						if (data1pos != -1) {data1 = resultkeep[data1pos];}
						if (data2pos != -1) {data2 = resultkeep[data2pos];}
					}
				}
				if (near >= datanear) {
					if (data1pos == -1) {
						fprintf(stdout,"\n");
					} else if (data2pos == -1) {
						fprintf(stdout,"\t\n");
					} else {
						fprintf(stdout,"\t\t\n");
					}
				} else {
					if (data1pos == -1) {
						fprintf(stdout,"%d\n",near);
					} else if (data2pos == -1) {
						fprintf(stdout,"%s\t%d\n", data1,near);
					} else {
						fprintf(stdout,"%s\t%s\t%d\n", data1, data2, near);
					}
				}
			} else {
				if (data2pos == -1) {
					fprintf(stdout,"\n");
				} else {
					fprintf(stdout,"\t\n");
				}
			}
		} else {
			if (datanear != -1) {
				near = end1-end2;
				if (near >= 0) {near = -1;}
				near2 = start2-start1-1;
				if (near2 >= 0) {near2 = -1;}
				if (near2 > near) {near = near2;}
				if (data1pos == -1) {
					fprintf(stdout,"%d\n",near);
				} else if (data2pos == -1) {
					fprintf(stdout,"%s\t%d\n", data1,near);
				} else {
					fprintf(stdout,"%s\t%s\t%d\n", data1, data2, near);
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
	}
	fclose(f1);
	fclose(f2);
	if (line1) {free(line1);}
	if (line2) {free(line2);}
	if (linekeep) {free(linekeep);}
	if (result1) {free(result1);}
	if (result2) {free(result2);}
	if (resultkeep) {free(resultkeep);}
	exit(EXIT_SUCCESS);
}

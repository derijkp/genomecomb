
#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
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

void output_resultline(FILE *f2,char **result1,DString *listr,int *sumr,int *grouppos,int groupnum,int *listpos,int listnum,int *sumpos,int sumnum) {
	char *sep = "";
	int i;
	if (groupnum) {
		for(i = 0 ; i < groupnum ; i++) {
			fprintf(f2,"%s%s",sep,result1[grouppos[i]]);
			sep = "\t";
		}
	}
	if (sumnum) {
		for(i = 0 ; i < sumnum ; i++) {
			fprintf(f2,"%s%d",sep,sumr[i]);
			sep = "\t";
		}
	}
	if (listnum) {
		for(i = 0 ; i < listnum ; i++) {
			fprintf(f2,"%s%s",sep,listr[i].string);
			sep = "\t";
		}
	}
	fprintf(f2,"\n");
}

int main(int argc, char *argv[]) {
	FILE *f1,*f2;
	int *grouppos=NULL, *listpos=NULL, *sumpos=NULL;
	int groupnum=0, listnum=0, sumnum=0;
	DString *listr=NULL;
	int *sumr=NULL;
	int *group1lens=NULL;
	char **result1=NULL,**result2=NULL,**tempresult;
	size_t len1=0,len2=0;
	char *line1 = NULL,*line2 = NULL,*templine = NULL;
	int error1,i,max=0,match,count=0,next=1000000;
	f1 = stdin; f2 = stdout;
	if (argc != 4) {
		fprintf(stderr,"Format is: groupby grouppos listpos sumpos\n");
		exit(EXIT_FAILURE);
	}
	/* fprintf(stdout,"\n-----\n%s,%s,%s\n-----\n",argv[1],argv[2],argv[3]); */
	parse_pos(argv[1],&grouppos,&groupnum);
	parse_pos(argv[2],&listpos,&listnum);
	parse_pos(argv[3],&sumpos,&sumnum);
	for(i = 0 ; i < groupnum ; i++) {if (grouppos[i] > max) {max = grouppos[i];}}
	for(i = 0 ; i < listnum ; i++) {if (listpos[i] > max) {max = listpos[i];}}
	for(i = 0 ; i < sumnum ; i++) {if (sumpos[i] > max) {max = sumpos[i];}}
	sumr = (int *)malloc(sumnum*sizeof(int));
	listr = (DString *)malloc(listnum*sizeof(DString));
	group1lens = (int *)malloc(groupnum*sizeof(int));
	result1 = (char **)malloc((max+1)*sizeof(char *));
	result2 = (char **)malloc((max+1)*sizeof(char *));
	/* ----- initialise from first ----- */
	error1 = get_tab(&line1,&len1,f1,max,result1);
	if (error1) {exit(EXIT_SUCCESS);}
	for(i = 0 ; i < groupnum ; i++) {
		group1lens[i] = strlen(result1[grouppos[i]]);
	}
	for(i = 0 ; i < listnum ; i++) {
		DStringInit(listr+i);
		DStringSet(listr+i,result1[listpos[i]]);
	}
	for(i = 0 ; i < sumnum ; i++) {
		sumr[i] = atoi(result1[sumpos[i]]);
	}
	/* ----- loop ----- */
	while (!get_tab(&line2,&len2,f1,max,result2)) {
		count++;
		if (count >= next) {
			fprintf(stderr,"%d\n",count);
			next += 1000000;
		}
		/* ----- check if matches ----- */
		match = 1;
		for(i = 0 ; i < groupnum ; i++) {
			if (strlen(result2[grouppos[i]]) != group1lens[i]) {
				match = 0;
				break;
			}
		}
		if (match == 1) {
			for(i = 0 ; i < groupnum ; i++) {
				if (strncmp(result2[grouppos[i]], result1[grouppos[i]], group1lens[i]) != 0) {
					match = 0;
					break;
				}
			}
		}
		/* ----- write and start new on non-match or add on match ----- */
		if (!match) {
			output_resultline(f2,result1,listr,sumr,grouppos,groupnum,listpos,listnum,sumpos,sumnum);
			for(i = 0 ; i < listnum ; i++) {
				DStringSet(listr+i,"");
			}
			for(i = 0 ; i < sumnum ; i++) {
				sumr[i] = 0;
			}
			/* put result2 in result1 */
			templine = line1; line1 = line2; line2 = templine;
			i = len1; len1 = len2; len2 = i;
			tempresult = result1; result1 = result2; result2 = tempresult;
			for(i = 0 ; i < groupnum ; i++) {
				group1lens[i] = strlen(result1[grouppos[i]]);
			}
			for(i = 0 ; i < listnum ; i++) {
				DStringSet(listr+i,result1[listpos[i]]);
			}
			for(i = 0 ; i < sumnum ; i++) {
				sumr[i] = atoi(result1[sumpos[i]]);
			}
		} else {
			for(i = 0 ; i < listnum ; i++) {
				DStringAppend(listr+i,",");
				DStringAppend(listr+i,result2[listpos[i]]);
			}
			for(i = 0 ; i < sumnum ; i++) {
				sumr[i] += atoi(result2[sumpos[i]]);
			}
		}
	}
	output_resultline(f2,result1,listr,sumr,grouppos,groupnum,listpos,listnum,sumpos,sumnum);
	exit(EXIT_SUCCESS);
}

/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>
#include "tools.h"

void output_resultline(FILE *f2,DStringArray *result1,DString *listr,long long int *sumr,int *grouppos,int groupnum,int *listpos,int listnum,int *sumpos,int sumnum,double *statsr,int statsnum) {
	char *sep = "";
	int i;
	if (groupnum) {
		for(i = 0 ; i < groupnum ; i++) {
			fprintf(f2,"%s%s",sep,result1->data[grouppos[i]].string);
			sep = "\t";
		}
	}
	if (sumnum) {
		for(i = 0 ; i < sumnum ; i++) {
			fprintf(f2,"%s%lld",sep,sumr[i]);
			sep = "\t";
		}
	}
	if (listnum) {
		for(i = 0 ; i < listnum ; i++) {
			fprintf(f2,"%s%s",sep,listr[i].string);
			sep = "\t";
		}
	}
	if (statsnum) {
		for(i = 0 ; i < statsnum ; i++) {
			fprintf(f2,"%s%.20g\t%.20g\t%d\t%.20g",sep,statsr[4*i],statsr[4*i+1],(int)statsr[4*i+2],statsr[4*i+3]);
			sep = "\t";
		}
	}
	fprintf(f2,"\n");
}

int main(int argc, char *argv[]) {
	FILE *f1,*f2;
	int *grouppos=NULL, *listpos=NULL, *sumpos=NULL, *statspos=NULL;
	int groupnum=0, listnum=0, sumnum=0, statsnum=0;
	DString *listr=NULL;
	long long int *sumr=NULL;
	double *statsr=NULL;
	int *group1lens=NULL;
	DString *line1 = NULL,*line2 = NULL,*templine = NULL;
	DStringArray *result1=NULL,*result2=NULL,*tempresult;
	size_t len1=0,len2=0;
	int error1,i,max=0,match,count=0;
	f1 = stdin; f2 = stdout;
	if (argc != 5) {
		fprintf(stderr,"Format is: groupby grouppos listpos sumpos statspos\n");
		exit(EXIT_FAILURE);
	}
	/* fprintf(stdout,"\n-----\n%s,%s,%s\n-----\n",argv[1],argv[2],argv[3]); */
	parse_pos(argv[1],&grouppos,&groupnum);
	parse_pos(argv[2],&listpos,&listnum);
	parse_pos(argv[3],&sumpos,&sumnum);
	parse_pos(argv[4],&statspos,&statsnum);
	for(i = 0 ; i < groupnum ; i++) {if (grouppos[i] > max) {max = grouppos[i];}}
	for(i = 0 ; i < listnum ; i++) {if (listpos[i] > max) {max = listpos[i];}}
	for(i = 0 ; i < sumnum ; i++) {if (sumpos[i] > max) {max = sumpos[i];}}
	for(i = 0 ; i < statsnum ; i++) {if (statspos[i] > max) {max = statspos[i];}}
	line1 = DStringNew(); line2=DStringNew();
	result1 = DStringArrayNew(max+2);
	result2 = DStringArrayNew(max+2);
	sumr = (long long int *)malloc(sumnum*sizeof(long long int));
	statsr = (double *)malloc(4*statsnum*sizeof(long long int));
	listr = (DString *)malloc(listnum*sizeof(DString));
	group1lens = (int *)malloc(groupnum*sizeof(int));
	/* ----- initialise from first ----- */
	error1 = DStringGetTab(line1,f1,max,result1,1,NULL);
	if (error1) {exit(EXIT_SUCCESS);}
	for(i = 0 ; i < groupnum ; i++) {
		group1lens[i] = result1->data[grouppos[i]].size;
	}
	for(i = 0 ; i < listnum ; i++) {
		DStringInit(listr+i);
		DStringCopy(listr+i,result1->data+listpos[i]);
	}
	for(i = 0 ; i < sumnum ; i++) {
		if (sumpos[i] == -1) {
			sumr[i] = 1;
		} else {
			sumr[i] = atoll(result1->data[sumpos[i]].string);
		}
	}
	for(i = 0 ; i < statsnum ; i++) {
		if (statspos[i] == -1) {
			statsr[4*i] = -1;
			statsr[4*i+1] = -1;
			statsr[4*i+2] = -1;
			statsr[4*i+3] = -1;
		} else {
			statsr[4*i] = atof(result1->data[statspos[i]].string);
			statsr[4*i+1] = statsr[4*i];
			statsr[4*i+2] = 1;
			statsr[4*i+3] = statsr[4*i];
		}
	}
	/* ----- loop ----- */
	while (!DStringGetTab(line2,f1,max,result2,1,NULL)) {
		count++;
		/* ----- check if matches ----- */
		match = 1;
		for(i = 0 ; i < groupnum ; i++) {
			if (result2->data[grouppos[i]].size != group1lens[i]) {
				match = 0;
				break;
			}
		}
		if (match == 1) {
			for(i = 0 ; i < groupnum ; i++) {
				if (strncmp(result2->data[grouppos[i]].string, result1->data[grouppos[i]].string, group1lens[i]) != 0) {
					match = 0;
					break;
				}
			}
		}
		/* ----- write and start new on non-match or add on match ----- */
		if (!match) {
			output_resultline(f2,result1,listr,sumr,grouppos,groupnum,listpos,listnum,sumpos,sumnum,statsr,statsnum);
			for(i = 0 ; i < listnum ; i++) {
				DStringSet(listr+i,"");
			}
			for(i = 0 ; i < sumnum ; i++) {
				sumr[i] = 0;
			}
			for(i = 0 ; i < statsnum ; i++) {
				statsr[4*i] = atof(result2->data[statspos[i]].string);
				statsr[4*i+1] = statsr[4*i];
				statsr[4*i+2] = 1;
				statsr[4*i+3] = statsr[4*i];
			}
			/* put result2 in result1 */
			templine = line1; line1 = line2; line2 = templine;
			i = len1; len1 = len2; len2 = i;
			tempresult = result1; result1 = result2; result2 = tempresult;
			for(i = 0 ; i < groupnum ; i++) {
				group1lens[i] = result1->data[grouppos[i]].size;
			}
			for(i = 0 ; i < listnum ; i++) {
				DStringCopy(listr+i,result1->data+listpos[i]);
			}
			for(i = 0 ; i < sumnum ; i++) {
				if (sumpos[i] == -1) {
					sumr[i] = 1;
				} else {
					sumr[i] = atoll(result1->data[sumpos[i]].string);
				}
			}
		} else {
			for(i = 0 ; i < listnum ; i++) {
				DStringAppendS(listr+i,",",1);
				DStringAppendS(listr+i,result2->data[listpos[i]].string,result2->data[listpos[i]].size);
			}
			for(i = 0 ; i < sumnum ; i++) {
				if (sumpos[i] == -1) {
					sumr[i] += 1;
				} else {
					sumr[i] += atoll(result2->data[sumpos[i]].string);
				}
			}
			for(i = 0 ; i < statsnum ; i++) {
				double newval = atof(result2->data[statspos[i]].string);
				if (newval < statsr[4*i]) statsr[4*i] = newval;
				statsr[4*i+1] += newval;
				statsr[4*i+2] += 1;
				if (newval > statsr[4*i+3]) statsr[4*i+3] = newval;
			}
		}
	}
	output_resultline(f2,result1,listr,sumr,grouppos,groupnum,listpos,listnum,sumpos,sumnum,statsr,statsnum);
	exit(EXIT_SUCCESS);
}

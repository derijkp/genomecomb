#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

int pos2bin(int start,int end) {
	int pre = 4681;
	start = start >> 14;
	end = end >> 14;
	while (1) {
		if (start == end) {
			return pre+start;
		}
		pre = pre >> 3;
		if (!pre) {return 0;}
		start = start >> 3;
		end = end >> 3;
	}
}

int main(int argc, char *argv[]) {
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	char chromosome[10],chromosomeb1[10],chromosomeb2[10],weight;
	char *chromosome1,*chromosome2;
	int flags,offsetInChr,gap1,gap2,gap3,mateRec,side,last,index=-1,cur=0;
	char strand1='?',strand2='?';
	int chr1=-1,start1=-1,end1=-1,side1=-1,weight1=-1;
	int chr2=-1,start2=-1,end2=-1,side2=-1,weight2=-1;
	int numl=0,numr=0,type,dist,min = 300,max = 600,bin,num,fnum=0;
	if (argc != 2) {
		fprintf(stderr,"Format is: map2sv num");
		exit(EXIT_FAILURE);
	}
	num = atoi(argv[1]);
	while ((read = getline(&line, &len, stdin)) != -1) {
		if (line[0] == '>') break;
	}
	while ((read = getline(&line, &len, stdin)) != -1) {
		sscanf(line,"%d\t%9s\t%d\t%d\t%d\t%d\t%c\t%d",&flags,chromosome,&offsetInChr,&gap1,&gap2,&gap3,&weight,&mateRec);
		last = flags & 0x01;
		side = flags & 0x02;
		if (side) {numr++;} else {numl++;}
		if (!cur) {
			chr1 = chromosomenum(chromosome);
			strcpy(chromosomeb1,chromosome);
			chromosome1 = chromosomeb1;
			start1 = offsetInChr;
			end1 = offsetInChr + 35 + gap1 + gap2 + gap3;
			if (flags & 0x04) {strand1 = '-';} else  {strand1 = '+';}
			side1 = side;
			weight1 = weight - 33;
			index = mateRec;
		} else if ((cur == index) || ((index == 0) && (side != side1))) {
			chr2 = chromosomenum(chromosome);
			strcpy(chromosomeb2,chromosome);
			chromosome2 = chromosomeb2;
			start2 = offsetInChr;
			end2 = offsetInChr + 35 + gap1 + gap2 + gap3;
			if (flags & 0x04) {strand2 = '-';} else  {strand2 = '+';}
			side2 = side;
			weight2 = weight - 33;
			index = -1;
		}
		cur++;
		if (last) {
			if (start2 == -1) {
				if (numl == -1) {numl = numr;side1=1;}
				bin = pos2bin(start1,end1);
				fprintf(stdout,"%s\t%d\t%c\t%d\t%d\t%d\t%d\t%c",
					chromosome1,bin,strand1,start1,end1,weight1,numl,'s');
				fprintf(stdout,"\t\t\t\t\t\t\t\t%d\t%d\t%d\n",
					num,fnum,(side1)?1:0);
			} else {
				if (end1 > start2) {
					int tempi;char tempc, *temps;
					tempi=chr1; chr1=chr2; chr2=tempi;
					temps=chromosome1; chromosome1=chromosome2; chromosome2=temps;
					tempi=start1; start1=start2; start2=tempi;
					tempi=end1; end1=end2; end2=tempi;
					tempi=side1; side1=side2; side2=tempi;
					tempi=weight1; weight1=weight2; weight2=tempi;
					tempc=strand1; strand1=strand2; strand2=tempc;
					tempi=numl; numl=numr; numr=numl;
				}
				dist = start2-end1;
				bin = pos2bin(start1,end2);
				if (chr1 != chr2) {
					type = 'c';
					dist = -1;
					bin = pos2bin(start1,end1);
				} else if (strand1 != strand2) {
					type = 'r';
				} else if (dist < min) {
					type = 'i';
				} else if (dist > max) {
					type = 'd';
				} else {
					type = 'n';
				}
				fprintf(stdout,"%s\t%d\t%c\t%d\t%d\t%d\t%d\t%c",
					chromosome1,bin,strand1,start1,end1,weight1,numl,type);
				fprintf(stdout,"\t%s\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
					chromosome2,strand2,start2,end2,weight2,numr,dist,num,fnum,(side1)?1:0);
				if (type == 'c') {
					bin = pos2bin(start2,end2);
					fprintf(stdout,"%s\t%d\t%c\t%d\t%d\t%d\t%d\t%c",
						chromosome2,bin,strand2,start2,end2,weight2,numr,type);
					fprintf(stdout,"\t%s\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
						chromosome1,strand1,start1,end1,weight1,numl,dist,num,fnum,(side2)?1:0);
				}
			}
			fnum++;
			cur = 0;
			index = -1;
			side1 = -1;
			numl = 0; numr = 0;
			start1 = -1; start2 = -1;
		}
	}
	if (line) {
		free(line);
	}
	exit(EXIT_SUCCESS);
}

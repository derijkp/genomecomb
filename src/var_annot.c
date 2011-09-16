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
#include "tools.h"

void connectalt(
	char *alt1, char *alt2, char *data
) {
	char *alt1keep,*alt2keep,*alt2move,*datakeep;
	int alt2num,found,i,pre = 0;
	alt1keep = alt1;
	while (1) {
		if (pre) {
			fprintf(stdout,",");
		} else {
			pre = 1;
		}
		while (*alt1 != ',' && *alt1 != '\0') {
			alt1++;
		}
		alt2keep = alt2;
		alt2move = alt2;
		alt2num = 0;
		found = 0;
		/*fprintf(stdout,"\n-----%.*s\n",alt1-alt1keep,alt1keep);*/
		while (1) {
			while (*alt2move != ',' && *alt2move != '\0') {
				alt2move++;
			}
			if ((alt1-alt1keep == alt2move-alt2keep) && (strncmp(alt1keep,alt2keep,alt1-alt1keep) == 0)) {
				found = 1;
				datakeep = data;
				while (alt2num) {
					if (*datakeep == '\0') break;
					if (*datakeep == ',') alt2num--;
					datakeep++;
				}
				if (!alt2num) {
					if (*datakeep == ',') datakeep++;
					i = 0;
					while (datakeep[i] != ',' && datakeep[i] != '\0') {i++;}
					if (i == 0) {
						fprintf(stdout,"-");
					} else {
						fprintf(stdout,"%.*s",i,datakeep);
					}
				}
				break;
			}
			if (*alt2move == '\0') break;
			alt2num++;
			alt2move++;
			alt2keep = alt2move;
		}
		if (!found) {
			fprintf(stdout,"-");
		}
		if (*alt1 == '\0') break;
		alt1++;
		alt1keep = alt1;
	}
}

int main(int argc, char *argv[]) {
	FILE *f1,*f2;
	DString *result1=NULL,*result2=NULL;
	DString *line1 = NULL,*line2 = NULL;
	char *chromosome1,*chromosome2;
	int chr1pos,start1pos,end1pos,type1pos,alt1pos,data1pos,max1;
	int chr2pos,start2pos,end2pos,type2pos,alt2pos,data2pos,max2;
	int nchr1=0,start1,end1;
	int nchr2=0,start2,end2;
	int error2,curchr=0,nextpos=0,sametype,cmp;
	if ((argc != 15)) {
		fprintf(stderr,"Format is: reg_annot file1 chrpos1 startpos1 endpos1 type1pos alt1pos file2 chrpos2 startpos2 endpos2 type2pos alt2pos data1pos data2pos");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64(argv[1],"r");
	chr1pos = atoi(argv[2]);
	start1pos = atoi(argv[3]);
	end1pos = atoi(argv[4]);
	type1pos = atoi(argv[5]);
	alt1pos = atoi(argv[6]);
	max1 = chr1pos;
	if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	if (type1pos > max1) {max1 = type1pos;} ; if (alt1pos > max1) {max1 = alt1pos;} ;
	line1 = DStringNew(); line2=DStringNew();
	result1 = DStringArrayNew(max1+1);
	f2 = fopen64(argv[7],"r");
	chr2pos = atoi(argv[8]);
	start2pos = atoi(argv[9]);
	end2pos = atoi(argv[10]);
	type2pos = atoi(argv[11]);
	alt2pos = atoi(argv[12]);
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	if (type2pos > max2) {max2 = type2pos;} ; if (alt2pos > max2) {max2 = alt2pos;} ;
	data1pos = atoi(argv[13]);
	data2pos = atoi(argv[14]);
	if (data1pos > max2) {max2 = data1pos;} ; if (data2pos > max2) {max2 = data2pos;} ;
	result2 = DStringArrayNew(max2+1);
	skip_header(f1,line1);
	skip_header(f2,line2);
	error2 = DStringGetTab(line2,f2,max2,result2);
	chromosome2 = result2[chr2pos].string;
	nchr2 = chromosomenum(chromosome2);
	sscanf(result2[start2pos].string,"%d",&start2);
	sscanf(result2[end2pos].string,"%d",&end2);
	while (!DStringGetTab(line1,f1,max1,result1)) {
		chromosome1 = result1[chr1pos].string;
		nchr1 = chromosomenum(chromosome1);
		sscanf(result1[start1pos].string,"%d",&start1);
		sscanf(result1[end1pos].string,"%d",&end1);
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
			nextpos += 50000000;
		}
		while (!error2) {
			sametype = 0;
			if (nchr2 == nchr1) {
				if (start2 == start1) {
					if (end2 == end1) {
						cmp = strcmp(result2[type2pos].string, result1[type1pos].string);
						if (cmp == 0) {
							sametype = 1; break;
						} else if (cmp > 0) break; 
					} else if (end2 > end1) break; 
				} else if (start2 > start1) break;
			} else if (nchr2 > nchr1) break;
			error2 = DStringGetTab(line2,f2,max2,result2);
			if (error2)  {break;}
			chromosome2 = result2[chr2pos].string;
			nchr2 = chromosomenum(chromosome2);
			sscanf(result2[start2pos].string,"%d",&start2);
			sscanf(result2[end2pos].string,"%d",&end2);
		}
		if (error2 || (nchr1 != nchr2) || (start2 != start1) || (end2 != end1) || !sametype) {
			if (data2pos == -1) {
				fprintf(stdout,"-\n");
			} else {
				fprintf(stdout,"-\t-\n");
			}
		} else {
			if (data1pos == -1) {
				fprintf(stdout,"1\n");
			} else if (data2pos == -1) {
				connectalt(result1[alt1pos].string,result2[alt2pos].string,result2[data1pos].string);
				fprintf(stdout,"\n");
				/*fprintf(stdout,"%s\n", result2[data1pos]);*/
			} else {
				connectalt(result1[alt1pos].string,result2[alt2pos].string,result2[data1pos].string);
				fprintf(stdout,"\t");
				connectalt(result1[alt1pos].string,result2[alt2pos].string,result2[data2pos].string);
				fprintf(stdout,"\n");
				/* fprintf(stdout,"%s\t%s\n", result2[data1pos].string, result2[data2pos].string); */
			}
		}
	}
	fclose(f1);
	fclose(f2);
	if (line1) {DStringDestroy(line1);}
	if (line2) {DStringDestroy(line2);}
	if (result1) {DStringArrayDestroy(result1);}
	if (result2) {DStringArrayDestroy(result2);}
	exit(EXIT_SUCCESS);
}

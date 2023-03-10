/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"
#include "debug.h"

void connectalt(
	char *alt1, char *alt2, char *data, char *notfound
) {
	char *alt1keep,*alt2keep,*alt2move,*datakeep,*prevdata;
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
				prevdata = data;
				while (alt2num) {
					if (*datakeep == '\0') break;
					if (*datakeep == ',') {
						alt2num--;
						prevdata = datakeep;
					}
					datakeep++;
				}
				if (*datakeep == '\0') {
					i = datakeep-prevdata;
					datakeep = prevdata;
				} else {
					if (*datakeep == ',') datakeep++;
					i = 0;
					while (datakeep[i] != ',' && datakeep[i] != '\0') {i++;}
				}
				if (i == 0) {
					fprintf(stdout,"%s",notfound);
				} else {
					fprintf(stdout,"%.*s",i,datakeep);
				}
				break;
			}
			if (*alt2move == '\0') break;
			alt2num++;
			alt2move++;
			alt2keep = alt2move;
		}
		if (!found) {
			fprintf(stdout,"%s",notfound);
		}
		if (*alt1 == '\0') break;
		alt1++;
		alt1keep = alt1;
	}
}

int main(int argc, char *argv[]) {
	FILE *f1,*f2;
	DStringArray *result1=NULL,*result2=NULL;
	DString *line1 = NULL,*line2 = NULL;
	DString *prevchromosome1 = DStringNew(), *prevchromosome2 = DStringNew();
	DString *prevtype1 = DStringNew(), *prevtype2 = DStringNew();
	DString *prevalt1 = DStringNew(), *prevalt2 = DStringNew();
	DString *chromosome1=NULL,*chromosome2=NULL,*type1 = NULL,*type2 = NULL,*alt1 = NULL,*alt2 = NULL;
	unsigned int numfields1,numfields2,numfields,pos1 = 0,pos2 = 0;
	char *notfound = "-", **data;
	int prevstart1 = -1,prevend1 = -1,prevstart2 = -1,prevend2 = -1;
	int chr1pos,start1pos,end1pos,type1pos,alt1pos,max1;
	int chr2pos,start2pos,end2pos,type2pos,alt2pos,max2;
	int comp;
	int datalen=0,*datapos=NULL;
	int start1,end1;
	int start2,end2;
	int error2,nextpos=0,sametype,cmp,i;
	if ((argc < 14)) {
		fprintf(stderr,"Format is: var_annot file1 chrpos1 startpos1 endpos1 type1pos alt1pos file2 chrpos2 startpos2 endpos2 type2pos alt2pos notfound datapos ...");
		exit(EXIT_FAILURE);
	}
	i = 1;
	f1 = fopen64_or_die(argv[i++],"r");
	chr1pos = atoi(argv[i++]);
	start1pos = atoi(argv[i++]);
	end1pos = atoi(argv[i++]);
	type1pos = atoi(argv[i++]);
	alt1pos = atoi(argv[i++]);
	max1 = chr1pos;
	if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	if (type1pos > max1) {max1 = type1pos;} ; if (alt1pos > max1) {max1 = alt1pos;} ;
	line1 = DStringNew(); line2=DStringNew();
	/* The following allocation is not destroyed at end as it may point to something else */
	/* This will leak mem, but as the prog is finished anyway ... */
	result1 = DStringArrayNew(max1+2);
	if (strlen(argv[i]) == 1 && argv[i][0] == '-') {
		f2 = stdin;
	} else {
		f2 = fopen64_or_die(argv[i],"r");
	}
	i++;
	chr2pos = atoi(argv[i++]);
	start2pos = atoi(argv[i++]);
	end2pos = atoi(argv[i++]);
	type2pos = atoi(argv[i++]);
	alt2pos = atoi(argv[i++]);
NODPRINT("var_annot %s %d %d %d %d %d %s %d %d %d %d %d ...",
	argv[1],chr1pos,start1pos,end1pos,type1pos,alt1pos,
	argv[7],chr2pos,start2pos,end2pos,type2pos,alt2pos
);
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	if (type2pos > max2) {max2 = type2pos;} ; if (alt2pos > max2) {max2 = alt2pos;} ;
	notfound = argv[i++];
	data = argv + i;
	datalen = argc-i;
	datapos = malloc(datalen*sizeof(int));
	for (i = 0 ; i < datalen ; i++) {
		datapos[i] = atoi(data[i]);
NODPRINT(" %d",datapos[i]);
	}
NODPRINT("\n");
	for (i = 0 ; i < datalen ; i++) {
		if (datapos[i] > max2) {max2 = datapos[i];}
	}
	result2 = DStringArrayNew(max2+2);
	skip_header(f1,line1,&numfields1,&pos1);
	skip_header(f2,line2,&numfields2,&pos2);
	error2 = DStringGetTab(line2,f2,max2,result2,1,&numfields); pos2++;
	if (!error2) {
		check_numfieldserror(numfields2,numfields,line2,argv[7],&pos2);
		chromosome2 = result2->data+chr2pos;
		type2 = result2->data+type2pos;
		alt2 = result2->data+alt2pos;
	}
	sscanf(result2->data[start2pos].string,"%d",&start2);
	sscanf(result2->data[end2pos].string,"%d",&end2);
NODPRINT("line2 %s,%d,%d %s",Loc_ChrString(chromosome2),start2,end2,line2->string)
	while (!DStringGetTab(line1,f1,max1,result1,1,&numfields)) {
		pos1++;
		check_numfieldserror(numfields,numfields1,line1,argv[1],&pos1);
		chromosome1 = result1->data+chr1pos;
		sscanf(result1->data[start1pos].string,"%d",&start1);
		sscanf(result1->data[end1pos].string,"%d",&end1);
		type1 = result1->data+type1pos;
		alt1 = result1->data+alt1pos;
NODPRINT("line1 (a=%3.3s) %s,%d,%d %s",type1->string,chromosome1->string,start1,end1,alt1->string)
/*
fprintf(stdout,"----- %d\t%s\t%d\t%d\n",1,Loc_ChrString(chromosome1),start1,end1);
fprintf(stdout,"--------- %d\t%s\t%d\t%d\n",2,Loc_ChrString(chromosome2),start2,end2);
*/
		checksort(prevchromosome1,&prevstart1,&prevend1,prevtype1,prevalt1,chromosome1,start1,end1,type1,alt1,argv[1],&nextpos,1);
/*
		if (start1 >= nextpos) {
			fprintf(stderr, "%s-%d\n",Loc_ChrString(chromosome1),start1);
			fflush(stderr);
			nextpos += 50000000;
		}
*/
		while (!error2) {
			sametype = 0;
NODPRINT("line2 %s,%d,%d %s %s",Loc_ChrString(prevchromosome2),prevstart2,prevend2,prevtype2->string,prevalt2->string)
			comp = DStringLocCompare(chromosome2,chromosome1);
			if (comp == 0) {
				if (start2 == start1) {
					if (end2 == end1) {
						cmp = DStringCompare(type2, type1);
						if (cmp == 0) {
							sametype = 1; break;
						} else if (cmp > 0) break; 
					} else if (end2 > end1) break; 
				} else if (start2 > start1) break;
			} else if (comp > 0) break;
			error2 = DStringGetTab(line2,f2,max2,result2,1,&numfields); pos2++;
			if (error2)  {break;} else {
				check_numfieldserror(numfields,numfields2,line2,argv[7],&pos2);
			}
			chromosome2 = result2->data+chr2pos;
			sscanf(result2->data[start2pos].string,"%d",&start2);
			sscanf(result2->data[end2pos].string,"%d",&end2);
			type2 = result2->data+type2pos;
			alt2 = result2->data+alt2pos;
			checksort(prevchromosome2,&prevstart2,&prevend2,prevtype2,prevalt2,chromosome2,start2,end2,type2,alt2,argv[7],&nextpos,1);
		}
		if (error2 || (DStringLocCompare(chromosome2,chromosome1) != 0) || (start2 != start1) || (end2 != end1) || !sametype) {
			if (!datalen) {
				fprintf(stdout,"%s\n",notfound);
			} else {
				for (i = 1 ; i < datalen ; i++) {
					fprintf(stdout,"%s\t",notfound);
				}
				fprintf(stdout,"%s\n",notfound);
			}
		} else {
			if (!datalen) {
				connectalt(alt1->string,alt2->string,"1",notfound);
				putc_unlocked('\n',stdout);
			} else {
				connectalt(alt1->string,alt2->string,result2->data[datapos[0]].string,notfound);
				for (i = 1 ; i < datalen ; i++) {
					fprintf(stdout,"\t");
					if (datapos[i] != -1) {
						connectalt(alt1->string,alt2->string,result2->data[datapos[i]].string,notfound);
					}
				}
				fprintf(stdout,"\n");
			}
		}
	}
	FCLOSE(f1);
	FCLOSE(f2);
	DStringDestroy(prevchromosome1); DStringDestroy(prevtype1); DStringDestroy(prevalt1);
	DStringDestroy(prevchromosome2); DStringDestroy(prevtype2); DStringDestroy(prevalt2);
	if (datapos) {free(datapos);}
	if (line1) {DStringDestroy(line1);}
	if (line2) {DStringDestroy(line2);}
	if (result1) {DStringArrayDestroy(result1);}
	if (result2) {DStringArrayDestroy(result2);}
	exit(EXIT_SUCCESS);
}

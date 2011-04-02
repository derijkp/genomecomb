#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	char chromosome[10],weight;
	int flags,offsetInChr,gap1,gap2,gap3,mateRec,end,side,strand,last,index=-1,cur=0,side1 = -1;
	
	if (argc != 1) {
		fprintf(stderr,"Format is: map2besthit");
		exit(EXIT_FAILURE);
	}
	while ((read = getline(&line, &len, stdin)) != -1) {
		if (line[0] == '>') break;
	}
	while ((read = getline(&line, &len, stdin)) != -1) {
		sscanf(line,"%d\t%9s\t%d\t%d\t%d\t%d\t%c\t%d",&flags,chromosome,&offsetInChr,&gap1,&gap2,&gap3,&weight,&mateRec);
		last = flags & 0x01;
		side = flags & 0x02;
		if (!cur) {
			strand = flags & 0x04;
			end = offsetInChr + 35 + gap1 + gap2 + gap3;
			fprintf(stdout,"%s\t%d\n", chromosome,(offsetInChr+end)/2);
			index = mateRec;
			side1 = side;
		} else if ((cur == index) || ((index == 0) && (side != side1))) {
			strand = flags & 0x04;
			end = offsetInChr + 35 + gap1 + gap2 + gap3;
			fprintf(stdout,"%s\t%d\n", chromosome,(offsetInChr+end)/2);
			index = -1;
		}
		cur++;
		if (last) {
			cur = 0;
			index = -1;
		}
	}
	if (line) {
		free(line);
	}
	exit(EXIT_SUCCESS);
}

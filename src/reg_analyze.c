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

void perror(const char *s);
void annotate (char *out1, char *in, int in_chrom, int in_begin, int in_end, char *db, int chrom, int begin, int end, int name, int score, int allel_freq);
char** split (char *str[], char delims[], int NUM);
void makeLine(FILE *output, int name, int score, int allel_freq, char *name_annot, char * score_annot, char * allel_annot);
void freeMemoryArr (char ** dArray, int NUM);

 

int main (int argc, char *argv[]) {
	
	if ((argc < 7)) {
		perror("Format is: directory inputfile, fileName, chrom_column, begin_column, end_column, databaseColumnList");
		exit(EXIT_FAILURE);
	}
	int input_chrom_column = atoi(argv[3]);
	int input_begin_column = atoi(argv[4]);
	int input_end_column = atoi(argv[5]);
	char *dbfile[20];
	int chrom[20];
	int begin[20];
	int end[20];
	int name[20];
	int score[20];
 	int allel_freq[20];

	//assigning the column data to variables
	int i = 6;
	int k = 0;
	while (i < argc-1) {
		if (argv[i][0] == '/') {
			dbfile[k] = argv[i];
			i++;
		} else {
			perror ("First argument in input_c from Aregio.tcl has to be a fileName");
			exit(EXIT_FAILURE);
		}
		while ((argv[i][0] != '/') && (argv[i] != NULL)) {
			chrom[k] = atoi( argv[i]);
			i++;
			begin[k] = atoi(argv[i]);
			i++;
			end[k] = atoi(argv[i]);
			i++;
			name[k] = atoi(argv[i]);
			i++;
			score[k] = atoi(argv[i]);
			i++;
			allel_freq[k] = atoi(argv[i]);
			i++;
			if (argv[i] == NULL) break;
		}
		k++;
	}
	
	// Open Databases 
	
	i = 0;
	while (i < k ) {
		// Vergelijk file met huidige db
		printf("Running input against %s ..... \n",dbfile[i]);
		fflush(stdout);
		annotate (argv[1], argv[2], input_chrom_column, input_begin_column, input_end_column, dbfile[i], chrom[i], begin[i], end[i], name[i], score[i], allel_freq[i]);
		i++;
	}

	exit(EXIT_SUCCESS);
}

void annotate (char *out1, char *in, int in_chrom, int in_begin, int in_end, char *db, int chrom_col, int begin, int end, int name, int score, int allel_freq) {
	FILE *input;
	input = fopen64_or_die(in, "r+");
	FILE *db_input;
	db_input = fopen64_or_die(db, "r");
	
	// Initiate

	int bytes_read;
	int bytes_read_db;
	size_t nbytes = 0;
	size_t nbytes_db = 0;
	char *line = NULL;
	char *db_line = NULL;
	int loop;
	char **db_column = 0;
	char **column = 0;
	int begin_db;
	int end_db;
	int snp;
	int chrom_db;	
	int chrom;
	FILE *output_id;
	char *out2 = NULL;
	char **line_arr = 0;
	char *annot_col_1 = NULL;
	char *annot_col_2 = NULL;
	char *annot_col_3 = NULL;
	int NUM_db;
	int NUM_in;
	
	//Reading headers

	bytes_read = getline (&line, &nbytes, input);
	//puts ("header:");
	//puts (line);
	bytes_read_db = getline (&db_line, &nbytes_db, db_input);
	//puts ("database header:");
	//puts (db_line);

	
	// Making header and file for output
	int NUM_arr = strlen(db);
	line_arr = split(&db, "./", NUM_arr);
	char *name_db = line_arr[4];
	out2 = ".tsv";
	char *outfile = malloc(strlen(out1) + strlen(out2) + strlen(name_db) + 3);
	
	if (outfile != NULL) {
		strcpy(outfile, out1);
		strcat(outfile, name_db);
		strcat(outfile, out2);
	} 

	output_id = fopen64_or_die(outfile, "w");

	annot_col_1 = "_name";
	annot_col_2 = "_score";
	annot_col_3 = "_allel_freq";
	char *header_1 = malloc(strlen(name_db) + strlen(annot_col_1) + 2);
	char *header_2 = malloc(strlen(name_db) + strlen(annot_col_2) + 2);
	char *header_3 = malloc(strlen(name_db) + strlen(annot_col_3) + 2);
	strcpy(header_2, name_db);
	strcat(header_2, annot_col_2);
	strcpy(header_1, name_db);
	strcat(header_1, annot_col_1);
	strcpy(header_3, name_db);
	strcat(header_3, annot_col_3);
	
	makeLine(output_id, name, score, allel_freq, header_1, header_2, header_3);	
	freeMemoryArr(line_arr, NUM_arr);
	free(header_1);
	free(header_2);
	free(header_3);

	//reading first real line from both input file as db_file
	int l = 1;
	int k = 1;
	
	free(line);
	free(db_line);
	line = NULL;
	db_line = NULL;
	nbytes = 0;
	nbytes_db = 0;

	bytes_read = getline (&line, &nbytes, input);					 
	bytes_read_db = getline (&db_line, &nbytes_db, db_input); 		
	loop = 3;
	printf("From input reading line %d  ..... \n",l);
	printf("From db reading line %d  ..... \n",k);
	
	while (!feof(input) && !feof(db_input)){
		
		// Print every 100000 lines
		if (l%100000==0) {	
			printf("From input reading line %d  ..... \n",l);
		}
		if (k%100000==0) {	
			printf("From db reading line %d  ..... \n",k);
		}
		fflush(stdout);

		//when in while loop, you don't have to do everything again,
		// only the line that has changed
		// Put everything in 2 if-loops

		if (loop > 1) {	
			//Splitting the line into columns
			
			NUM_db = strlen(db_line);
			db_column = split(&db_line, "\t\n", NUM_db);
			//printf( "DB_Column 3 is %s \n", db_column[2] );
			// turning some chars into int for later use
			// Only dealing with snp for now
			// Otherwise make in_end also int
			// And re-think the annotate if annotate-loop
			begin_db = atoi(db_column[begin]);
			end_db = atoi(db_column[end]);
			chrom_db = atoi(db_column[chrom_col]);

			if (loop == 3) { 
				loop = 1; 
			}
		}
		if (loop == 1) {
			NUM_in = strlen(line);
			column = split(&line, "\t\n", NUM_in);
			snp = atoi(column[in_begin]);
			chrom = atoi(column[in_chrom]);

		}


		// Check the chromosomes
		// Make sure that there are no X,Y,M chromosomes but chrom 23, 24 and 25 (by Asort.tcl)
		// Otherwise you will lose data using this
				
		if (chrom < chrom_db) {		
			makeLine(output_id, name, score, allel_freq, "", "", "");	
			free(line);
			line = NULL;
			nbytes = 0;
			bytes_read = getline (&line, &nbytes, input);
			l++;
			loop = 1;
			freeMemoryArr(column, NUM_in);
		} else if (chrom > chrom_db) {
			free(db_line);
			db_line = NULL;
			nbytes_db = 0;
			bytes_read_db = getline (&db_line, &nbytes_db, db_input);
			k++;
			loop = 2;
			freeMemoryArr(db_column, NUM_db);
		} else {
			// Annotate
			if (snp < end_db) {
				if (snp < begin_db)		{
					makeLine(output_id, name, score, allel_freq, "", "", "");
					free(line);
					line = NULL;
					nbytes = 0;
					bytes_read = getline (&line, &nbytes, input);
					l++;
					loop = 1;
					freeMemoryArr(column, NUM_in);
				} else {
					//annoteer
					makeLine(output_id, name, score, allel_freq, db_column[name], db_column[score], db_column[allel_freq]);
					free(line);
					line = NULL;
					nbytes = 0;
					bytes_read = getline (&line, &nbytes, input);
					l++;
					loop = 1;
					freeMemoryArr(column, NUM_in);
				}
			} else {
				free(db_line);
				db_line = NULL;
				nbytes_db = 0;
				bytes_read_db = getline (&db_line, &nbytes_db, db_input);
				k++;
				loop = 2;
				freeMemoryArr(db_column, NUM_db);
			}
		}
	}

	if (loop == 2) {
		freeMemoryArr(column, NUM_in);
	} else if (loop == 1) {
		freeMemoryArr(db_column, NUM_db);
	}

	free(db_line);
	free(line);
	printf("From input read last line %d  ..... \n",l);
	printf("From db read last line %d  ..... \n",k);
	free(outfile);
	fclose(output_id);
	fclose(input);
	fclose(db_input);
}

void makeLine(FILE *output, int name, int score, int allel_freq, char *name_annot, char * score_annot, char * allel_annot) {
	if (allel_freq == 100) {
		if (name == 100 || score == 100 ) {
			if (name == 100) {
				fprintf(output, "%s\n", score_annot);
			} else if (score == 100) {
				fprintf(output, "%s\n", name_annot);
			} 
		} else {
			fprintf (output, "%s\t%s\n", name_annot, score_annot);
		}
	} else if (allel_freq != 100) {
		fprintf(output, "%s\n", allel_annot);
	}
	
}

char** split ( char *str[], char delims[], int NUM) {
	//char delims[] = "\t";

	char *part;
	char **result;
	int i;
	//Initiation array fo pointers to array
	result = (char**) malloc(NUM*sizeof(char*));    
	for(i = 0; i < NUM; i++) {
  		result[i] = (char*)malloc(NUM* sizeof(char));
	}
	if (result == NULL) {
 		printf("Memory not allocated.\n");
  		exit(EXIT_FAILURE);
 	}
	
	//Putting the different columns into result array
	i = 0;
	part = strtok (*str, delims);
	while( part != NULL ) {
		int j=0;
		while (j <= strlen(part)) {
			result[i][j] = (char) part[j];
			j++;
		}
		//printf( "result %d is %s \n", i, result[i] );
		part = strtok(NULL, delims);
		i++; 
	}
	return result;
}

void freeMemoryArr (char ** dArray, int NUM) {
	int i = 0;
	while (i < NUM) {
		free(dArray[i]);
		i++;
	}
	free(dArray);
	dArray = NULL; 
}


	

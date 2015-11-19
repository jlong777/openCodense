/*
 ReadMicroarrayData.c                                01/2008 jlong@jimlong.org

 infileName format: a header line followed by data, i.e:
 probeID        Array1Dye1      Array2Dye1      Array3Dye1      Array4Dye1 ...etc
 1415670_at     9.172966        9.236673        8.789873        8.698874   ...etc
 1415671_at     10.665967       10.677811       10.463303       10.464032  ...etc
    .               .               .               .               .
    .               .               .               .               .
   etc             etc             etc             etc             etc
 
 Copyright (C) 2008 James Long
 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CODENSE_H
  #include "codense.h"
  #define CODENSE_H
#endif

#define BUFSZ 2048

microarrayDataSet * ReadMicroarrayData(char *infileName)
{
  int i, j, cols, rows;
  char buf[BUFSZ], *tok;
  
  microarrayDataSet *maDataSet;
  FILE *infile;
  
  infile = fopen(infileName, "r");
  if(infile==NULL)
  {
    fprintf(stderr, "ReadMicroarrayData: Error opening %s, exiting...\n", infileName);
    return NULL;
  }
  
  /* count rows and columns */
  fgets(buf, BUFSZ, infile);
  tok = strtok(buf, " \t,"); /* delimiters are spaces, tabs, and commas */
  cols = -1; /* don't count the first col */
  do
  {
    cols++;
    tok = strtok(NULL, " \t,");
  } while(tok);

  rows = 0; /* don't count the first row */
  while(fgets(buf, BUFSZ, infile))
    rows++;

  rewind(infile);
  
  maDataSet = (microarrayDataSet *) malloc(sizeof(microarrayDataSet));
  if(maDataSet==NULL)
  {
    fprintf(stderr, "ReadMicroarrayData: Error allocating microarrayDataSet, exiting...\n");
    return NULL;
  }
  
  maDataSet->numCols = cols;
  maDataSet->numRows = rows;
  
  /* labels are the probeIDs */
  maDataSet->labels = (char **) malloc(rows*sizeof(char *));
  if((maDataSet->labels)==NULL)
  {
    fprintf(stderr, "ReadMicroarrayData: Error allocating label array, exiting...\n");
    return 0;
  }
  
  /* return a pointer to the data */
  maDataSet->data = (double *) malloc(cols*rows*sizeof(double));
  if((maDataSet->data)==NULL)
  {
    fprintf(stderr, "ReadMicroarrayData: Error allocating data array, exiting...\n");
    return 0;
  }
  
  i = 0;
  fgets(buf, BUFSZ, infile);   /* discard the first row */
  while(fgets(buf, BUFSZ, infile))
  {
    /* first column is label */
    tok = strtok(buf, " \t,"); 
    maDataSet->labels[i] = (char *) malloc(strlen(tok)+1);
    if((maDataSet->labels[i])==NULL)
    {
      fprintf(stderr, "ReadMicroarrayData: Error allocating label, exiting...\n");
      return 0;
    }
    strcpy(maDataSet->labels[i], tok);
    
    /* data */
    j = 0;
    do
    {
      tok = strtok(NULL, " \t,");
      if(tok) maDataSet->data[i*cols + j++] = atof(tok);
    } while(tok);
    
    i++;
  }
  
  fclose(infile);
  
  return maDataSet;
}

void DeallocateMicroarrayDataSet(microarrayDataSet *maDataSet)
{
  int i;
  
  free(maDataSet->data);
  
  for(i=0; i<maDataSet->numRows; i++)
    free(maDataSet->labels[i]);
    
  free(maDataSet->labels);
  free(maDataSet);         /**recent**/
}

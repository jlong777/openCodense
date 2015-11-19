/*
 codense.c                                     04/2015 jlong@jimlong.org

 1) read in a series of microarray data
 2) use pearson's/mutual-information techniques to generate a graph G_i for each
    microarray set, i.e. a microarrayGraphSet
 3) create a summaryGraph from the microarrayGraphSet, and create list of dense 
    subgraphs from the summaryGraph (summaryGraphDenseSubGraphSet)
 4) generate edge support matrix for each entry in summaryGraphDenseSubGraphSet
    note - steps 2 & 4 use pearson's/mutual-information as a correlation measure
 5) generate second order graphs S_j for each edge support matrix
 6) mine dense subgraphs sub_k(S_j) of the second order graphs and
 7) convert dense subgraphs into coherent dense subgraphs of the original G_i
 
 The above is an implementation of CODENSE:
 
 Hu H, Yan X, Huang Y, Han J, Zhou XJ: "Mining coherent dense subgraphs across
 massive biological networks for functional discovery."
 Bioinformatics (ISMB 2005), Vol. 21 Suppl. 1 2005, pages 213-221
 
 with the addition of some information theory and z-score ideas to infer edges
 in a graph repesenting correlation between mRNAs.
 
 this code takes as an argument a file that contains a list of files, 
 one file name per line, each one a data set consisting of a header line 
 followed by data, i.e:
 
 probeID        Array1Dye1      Array2Dye1      Array3Dye1      Array4Dye1 ...etc
 1415670_at     9.172966        9.236673        8.789873        8.698874   ...etc
 1415671_at     10.665967       10.677811       10.463303       10.464032  ...etc
    .               .               .               .               .
    .               .               .               .               .
   etc             etc             etc             etc             etc


 Copyright (C) 2015 James Long

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
 
#include "codense.h"
#include <unistd.h>

#define BUFSZ 512
#define VERSION "1.0"

int usage()
{
  printf("\nAn implementation of the CODENSE algortihm:\n");
  printf("Hu H, Yan X, Huang Y, Han J, Zhou XJ: \"Mining coherent dense subgraphs\n");
  printf("across massive biological networks for functional discovery.\"\n");
  printf("Bioinformatics (ISMB 2005), Vol. 21 Suppl. 1 2005, pages 213-221\n\n");
  printf("usage: codense [options] <input text file>\n");
  printf("                -h -help\n");
  printf("                -V (version)\n");
  printf("This code takes an <input text file> argument, containing a list of files,\n");
  printf("one file name (with path) per line, each one a data set consisting of a\n");
  printf("header line followed by mRNA expression data, i.e:\n\n");
  printf("probeID        Array1Dye1      Array2Dye1      Array3Dye1      Array4Dye1 ...etc\n");
  printf("1415670_at     9.172966        9.236673        8.789873        8.698874   ...etc\n");
  printf("1415671_at     10.665967       10.677811       10.463303       10.464032  ...etc\n");
  printf("   .               .               .               .               .\n");
  printf("   .               .               .               .               .\n");
  printf("  etc             etc             etc             etc             etc\n\n");
  printf("See the README file.\n\n");
  
  return 0;
}

int main(int argc, char *argv[])
{
  int i, option, paramsPrinted=0, supp;
  double SUPP;
  gnGraph **coherentDenseSubGraphs;
  
  /* pipeline parameters
  DO_MI          mutual information flag, either 0 or 1 for all runs
  DO_Z           z-score flag, used for PC & MI, 0 or 1 for all runs
  CO_DENSITY     minimum density for a dense subgraph in coherent dense subgraphs
  SG_DENSITY     minimum density for a dense subgraph in summary graph dense subgraphs
  SO_DENSITY     minimum density for a dense subgraph in second order dense subgraphs
  GENE_THRESHOLD gene threshold: average concentration a gene has to exceed to
                 be considered 'on' in Pearson's correlation
  PC_CUTOFF      minimum PC score for an edge to be included in a 1st-order graph:
                 for simple Pearson's, values >= PC are an edge
                 for DOZ, values >= PC and in top z-score slice are an edge
  SOG_CUTOFF     minimum similarity score for an edge to be included in a 
                 2nd-order graph for GenerateSecondOrderGraph()
  SUPPORT        an edge must exist in at least 'support' of the microarray 
                 data sets in order to be included in the summary graph 
                   
  ZPC_SLICE      z-score cutoff for PC
  ZMI_SLICE      z-score cutoff for MI
  
  mutual information
  DIM            joint probability matrix Dimension, mod 32 = 0
  MI_CUTOFF,     score that MI must be >= for an edge to be considered:
                 for DO_Z, values >= MI_CUTOFF and in top z-score slice are an edge
                 for simple mutual information, values >= MI_CUTOFF are an edge
  NEG_CUTOFF,    1) if using Mutual Information and only positive Pearson's 
                    Correlation to infer an edge, it must be above this value, 
                    i.e. no strong negative Pearson's should infer an edge, 
                    even if MI is high (the usual mode) 
                 2) if using negative Pearson's Correlation to infer an edge, it 
                    must be below this value (future experimental feature, only 
                    do positive PC for now)
  PEAK           central max of kernel density function
  SFACT          factor by which the sum of one gene's expression levels must
                 be less than the other for an MI calculation to proceed
  */

  FILE *fout;
  
  /* options parsing */
  while((option = getopt(argc, argv, "hV")) > 0)
  {
    switch(option)
    {
      case 'h':
        usage();
        return 0;
     
      case 'V':
        printf("ver: %s\n", VERSION);
        return 0;
       
      default: break;
    }
  }
  
  if(argc == 2)
  {
    /* start with SUPP high, and ratchet down until modules are detected */
    for(supp=9; supp>0; supp--)
    {
      SUPP = (double)supp * 0.1;
  
      fout = fopen("parameters", "w");
      if(fout==NULL)
      {
        fprintf(stderr, "codense: Error opening parameters file, exiting...\n");
        fflush(NULL);
        return 1;
      }
      
      /* set the pipeline parameters */
      fprintf(fout, "#define DO_MI             1\n");
      fprintf(fout, "#define DO_Z              1\n");
      fprintf(fout, "#define CO_DENSITY        0.7\n");
      fprintf(fout, "#define SG_DENSITY        0.7\n");
      fprintf(fout, "#define SO_DENSITY        0.7\n");
      fprintf(fout, "#define GENE_THRESHOLD    0.005\n");
      fprintf(fout, "#define PC_CUTOFF         0.7\n");
      fprintf(fout, "#define SOG_CUTOFF        0.5\n");
      fprintf(fout, "#define SUPPORT           %3.1lf\n", SUPP);
      fprintf(fout, "#define PLUS_CORREL       1\n");
      fprintf(fout, "#define COMINNODES        4\n");
      fprintf(fout, "#define SGMINNODES        4\n");
      fprintf(fout, "#define SOMINNODES        4\n");
      fprintf(fout, "#define ZPC_SLICE         0.015\n");
      fprintf(fout, "#define ZMI_SLICE         0.015\n");

      /* pipeline mutual information parameters */
      fprintf(fout, "#define DIM               512\n");
      fprintf(fout, "#define MI_CUTOFF         3.0\n");
      fprintf(fout, "#define NEG_CUTOFF        -0.25\n");
      fprintf(fout, "#define PEAK              256.0\n");
      fprintf(fout, "#define SFACT             10.0\n");
      fclose(fout);
      
      if(! paramsPrinted)
      {
        printf("\nparameters file:\n");
        system("cat parameters");
        paramsPrinted = 1;
      }
      
      printf("\ncomputing for SUPPORT = %3.1lf:\n", SUPP);
      
      coherentDenseSubGraphs = CodenseMI(argv[1]);
      if(coherentDenseSubGraphs != NULL)
      {
        i = 0;
        printf("coherentDenseSubGraphs returned:\n\n");
        while(coherentDenseSubGraphs[i])
          PrintGraph(coherentDenseSubGraphs[i++]);
	
        /* free coherentDenseSubGraphs */
        i = 0;
        while(coherentDenseSubGraphs[i])
          FreeGraph(coherentDenseSubGraphs[i++]);
	
        free(coherentDenseSubGraphs);
        printf("\n\n");
      }
      else
        printf("no coherentDenseSubGraphs returned\n\n");
    }
  }
  else usage();
  
  return 0;
}

/*
 CodenseMI.c                                          01/2008 jlong@jimlong.org

 1) read in a parameters file
 2) read in a series of microarray data
 3) use pearson's/mutual-information techniques to generate a graph G_i for 
    each microarray set, i.e. a microarrayGraphSet
 4) create a summaryGraph from the microarrayGraphSet, and create a list of 
    dense subgraphs from the summaryGraph (summaryGraphDenseSubGraphSet)
 5) form an edge support matrix for each entry in summaryGraphDenseSubGraphSet
 6) generate second order graphs S_j for each edge support matrix
 7) mine dense subgraphs sub_k(S_j) of the second order graphs and
 8) convert dense subgraphs into coherent dense subgraphs of the original G_i
 
 The above is an implementation of CODENSE:
 
 Hu H, Yan X, Huang Y, Han J, Zhou XJ: "Mining coherent dense subgraphs across
 massive biological networks for functional discovery."
 Bioinformatics (ISMB 2005), Vol. 21 Suppl. 1 2005, pages 213-221
 
 with the addition of some information theory and z-score ideas to infer edges
 in a graph repesenting correlation between mRNAs.
 
 this code takes as an argument a file that contains a list of files, one file
 name per line, each one a data set; see ReadMicroarrayData.c for format
 
 returns a list of coherent dense subgraphs.
 
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
#include <unistd.h>

#define BUFSZ 512
#define NUMEDGES 210 /* maximum number of edges in a graph for ODES,
                        a clique of 21 nodes (21*20/2 = 210)
                     */

/* debugging */
//#define DELAY
//#define PROGRESS
//#define TRACE

int SubGraphEquality(gnGraph *s01, gnGraph *s02);

gnGraph ** CodenseMI(char * fileList)
{
  int i, j, k, m, n, dim, doMI, doZ, listLen=0, coMinNodes, numEdges,
      sgMinNodes, soMinNodes, numberOfDataSets, savNumDataSets, plusCorrelateOnly;
  char buf[BUFSZ], *p;
  double coDensity, geneThreshold, sgDensity, soDensity, miCutOff, negCutOff, 
         pcCutOff,  sogCutOff, sumfactor, peak, support, zmiSlice, zpcSlice;
  
  edgeSupportMatrix *esMatrix=NULL;                   /*  */
  gnGraph **coherentDenseSubGraphs=NULL,              /*  */
          **coherentDenseSubGraphsToReturn=NULL,      /*  */
          **microarrayGraphSet=NULL,                  /*  */
           *secondOrderGraph=NULL,                    /*  */
          **secondOrderGraphDenseSubGraphSet=NULL,    /*  */
           *summaryGraph=NULL,                        /*  */
          **summaryGraphDenseSubGraphSet=NULL,        /*  */
           *v2e;
  microarrayDataSet *maDataSet;
  FILE *infile, *params;
  
  params = fopen("parameters", "r");
  if(params==NULL)
  {
    fprintf(stderr, "CodenseMI: Error opening parameters file, exiting...\n");
    return NULL;
  }
  
  /* set defaults, these are modified at runtime by reading parameters file */
  dim = doMI = doZ = 0;
  coDensity  = sgDensity  = soDensity = 0.7;
  coMinNodes = sgMinNodes = soMinNodes = 3;
  geneThreshold = 0.0001;
  miCutOff = pcCutOff = sogCutOff = 0.5;
  negCutOff = -0.5;
  peak = 0.0;
  plusCorrelateOnly = 1; /* this one can not be modified by parameters file, may change in the future */
  sumfactor = 10.0;
  support = 0.5;
  zmiSlice = zpcSlice = 0.000001;
  
  while(fgets(buf, BUFSZ, params))
  {
    if(strstr(buf, "#define"))
    {
      p = strtok(buf,  " \t"); /* 1st token */
      p = strtok(NULL, " \t"); /* 2nd token */
      
      if (!strcmp(p, "DO_MI"))
      {
        p = strtok(NULL, " \t");
        doMI = (int) strtol(p, NULL, 10);
      }
      else if (!strcmp(p, "DO_Z"))
      {
        p = strtok(NULL, " \t");
        doZ = (int) strtol(p, NULL, 10);
      }
      else if(!strcmp(p, "CO_DENSITY"))
      {
        p = strtok(NULL, " \t");
        coDensity = strtod(p, NULL);
      }
      else if(!strcmp(p, "SG_DENSITY"))
      {
        p = strtok(NULL, " \t");
        sgDensity = strtod(p, NULL);
      }
      else if(!strcmp(p, "SO_DENSITY"))
      {
        p = strtok(NULL, " \t");
        soDensity = strtod(p, NULL);
      }
      else if(!strcmp(p, "GENE_THRESHOLD"))
      {
        p = strtok(NULL, " \t");
        geneThreshold = strtod(p, NULL);
      }
      else if(!strcmp(p, "NEG_CUTOFF"))
      {
        p = strtok(NULL, " \t");
        negCutOff = strtod(p, NULL);
        
        if(negCutOff<-1.0 || negCutOff>1.0)
        {
          fprintf(stderr, "CodenseMI: negCutOff out of bounds, exiting...\n");
          return NULL;
        }
      }
      else if(!strcmp(p, "PC_CUTOFF"))
      {
        p = strtok(NULL, " \t");
        pcCutOff = strtod(p, NULL);
      }
      else if(!strcmp(p, "SOG_CUTOFF"))
      {
        p = strtok(NULL, " \t");
        sogCutOff = strtod(p, NULL);
      }
      else if(!strcmp(p, "SUPPORT"))
      {
        p = strtok(NULL, " \t");
        support = strtod(p, NULL);
      }
      else if (!strcmp(p, "COMINNODES"))
      {
        p = strtok(NULL, " \t");
        coMinNodes = (int) strtol(p, NULL, 10);
      }
      else if (!strcmp(p, "SGMINNODES"))
      {
        p = strtok(NULL, " \t");
        sgMinNodes = (int) strtol(p, NULL, 10);
      }
      else if (!strcmp(p, "SOMINNODES"))
      {
        p = strtok(NULL, " \t");
        soMinNodes = (int) strtol(p, NULL, 10);
      }
      else if (!strcmp(p, "DIM"))
      {
        p = strtok(NULL, " \t");
        dim = (int) strtol(p, NULL, 10);
      }
      if(!strcmp(p, "MI_CUTOFF"))
      {
        p = strtok(NULL, " \t");
        miCutOff = strtod(p, NULL);
      }
      else if (!strcmp(p, "PEAK"))
      {
        p = strtok(NULL, " \t");
        peak = strtod(p, NULL);
      }
      else if (!strcmp(p, "SFACT"))
      {
        p = strtok(NULL, " \t");
        sumfactor = strtod(p, NULL);
      }
      else if(!strcmp(p, "ZPC_SLICE"))
      {
        p = strtok(NULL, " \t");
        zpcSlice = strtod(p, NULL);
      }
      else if(!strcmp(p, "ZMI_SLICE"))
      {
        p = strtok(NULL, " \t");
        zmiSlice = strtod(p, NULL);
      }
    }
  }
  
#ifdef TRACE
  /* print parameters */
  printf("#define DO_MI             %d\n",  doMI);
  printf("#define DO_Z              %d\n",  doZ);
  printf("#define CO_DENSITY        %lf\n", coDensity);
  printf("#define SG_DENSITY        %lf\n", sgDensity);
  printf("#define SO_DENSITY        %lf\n", soDensity);
  printf("#define GENE_THRESHOLD    %lf\n", geneThreshold);
  printf("#define NEG_CUTOFF        %lf\n", negCutOff);
  printf("#define PC_CUTOFF         %lf\n", pcCutOff);
  printf("#define SOG_CUTOFF        %lf\n", sogCutOff);
  printf("#define SUPPORT           %lf\n", support);
  printf("#define PLUS_CORREL       1\n");
  printf("#define COMINNODES        %d\n",  coMinNodes);
  printf("#define SGMINNODES        %d\n",  sgMinNodes);
  printf("#define SOMINNODES        %d\n",  soMinNodes);
  
  if(doMI)
  {
    printf("#define DIM               %d\n",  dim);
    printf("#define MI_CUTOFF         %lf\n", miCutOff);
    printf("#define PEAK              %lf\n", peak);
    printf("#define SFACT             %lf\n", sumfactor);
  }
  
  if(doZ)
  {
    printf(         "#define ZPC_SLICE         %10.8lf\n", zpcSlice);
    if(doMI) printf("#define ZMI_SLICE         %10.8lf\n", zmiSlice);
  }
#endif
  
  fclose(params);
  
  infile = fopen(fileList, "r");
  if(infile==NULL)
  {
    fprintf(stderr, "CodenseMI: Error opening %s, exiting...\n", fileList);
    return NULL;
  }
  
  /* count number (numberOfDataSets) of experiments in microarrayGraphSet */
  numberOfDataSets = 0;
  while(fgets(buf, BUFSZ, infile))
  {
    if(buf[0] == '\n') continue;
    numberOfDataSets++;          /* upper bound */
  }
  rewind(infile);
  
#ifdef TRACE
  printf("Number of Data Sets = %d\n\n", numberOfDataSets);
#endif
  
#ifdef DELAY
  fflush(NULL);
  sleep(3);
#endif

  microarrayGraphSet = (gnGraph **) malloc(numberOfDataSets*sizeof(gnGraph *));
  if(microarrayGraphSet==NULL)
  {
    fprintf(stderr, "CodenseMI: malloc error, exiting...\n");
    return NULL;
  }
  
  /* init */
  for(i=0; i<numberOfDataSets; i++)
    microarrayGraphSet[i] = NULL;
    
  savNumDataSets = numberOfDataSets;

#ifdef PROGRESS
  printf("Creating Microarray Data Set Graph\n"); 
  system("date");
  printf("\n");
  fflush(NULL);
#endif
  
  i = 0;
  while(fgets(buf, BUFSZ, infile))
  {
    /* scrub */
    for(j=0; j<strlen(buf); j++)
      if(buf[j]=='\n' || buf[j]=='\r')
        buf[j] = '\0';

    if(buf[0])
    {
      maDataSet = ReadMicroarrayData(buf);                                 /* 1 */
      if(maDataSet)
      {
        microarrayGraphSet[i++] = CreateMicroarrayDataSetGraph(maDataSet,  /* 2 */
                                                               zpcSlice,      pcCutOff,
                                                               zmiSlice,      miCutOff,
                                                               negCutOff,     peak,
                                                               geneThreshold, sumfactor,
                                                               dim,  doMI,  0,  doZ,
                                                               plusCorrelateOnly);
#ifdef TRACE
        printf("microarrayGraphSet[%d]: \n", i-1);
        PrintGraph(microarrayGraphSet[i-1]);
#ifdef DELAY 
        sleep(3);
#endif
#endif
        DeallocateMicroarrayDataSet(maDataSet);
      }
      else
        numberOfDataSets--;
    }
    else
      numberOfDataSets--;
  }

  fclose(infile);
  
  if(numberOfDataSets==0)
  {
    fprintf(stderr, "CodenseMI: no datasets to work on, returning...\n\n");
    return NULL;
  }
  
#ifdef PROGRESS
  printf("Creating Summary Graph\n");
  system("date");
  printf("\n");
  fflush(NULL);
#endif

  summaryGraph = CreateSummaryGraph(microarrayGraphSet,                    /* 3 */
                                    numberOfDataSets,
                                    support);
#ifdef TRACE
  printf("summary graph: \n");
  PrintGraph(summaryGraph);
#endif
  numEdges = NumGraphEdges(summaryGraph);
  if(numEdges > NUMEDGES)
  {
#ifdef TRACE
    printf("\nsummaryGraph has %d edges, max for ODES is %d, returning...\n", numEdges, NUMEDGES);
#endif
    FreeGraph(summaryGraph);
    return NULL;
  }
#ifdef DELAY
  sleep(3);
#endif
  
#ifdef PROGRESS
  printf("Creating Summary Graph Dense SubGraph Set\n"); 
  system("date");
  printf("\n");
  fflush(NULL);
#endif

  summaryGraphDenseSubGraphSet = ODES(summaryGraph, sgDensity, sgMinNodes);/* 3 */
  if(summaryGraphDenseSubGraphSet==NULL || summaryGraphDenseSubGraphSet[0]==NULL)
  {
#ifdef TRACE
    printf("CodenseMI: no dense subgraphs found in summaryGraph, returning...\n\n");
#endif
    FreeGraph(summaryGraph);
    return NULL;
  }
  
  FreeGraph(summaryGraph);
  
  i = 0;
  while(summaryGraphDenseSubGraphSet[i])
  {
#ifdef TRACE
  printf("summaryGraphDenseSubGraphSet[%d]: \n", i);
  PrintGraph(summaryGraphDenseSubGraphSet[i]);
#ifdef DELAY
  sleep(3);
#endif
#endif

#ifdef PROGRESS
    printf("Generating Edge Support Matrix\n"); 
    system("date");
    printf("\n");
    fflush(NULL);
#endif
    /* 4 & 5 */
    esMatrix = GenerateEdgeSupportMatrix(summaryGraphDenseSubGraphSet[i],
                                         microarrayGraphSet,
                                         numberOfDataSets);
#ifdef PROGRESS
    printf("Generating Second Order Graph\n"); 
    system("date");
    printf("\n");fflush(NULL);
#endif

    /* only esMatrix, sogCutOff, and hard coded values matter, other variables
     * not used here since doMI = doZ = 0, & doS = plusCorrelateOnly = 1, see 
     * source code
     */
    secondOrderGraph = GenerateSecondOrderGraph(esMatrix,  zpcSlice,
                                                sogCutOff, zmiSlice, 
                                                miCutOff,  negCutOff,
                                                peak, 0.0, 0.0, dim, 0, 1, 0, 1);

#ifdef TRACE
    printf("Edge Support Matrix:\n\n");
    k = 0;
    while(esMatrix->esvList[k] != NULL)
    {
      printf("%s  \t", esMatrix->esvList[k]->edgeName);
  
      for(j=0; j<esMatrix->vecLen; j++)
        if(esMatrix->esvList[k]->supportVector[j] > 0.0)
          printf("%2.1lf  ", esMatrix->esvList[k]->supportVector[j]);
        else
          printf("---  ");
  
      printf("\n");
    
      k++;
    }
    
    if(esMatrix != NULL) 
      FreeEdgeSupportMatrix(esMatrix);
      
    printf("\nsecondOrderGraph:\n"); 
    PrintGraph(secondOrderGraph);
#endif
    numEdges = NumGraphEdges(secondOrderGraph);
    if(numEdges > NUMEDGES)
    {
#ifdef TRACE
      printf("\nsecondOrderGraph has %d edges, max for ODES is %d, moving on to next summaryGraphDenseSubGraphSet...\n", numEdges, NUMEDGES);
#endif
      FreeGraph(secondOrderGraph);
      FreeGraph(summaryGraphDenseSubGraphSet[i]);
      i++;
      continue;
    }
#ifdef DELAY
    sleep(3);
#endif

#ifdef PROGRESS
    printf("Generating Second Order Graph Dense SubGraph Set\n"); 
    system("date");
    printf("\n");
    fflush(NULL);
#endif
    /* 6 */
    secondOrderGraphDenseSubGraphSet = ODES(secondOrderGraph, soDensity, soMinNodes);
    
    FreeGraph(secondOrderGraph);
      
    j = 0;
    while(secondOrderGraphDenseSubGraphSet[j])
    {
#ifdef TRACE
      printf("\nsecondOrderGraphDenseSubGraphSet[%d] = \n", j);
      PrintGraph(secondOrderGraphDenseSubGraphSet[j]);
#ifdef DELAY
      sleep(3);
#endif
#endif

#ifdef PROGRESS
      printf("Generating Coherent Dense Sub Graphs\n"); 
      system("date");
      printf("\n");
      fflush(NULL);
#endif
      /* 7 */
      v2e = Vertex2Edge(secondOrderGraphDenseSubGraphSet[j]);
      coherentDenseSubGraphs = ODES(v2e, coDensity, coMinNodes);
      
      FreeGraph(v2e);

      /* what is the length of the list */
      k = 0;
      while(coherentDenseSubGraphs[k]) k++;

#ifdef PROGRESS
      printf("There was/were %d coherentDenseSubGraphs identified, checking for uniqueness...\n", k);
      fflush(NULL);
#endif

      if(k)
      {
        if(coherentDenseSubGraphsToReturn != NULL)
        {
          /* check for duplicates against coherentDenseSubGraphsToReturn */
          m = 0;
          while(coherentDenseSubGraphsToReturn[m] != NULL)
          {
            for(n=0; n<k; n++)
              if(SubGraphEquality(coherentDenseSubGraphs[n], coherentDenseSubGraphsToReturn[m]))
                break;
          
            if(n<k)
            {
              FreeGraph(coherentDenseSubGraphs[n]);
              for(;n<k; n++) /* value at k is NULL */
                coherentDenseSubGraphs[n] = coherentDenseSubGraphs[n+1];
          
              k--;
            }
            
            m++;
          }
        }
      }
      
#ifdef PROGRESS
      printf("There are %d new unique coherentDenseSubGraphs\n", k);
      fflush(NULL);
#endif

      if(k)
      {
        /* reallocate coherentDenseSubGraphsToReturn */
        coherentDenseSubGraphsToReturn = (gnGraph **) 
                                         realloc(coherentDenseSubGraphsToReturn, 
                                                (listLen+k+1)*sizeof(gnGraph *));
        if(coherentDenseSubGraphsToReturn==NULL)
        {
          fprintf(stderr, "CodenseMI: realloc error, returning...\n");
          return NULL;
        }
      
        for(m=0; m<k; m++)
          coherentDenseSubGraphsToReturn[listLen+m] = coherentDenseSubGraphs[m];
      
        coherentDenseSubGraphsToReturn[listLen+k] = NULL;
        listLen += k;
      }
      
      FreeGraph(secondOrderGraphDenseSubGraphSet[j]);                    /**recent**/
      j++;
    }
    
    FreeGraph(summaryGraphDenseSubGraphSet[i]);                          /**recent**/
    
    i++;
  } /* end while(summaryGraphDenseSubGraphSet[i]) */
  
  for(i=0; i<savNumDataSets; i++)
    FreeGraph(microarrayGraphSet[i]);
  
  free(microarrayGraphSet);
  free(secondOrderGraphDenseSubGraphSet);                                /**recent**/
  free(summaryGraphDenseSubGraphSet);                                    /**recent**/
    
#ifdef PROGRESS
  if(coherentDenseSubGraphsToReturn != NULL)
  {
    printf("Returned coherentDenseSubGraphs:\n");
    i = 0;
    while(coherentDenseSubGraphsToReturn[i])
      PrintGraph(coherentDenseSubGraphsToReturn[i++]);
  }
#endif

  return coherentDenseSubGraphsToReturn;
}


/* is subgraph S1 = subgraph S2? */
int SubGraphEquality(gnGraph *s01, gnGraph *s02)
{
  int i, j, nodeIsInGraph;

  /* s01 != s02 if it does not contain same # nodes as s02 */
  if(s01->numNodes != s02->numNodes)
    return 0;

  /*
   s01 != s02 if it contains a node not in s02,
   no need to worry about edges for our problem
  */
  for(i=0; i<s01->numNodes; i++)
  {
    nodeIsInGraph = 0;
    for(j=0; j<s02->numNodes; j++)
      if(!strcasecmp(s01->gnList[i]->label, s02->gnList[j]->label))
        {nodeIsInGraph = 1; break;}

    if(!nodeIsInGraph)
      return 0; /* is not equal */
  }

  return 1;     /* is equal */
}

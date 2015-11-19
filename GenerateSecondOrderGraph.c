/*
 GenerateSecondOrderGraph.c                      12/2009 jlong@jimlong.org
 
 given an edgeSupportMatrix, construct a graph where edges in the matrix are 
 vertices, and an edge exists between two vertices if they are correlated in 
 the matrix (but see definition in Zhou paper below). 
 We did this already in CreateMicroarrayDataSetGraph(), so can just turn the 
 edgeSupportMatrix into a microarrayDataSet, and call this routine,
 GenerateSecondOrderGraph(), with doMI = doZ = 0, doS = plusCorrelateOnly = 1, 
 which will call CreateMicroarrayDataSetGraph() with geneThreshold = 0.0, where
 Similarity (doS) is used in place of PC to determine an edge, and pcCutOff is 
 the fraction to meet or exceed in order for an edge to be declared. 
 
 Definition in Zhou paper (paraphrased):
 Given a relation graph dataset D={G1,G2,...,Gn}, the Second Order Graph (SOG) 
 is an unweighted graph, where the vertex set of SOG is the edge set of G, and 
 an edge connects two vertices u and v in SOG if the "similarity" between the 
 Edge Support Vectors of u and v is greater than a threshold.
 
 
 Note that the edges are unweighted, but using CreateMicroarrayDataSetGraph() 
 gives them a weight, so before returning the graph, we change the edge weights
 to 1.0 in case the weights are ever used for anything.
 
 Copyright (C) 2009 James Long
 
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

gnGraph * GenerateSecondOrderGraph(const edgeSupportMatrix *esMatrix,
                                   double zpcSlice,      double pcCutOff,
                                   double zmiSlice,      double miCutOff,
                                   double negCutOff,     double peak,
                                   double geneThreshold, double sumfactor,
                                   int dim, int doMI, int doS, int doZ,
                                   int plusCorrelateOnly)
{
  int i, j, k;
  gnGraph *g;
  microarrayDataSet *maDataSet;

  if(esMatrix == NULL)
    return NULL;
  
  maDataSet = (microarrayDataSet *) malloc(sizeof(microarrayDataSet));
  if(maDataSet==NULL)
  {
    fprintf(stderr, "GenerateSecondOrderGraph: Error allocating maDataSet, exiting...\n");
    return NULL;
  }
  
  maDataSet->numCols = esMatrix->vecLen;
  
  i = 0;
  while(esMatrix->esvList[i]!=NULL) i++;
  
  maDataSet->numRows = i;
  maDataSet->data   = (double *) malloc(i*maDataSet->numCols*sizeof(double));
  maDataSet->labels = (char **)  malloc(i*sizeof(char *));
  if(maDataSet->labels==NULL)
  {
    fprintf(stderr, "GenerateSecondOrderGraph: Error allocating maDataSet->labels, exiting...\n");
    return NULL;
  }
  
  k = 0;
  for(i=0; i<maDataSet->numRows; i++)
  {
    maDataSet->labels[i] = (char *) malloc(strlen((esMatrix->esvList[i])->edgeName)+1);
    if(maDataSet->labels[i]==NULL)
    {
      fprintf(stderr, "GenerateSecondOrderGraph: Error allocating maDataSet->labels[%d], exiting...\n", i);
      return NULL;
    }
    
    strcpy(maDataSet->labels[i], (esMatrix->esvList[i])->edgeName);
    
    for(j=0; j<maDataSet->numCols; j++)
    {
      maDataSet->data[k++] = (esMatrix->esvList[i])->supportVector[j];
    }
  }
  
  g = CreateMicroarrayDataSetGraph(maDataSet, zpcSlice, pcCutOff, 
                                   zmiSlice,  miCutOff, negCutOff, 
                                   peak, 0.0, 0.0, dim, doMI, doS, 
                                   doZ, plusCorrelateOnly);

  /* change edge weights to 1.0 */
  for(i=0; i<g->numNodes; i++)
    for(j=0; j<g->gnList[i]->numEdges; j++)
      g->gnList[i]->edgeWeights[j] = 1.0;
      
  DeallocateMicroarrayDataSet(maDataSet);
  
  return g;
}


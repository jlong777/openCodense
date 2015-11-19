/*
  GenerateEdgeSupportMatrix.c                        12/2009 jlong@jimlong.org
  
  given a subGraph and a graphSet, generate an edgeSupportMatrix, i.e.
  a list consisting of every edge in subGraph along with a vector for 
  each edge enumerating its weight in each graphSet
 
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

edgeSupportMatrix * GenerateEdgeSupportMatrix(gnGraph *subGraph,
                                               gnGraph **graphSet,
                                               int numberOfDataSets)
{
  int i, j, k, m, len, numEdgeSupportVectors=0;
  char *head, *tail;
  void *vp;
  graphNode *p, *q;
  edgeSupportMatrix *esMatrix;
  edgeSupportVector **esvList;
  
  esMatrix = (edgeSupportMatrix *)  malloc(sizeof(edgeSupportMatrix));
  esvList  = (edgeSupportVector **) malloc(sizeof(edgeSupportVector *));
  if(esMatrix==NULL || esvList==NULL)
  {
    fprintf(stderr, "GenerateEdgeSupportMatrix: malloc error, returning...\n");
    return NULL;
  }
  esvList[numEdgeSupportVectors] = NULL;
  
  /* for fast searching, sort each graph in graphSet on labels */
  for(i=0; i<numberOfDataSets; i++)
    qsort((void *) graphSet[i]->gnList, graphSet[i]->numNodes, sizeof(graphNode *), compGN);
  
  /*
   Walk through the edges of subGraph. Do this by visiting each node, and 
   looking only at edges to nodes lexicographically greater than the node 
   we are at. Then build the edgeSupportVector.
  */
  numEdgeSupportVectors = 0;
  for(i=0; i<subGraph->numNodes; i++)
  {
    p    = subGraph->gnList[i];
    head = p->label;
    
    for(j=0; j<p->numEdges; j++)
    {
      tail = p->edges[j];
      if(strcasecmp(head, tail)<0)
      {
        numEdgeSupportVectors++;
	esvList = (edgeSupportVector **) realloc(esvList, (numEdgeSupportVectors+1)*
	                                         sizeof(edgeSupportVector *));
	if(esvList==NULL)
        {
          fprintf(stderr, "GenerateEdgeSupportMatrix: realloc error, returning...\n");
          return NULL;
        }
	esvList[numEdgeSupportVectors] = NULL;
	
	esvList[numEdgeSupportVectors-1] = (edgeSupportVector *)
	                                   malloc(sizeof(edgeSupportVector));
	if(esvList[numEdgeSupportVectors-1]==NULL)
        {
	  fprintf(stderr, "GenerateEdgeSupportMatrix: malloc error, returning...\n");
          return NULL;
	}
        
        /* create an edgeName consisting of head<-->tail */
	len = strlen(head) + strlen(tail) + 5;
	esvList[numEdgeSupportVectors-1]->edgeName = (char *) malloc(len);
        if(esvList[numEdgeSupportVectors-1]->edgeName == NULL)
        {
          fprintf(stderr, "GenerateEdgeSupportMatrix: malloc error, returning...\n");
          return NULL;
        }
        
	strcpy(esvList[numEdgeSupportVectors-1]->edgeName, head);
	strcat(esvList[numEdgeSupportVectors-1]->edgeName, "<-->");
	strcat(esvList[numEdgeSupportVectors-1]->edgeName, tail);

	if(numberOfDataSets)
	{				   
          esvList[numEdgeSupportVectors-1]->supportVector = (double *) 
                                                            malloc(numberOfDataSets*
                                                                   sizeof(double));
          if(esvList[numEdgeSupportVectors-1]->supportVector==NULL)
          {
	    fprintf(stderr, "GenerateEdgeSupportMatrix: malloc error, returning...\n");
            return NULL;
	  }
	}

        /*
	 now compute the support for this edge in the graphSet
	 i.e, generate the supportVector
	*/
	for(k=0; k<numberOfDataSets; k++)
	{
          /* find the graphNode whose label is head */
          vp = (void *)(graphSet[k]->gnList);
          vp = bsearch((void *)head, vp, graphSet[k]->numNodes,
                        sizeof(graphNode *), compStrGN);
          q = *(graphNode **)vp;
	  
	  /* now look to see if q has an edge to tail */
	  esvList[numEdgeSupportVectors-1]->supportVector[k] = 0.0;
	  for(m=0; m<q->numEdges; m++)
	    if(!strcasecmp(tail, q->edges[m]))
            {
	      esvList[numEdgeSupportVectors-1]->supportVector[k] = q->edgeWeights[m];
	      break;
            }
	}
      }
    }
  }
  
  esMatrix->esvList = esvList;
  esMatrix->vecLen  = numberOfDataSets;
  
  return esMatrix;
}

void FreeEdgeSupportMatrix(edgeSupportMatrix *esMatrix)
{
  int i=0;
  
  while(esMatrix->esvList[i] != NULL)
  {
    if(esMatrix->esvList[i]->edgeName != NULL)
      free(esMatrix->esvList[i]->edgeName);
    
    if(esMatrix->esvList[i]->supportVector != NULL)
      free(esMatrix->esvList[i]->supportVector);
      
    free(esMatrix->esvList[i]);
  
    i++;
  }
  
  free(esMatrix->esvList);
  free(esMatrix);
}

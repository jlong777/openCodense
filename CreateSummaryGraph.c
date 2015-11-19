/*
 CreateSummaryGraph.c                            02/2008 jlong@jimlong.org
 
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

//#define DEBUG_CSG

/* returns a graph sorted on labels */
#define BUFLEN 1024
gnGraph * CreateSummaryGraph(gnGraph **graphSet, 
                             int numberOfDataSets, 
                             double support)
{
  long i, j, k, edgeBufLen=0, numEdgeBufEdges, numLabels, numUniq, sum;
  char **edgeBuf, **uniqBuf, **uniqLabels;
  void *hit;
  double *edgeWghtBuf, sumOfWeights, threshold;
  gnGraph *sg;
  graphNode *tmp;
  
  threshold = support * numberOfDataSets;

#ifdef DEBUG_CSG
  printf("threshold = %lf = support = %lf * numberOfDataSets = %d\n",
          threshold, support, numberOfDataSets);
#endif

  /* temporary buffers to hold summary graph edges and weights */
  edgeBuf     = (char **)  malloc(BUFLEN*sizeof(char *));
  edgeWghtBuf = (double *) malloc(BUFLEN*sizeof(double));
  if(edgeBuf==NULL || edgeWghtBuf==NULL)
  {
    fprintf(stderr, "CreateSummaryGraph: malloc error, returning...\n");
    return NULL;
  }
  edgeBufLen = BUFLEN;

  /* count all labels */
  sum = 0;
  for(i=0; i<numberOfDataSets; i++)
    sum += graphSet[i]->numNodes;
    
  uniqLabels = (char **) malloc(sum*sizeof(char *));
  if(uniqLabels==NULL)
  {
    fprintf(stderr, "CreateSummaryGraph: malloc error, returning...\n");
    return NULL;
  }
  
  /* get all the labels from the datasets */
  k = 0;
  for(i=0; i<numberOfDataSets; i++)
    for(j=0; j<graphSet[i]->numNodes; j++)
      if(graphSet[i]->gnList)
        uniqLabels[k++] = graphSet[i]->gnList[j]->label;
  numLabels = k;
  
  /* uniq the list of all vertex labels */
  qsort((void *) uniqLabels, numLabels, sizeof(char *), compString);
  j = 0;
  for(i=j+1; i<numLabels; i++)
  {
    while(i<numLabels && !strcasecmp(uniqLabels[j], uniqLabels[i])) i++;
    if(i<numLabels) uniqLabels[++j] = uniqLabels[i];
  }
  numLabels = j+1;
  
  /* sort each gnList on labels so it can be searched with bsearch */
  for(i=0; i<numberOfDataSets; i++)
  {
#ifdef DEBUG_CSG
    printf("unsorted graph\n==============\n"); PrintGraph(graphSet[i]);
#endif
    qsort((void *) graphSet[i]->gnList,
                   graphSet[i]->numNodes,
                   sizeof(graphNode *), compGN);
#ifdef DEBUG_CSG
    printf("sorted graph\n============\n"); PrintGraph(graphSet[i]);
#endif
  }

  /* the summary graph */
  sg = (gnGraph *) malloc(sizeof(gnGraph));
  if(sg==NULL)
  {
    fprintf(stderr, "CreateSummaryGraph: malloc error, returning...\n");
    return NULL;
  }
  
  /* init graph with vertices */
  sg->numNodes  = numLabels;
  sg->gnList = (graphNode **) malloc(numLabels * sizeof(graphNode *));
  if((sg->gnList)==NULL)
  {
    fprintf(stderr, "CreateSummaryGraph: malloc error, returning...\n");
    return NULL;
  }
  
  for(i=0; i<numLabels; i++)
  {
    sg->gnList[i] = (graphNode *) malloc(sizeof(graphNode));
    if((sg->gnList[i])==NULL)
    {
      fprintf(stderr, "CreateSummaryGraph: malloc error, returning...\n");
      return NULL;
    }
    
    sg->gnList[i]->label = (char *) malloc(strlen(uniqLabels[i])+1);
    if((sg->gnList[i]->label)==NULL)
    {
      fprintf(stderr, "CreateSummaryGraph: malloc error, returning...\n");
      return NULL;
    }

    strcpy(sg->gnList[i]->label, uniqLabels[i]);
    sg->gnList[i]->numEdges = 0;
    
    sg->gnList[i]->edges       = NULL;
    sg->gnList[i]->edgeWeights = NULL;
  }

  /*
   We now have each graph of graphSet and the summary graph sorted on labels
   (since uniqLabels is sorted, summary graph is sorted).
   Add edges to the summary graph only if they occur enough times across
   the graphSet, i.e. the sum of the edge weights for an edge exceeds 
   threshold = support * numberOfDataSets
   So for every uniqLabel:
   1) look for that label in every graph of graphset
   2) collect all edges attached to it in every graph
   3) sum the edge weights for each edge.
  */
  for(i=0; i<numLabels; i++)
  {
    edgeBuf[0] = NULL;
    numEdgeBufEdges = 0;
    for(j=0; j<numberOfDataSets; j++)
    {
      hit = (void *)(graphSet[j]->gnList);
      hit = bsearch((void *)uniqLabels[i], hit, graphSet[j]->numNodes,    /* 1 */
                     sizeof(graphNode *), compStrGN);
      
      if(hit) /* collect edges attached to this vertex */
      {
        tmp = *(graphNode **)hit;
        for(k=0; k<tmp->numEdges; k++)
        {
          if(edgeBufLen < numEdgeBufEdges + tmp->numEdges)
          {
            edgeBuf     = (char **)  realloc(edgeBuf, 
                                            (numEdgeBufEdges+tmp->numEdges+BUFLEN) *
                                             sizeof(char *));
            edgeWghtBuf = (double *) realloc(edgeWghtBuf, 
                                            (numEdgeBufEdges+tmp->numEdges+BUFLEN) *
                                             sizeof(double));
            if(edgeBuf==NULL || edgeWghtBuf==NULL)
            {
              fprintf(stderr, "CreateSummaryGraph: malloc error, returning...\n");
              return NULL;
            }
            edgeBufLen = numEdgeBufEdges + tmp->numEdges + BUFLEN;
          }
          
          edgeBuf[numEdgeBufEdges]       = tmp->edges[k];                 /* 2 */
          edgeWghtBuf[numEdgeBufEdges++] = tmp->edgeWeights[k];
        }
      }
    }
    
    /*
     edgeBuf contains pointers to every vertex in graphSet that is connected to 
     the same label as the current vertex in the summary graph. Add those edges
     to the summary graph that exceed 'threshold', using as their weight the 
     average.
    */
    
    if(edgeBuf[0] != NULL) /* edges exist */
    { 
      /* uniq the list of vertices in edgeBuf */
      uniqBuf = (char **)  malloc(numEdgeBufEdges*sizeof(char *));
      if(uniqBuf==NULL)
      {
        fprintf(stderr, "CreateSummaryGraph: malloc error, returning...\n");
        return NULL;
      }
      
      for(k=0; k<numEdgeBufEdges; k++)
        uniqBuf[k] = edgeBuf[k];
        
      qsort((void *) uniqBuf, numEdgeBufEdges, sizeof(char *), compString);
      
      j = 0;
      for(k=1; k<numEdgeBufEdges; k++)
      {
        while(k<numEdgeBufEdges && !strcasecmp(uniqBuf[j], uniqBuf[k])) k++;
        if(k<numEdgeBufEdges) uniqBuf[++j] = uniqBuf[k];
      }
      numUniq = j+1;
      
      /* count support for each edge */
      for(j=0; j<numUniq; j++)
      {
        sumOfWeights = 0.0;
        
        for(k=0; k<numEdgeBufEdges; k++)
        {
          if(!strcasecmp(uniqBuf[j], edgeBuf[k]))
            sumOfWeights += edgeWghtBuf[k];
        }
        
        /* add edge to vertex uniqBuf[j] from vertex uniqLabels[i] */
        if(sumOfWeights >= threshold)
        {
#ifdef DEBUG_CSG
          printf("edge to %s will be added to sg\n", uniqBuf[j]);
#endif
          /* go to node uniqLabels[i] in sg */
          hit = (void *)(sg->gnList);
          hit = bsearch((void *)uniqLabels[i], hit, sg->numNodes,
                         sizeof(graphNode *), compStrGN);
          tmp = *(graphNode **)hit;
          
          /* allocate space */
          if(tmp->numEdges == 0)
          {
            tmp->edges       = (char **)  malloc(sizeof(char *));
            tmp->edgeWeights = (double *) malloc(sizeof(double));
            if((tmp->edges)==NULL || (tmp->edgeWeights)==NULL)
            {
              fprintf(stderr, "CreateSummaryGraph: malloc error, returning...\n");
              return NULL;
            }
          }
          else
          {
            tmp->edges       = (char **)  realloc(tmp->edges, 
                                                 (tmp->numEdges + 1) *
                                                  sizeof(char *));
            tmp->edgeWeights = (double *) realloc(tmp->edgeWeights, 
                                                 (tmp->numEdges + 1) *
                                                  sizeof(double));
            if((tmp->edges)==NULL || (tmp->edgeWeights)==NULL)
            {
              fprintf(stderr, "CreateSummaryGraph: malloc error, returning...\n");
              return NULL;
            }
          }
#ifdef DEBUG_CSG
          printf("adding %s to %s\n", uniqBuf[j], tmp->label);
#endif
          /* finally add the edge */
          tmp->edges[tmp->numEdges] = (char *) malloc(strlen(uniqBuf[j])+1);
          if((tmp->edges[tmp->numEdges])==NULL)
          {
            fprintf(stderr, "CreateSummaryGraph: malloc error, returning...\n");
            return NULL;
          }
          strcpy(tmp->edges[tmp->numEdges], uniqBuf[j]);
          
          tmp->edgeWeights[tmp->numEdges] = sumOfWeights/numberOfDataSets;
#ifdef DEBUG_CSG
          printf("added %s with weight %lf\n", tmp->edges[tmp->numEdges],
                                               tmp->edgeWeights[tmp->numEdges]);
#endif
          tmp->numEdges++;
        }
      }
      
      free(uniqBuf);
    }
  }
  
#ifdef DEBUG_CSG
  PrintGraph(sg);
#endif

  free(edgeBuf);
  free(edgeWghtBuf);
  free(uniqLabels);

  return sg;
}

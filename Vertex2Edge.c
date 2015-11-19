/*
  Vertex2Edge.c                                12/2009 jlong@jimlong.org
  
  Each dense subgraph DS of a Second Order Graph is a graph whose vertices
  are labeled "L1<-->L2", where L1 and L2 are vertices in a graph G that
  are adjacent. We will call G the Vertex-to-Edge graph of DS.
  For CODENSE, every coherent dense subgraph is a dense subgraph in some
  Vertex-to-Edge graph G of some dense subgraph D in a Second Order Graph
  of the Summary Graph.
  
  Definition in Zhou paper:
  Given a relation graph dataset D={G1,G2,...,Gn}, a subgraph SG is coherent
  if all the edges of SG have support higher than k and the Second Order Graph
  of SG is dense.
 
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

//#define DEBUG

gnGraph * Vertex2Edge(gnGraph *ds)
{
  int i, j, found_L1, found_L2, len, numEdges, numNodes;
  char *L1, *L2;
  gnGraph *g;
  
  g = (gnGraph *) malloc(sizeof(gnGraph));
  if(g==NULL)
  {
    fprintf(stderr, "Vertex2Edge: malloc error, returning...\n");
    return NULL;
  }
  g->gnList   = NULL;
  g->numNodes = 0;
  
  
  /* get the vertices */
  for(i=0; i<ds->numNodes; i++)
  {
    /* get L1 and L2 */
    len = 0;
    while(ds->gnList[i]->label[len] != '<') len++;
    
    L1 = (char *) malloc(len+1);
    L2 = (char *) malloc(strlen(ds->gnList[i]->label)-len-3);
    if(L1==NULL || L2==NULL)
    {
      fprintf(stderr, "Vertex2Edge: malloc error, returning...\n");
      return NULL;
    }
    
    strncpy(L1, ds->gnList[i]->label, len);
    L1[len] = '\0';
    strcpy (L2, ds->gnList[i]->label + (len+4));
    
#ifdef DEBUG
  printf("L1 = %s  L2 = %s\n", L1, L2);
#endif

    /* now march through the g nodes to see if we have these two */
    found_L1 = found_L2 = 0;
    for(j=0; j<g->numNodes; j++)
    {
      if(!strcmp(g->gnList[j]->label, L1)) found_L1 = 1;
      if(!strcmp(g->gnList[j]->label, L2)) found_L2 = 1;
      if(found_L1 && found_L2) break;
    }

    if(!found_L1) /* not there, so create */
    {
#ifdef DEBUG
  printf("making node %s\n", L1);
#endif
      numNodes = ++(g->numNodes);
      
      g->gnList = (graphNode **) realloc(g->gnList, numNodes*sizeof(graphNode *));
      if(g->gnList==NULL)
      {
        fprintf(stderr, "Vertex2Edge: realloc error, returning...\n");
        return NULL;
      }
      
      g->gnList[numNodes-1] = (graphNode *) malloc(sizeof(graphNode));
      if(g->gnList[numNodes-1]==NULL)
      {
        fprintf(stderr, "Vertex2Edge: malloc error, returning...\n");
        return NULL;
      }
      
      g->gnList[numNodes-1]->label = (char *) malloc(len+1);
      if(g->gnList[numNodes-1]->label==NULL)
      {
        fprintf(stderr, "Vertex2Edge: malloc error, returning...\n");
        return NULL;
      }
      
      strcpy(g->gnList[numNodes-1]->label, L1);
      
      g->gnList[numNodes-1]->numEdges    = 0;
      g->gnList[numNodes-1]->edges       = NULL;
      g->gnList[numNodes-1]->edgeWeights = NULL;
    }
    
    /* make an edge from L1 to L2 */
    for(j=0; j<g->numNodes; j++)
      if(!strcmp(g->gnList[j]->label, L1))
        break;
    
    numEdges = ++(g->gnList[j]->numEdges);
      
    g->gnList[j]->edges       = (char **)  realloc(g->gnList[j]->edges,
                                                   numEdges*sizeof(char *));
    g->gnList[j]->edgeWeights = (double *) realloc(g->gnList[j]->edgeWeights,
                                                   numEdges*sizeof(double));
    if(g->gnList[j]->edges==NULL || g->gnList[j]->edgeWeights==NULL)
    {
      fprintf(stderr, "Vertex2Edge: realloc error, returning...\n");
      return NULL;
    }
      
    g->gnList[j]->edges[numEdges-1] = (char *) malloc(strlen(L2)+1);
    if(g->gnList[j]->edges[numEdges-1]==NULL)
    {
      fprintf(stderr, "Vertex2Edge: malloc error, returning...\n");
      return NULL;
    }
    
    /* edge and weight */
    strcpy(g->gnList[j]->edges[numEdges-1], L2);
    g->gnList[j]->edgeWeights[numEdges-1] = 1.0;


    if(!found_L2) /* not there, so create */
    {
#ifdef DEBUG
  printf("making node %s\n", L2);
#endif
      numNodes = ++(g->numNodes);
      
      g->gnList = (graphNode **) realloc(g->gnList, numNodes*sizeof(graphNode *));
      if(g->gnList==NULL)
      {
        fprintf(stderr, "Vertex2Edge: realloc error, returning...\n");
        return NULL;
      }
      
      g->gnList[numNodes-1] = (graphNode *) malloc(sizeof(graphNode));
      if(g->gnList[numNodes-1]==NULL)
      {
        fprintf(stderr, "Vertex2Edge: malloc error, returning...\n");
        return NULL;
      }
      
      g->gnList[numNodes-1]->label = (char *) malloc(strlen(L2)+1);
      if(g->gnList[numNodes-1]->label==NULL)
      {
        fprintf(stderr, "Vertex2Edge: malloc error, returning...\n");
        return NULL;
      }
      
      strcpy(g->gnList[numNodes-1]->label, L2);
      
      g->gnList[numNodes-1]->numEdges    = 0;
      g->gnList[numNodes-1]->edges       = NULL;
      g->gnList[numNodes-1]->edgeWeights = NULL;
    }
    
    /* find L2, and make an edge from L2 to L1 */
    for(j=0; j<g->numNodes; j++)
      if(!strcmp(g->gnList[j]->label, L2))
        break;
        
    numEdges = ++(g->gnList[j]->numEdges);
      
    g->gnList[j]->edges       = (char **)  realloc(g->gnList[j]->edges,
                                                   numEdges*sizeof(char *));
    g->gnList[j]->edgeWeights = (double *) realloc(g->gnList[j]->edgeWeights,
                                                   numEdges*sizeof(double));
    if(g->gnList[j]->edges==NULL || g->gnList[j]->edgeWeights==NULL)
    {
      fprintf(stderr, "Vertex2Edge: realloc error, returning...\n");
      return NULL;
    }
      
    g->gnList[j]->edges[numEdges-1] = (char *) malloc(strlen(L1)+1);
    if(g->gnList[j]->edges[numEdges-1]==NULL)
    {
      fprintf(stderr, "Vertex2Edge: malloc error, returning...\n");
      return NULL;
    }
    
    /* edge and weight */
    strcpy(g->gnList[j]->edges[numEdges-1], L1);
    g->gnList[j]->edgeWeights[numEdges-1] = 1.0;
    
    free(L1);
    free(L2);
  }
  
  return g;
}

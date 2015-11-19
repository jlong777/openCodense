/*
 ODES_master.c                                           jlong@jimlong.org
 
 A breadth-first algorithm to find dense subgraphs with density >= 1/2
 
 Copyright (C) 2010 James Long
 
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
 
 Experimental:
 This exact algorithm can be combined with a faster heuristic routine by passing
 it a NULL terminated list of edges found by the heuristic routine. This list 
 will be excluded from the list of edges that initialize the algorithm.
 
 The format of each list entry is "label<whitespace>label". 
 
 Using this feature runs the risk of not finding dense subgraphs that overlap 
 the set of excluded edges.
*/

#ifndef CODENSE_H
  #include "codense.h"
  #define CODENSE_H
#endif

#define BUFSZ 1024
#define MAX_SUBG_VERTS 100
#define REALLOCLEN     100
//#define DEBUG_INSERTNODE
//#define DEBUG_THREAD
//#define UNPROVEN_CONJECTURE

int InsertNodeIntoSortedGraph(gintNode *node, gintGraph **subGA);
int * Limbs(gintGraph *inG, gintGraph *subG);
void FreeIntGraph(gintGraph *gp);
void FreeSubGraphChain(subGraphChain *sgc);
void PrintIntGraph(gintGraph *gp);
void ThreadFunc(thread_args *ta);

int lenL[NUM_ODE_THREADS], product[NUM_ODE_THREADS], rLen[NUM_ODE_THREADS], 
    totNumSubGraphs[NUM_ODE_THREADS];
int **hits[NUM_ODE_THREADS], subGnodes[NUM_ODE_THREADS][MAX_SUBG_VERTS];
char **labelMap;
gintGraph *gintG, ** potentialGraphsToReturn[NUM_ODE_THREADS];
gnGraph **denseSubGraphsToReturn[NUM_ODE_THREADS];
subGraphChain *L[NUM_ODE_THREADS], *masterCurr, *sgcCurr[NUM_ODE_THREADS];

#ifndef EXCLUDED_EDGES /* experimental, set in codense.h */
gnGraph ** ODES(gnGraph *inG, 
                double desiredDensity, 
                int minNodes)
#else
gnGraph ** ODES(gnGraph *inG, 
                double desiredDensity, 
                int minNodes,
                char ** excludedEdges)
#endif
{
  int i, j, k, returnLen;
  void *vp;
  gintNode *node;
  gnGraph **denseSubGraphs=NULL;
  subGraphChain *sgcPt;
  
  pthread_t   pt[NUM_ODE_THREADS];
  thread_args ta[NUM_ODE_THREADS];
  
#ifdef EXCLUDED_EDGES
  int lenExcludedIntEdges=0, **excludedIntEdges=NULL, tmp[2];
  char buf[BUFSZ], *tok_0, *tok_1;
#endif

#ifdef UNPROVEN_CONJECTURE
  int minDegree;
  
  /*
   an initial edge must connect vertices that have at least the average degree 
   required for a dense subgraph with the minimum number of vertices.
  */
  minDegree = (int)(desiredDensity * (double)(minNodes-1));
#endif
  
  desiredDensity  = desiredDensity - 0.000000001;

  if(minNodes < 2)
  {
    printf("ODES: minNodes must be > 1, and should be > 3, exiting...\n");
    return NULL;
  }
  
  if(inG==NULL)
  {
    fprintf(stderr, "ODES: input graph is NULL, returning...\n");
    return NULL;
  }
  
  /* sort inG on labels so labelMap can be searched */
  qsort((void *)(inG->gnList), inG->numNodes, sizeof(graphNode *), compGN);
  
  /*
   inG represents nodes as character strings, which are time consuming for
   a searching algorithm due to string comparisons. So we will convert the 
   strings to integers, and store the mapping in labelMap.
  */
  gintG = (gintGraph *) malloc(sizeof(gintGraph));
  if(gintG==NULL)
  {
    fprintf(stderr, "ODES: malloc error, returning...\n");
    return NULL;
  }
  gintG->numNodes = inG->numNodes;
  
  gintG->nodeList = (gintNode **) malloc(gintG->numNodes*sizeof(gintGraph *));
  if(gintG->nodeList==NULL)
  {
    fprintf(stderr, "ODES: malloc error, returning...\n");
    return NULL;
  }
  
  labelMap = (char **) malloc(gintG->numNodes*sizeof(char *));
  if(labelMap==NULL)
  {
    fprintf(stderr, "ODES: malloc error, returning...\n");
    return NULL;
  }
  
  for(i=0; i<gintG->numNodes; i++)
  {
    gintG->nodeList[i] = (gintNode *) malloc(sizeof(gintNode));
    if(gintG->nodeList[i]==NULL)
    {
      fprintf(stderr, "ODES: malloc error, returning...\n");
      return NULL;
    }
    
    gintG->nodeList[i]->label       = i;
    gintG->nodeList[i]->edgeWeights = inG->gnList[i]->edgeWeights;
    gintG->nodeList[i]->numEdges    = inG->gnList[i]->numEdges;
    
    gintG->nodeList[i]->edges = (int *) malloc(gintG->nodeList[i]->numEdges*sizeof(int));
    if(gintG->nodeList[i]->edges==NULL)
    {
      fprintf(stderr, "ODES: malloc error, returning...\n");
      return NULL;
    }
    
    labelMap[i] = inG->gnList[i]->label;
  }
  
#ifdef DEBUG_THREAD
  for(i=0; i<inG->numNodes; i++) printf("labelMap[%d] = %s\n", i, labelMap[i]);
#endif

  /* search the labelMap to build the edge lists */
  for(i=0; i<inG->numNodes; i++)
  {
    for(j=0; j<inG->gnList[i]->numEdges; j++)
    {
      /* find edge in labelMap */
      vp = (void *)labelMap;
      vp = bsearch((void *)inG->gnList[i]->edges[j], vp, inG->numNodes,
                   sizeof(char *), compStrLM);
      /*
       (*(char **)vp) is = inG->gnList[i]->edges[j] 
       (char **)vp - labelMap is the index in the labelMap
      */
      gintG->nodeList[i]->edges[j] = (char **)vp - labelMap;
    }
  }
  
#ifdef EXCLUDED_EDGES
  /*
   search the labelMap to convert labels in excludedEdges to ints, which are
   then placed into an int array with the lowest one first
  */
  if(excludedEdges != NULL)
  {
    i = lenExcludedIntEdges = 0;
    while(excludedEdges[i] != NULL)
    {
      lenExcludedIntEdges++;
      excludedIntEdges = (int **) realloc(excludedIntEdges, (i+1)*sizeof(int *));
      if(excludedIntEdges==NULL)
      {
        fprintf(stderr, "ODES: realloc error, returning...\n");
        return NULL;
      }
      
      excludedIntEdges[i] = (int *) malloc(2*sizeof(int));
      if(excludedIntEdges[i]==NULL)
      {
        fprintf(stderr, "ODES: malloc error, returning...\n");
        return NULL;
      }

      strncpy(buf, excludedEdges[i], BUFSZ);
      tok_0 = strtok(buf,  " \t"); /* spaces, tabs */
      tok_1 = strtok(NULL, " \t");
      
#ifdef DEBUG_THREAD 
      printf("tok_0 = %s  tok_1 = %s\n", tok_0, tok_1);
#endif

      /* convert to ints via LabelMap */
      vp = (void *)labelMap;
      vp = bsearch((void *)tok_0, vp, inG->numNodes, sizeof(char *), compStrLM);
      excludedIntEdges[i][0] = (char **)vp - labelMap;
      
      vp = (void *)labelMap;
      vp = bsearch((void *)tok_1, vp, inG->numNodes, sizeof(char *), compStrLM);
      excludedIntEdges[i][1] = (char **)vp - labelMap;
      
      if(excludedIntEdges[i][0] > excludedIntEdges[i][1])
      {
                             j = excludedIntEdges[i][0]; 
        excludedIntEdges[i][0] = excludedIntEdges[i][1]; 
        excludedIntEdges[i][1] = j;
      }

      i++;
    }
    
#ifdef DEBUG_THREAD 
    printf("before:\n");
    for(j=0; j<i; j++)
      printf("%d   %d\n", excludedIntEdges[j][0], excludedIntEdges[j][1]);
#endif

    qsort((void *)excludedIntEdges, lenExcludedIntEdges, sizeof(int *), compIntA);
    
#ifdef DEBUG_THREAD 
    printf("after:\n");
    for(j=0; j<i; j++)
      printf("%d   %d\n", excludedIntEdges[j][0], excludedIntEdges[j][1]);
#endif
  }
#endif
  
  
  /* init */
  for(i=0; i<NUM_ODE_THREADS; i++)
  {
    L[i] = (subGraphChain *) malloc(sizeof(subGraphChain));
    if(L[i]==NULL)
    {
      fprintf(stderr, "ODES: malloc error, returning...\n");
      return NULL;
    }
    L[i]->subG = NULL;
    L[i]->lmbs = NULL;
    L[i]->next = NULL;
    
    denseSubGraphsToReturn[i]  = NULL;
    potentialGraphsToReturn[i] = NULL;
    totNumSubGraphs[i]         = 0;
    ta[i].idx                  = i;
    ta[i].minNodes             = minNodes;
    ta[i].numThreads           = NUM_ODE_THREADS;
    ta[i].desiredDensity       = desiredDensity;
  }
  
  /* spawn all but last thread */
  for(i=0; i<NUM_ODE_THREADS-1; i++)
    pthread_create(&pt[i], NULL, (void *(*)(void*))ThreadFunc, (void *)(ta+i));
  
  /*
   Initialize list L[0] with sorted dense subgraphs consisting of one edge.
   Also initialize a list of limbs for each subG in L[0].
  */
  masterCurr = L[0];
  for(i=0; i<gintG->numNodes; i++)
  {
    node = gintG->nodeList[i];
    
#ifdef UNPROVEN_CONJECTURE
    if(node->numEdges < minDegree) continue;
#endif

    for(j=0; j<node->numEdges; j++)
    {
      if(node->edges[j] < node->label) continue; /* pairs are sorted */
      
#ifdef EXCLUDED_EDGES
      /* check if this edge "label - edges[j]" is an excluded edge */
      if(excludedEdges != NULL)
      {
        tmp[0] = node->label;
        tmp[1] = node->edges[j];
        
        vp = (void *)excludedIntEdges;
        vp = bsearch((void *)tmp, vp, lenExcludedIntEdges, sizeof(int *), compIntA2);
      
        if(vp != NULL) continue;
      }
#endif
      
      sgcPt = (subGraphChain *) malloc(sizeof(subGraphChain));
      if(sgcPt==NULL)
      {
        fprintf(stderr, "ODES: malloc error, returning...\n");
        return NULL;
      }
      sgcPt->subG = NULL;
      sgcPt->lmbs = NULL;
      sgcPt->next = NULL;
      
      InsertNodeIntoSortedGraph(node,  &(sgcPt->subG));
      
      /* find node in gintG whose label is node->edges[j] */
      vp = (void *)(gintG->nodeList);
      vp = bsearch((void *)&(node->edges[j]), vp, gintG->numNodes,
                   sizeof(gintNode *), compIntGintN);

      /* (*(gintNode **)vp)->label is = node->edges[j] */
#ifdef UNPROVEN_CONJECTURE
      if((*(gintNode **)vp)->numEdges < minDegree) continue;
#endif
      InsertNodeIntoSortedGraph(*(gintNode **)vp, &(sgcPt->subG));

      sgcPt->lmbs = Limbs(gintG, sgcPt->subG);
      
#ifdef DEBUG_THREAD
      printf("master gave ");fflush(NULL);
      for(k=0; k<(sgcPt->subG)->numNodes; k++)
        printf("%d ", sgcPt->subG->nodeList[k]->label);
      printf("to 0\n");
      fflush(NULL);
#endif

      /* place sgcPt in thread[0] buffer */
      masterCurr->next = sgcPt;
      masterCurr       = sgcPt;
    }
  }

  /* place NULL subgraph in buffer */
  sgcPt = (subGraphChain *) malloc(sizeof(subGraphChain));
  if(sgcPt==NULL)
  {
    fprintf(stderr, "ODES: malloc error, returning...\n");
    return NULL;
  }
  sgcPt->subG = NULL;
  sgcPt->lmbs = NULL;
  sgcPt->next = NULL;
  
  masterCurr->next = sgcPt;
  masterCurr       = sgcPt;
  
  /* spawn last thread (needs to know masterCurr) */
  pthread_create(&pt[NUM_ODE_THREADS-1], NULL, (void *(*)(void*))ThreadFunc, 
                 (void *)(ta+NUM_ODE_THREADS-1));
  
#ifdef DEBUG_THREAD
  printf("master waiting on threads\n");
  fflush(NULL);
#endif

  /* wait on threads */
  for(i=0; i<NUM_ODE_THREADS; i++)
    pthread_join(pt[i], NULL);
  
  /* gather the results */
  returnLen = 0;
  for(i=0; i<NUM_ODE_THREADS; i++)
  {
    j = 0;
    while(denseSubGraphsToReturn[i][j++] != NULL);
    
    returnLen += j-1;
  }

  denseSubGraphs = (gnGraph **) realloc(denseSubGraphs, 
                                        (returnLen+1)*sizeof(gnGraph *));
  if(denseSubGraphs==NULL)
  {
    fprintf(stderr, "ODES: realloc error, returning...\n");
    return NULL;
  }
  
  k = 0;
  for(i=0; i<NUM_ODE_THREADS; i++)
  { 
    j = 0;
    while(denseSubGraphsToReturn[i][j] != NULL)
      denseSubGraphs[k++] = denseSubGraphsToReturn[i][j++];
  
    free(denseSubGraphsToReturn[i]);
  }
  
  denseSubGraphs[returnLen] = NULL;
  
  /* clean up */
  /* can't do a FreeIntGraph(gintG), because some things point to inG, so: */
  for(i=0; i<gintG->numNodes; i++)
  {
    if(gintG->nodeList[i]->edges != NULL)
      free(gintG->nodeList[i]->edges);
      
    free(gintG->nodeList[i]);
  }
  free(gintG->nodeList);
  free(gintG);
  free(labelMap);
  
#ifdef EXCLUDED_EDGES
  i = 0;
  while(excludedEdges[i] != NULL)
    free(excludedEdges[i++]);

  free(excludedEdges);
#endif
  
  for(i=0; i<NUM_ODE_THREADS; i++)
    FreeSubGraphChain(sgcCurr[i]);
  
  return denseSubGraphs;
}

/*
 inserts a node into a sorted graph and its edges that connect to the graph 
*/
int InsertNodeIntoSortedGraph(gintNode *node, gintGraph **subGA)
{
  int i, j, newNodeIndex;
  gintGraph *subG;
  gintNode *n, *p;

#ifdef DEBUG_INSERTNODE
  printf("InsertNodeIntoSortedGraph: start *********************************************\n");
  printf("InsertNodeIntoSortedGraph: adding node %d to subG\n", node->label);
  fflush(NULL);
#endif

  /* allocate space for node */
  if(*subGA==NULL)
  {
    *subGA = (gintGraph *) malloc(sizeof(gintGraph));
    if(*subGA==NULL)
    {
      fprintf(stderr, "InsertNodeIntoSortedGraph: malloc error, returning...\n");
      return 1;
    }
    
    (*subGA)->nodeList = (gintNode **) malloc(sizeof(gintNode *));
    if((*subGA)->nodeList==NULL)
    {
      fprintf(stderr, "InsertNodeIntoSortedGraph: malloc error, returning...\n");
      return 1;
    }
    
    (*subGA)->numNodes = 0;
  }
  else
  {
    (*subGA)->nodeList = (gintNode **) realloc((*subGA)->nodeList, 
                                              ((*subGA)->numNodes+1) *
                                               sizeof(gintNode *));
    if((*subGA)->nodeList==NULL)
    {
      fprintf(stderr, "InsertNodeIntoSortedGraph: realloc error, returning...\n");
      return 1;
    }
  }
  
  subG = *subGA;
  subG->numNodes += 1;

  /* merge node->label into sorted subG->nodeList */
  for(newNodeIndex=0; newNodeIndex<subG->numNodes-1; newNodeIndex++)
    if(node->label < subG->nodeList[newNodeIndex]->label)
      break;
  
  /* shift elements of subG->nodeList and insert node */
  for(i=subG->numNodes-1; i>newNodeIndex; i--)
    subG->nodeList[i] = subG->nodeList[i-1];

  subG->nodeList[newNodeIndex] = (gintNode *) malloc(sizeof(gintNode));
  if(subG->nodeList[newNodeIndex]==NULL)
  {
    fprintf(stderr, "InsertNodeIntoSortedGraph: malloc error, returning...\n");
    return 1;
  }

  n = subG->nodeList[newNodeIndex]; /* the new node for the graph */
  
  n->label       = node->label;
  n->edges       = NULL;
  n->edgeWeights = NULL;
  n->numEdges    = 0;
  
#ifdef DEBUG_INSERTNODE
  printf("InsertNodeIntoSortedGraph: subG->numNodes = %d\n", subG->numNodes); 
  PrintIntGraph(subG);
#endif
  
  /*
   An edge of "node" that connects to another node in subG
   is added to the "edges" list of the node it connects to
  */
  for(i=0; i<node->numEdges; i++)
  {
    for(j=0; j<subG->numNodes; j++)
    {
      if(j==newNodeIndex) continue; /* don't look at ourself */
      
      p = subG->nodeList[j]; /* an existing node in subG */
      if(node->edges[i] == p->label)
      {
        if(p->edges==NULL)
        {
          p->edges       = (int *)    malloc(sizeof(int));
          p->edgeWeights = (double *) malloc(sizeof(double));
          p->numEdges = 0;
        }
        else
        {
          p->edges       = (int *)    realloc(p->edges,       (p->numEdges+1)*sizeof(int));
          p->edgeWeights = (double *) realloc(p->edgeWeights, (p->numEdges+1)*sizeof(double));
        }
        
        if(p->edges==NULL || p->edgeWeights==NULL)
        {
          fprintf(stderr, "InsertNodeIntoSortedGraph: malloc error, returning...\n");
          return 1;
        }
        
        p->edges[p->numEdges]       = node->label; /* adds edge to p */
        p->edgeWeights[p->numEdges] = node->edgeWeights[i];
        p->numEdges++;
        
        /* add edge to node */
        if(n->edges==NULL)
        {
          n->edges       = (int *)    malloc(sizeof(int));
          n->edgeWeights = (double *) malloc(sizeof(double));
          n->numEdges = 0;
        }
        else
        {
          n->edges       = (int *)    realloc(n->edges,       (n->numEdges+1)*sizeof(int));
          n->edgeWeights = (double *) realloc(n->edgeWeights, (n->numEdges+1)*sizeof(double));
        }
        
        if(n->edges==NULL || n->edgeWeights==NULL)
        {
          fprintf(stderr, "InsertNodeIntoSortedGraph: malloc error, returning...\n");
          return 1;
        }
        
        n->edges[n->numEdges]       = p->label; /* adds edge to n */
        n->edgeWeights[n->numEdges] = node->edgeWeights[i];
        n->numEdges++;
        
        break;
      }
    }
  }
  
#ifdef DEBUG_INSERTNODE
  PrintIntGraph(subG);
  printf("InsertNodeIntoSortedGraph: end ***********************************************\n");
  fflush(NULL);
#endif

  return 0;
}

/*
 find limbs of a 2-node graph subG:
   1) for each of the 2 nodes of subG
   2)   find the node in gintG
   3)   for each edge of the node in gintG
   4)     add to list those not internal to subG and not already in list,
          keeping the list sorted, with its length in limbs[0], i.e.
          if there are no edges in limbs, limbs[0] = 0, even though
          the length of limbs is 1.
*/
int * Limbs(gintGraph *gintG, gintGraph *subG)
{
  int i, j, k, add, inSubG;
  int *limbs=NULL;
  void *vp;
  gintNode *n, **sl;
  
  if(!subG->numNodes)
  {
    fprintf(stderr, "Limbs: subgraph has no nodes, returning NULL...\n");
    return NULL;
  }
  
  sl = subG->nodeList;

  limbs = (int *) malloc(sizeof(int));
  if(limbs==NULL)
  {
    fprintf(stderr, "Limbs: malloc error, returning...\n");
    return NULL;
  }
  limbs[0] = 0;

  /* 1) */
  for(i=0; i<2; i++)
  {
    /* 2) */
    vp = (void *)(gintG->nodeList);
    vp = bsearch((void *)&(sl[i]->label), vp, gintG->numNodes,
                  sizeof(gintNode *), compIntGintN);
    
    n = *(gintNode **)vp; /* n->label is = sl[i]->label in inG */

    /* 3) */
    for(j=0; j<n->numEdges; j++)
    {
      inSubG = 0;
      for(k=0; k<2; k++)
      {
        if(n->edges[j]==sl[k]->label)
        {
          inSubG = 1;
          break;
        }
      }

      /* 4) */
      if(!inSubG)
      {
        add = 1;
        if(limbs[0])
        {
          vp = (void *)(limbs+1);
          vp = bsearch((void *)&(n->edges[j]), vp, limbs[0],
                       sizeof(int), compInt);
          if(vp!=NULL) add = 0;
        }
        
        if(add)
        {
          limbs[0]++;
          limbs = (int *) realloc(limbs, (limbs[0]+1)*sizeof(int));
          if(limbs==NULL)
          {
            fprintf(stderr, "Limbs: realloc error, returning...\n");
            return NULL;
          }

          /* keep limbs sorted */
          k = limbs[0];
          while(n->edges[j] < limbs[k-1] && k>1)
          {
            limbs[k] = limbs[k-1];
            k--;
          }
          limbs[k] = n->edges[j];
        }
      }
    }
  }

  return limbs;
}

/*
 a binary search to determine if subG with "node" added is already in a 
 sorted list of graphs L, returning -1 if it is there, or the index of 
 where it should be inserted if it is not.
*/
int BinarySearchOfSortedGraphList(int *subGnodes, int numNodes, int **L, int lenL)
{
  int i, comp, found, hi, isGreaterThan, isLessThan, lo, mid;
  
  found = lo = 0;
  hi = lenL;
  while(lo < hi)
  {
    mid = lo + ((hi-lo)/2);
    i = isGreaterThan = isLessThan = 0;
    
    /* compare portion of subG with nodes less than label */
    while(i < numNodes)
    {
      comp = L[mid][i] - subGnodes[i];
      
      if(comp < 0)
      {
        isLessThan = 1;
        break;
      }
      else if(comp > 0) 
      {
        isGreaterThan = 1;
        break; 
      }
      
      i++;
    }
    
    if(!isLessThan && !isGreaterThan) /* must be equal now, so found it */
      return -1;
    
    if(isLessThan)
      lo = mid + 1;
    else
      hi = mid;
  }

  /* lo == hi , are they equal at lo? */
  if(lo < lenL) /* lo == lenL is insertion at end */
  {
    found = 1;
    i = 0;
    while(i < numNodes)
    {
      if(L[mid][i] != subGnodes[i])
      {
        found = 0; 
        break;
      }
    
      i++;
    }
  }
  
  if(found)
    return -1;
  else
    return lo;
}

/*
 make a copy of limbs, removing the entry = node->label, and adding
 entries from the edges of the node that are not already in the subG,
 returning a sorted list
*/
int * AddNode2Limbs(int *limbs, gintNode *node, gintGraph *subG)
{
  int i, j, inLimbs, inSubG, matchIndex;
  int *p;
  void *vp;
  
  /* make a copy of limbs, removing node->label */
  matchIndex = -1;
  for(i=1; i<limbs[0]+1; i++)
  {
    if(node->label == limbs[i])
    {
      matchIndex = i;
      break;
    }
  }
  
  p = (int *) malloc(limbs[0]*sizeof(int));
  if(p==NULL)
  {
    fprintf(stderr, "AddNode2Limbs: malloc error, returning...\n");
    return NULL;
  }
  
  j = 1;
  for(i=1; i<limbs[0]; i++)
  {
    if(i==matchIndex) j++; /* don't copy node->label */

    p[i] = limbs[j++];
  }
  p[0] = limbs[0] - 1;
  
  /* add entries from the edges of node that are not already in subG */
  for(i=0; i<node->numEdges; i++)
  {
    inSubG = 0;

    /* is node->edges[i] in subG node list? */
    vp = (void *)(subG->nodeList);
    vp = bsearch((void *)&(node->edges[i]), vp, subG->numNodes, 
                 sizeof(gintNode *), compIntGintN);
    if(vp!=NULL) inSubG = 1;

    /* make sure the edge is not already in limbs */
    if(!inSubG)
    {
      inLimbs = 0;
      if(limbs[0])
      {
        vp = (void *)(limbs+1);
        vp = bsearch((void *)&(node->edges[i]), vp, limbs[0],
                     sizeof(int), compInt);
        if(vp!=NULL) inLimbs = 1;
      }

      if(!inLimbs)
      {
        p[0]++;   /* the number of elements, not counting itself */
        p = (int *) realloc(p, (p[0]+1)*sizeof(int));
        if(p==NULL)
        {
          fprintf(stderr, "AddNode2Limbs: realloc error, returning...\n");
          return NULL;
        }
    
        p[p[0]] = node->edges[i];
      }
    }
  }
  
  if(limbs[0])
    qsort((void *)(p+1), p[0], sizeof(int), compInt);
  
  return p;
}

/* mainly for debugging */
void PrintIntGraph(gintGraph *gp)
{
  int i, j, len;
  gintNode *p;
  
  printf("    Node        Edges                   Weights\n");
  printf("    ===========================================\n");
  
  for(i=0; i<gp->numNodes; i++)
  {
    p = gp->nodeList[i];
  
    printf("%3d) %d\t\t", i, p->label); 
  
    len = 0;
    if(p->edges!=NULL)
    {
      for(j=0; j<p->numEdges; j++)
      {
        printf("%d ", p->edges[j]);
        fflush(stdout);
        
        if(p->edges[j] < 10)
          len += 2;
        else if(p->edges[j] < 100)
          len += 3;
        else
          len += 4;
      }
      
      if(len < 8)
        printf("\t\t\t");
      else if(len < 16)
        printf("\t\t");
      else
        printf("\t");
    
      if(p->edgeWeights!=NULL)
        for(j=0; j<p->numEdges; j++)
          printf("%3.1lf ", p->edgeWeights[j]);
          
    }
    printf("\n");
  }
  
  printf("    ===========================================\n\n");
  fflush(stdout);
}

/* make a complete copy of a graph */ 
gintGraph * CopyIntGraph(gintGraph *graph)
{
  int i, j;
  gintGraph *copy;
  gintNode *cn, *gn;
  
  copy = (gintGraph *) malloc(sizeof(gintGraph));
  if(copy==NULL)
  {
    fprintf(stderr, "CopyIntGraph: malloc error, returning...\n");
    return NULL;
  }
  copy->numNodes = graph->numNodes;
  
  copy->nodeList = (gintNode **) malloc(copy->numNodes * sizeof(gintNode *));
  if(copy->nodeList==NULL)
  {
    fprintf(stderr, "CopyIntGraph: malloc error, returning...\n");
    return NULL;
  }
  
  /* copy every node of graph */
  for(i=0; i<graph->numNodes; i++)
  {
    gn = graph->nodeList[i];
    
    cn = copy->nodeList[i] = (gintNode *) malloc(sizeof(gintNode));
    if(cn==NULL)
    {
      fprintf(stderr, "CopyIntGraph: malloc error, returning...\n");
      return NULL;
    }
    
    cn->label = gn->label;

    if(gn->numEdges)
    {
      cn->edges       = (int *)    malloc(gn->numEdges * sizeof(int));
      cn->edgeWeights = (double *) malloc(gn->numEdges * sizeof(double));
      if(cn->edges==NULL || cn->edgeWeights==NULL)
      {
        fprintf(stderr, "CopyIntGraph: malloc error, returning...\n");
        return NULL;
      }
      cn->numEdges = gn->numEdges;
        
      for(j=0; j<cn->numEdges; j++)
      {
        cn->edges[j]       = gn->edges[j];
        cn->edgeWeights[j] = gn->edgeWeights[j];
      }  
    }
    else
    {
      cn->edges = NULL;
      cn->edgeWeights = NULL;
      cn->numEdges = 0;
    }
  }
  
  return copy;
}

/* free everything that gp points to, and gp itself */
void FreeIntGraph(gintGraph *gp)
{
  int i;
  gintNode *p;
  
  if(gp==NULL) return;

  for(i=0; i<gp->numNodes; i++)
  {
    p = gp->nodeList[i];
    
    if(p->edges != NULL)
      free(p->edges);
      
    if(p->edgeWeights != NULL)
      free(p->edgeWeights);
    
    free(p);
  }
  
  free(gp->nodeList);
  free(gp);
  return;
}

void FreeSubGraphChain(subGraphChain *sgc)
{
  if(sgc==NULL) return;
  if(sgc->subG != NULL) FreeIntGraph(sgc->subG);
  if(sgc->lmbs != NULL) free        (sgc->lmbs);
  free(sgc);
}

void GIntGraphs2GnGraphs(int idx)
{
  int i, j, k;
  gnGraph *gnG;

  /* convert surviving subgraphs back into gnGraphs */
  denseSubGraphsToReturn[idx] = (gnGraph **) malloc((totNumSubGraphs[idx]+1)*sizeof(gnGraph *));
  if(denseSubGraphsToReturn[idx]==NULL)
  {
    fprintf(stderr, "GIntGraphs2GnGraphs: malloc error, thread exiting...\n");
    fflush(stderr);
    return;
  }

  for(i=0; i<totNumSubGraphs[idx]; i++)
  {
    gnG = (gnGraph *) malloc(sizeof(gnGraph));
    if(gnG==NULL)
    {
      fprintf(stderr, "GIntGraphs2GnGraphs: malloc error, thread exiting...\n");
      fflush(stderr);
      return;
    }
    gnG->numNodes = potentialGraphsToReturn[idx][i]->numNodes;
  
    gnG->gnList = (graphNode **) malloc(gnG->numNodes*sizeof(gnGraph *));
    if(gnG->gnList==NULL)
    {
      fprintf(stderr, "GIntGraphs2GnGraphs: malloc error, thread exiting...\n");
      fflush(stderr);
      return;
    }
  
    for(j=0; j<gnG->numNodes; j++)
    {
      gnG->gnList[j] = (graphNode *) malloc(sizeof(graphNode));
      if(gnG->gnList[j]==NULL)
      {
        fprintf(stderr, "GIntGraphs2GnGraphs: malloc error, thread exiting...\n");
        fflush(stderr);
        return;
      }
  
      gnG->gnList[j]->label = (char *) malloc(strlen(labelMap[potentialGraphsToReturn[idx][i]->nodeList[j]->label])+1);
      if(gnG->gnList[j]->label==NULL)
      {
        fprintf(stderr, "GIntGraphs2GnGraphs: malloc error, thread exiting...\n");
        fflush(stderr);
        return;
      }
  
      strcpy(gnG->gnList[j]->label, labelMap[potentialGraphsToReturn[idx][i]->nodeList[j]->label]);
      gnG->gnList[j]->numEdges = potentialGraphsToReturn[idx][i]->nodeList[j]->numEdges;

      gnG->gnList[j]->edges       = (char **) malloc(gnG->gnList[j]->numEdges*sizeof(char *));
      gnG->gnList[j]->edgeWeights = (double *)malloc(gnG->gnList[j]->numEdges*sizeof(double));
      if(gnG->gnList[j]->edges==NULL || gnG->gnList[j]->edgeWeights==NULL)
      {
        fprintf(stderr, "GIntGraphs2GnGraphs: malloc error, thread exiting...\n");
        fflush(stderr);
        return;
      }
  
      for(k=0; k<gnG->gnList[j]->numEdges; k++)
      {
        gnG->gnList[j]->edges[k] = (char *) malloc(strlen(labelMap[potentialGraphsToReturn[idx][i]->nodeList[j]->edges[k]])+1);
        if(gnG->gnList[j]->edges[k]==NULL)
        {
          fprintf(stderr, "GIntGraphs2GnGraphs: malloc error, thread exiting...\n");
          fflush(stderr);
          return;
        }
  
        strcpy(gnG->gnList[j]->edges[k], labelMap[potentialGraphsToReturn[idx][i]->nodeList[j]->edges[k]]);
        gnG->gnList[j]->edgeWeights[k] = potentialGraphsToReturn[idx][i]->nodeList[j]->edgeWeights[k];
      }
    }
  
    FreeIntGraph(potentialGraphsToReturn[idx][i]);
    denseSubGraphsToReturn[idx][i] = gnG;
  }
  
  free(potentialGraphsToReturn[idx]);
  denseSubGraphsToReturn[idx][i] = NULL;
}

int ProcessSubGraph(int idx, gintGraph *subG, int *limbs, int minNodes,
                    double desiredDensity)
{
  int i, j, k, canBeExtended, existingEdges, insertionIndex, **ipt, numNewEdges;
  void *vp;
  double den;
  register int tmp;
  gintNode *node;
  subGraphChain *sgcPt;

  if(limbs==NULL)/* list of nodes outside of subG connected to
                    nodes inside subG, limbs[0] contains count */
  {
    fprintf(stderr, "ProcessSubGraph: limbs error, exiting...\n");
    fflush(stderr);
    return 1;
  }

  /* count existing edges in subG */
  existingEdges = 0;
  for(j=0; j<subG->numNodes; j++)
    existingEdges += subG->nodeList[j]->numEdges;
  existingEdges /= 2; /* each edge was counted twice */

  canBeExtended = 0;
  for(k=1; k<limbs[0]+1; k++)
  {
    /* find "node" in gintG whose label is limbs[k] */
    vp = (void *)(gintG->nodeList);
    vp = bsearch((void *)&(limbs[k]), vp, gintG->numNodes,
                 sizeof(gintNode *), compIntGintN);
    node = *(gintNode **)vp; /* node->label is = limbs[k] */
    
    /* if an edge of "node" goes into subG, it is a new edge */
    numNewEdges = 0;
    for(j=0; j<node->numEdges; j++)
    {
      /* is node->edges[j] in subG node list? */
      vp = (void *)(subG->nodeList);
      vp = bsearch((void *)&(node->edges[j]), vp, subG->numNodes, 
                   sizeof(gintNode *), compIntGintN);
      if(vp!=NULL)
        numNewEdges++;
    }
  
    /*
     compute density with limbs[k] added = (existingEdges + 
     number of edges of limbs[k] that connect to nodes in subG) / 
     (((subG->numNodes + 1) * (subG->numNodes))/2)
    */
    den = (double)(existingEdges + numNewEdges) /
          (double)(((subG->numNodes + 1) * (subG->numNodes))/2.0);

    if(den >= desiredDensity)
    {
      canBeExtended = 1;
      
      /*
       A binary search to determine if the subG with "node" was already 
       added to L[1]; if not, add it.
       I used to do the binary search on L[1], but that won't work for the 
       threaded version, because I want to pass what would go into L[1] on
       to the next thread. 
      */
      i   = 0;
      tmp = 1;
      for(j=0; j<subG->numNodes; j++)
      {
        if(node->label < subG->nodeList[j]->label && tmp)
        {
          subGnodes[idx][i++] = node->label;
          tmp = 0;
          j--;
        }
        else
          subGnodes[idx][i++] = subG->nodeList[j]->label;
      }
      
      if(tmp)
        subGnodes[idx][i] = node->label;
      
      insertionIndex = BinarySearchOfSortedGraphList(subGnodes[idx], 
                                                     subG->numNodes+1,
                                                     hits[idx], lenL[idx]);

      if(insertionIndex == -1) /* it's already there */
       continue;
  
      /* it's not there, so add to buffer L[(idx+1)%NUM_ODE_THREADS] */
      lenL[idx]++;
      product[idx] = 1;
      if(lenL[idx] > rLen[idx])
      {
        rLen[idx] += REALLOCLEN;
        ipt = (int **) realloc(hits[idx], rLen[idx]*sizeof(int *));
        if(ipt==NULL)
        {
          fprintf(stderr, "ProcessSubGraph: realloc error, exiting...\n");
          fflush(stderr);
          return 1;
        }
	else
	  hits[idx] = ipt;
      }
     
      sgcPt = (subGraphChain *) malloc(sizeof(subGraphChain));
      if(sgcPt==NULL)
      {
        fprintf(stderr, "ProcessSubGraph: malloc error, exiting...\n");
        fflush(stderr);
        return 1;
      }
      sgcPt->subG = NULL;
      sgcPt->lmbs = NULL;
      sgcPt->next = NULL;
     
      /* keep hits ordered */
      for(j=lenL[idx]-1; j>insertionIndex; j--)
        hits[idx][j] = hits[idx][j-1];
     
      hits[idx][insertionIndex] = (int *) malloc((subG->numNodes+1)*sizeof(int));
      if(hits[idx][insertionIndex]==NULL)
      {
        fprintf(stderr, "ProcessSubGraph: malloc error, exiting...\n");
        fflush(stderr);
        return 1;
      }
     
      memcpy((void *)hits[idx][insertionIndex], (void *)subGnodes[idx], 
             (subG->numNodes+1)*sizeof(int));
            
      sgcPt->subG = CopyIntGraph(subG);
      InsertNodeIntoSortedGraph(node, &(sgcPt->subG));
     
      sgcPt->lmbs = AddNode2Limbs(limbs, node, subG);

      /* place sgcPt in thread[(idx+1)%NUM_ODE_THREADS] buffer */
      sgcCurr[idx]->next = sgcPt;
      sgcCurr[idx]       = sgcPt;
    }
  }

  /*
   if subG was not extended, check if it has sufficient 
   nodes to be added to potentialGraphsToReturn
  */
  if(!canBeExtended && subG->numNodes >= minNodes)
  {
    totNumSubGraphs[idx]++;
    potentialGraphsToReturn[idx] = (gintGraph **) realloc(potentialGraphsToReturn[idx],
                                                          (totNumSubGraphs[idx]+1) * 
                                                          sizeof(gintGraph *));
    if(potentialGraphsToReturn[idx]==NULL)
    {
      fprintf(stderr, "ProcessSubGraph: realloc error, exiting...\n");
      fflush(stderr);
      return 1;
    }
    
    potentialGraphsToReturn[idx][totNumSubGraphs[idx]-1] = CopyIntGraph(subG);
    potentialGraphsToReturn[idx][totNumSubGraphs[idx]]   = NULL;
  }
  
  return 0;
}

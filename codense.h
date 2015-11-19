/*
 codense.h                                          01/2008 jlong@jimlong.org
 
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

#define CODENSE_H
#define NUM_CMD_THREADS    8 /* number of omp threads for CreateMicroarrayDataSetGraph */
#define NUM_ODE_THREADS    8 /* number of p-threads for ODES, master and slave */

#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <search.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

typedef struct microarrayDataSet_jlong
{
  int numCols, numRows;
  char **labels; /* each label is unique */
  double *data;  /* really a data table  */
} microarrayDataSet;

typedef struct graphNode_jlong
{
  char *label;   /* from microarrayDataSet */
  char **edges;  /* array of graphNode labels that the node connects to */
  double *edgeWeights;
  int numEdges;
} graphNode;

typedef struct gnGraph_jlong
{
  graphNode **gnList;
  int numNodes;
} gnGraph;

typedef struct edgeSupportVector_jlong
{
  char *edgeName;
  double *supportVector;
} edgeSupportVector;

typedef struct edgeSupportMatrix_jlong
{
  edgeSupportVector **esvList; /* NULL terminated */
  int vecLen;                  /* length of supportVectors */
} edgeSupportMatrix;

/* structures for more efficient searching */
typedef struct gintNode_jlong
{
  int label;
  int *edges;
  double *edgeWeights;
  int numEdges;
} gintNode;

typedef struct gintGraph_jlong
{
  gintNode **nodeList;
  int numNodes;
} gintGraph;

typedef struct subGraphChain_jlong subGraphChain;
struct subGraphChain_jlong
{
  gintGraph *subG;
  int *lmbs;
  subGraphChain *next;
};

typedef struct thread_args_jlong
{
  int idx;
  int minNodes;
  int numThreads;
  double desiredDensity;
} thread_args;


/* codense routines */
gnGraph ** CodenseMI(char * fileList);
gnGraph * CreateMicroarrayDataSetGraph(const microarrayDataSet *dataSet,
                                       double zpcSlice,      double pcCutOff, 
                                       double zmiSlice,      double miCutOff,
                                       double negCutOff,     double peak,
                                       double geneThreshold, double sumfactor,
                                       int dim, int doMI,    int doS, int doZ,
                                       int plusCorrelateOnly);
gnGraph * CreateSummaryGraph(gnGraph **graphSet,
                             int numberOfDataSets,
                             double support);
void DeallocateMicroarrayDataSet(microarrayDataSet *);
void FreeEdgeSupportMatrix(edgeSupportMatrix *esMatrix);
edgeSupportMatrix * GenerateEdgeSupportMatrix(gnGraph *subGraph,
                                               gnGraph **graphSet,
                                               int numberOfDataSets);
gnGraph * GenerateSecondOrderGraph(const edgeSupportMatrix *esMatrix,
                                   double zpcSlice,      double pcCutOff, 
                                   double zmiSlice,      double miCutOff,
                                   double negCutOff,     double peak,
                                   double geneThreshold, double sumfactor,
                                   int dim, int doMI,    int doS, int doZ,
                                   int plusCorrelateOnly);
gnGraph * Vertex2Edge(gnGraph *ds);
microarrayDataSet * ReadMicroarrayData(char *infileName);

/* utility functions defined in util.c */
int compGN      (const void *p1, const void *p2);
int compInt     (const void *p1, const void *p2);
int compIntGintN(const void *p1, const void *p2);
int compString  (const void *p1, const void *p2);
int compStrGN   (const void *p1, const void *p2);
int compStrLM   (const void *p1, const void *p2);
double Density(gnGraph *gp);
void FreeGraph(gnGraph *gp);
int NumGraphEdges(gnGraph *gp);
int PrintGraph(gnGraph *gp);
gnGraph * ReadGraph(char *filename);
gnGraph * CopyGraph(gnGraph *graph);

/* dense subgraph algorithm */
gnGraph ** ODES(gnGraph *bigGraph, 
                double density, 
                int minNodes);

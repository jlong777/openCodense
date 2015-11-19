/*
 CreateMicroarrayDataSetGraph.c                   01/2008 jlong@jimlong.org
 
 input a microarrayDataSet and return a graph, where nodes in the graph
 represent genes, and edges represent significant expression correlation
 between genes.
 
 Algorithm:
 for each gene
   compute Pearson's correlation between the gene and all others
     "significant" correlation constitutes an edge
   compute mutual information (MI) between the gene and all others
     "significant" MI constitutes an edge
     
 note - "significance" for a correlation is computed in a manner similar to
 the CLR algorithm, using z-scores for each gene's correlation value.
 
 If we are interested in only "+" correlations, then only "+" Pearson's are
 used, and only those MI for which the Pearson's is not strongly "-". The 
 z-score computation for Pearson's in this case is just 
 z_i + z_j > max(z_i) - ZPC_SLICE * (max(z_i) - min(z_i)) + 
             max(z_j) - ZPC_SLICE * (max(z_j) - min(z_j))
 where z_i is the z-score for gene i computed from its set of correlation
 scores, and similarly for the z_j of gene j.
 
 For "+" and "-" correlations, the z-score computation for Pearson's is
 sqrt(z_i^2 + z_j^2) > max(z_i) - ZPC_SLICE * (max(z_i) - min(z_i)) + 
                       max(z_j) - ZPC_SLICE * (max(z_j) - min(z_j))
 
 Since MI is always positive, it's z-score computation is always 
 z_i + z_j > max(z_i) - ZMI_SLICE * (max(z_i) - min(z_i)) + 
             max(z_j) - ZMI_SLICE * (max(z_j) - min(z_j))
 
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

#define USE_OPENMP

#ifdef USE_OPENMP
#include <alloca.h>
#endif

//#define FAST_LOG2 /* only effective if doMI != 0; recommended only for testing */

#define FUZ 7 /* must be odd, and change kdf matrix accordingly */
/* 
  kernel density function that will estimate each observation, i.e. 
  observations are "fuzzy", & center value will be supplied at runtime
*/
double kdf[FUZ][FUZ]=  /* make observations "fuzzy" */
{
  {16.0,  32.0,  48.0,  64.0,  48.0,  32.0,  16.0},
  {32.0,  64.0,  96.0, 128.0,  96.0,  64.0,  32.0},
  {48.0,  96.0, 192.0, 256.0, 192.0,  96.0,  48.0},
  {64.0, 128.0, 256.0,   0.0, 256.0, 128.0,  64.0},
  {48.0,  96.0, 192.0, 256.0, 192.0,  96.0,  48.0},
  {32.0,  64.0,  96.0, 128.0,  96.0,  64.0,  32.0},
  {16.0,  32.0,  48.0,  64.0,  48.0,  32.0,  16.0}
};
  
/* stats */
double Mean(const double *data, int num);
double Sigma(const double *data, int num); /* standard deviation */
double Similarity(const double *x, const double *y, int len);
double PearsonsCorrelation(const double *x, const double *y, 
                           const double geneThreshold, int len);
double MutualInformation(const double *x, const double *y, 
                               double sumfactor, int dim, int len);
#ifdef USE_OPENMP
double MutualInformationOMP(const double *x, const double *y, 
                                  double sumfactor, int dim, int len);
#endif
gnGraph * CreateMicroarrayDataSetGraph(const microarrayDataSet *dataSet,
                                       double zpcSlice,      double pcCutOff,
                                       double zmiSlice,      double miCutOff,
                                       double negCutOff,     double peak,
                                       double geneThreshold, double sumfactor, 
                                       int dim,  int doMI,   int doS, int doZ,
                                       int plusCorrelateOnly)
{
  int i, j, edgePC, gene, numEdges=0, numGenes, row, col;
  char **cpt;
  double *dpt, *maxZPC=NULL, *minZPC=NULL, pc, *rho=NULL, *rhoMean=NULL, *rhoSigma=NULL;
  graphNode **nodeList;
  gnGraph *graphPt;
  
  /* for MI */
  int edgeMI;
  double *maxZMI=NULL, *minZMI=NULL, mi, *ti=NULL, *tiMean=NULL, *tiSigma=NULL;
  
  if(doS && doZ) 
  {
    fprintf(stderr, "CreateMicroarrayDataSetGraph: doS & doZ cannot both be set, returning...\n");
    return NULL;
  }
  
  if(doS && doMI) 
  {
    fprintf(stderr, "CreateMicroarrayDataSetGraph: doS & doMI cannot both be set, returning...\n");
    return NULL;
  }
  
  if(doMI && (dim<=0 || peak<=0.0))
  {
    fprintf(stderr, "CreateMicroarrayDataSetGraph: initialization error, returning...\n");
    return NULL;
  }

  /* the graph */
  graphPt = (gnGraph *) malloc(sizeof(gnGraph));
  if(graphPt==NULL)
  {
    fprintf(stderr, "CreateMicroarrayDataSetGraph: malloc error, returning...\n");
    return NULL;
  }
  
  col = dataSet->numCols;
  numGenes = row = dataSet->numRows;

  if(doMI)
    kdf[FUZ/2][FUZ/2] = peak; /* regulates bandwidth of KDF */
  
  /* init graph with nodes, i.e. one graphNode per gene */
  graphPt->numNodes  = numGenes;
  graphPt->gnList = (graphNode **) malloc(numGenes * sizeof(graphNode *));
  if((graphPt->gnList)==NULL)
  {
    fprintf(stderr, "CreateMicroarrayDataSetGraph: malloc error, returning...\n");
    return NULL;
  }
  nodeList = graphPt->gnList;
  
  /*
   instantiate graphNodes with their nodes, the rest of 
   this code will determine the edges by various means
  */
  for(i=0; i<numGenes; i++)
  {
    nodeList[i] = (graphNode *) malloc(sizeof(graphNode));
    if(nodeList[i]==NULL)
    {
      fprintf(stderr, "CreateMicroarrayDataSetGraph: malloc error, returning...\n");
      return NULL;
    }
    
    nodeList[i]->label = (char *) malloc(strlen(dataSet->labels[i])+1);
    if((nodeList[i]->label)==NULL)
    {
      fprintf(stderr, "CreateMicroarrayDataSetGraph: malloc error, returning...\n");
      return NULL;
    }
    
    strcpy(nodeList[i]->label, dataSet->labels[i]);
    nodeList[i]->numEdges    = 0;
    nodeList[i]->edges       = NULL;
    nodeList[i]->edgeWeights = NULL;
  }
  
  /* distribution of correlation scores */
  if(doZ)
  {
    /*
     compute mean and standard deviation for the 
     set of Pearson's correlation for each gene 
    */
    maxZPC   = (double *) malloc(numGenes * sizeof(double));
    minZPC   = (double *) malloc(numGenes * sizeof(double));
    rhoMean  = (double *) malloc(numGenes * sizeof(double));
    rhoSigma = (double *) malloc(numGenes * sizeof(double));
#ifdef USE_OPENMP
    if(maxZPC==NULL || minZPC==NULL || rhoMean==NULL || rhoSigma==NULL)
#else
    rho      = (double *) malloc(numGenes * sizeof(double));
    if(maxZPC==NULL || minZPC==NULL || rho==NULL || rhoMean==NULL || rhoSigma==NULL)
#endif
    {
      fprintf(stderr, "CreateMicroarrayDataSetGraph: malloc error, returning...\n");
      return NULL;
    }
    
#ifdef USE_OPENMP
    #pragma omp parallel for private(i,j,rho) num_threads(NUM_CMD_THREADS)
#endif
    for(gene=0; gene<numGenes; gene++)
    {
      i = 0;
#ifdef USE_OPENMP
      if(!i)
      {
        /* 16,000 genes is only 128KB, which should fit in any modern stack */
        rho = (double *) alloca(numGenes * sizeof(double)); /* since allocated in stack, freed on function return */
        if(rho==NULL) fprintf(stderr, "CreateMicroarrayDataSetGraph: rho alloca error in thread, do manual abort!!!\n");
      }
#endif
      maxZPC[gene] = 0.0;
      minZPC[gene] = 99999.9;
      for(j=0; j<numGenes; j++)
      {
        if(j == gene)
          continue;
        else
	{
          rho[i] = PearsonsCorrelation((dataSet->data)+(gene*col),
                                       (dataSet->data)+   (j*col), 
                                       geneThreshold, col);
          if(rho[i] > maxZPC[gene])
	    maxZPC[gene] = rho[i];  /* max PC */
	    
          if(rho[i] < minZPC[gene])
	    minZPC[gene] = rho[i];  /* min PC */
	    
	  i++;
	}
      }
    
      rhoMean[gene]  = Mean (rho, numGenes-1);
      rhoSigma[gene] = Sigma(rho, numGenes-1);
      
      maxZPC[gene] = (maxZPC[gene]-rhoMean[gene])/rhoSigma[gene];  /* max PC z-score */
      minZPC[gene] = (minZPC[gene]-rhoMean[gene])/rhoSigma[gene];  /* min PC z-score */
    }
  }
  
  if(doMI && doZ)
  {
    /*
     compute mean and standard deviation for the 
     set of MutualInformation scores for each gene 
    */
    maxZMI   = (double *) malloc(numGenes * sizeof(double));
    minZMI   = (double *) malloc(numGenes * sizeof(double));
    tiMean   = (double *) malloc(numGenes * sizeof(double));
    tiSigma  = (double *) malloc(numGenes * sizeof(double));
#ifdef USE_OPENMP
    if(maxZMI==NULL || minZMI==NULL || tiMean==NULL  || tiSigma==NULL)
#else
    ti       = (double *) malloc(numGenes * sizeof(double));
    if(maxZMI==NULL || minZMI==NULL || ti==NULL || tiMean==NULL  || tiSigma==NULL)
#endif
    {
      fprintf(stderr, "CreateMicroarrayDataSetGraph: malloc error, returning...\n");
      return NULL;
    }
    
    /* distribution of MutualInformation scores */
#ifdef USE_OPENMP
    #pragma omp parallel for private(i,j,ti) num_threads(NUM_CMD_THREADS)
#endif
    for(gene=0; gene<numGenes; gene++)
    {
      i = 0;
#ifdef USE_OPENMP
      if(!i)
      {
        /* 16,000 genes is only 128KB, which should fit in any modern stack */
        ti = (double *) alloca(numGenes * sizeof(double)); /* since allocated in stack, freed on function return */
        if(ti==NULL) fprintf(stderr, "CreateMicroarrayDataSetGraph: ti alloca error in thread, do manual abort!!!\n");
      }
#endif
      maxZMI[gene] = 0.0;
      minZMI[gene] = 99999.9;
      for(j=0; j<numGenes; j++)
      {
        if(j == gene)
          continue;
        else
	{
#ifdef USE_OPENMP
          ti[i] = MutualInformationOMP((dataSet->data)+(gene*col),
                                       (dataSet->data)+   (j*col), 
                                        sumfactor,      dim, col);
/* test omp vs non-omp scores */
/*double foo = MutualInformation((dataSet->data)+(gene*col), (dataSet->data)+(j*col), sumfactor, dim, col);
if(fabs(foo - ti[i]) > .000001)
{
  fprintf(stderr, "mi = %g  miomp = %g\n", foo, ti[i]);
  fprintf(stderr, "mutual information mismatch, exiting...\n"); 
  exit(1);
}*/
#else
          ti[i] = MutualInformation((dataSet->data)+(gene*col),
                                    (dataSet->data)+   (j*col), 
                                     sumfactor,      dim, col);
#endif
          if(ti[i] > maxZMI[gene])
	    maxZMI[gene] = ti[i];  /* max MI, will be max z-score below */
	    
          if(ti[i] < minZMI[gene])
	    minZMI[gene] = ti[i];  /* min MI, will be min z-score below */
            
	  i++;
	}
      }
    
      tiMean[gene]  = Mean (ti, numGenes-1);
      tiSigma[gene] = Sigma(ti, numGenes-1);
      
      maxZMI[gene] = (maxZMI[gene]-tiMean[gene])/tiSigma[gene];  /* max MI z-score */
      minZMI[gene] = (minZMI[gene]-tiMean[gene])/tiSigma[gene];  /* min MI z-score */
    }
  }

  /*
     now look for edge-declaring correlations, i.e. 
     1) those with high correlation (>= cutoff score), or
     2) those with high correlation (>= cutoff score), and
        whose combined z-scores are significant (with or without MI) 
  */
  for(gene=0; gene<numGenes; gene++)
  {
    for(i=gene+1; i<numGenes; i++)
    {
      edgeMI = edgePC = 0;
      
      if(doS)
      {
        pc = Similarity((dataSet->data)+(gene*col),
                        (dataSet->data)+   (i*col), 
                         col);
      }
      else
      {
        /* Pearson's correlation edge? */
        pc = PearsonsCorrelation((dataSet->data)+(gene*col),
                                 (dataSet->data)+   (i*col), 
                                  geneThreshold, col);
      }
      
      if(pc>1.000001 || pc<-1.000001) fprintf(stderr, "CreateMicroarrayDataSetGraph: WARNING, pearsons out-of-bounds: %lf\n", pc);

      if(plusCorrelateOnly) /* positive correlations only */
      {
        if(doZ)
        {
          if(isnan(rhoSigma[gene]) || isnan(rhoSigma[i])) edgePC = 0; /* stays zero */
          else if((pc>=pcCutOff)                   && 
                 ((pc-rhoMean[gene])/rhoSigma[gene] +                 /* combined z-scores */
                  (pc-rhoMean[i]   )/rhoSigma[i] > maxZPC[gene] - zpcSlice*(maxZPC[gene] - minZPC[gene]) + 
	                                           maxZPC[i]    - zpcSlice*(maxZPC[i]    - minZPC[i])))
            edgePC = 1;
        }
        else if(pc>=pcCutOff) edgePC = 1;
      }
      else /* negative correlations also */
      {
        if(doZ)
        {
          if(isnan(rhoSigma[gene]) || isnan(rhoSigma[i])) edgePC = 0; /* stays zero */
          else if((pc>=pcCutOff || pc<=negCutOff)          && 
                  (sqrt(((pc-rhoMean[gene])/rhoSigma[gene]) *
                        ((pc-rhoMean[gene])/rhoSigma[gene]) + 
                        ((pc-rhoMean[i]   )/rhoSigma[i])    *
                        ((pc-rhoMean[i]   )/rhoSigma[i])) > maxZPC[gene] - zpcSlice*(maxZPC[gene] - minZPC[gene]) + 
	                                                    maxZPC[i]    - zpcSlice*(maxZPC[i]    - minZPC[i])))
            edgePC = 1;
        }
        else if(pc>=pcCutOff || pc<=negCutOff) edgePC = 1;
      }
/*if(edgePC)
{
  printf("pc edge declared between %d and %d, pc = %lf\n", gene, i, pc);
  printf("rhoMean[%d] = %lf  rhoSigma[%d] = %lf\n", gene, rhoMean[gene], gene, rhoSigma[gene]);
  printf("rhoMean[%d] = %lf  rhoSigma[%d] = %lf\n",    i, rhoMean[i],       i, rhoSigma[i]);
}*/
      /* no need to doMI if an edge is already declared, or if pc is low */
      if(doMI && !edgePC && (pc>0.6))
      {
        /* mutual information edge? */
#ifdef USE_OPENMP
        mi = MutualInformationOMP((dataSet->data)+(gene*col), 
                                  (dataSet->data)+   (i*col),
                                   sumfactor,      dim, col);
/* test omp vs non-omp scores */
/*if(fabs(mi - MutualInformation((dataSet->data)+(gene*col), (dataSet->data)+(i*col), sumfactor, dim, col)) > .000001)
{
  fprintf(stderr, "mi = %10.9lf   miomp = %10.9lf\n", MutualInformation((dataSet->data)+(gene*col), (dataSet->data)+(i*col), sumfactor, dim, col), mi);
  fprintf(stderr, "mutual information mismatch, exiting...\n"); 
  exit(1);
}*/
#else
        mi = MutualInformation((dataSet->data)+(gene*col), 
                               (dataSet->data)+   (i*col),
                                sumfactor,      dim, col);
#endif
        if(mi<0.0) fprintf(stderr, "CreateMicroarrayDataSetGraph: WARNING, mi out-of-bounds: %lf\n", mi);
      
        if(plusCorrelateOnly) /* no strong negative Pearson's */
        {
          if(pc>negCutOff)
          {
            if(doZ)
            {
              if(isnan(tiSigma[gene]) || isnan(tiSigma[i])) edgeMI = 0; /* stays zero */
              else if((mi>=miCutOff)                 && 
                     ((mi-tiMean[gene])/tiSigma[gene] + 
                      (mi-tiMean[i])/tiSigma[i] > maxZMI[gene] - zmiSlice*(maxZMI[gene] - minZMI[gene]) + 
	                                          maxZMI[i]    - zmiSlice*(maxZMI[i]    - minZMI[i])))
                edgeMI = 1;
            }
            else if(mi>=miCutOff) edgeMI = 1;
          }
        }
        else /* any strong MI is OK */
        {
          if(doZ)
          {
            if(isnan(tiSigma[gene]) || isnan(tiSigma[i])) edgeMI = 0; /* stays zero */
            else if((mi>=miCutOff)                 && 
                   ((mi-tiMean[gene])/tiSigma[gene] + 
                    (mi-tiMean[i])/tiSigma[i] > maxZMI[gene] - zmiSlice*(maxZMI[gene] - minZMI[gene]) + 
	                                        maxZMI[i]    - zmiSlice*(maxZMI[i]    - minZMI[i])))
              edgeMI = 1;
          }
          else if(mi>=miCutOff) edgeMI = 1;
        }
      }
      
      if(edgeMI || edgePC)
      {
        /* instantiate an edge between gene and i */
        if(nodeList[gene]->numEdges == 0)
        {
          nodeList[gene]->edges       = (char **)  malloc(sizeof(char *));
          nodeList[gene]->edgeWeights = (double *) malloc(sizeof(double));
          if((nodeList[gene]->edges)==NULL || 
             (nodeList[gene]->edgeWeights)==NULL)
          {
            fprintf(stderr, "CreateMicroarrayDataSetGraph: malloc error, returning...\n");
            return NULL;
          }
        }
        else
        {
          cpt =  (char **) realloc(nodeList[gene]->edges,       (nodeList[gene]->numEdges + 1) *
                                                                 sizeof(char *));
          dpt = (double *) realloc(nodeList[gene]->edgeWeights, (nodeList[gene]->numEdges + 1) *
                                                                 sizeof(double));
          if(cpt==NULL || dpt==NULL)
          {
            fprintf(stderr, "CreateMicroarrayDataSetGraph: realloc error, returning...\n");
            return NULL;
          }
	  else
	  {
	     nodeList[gene]->edges       = cpt;
	     nodeList[gene]->edgeWeights = dpt;
	  }
        }

        nodeList[gene]->edges[nodeList[gene]->numEdges] = (char *)malloc(strlen(nodeList[i]->label)+1);
        if((nodeList[gene]->edges[nodeList[gene]->numEdges])==NULL)
        {
          fprintf(stderr, "CreateMicroarrayDataSetGraph: malloc error, returning...\n");
          return NULL;
        }
        strcpy(nodeList[gene]->edges[nodeList[gene]->numEdges], 
               nodeList[i]->label);

        /* if an edge is inferred, set weight to 1.0 */
        nodeList[gene]->edgeWeights[nodeList[gene]->numEdges] = 1.0;
        
        /*
         future thought:
         if edge was inferred only by PC, use pc for edge weight
         if edge was inferred only by MI, need to set to whatever between 0 & 1
         
         i.e.
         if(edgePC)
           nodeList[gene]->edgeWeights[nodeList[gene]->numEdges] = pc;
         else if(edgeMI)
           nodeList[gene]->edgeWeights[nodeList[gene]->numEdges] = whatever;
        */
        
        nodeList[gene]->numEdges++;

        /* instantiate an edge between i and gene */
        if(nodeList[i]->numEdges == 0)
        {
          nodeList[i]->edges       = (char **)  malloc(sizeof(char *));
          nodeList[i]->edgeWeights = (double *) malloc(sizeof(double));
          if((nodeList[i]->edges)==NULL || (nodeList[i]->edgeWeights)==NULL)
          {
            fprintf(stderr, "CreateMicroarrayDataSetGraph: malloc error, returning...\n");
            return NULL;
          }
        }
        else
        {
          nodeList[i]->edges       = (char **)  realloc(nodeList[i]->edges,
                                                       (nodeList[i]->numEdges + 1) *
                                                       sizeof(char *));
          nodeList[i]->edgeWeights = (double *) realloc(nodeList[i]->edgeWeights,
                                                       (nodeList[i]->numEdges + 1) *
                                                       sizeof(double));
          if(nodeList[i]->edges==NULL || nodeList[i]->edgeWeights==NULL)
          {
            fprintf(stderr, "CreateMicroarrayDataSetGraph: realloc error, returning...\n");
            return NULL;
          }
        }

        nodeList[i]->edges[nodeList[i]->numEdges] = (char *)malloc(strlen(nodeList[gene]->label)+1);
        if((nodeList[i]->edges[nodeList[i]->numEdges])==NULL)
        {
          fprintf(stderr, "CreateMicroarrayDataSetGraph: malloc error, returning...\n");
          return NULL;
        }
        strcpy(nodeList[i]->edges[nodeList[i]->numEdges], 
               nodeList[gene]->label);

        /* if an edge is inferred, set weight to 1.0 */
        nodeList[i]->edgeWeights[nodeList[i]->numEdges] = 1.0;
        
        /*
         future thought:
         if edge was inferred only by PC, use pc for edge weight
         if edge was inferred only by MI, need to set to whatever between 0 & 1
         
         i.e.
         if(edgePC)
           nodeList[i]->edgeWeights[nodeList[i]->numEdges] = pc;
         else if(edgeMI)
           nodeList[i]->edgeWeights[nodeList[i]->numEdges] = whatever;
        */
        
        nodeList[i]->numEdges++;
        
        edgeMI = edgePC = 0;
        
        numEdges++;
      }
    }
  }
  
  if(doZ)
  {
    free(maxZPC);
    free(minZPC);
#ifndef USE_OPENMP
    free(rho);
#endif
    free(rhoMean);
    free(rhoSigma);
  
    if(doMI)
    {
      free(maxZMI);
      free(minZMI);
#ifndef USE_OPENMP
      free(ti);
#endif
      free(tiMean);
      free(tiSigma);
    }
  }
/*  
printf("numEdges  %d\n", numEdges);
PrintGraph(graphPt);
*/
  return graphPt;
}

double Mean(const double *data, int num)
{
  int i;
  double mean=0.0;
  
  for(i=0; i<num; i++)
    mean += data[i];
  
  return mean/(double)num;
}

/* standard deviation */
double Sigma(const double *data, int num)
{
  int i;
  double mean=0.0, sum=0.0, tmp;
  
  for(i=0; i<num; i++)
    mean += data[i];
  
  mean /= (double) num;
  
  for(i=0; i<num; i++)
  {
    tmp = data[i] - mean;
    sum += tmp*tmp;
  }
  
  return sqrt(sum/(double)num);
}

/* from http://davidmlane.com/hyperstat/A51911.html */
double PearsonsCorrelation(const double *x, const double *y, 
                           const double geneThreshold, int len)
{
  int i, numSmallX=0, numSmallY=0, same;
  double dblBuf, sumx=0.0, sumxx=0.0, sumxy=0.0, sumy=0.0, sumyy=0.0, xx, xy, yy;
  
  for(i=0; i<len; i++)
  {
    sumx  += x[i];
    sumxx += x[i] * x[i];
    sumxy += x[i] * y[i];
    sumy  += y[i];
    sumyy += y[i] * y[i];
  }

  // old criteria
  /* at least one gene is considered off too much of the time */
  //if((sumx < len*geneThreshold) || (sumy < len*geneThreshold)) return 0.0;
  
  // new criteria
  /*
   one sum is << than the other, or both are < 10^-12 (considered noise),
   or > len/2 of the entries for one gene are < 10^-13
  */
  if(((sumx < geneThreshold * sumy) || (sumy < geneThreshold * sumx)) ||
     ((sumx < 0.000000000001)       && (sumy < 0.000000000001)))
    return 0.0;

  /* < 10^-13 */
  for(i=0; i<len; i++)
  {
    if(x[i] < 0.0000000000001) numSmallX++;
    if(y[i] < 0.0000000000001) numSmallY++;
  }
  
  if((numSmallX > len/2) || (numSmallY > len/2)) return 0.0;
  
  /* 
   if x[i] = y[i] for each i, and sumx = sumxx = len, we'll get isnan(), 
   but should return 1.0 anytime x[i] = y[i] for each i. geneThreshold
   criteria above will eliminate cases where x[i] = y[i] = 0 for each i
  */
  if(sumx == sumy)
  {
    same = 1;
    for(i=0; i<len; i++)
      if(x[i] != y[i])
      {
        same = 0;
        break;
      }
    if(same) return 1.0;
  }
  
  xx = sumxx - (sumx * sumx)/len;
  xy = sumxy - (sumx * sumy)/len;
  yy = sumyy - (sumy * sumy)/len;
   
  dblBuf = xy / sqrt(xx * yy);
  
  if(isnan(dblBuf))      return 0.0;
  else if(isinf(dblBuf)) return 0.0;
  else                   return dblBuf;
}

#ifdef FAST_LOG2
/* fast log2, thanks to http://www.musicdsp.org/showone.php?id=63 

Fast log2
References : Posted by Laurent de SorasCode :
inline float fast_log2 (float val)
{
   assert (val > 0);

   int * const  exp_ptr = reinterpret_cast <int *> (&val);
   int          x = *exp_ptr;
   const int    log_2 = ((x >> 23) & 255) - 128;
   x &= ~(255 << 23);
   x += 127 << 23;
   *exp_ptr = x;

   return (val + log_2);
}
Comments
from : tobybear[AT]web[DOT]de
comment : And here is some native Delphi/Pascal code that does the same thing: function fast_log2(val:single):single; var log2,x:longint; begin x:=longint((@val)^); log2:=((x shr 23) and 255)-128; x:=x and (not(255 shl 23)); x:=x+127 shl 23; result:=single((@x)^)+log2; end; Cheers Toby www.tobybear.de

from : henry[AT]cordylus[DOT]de
comment : instead of using this pointer casting expressions one can also use a enum like this: enum FloatInt { float f;

from : henry[AT]cordylus[DOT]de
comment : instead of using this pointer casting expressions one can also use a enum like this: enum FloatInt { float f; int l; } p; and then access the data with: p.f = x; p.l >>= 23; Greetings, Henry

from : henry[AT]cordylus[DOT]de
comment : Sorry : didnt mean enum, ment UNION !!!

from : Laurent de Soras
comment : More precision can be obtained by adding the following line just before the return() : val = map_lin_2_exp (val, 1.0f / 2); Below is the function (everything is constant, so most operations should be done at compile time) : inline float map_lin_2_exp (float val, float k) { const float a = (k - 1) / (k + 1); const float b = (4 - 2*k) / (k + 1); // 1 - 3*a const float c = 2*a; val = (a * val + b) * val + c; return (val); } You can do the mapping you want for the range [1;2] -> [1;2] to approximate the function log(x)/log(2).

from : Laurent de Soras
comment : Sorry I meant log(x)/log(2) + 1
*/

float fastlog2(float f)
{
  int i, log_2, *pe;
  
  pe    = (int *) &f;
  i     = *pe;
  log_2 = ((i >> 23) & 0xFF) - 128;

  i &= ~(0xFF << 23);
  i +=   0x7F << 23;

  *pe = i;

  /*
   the following is for k = 0.455 from
   
   float map_lin_2_exp (float f, float k) 
   { 
     const float a = (k - 1) / (k + 1); 
     const float b = (4 - 2*k) / (k + 1); 
     const float c = 2*a; 
  
     f = (a * f + b) * f + c; 
  
     return (val); 
   }
   
   
  */
  f = (-0.374570447 * f + 2.12371134) * f - 0.749140893; /* k = 0.455 */
  //f = (-0.33333333333333 * f + 2) * f - 0.66666666666666; /* k = 0.5 */

  return (f + (float)log_2);
}
#endif

double MutualInformation(const double *x, const double *y, double sumfactor, int dim, int len)
{
  int i, j, k, row, col, numSmallX=0, numSmallY=0, trow, tcol;
  double ** jProbMatrix;
  double max, mi, min, minlog, minsum, minsumlog, sumx, sumy, sumxy, tsumx, tsumy;
  register double offset, tmp;

  if(dim%32)
  {
    fprintf(stderr, "MutualInformation: dim = %d must be zero mod 32, returning...\n", dim);
    return -1.0;
  }

  sumx = sumy = 0.0;
  for(i=0; i<len; i++){sumx += x[i]; sumy += y[i];}
  
  /*
   don't bother if one sum is > sumfactor times the other, 
   or both are < 10^-12 (considered noise),
   or > len/2 of the entries for one gene are < 10^-12
  */
  if(((sumx > sumfactor * sumy) || (sumy > sumfactor * sumx))  ||
     ((sumx < 0.000000000001)   && (sumy < 0.000000000001)))
    return 0.0;
    
  for(i=0; i<len; i++)
  {
    if(x[i] < 0.000000000001) numSmallX++;
    if(y[i] < 0.000000000001) numSmallY++;
  }
  
  if((numSmallX > len/2) || (numSmallY > len/2)) return 0.0;
  
  /* init */
  jProbMatrix = (double **) alloca(dim * sizeof(double*));
  if(jProbMatrix==NULL)
  {
    fprintf(stderr, "MutualInformation: FATAL malloc error, exiting...\n");
    exit(1);
  }
  
  tmp = 0.01/(double)(dim*dim);
  for(i=0; i<dim; i++)
  {
    jProbMatrix[i] = (double *) malloc(dim * sizeof(double));
    if(jProbMatrix[i]==NULL)
    {
      fprintf(stderr, "MutualInformation: FATAL malloc error, exiting...\n");
      exit(1);
    }
    
    for(j=0; j<dim; j++)
      jProbMatrix[i][j] = tmp;
  }

  /* 
  compute row and col factors for joint probability matrix, 
  so that all kdfs will fit into jProbMatrix w/o spilling over,
  and leaving a row/col border around the perimeter of jProbMatrix
  */
  max = 0.0; min = 9999.9;
  for(i=0; i<len; i++)
  {
    if(x[i] > max) max = x[i];
    if(y[i] > max) max = y[i];
    if(x[i] < min) min = x[i];
    if(y[i] < min) min = y[i];
  }

  tmp = (max - min)/(double)(dim-FUZ-2);
  offset = (int)floor(min/tmp) - 1;

  /* build the joint probability matrix */
  for(i=0; i<len; i++)
  { 
    row = (int)floor(x[i]/tmp) - offset + FUZ/2; /* make room for kdf */
    col = (int)floor(y[i]/tmp) - offset + FUZ/2;

    for(j=0; j<FUZ; j++)
    {
      trow = row-FUZ/2 + j;
      if(trow<0 || trow>=dim) continue;
      for(k=0; k<FUZ; k++)
      {
        tcol = col-FUZ/2 + k;
        if(tcol<0 || tcol>=dim) continue;
        
        jProbMatrix[trow][tcol] +=  kdf[j][k];
      }
    }
  }
  
  /* normalize, find min entry, and min row sum */
  tmp = 0.0;
  for(i=0; i<dim; i++)
    for(j=0; j<dim; j++)
      tmp += jProbMatrix[i][j];
      
  min = 1.0;
  for(i=0; i<dim; i++)
  {
    for(j=0; j<dim; j++)
    {
      jProbMatrix[i][j] /= tmp;
      
      if(jProbMatrix[i][j] < min)
        min = jProbMatrix[i][j];
    }
  }
  minlog    = min*log2(min);
  minsum    = dim*min;
  minsumlog = minsum*log2(minsum);
  
  /* print jProbMatrix */
/*
  for(i=0; i<dim; i++)
  {
    for(j=0; j<dim; j++)
      printf("%e  ", jProbMatrix[i][j]);
    printf("\n");
  }
  fflush(NULL);
*/

  /* compute H(X), H(Y), and H(X,Y) */
  sumx = sumy = sumxy = 0.0;
  for(i=0; i<dim; i++)
  {
    tsumx = tsumy = 0.0;
    for(j=0; j<dim; j++)
    {
      tmp    = jProbMatrix[i][j];
      tsumx += tmp;
      tsumy += jProbMatrix[j][i];

      /* H(X,Y) summation */

#ifdef FAST_LOG2
      sumxy += tmp*fastlog2((float)tmp);
#else
      sumxy += (tmp == min)? minlog: tmp*log2(tmp);
#endif
    }
    
    /* H(X), H(Y) summations */
#ifdef FAST_LOG2
    sumx += tsumx*fastlog2((float)tsumx);
    sumy += tsumy*fastlog2((float)tsumy);
#else
    sumx += (tsumx == minsum)? minsumlog: tsumx*log2(tsumx);
    sumy += (tsumy == minsum)? minsumlog: tsumy*log2(tsumy);
#endif
  }
  
  for(i=0; i<dim; i++)
    free(jProbMatrix[i]);

  sumx  = -sumx;  /* H(X) */
  sumy  = -sumy;  /* H(Y) */
  sumxy = -sumxy; /* H(X,Y) */
  
  mi = sumx + sumy - sumxy; /* mutual information */

  if(isnan(mi))      return 0.0;
  else if(isinf(mi)) return 0.0;
  else               return mi;
}

#ifdef USE_OPENMP
double MutualInformationOMP(const double *x, const double *y, double sumfactor, int dim, int len)
{
  int i, j, k, row, col, numSmallX=0, numSmallY=0, trow, tcol;
  double ** jProbMatrix;
  double max, mi, min, minlog, sumx, sumy, sumxy, tsumx, tsumy;
  register double offset, tmp;
  double *sumxi, *sumyi, *sumxyi;
  
  if(dim%32)
  {
    fprintf(stderr, "MutualInformationOMP: dim = %d must be zero mod 32, returning...\n", dim);
    return -1.0;
  }

  sumx = sumy = 0.0;
  for(i=0; i<len; i++){sumx += x[i]; sumy += y[i];}
  
  /*
   don't bother if one sum is > sumfactor times the other, 
   or both are < 10^-12 (considered noise),
   or > len/2 of the entries for one gene are < 10^-12
  */
  if(((sumx > sumfactor * sumy) || (sumy > sumfactor * sumx))  ||
     ((sumx < 0.000000000001)   && (sumy < 0.000000000001)))
    return 0.0;

  for(i=0; i<len; i++)
  {
    if(x[i] < 0.000000000001) numSmallX++;
    if(y[i] < 0.000000000001) numSmallY++;
  }
  
  if((numSmallX > len/2) || (numSmallY > len/2)) return 0.0;
  
  /* init */
  jProbMatrix = (double **) alloca(dim * sizeof(double*));
  sumxi       = (double *)  alloca(dim * sizeof(double));
  sumyi       = (double *)  alloca(dim * sizeof(double));
  sumxyi      = (double *)  alloca(dim * sizeof(double));
  if(jProbMatrix==NULL || sumxi==NULL || sumyi==NULL || sumxyi==NULL)
  {
    fprintf(stderr, "MutualInformationOMP: alloca error in thread, do manual abort!!!\n");
    return -1.0;
  }
  
  tmp = 0.01/(double)(dim*dim);
  for(i=0; i<dim; i++)
  {
    jProbMatrix[i] = (double *) malloc(dim * sizeof(double));
    if(jProbMatrix[i]==NULL)
    {
      fprintf(stderr, "MutualInformationOMP: malloc error in thread, do manual abort!!!\n");
      for(j=0; j<i; j++)
        free(jProbMatrix[j]);

      return -1.0;
    }
    
    for(j=0; j<dim; j++)
      jProbMatrix[i][j] = tmp;
  }
          
  /* 
  compute row and col factors for joint probability matrix, 
  so that all kdfs will fit into jProbMatrix w/o spilling over,
  and leaving a row/col border around the perimeter of jProbMatrix
  */
  max = 0.0; min = 9999.9;
  for(i=0; i<len; i++)
  {
    if(x[i] > max) max = x[i];
    if(y[i] > max) max = y[i];
    if(x[i] < min) min = x[i];
    if(y[i] < min) min = y[i];
  }

  tmp = (max - min)/(double)(dim-FUZ-2);
  offset = (int)floor(min/tmp) - 1;

  /* build the joint probability matrix */
  for(i=0; i<len; i++)
  { 
    row = (int)floor(x[i]/tmp) - offset + FUZ/2; /* make room for kdf */
    col = (int)floor(y[i]/tmp) - offset + FUZ/2;

    for(j=0; j<FUZ; j++)
    {
      trow = row-FUZ/2 + j;
      if(trow<0 || trow>=dim) continue;
      for(k=0; k<FUZ; k++)
      {
        tcol = col-FUZ/2 + k;
        if(tcol<0 || tcol>=dim) continue;
        
        jProbMatrix[trow][tcol] +=  kdf[j][k];
      }
    }
  }
  
  /* normalize, find min entry, and min row sum */
  tmp = 0.0;
  for(i=0; i<dim; i++)
    for(j=0; j<dim; j++)
      tmp += jProbMatrix[i][j];
      
  min = 1.0;
  for(i=0; i<dim; i++)
  {
    for(j=0; j<dim; j++)
    {
      jProbMatrix[i][j] /= tmp;
      
      if(jProbMatrix[i][j] < min)
        min = jProbMatrix[i][j];
    }
  }
  minlog = min*log2(min);
  
  
  /* print jProbMatrix */
/*
  for(i=0; i<dim; i++)
  {
    for(j=0; j<dim; j++)
      printf("%13.12lf  ", jProbMatrix[i][j]);
    printf("\n");
  }
*/

  /* compute H(X), H(Y), and H(X,Y) */
  sumx = sumy = sumxy = 0.0;
  for(i=0; i<dim; i++)
  {
    sumxi[i]  = 0.0;
    sumyi[i]  = 0.0;
    sumxyi[i] = 0.0;
  }
    
  #pragma omp parallel for private(j,tmp,tsumx,tsumy) num_threads(NUM_CMD_THREADS)
  for(i=0; i<dim; i++)
  {
    tsumx = tsumy = 0.0;
    for(j=0; j<dim; j++)
    {
      tmp    = jProbMatrix[i][j];
      tsumx += tmp;
      tsumy += jProbMatrix[j][i];

      /* H(X,Y) summation */
#ifdef FAST_LOG2
      sumxyi[i] += tmp*fastlog2((float)tmp);
#else
      sumxyi[i] += (tmp == min)? minlog: tmp*log2(tmp);
#endif
    }
    
    /* H(X), H(Y) summations */
#ifdef FAST_LOG2
    sumxi[i] += tsumx*fastlog2((float)tsumx);
    sumyi[i] += tsumy*fastlog2((float)tsumy);
#else
    sumxi[i] += tsumx*log2(tsumx);
    sumyi[i] += tsumy*log2(tsumy);
#endif
  }
  
  for(i=0; i<dim; i++)
    free(jProbMatrix[i]);

  for(i=0; i<dim; i++)
  {
    sumx  += sumxi[i];
    sumy  += sumyi[i];
    sumxy += sumxyi[i];
  }
  
  sumx  = -sumx;  /* H(X) */
  sumy  = -sumy;  /* H(Y) */
  sumxy = -sumxy; /* H(X,Y) */
  
  mi = sumx + sumy - sumxy; /* mutual information */

  if(isnan(mi))      return 0.0;
  else if(isinf(mi)) return 0.0;
  else               return mi;
}
#endif

/* What fraction of array values are the same in x & y? */
double Similarity(const double *x, const double *y, int len)
{
  int i, numeq=0;
  double xy, yx;
  
  if(len < 1) return 0.0;
  
  for(i=0; i<len; i++)
  {
    xy = x[i] - y[i];
    yx = y[i] - x[i];
    if(xy<0.00001 && yx<0.00001) numeq++;
  }
  
  return((double)numeq/(double)len);
}


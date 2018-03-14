#define REAL double

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <unistd.h>
#include "triangle.h"
#include "tristruct.h"
#include "util.h"
#include "base.h"
#include "mex.h"

#define check_size(p, m, n) \
	(p && (mxGetN(p) == n) && (mxGetM(p) == m))
	
void import_data(char *, int nrhs, const mxArray *prhs[]);

void export_data(char *, int nlhs, mxArray *plhs[]);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  char Temp[20] = "TRIXXXXXX";
  char format[300];

  mkstemp(Temp);
 
  debug_begin("triangle_debug.txt");
printf("Coucou1");
  import_data(Temp, nrhs, prhs);
printf("Coucou1");

  sprintf(format, "%s/triangle -rQea %s", BASE, Temp); 
  debug_printf("<%s>\n", format);
print(format);
//  system(format);


  export_data(Temp, nlhs, plhs);

  debug_end();
}

void import_data(char *Temp, int nrhs, const mxArray *prhs[])
{
  int i, nPoints = 0, nTriangles = 0;
  double *points = NULL, *triangles = NULL, 
         *areaTriangles = NULL;
  double *labelPoints = NULL, *labelTriangles = NULL;
  char fname[50];
  FILE *f;
  const mxArray *mesh, *p;

  if (nrhs != 2)
    mexErrMsgTxt("usage refined_mesh = refine(mesh, areas)");
 
  mesh = prhs[0];
  if (!mxIsStruct(mesh))
    mexErrMsgTxt("first parameter is not a structure");
  
  p = mxGetField(mesh, 0, "vertices");
  if (!p) 
    mexErrMsgTxt("no vertices array in mesh structure");

  if (mxGetN(p) != 2)
    mexErrMsgTxt("vertices must be in R^2");
  nPoints = mxGetM(p);
  points = mxGetPr(p);

  debug_printf("\n%d points\n", nPoints);
  for (i=0; i<nPoints; i++)
    debug_printf("%5d %10.4g %10.4g\n", i, points[i], points[nPoints+i]);

  p = mxGetField(mesh, 0, "triangles");
  if (!p) 
    mexErrMsgTxt("no triangles array in mesh structure");

  if (mxGetN(p) != 3)
    mexErrMsgTxt("triangles must have 3 and only 3 vertices");
  nTriangles = mxGetM(p);
  triangles = mxGetPr(p);
     
  p = mxGetField(mesh, 0, "lab_vertices");
  if (p && !check_size(p, nPoints, 1)) 
     mexErrMsgTxt("vertex labels must be a "
                   "vector whose size is the number of vertices");
  labelPoints = mxGetPr(p);

  p = mxGetField(mesh, 0, "lab_triangles");
  if (p && !check_size(p, nTriangles, 1)) 
     mexErrMsgTxt("triangle labels must be a "
                   "vector whose size is the number of triangles");
  labelTriangles = mxGetPr(p);

  if ((mxGetN(prhs[1]) != 1) && (mxGetM(prhs[1]) != nTriangles))
    mexErrMsgTxt("areas (2th parameter) is a vector whose length "
		 "must be the number of triangles");
  areaTriangles = mxGetPr(prhs[1]);

  debug_printf("\n%d triangles\n", nPoints);
  for (i=0; i<nTriangles; i++)
    debug_printf("%5d %5d %5d %5d area %10.4g\n", i, 
		                           (int) triangles[i], 
		                           (int) triangles[2*nTriangles+i],
		                           (int) triangles[nTriangles+i], 
					   areaTriangles[i]);


  sprintf(fname, "%s.node", Temp);
  f = fopen(fname, "w");
  fprintf(f, "%d 2 0 %d\n", nPoints, labelPoints ? 1 : 0);
  for (i=0; i<nPoints; i++) {
    fprintf(f, "%d %g %g ", i+1, points[i], points[i+nPoints]);
    fprintf(f, "%d", (int) labelPoints[i]);
    fprintf(f, "\n");
    }
  fclose(f);

  sprintf(fname, "%s.ele", Temp); 
  f = fopen(fname, "w");
  fprintf(f, "%d 3 %d\n", nTriangles, labelTriangles ? 1 : 0);
  for (i=0; i<nTriangles; i++) {
    fprintf(f, "%d %d %d %d", i+1, 
            (int) triangles[i], 
            (int) triangles[i+nTriangles], 
            (int) triangles[i+2*nTriangles]);
    fprintf(f, " %d", (int) labelTriangles[i]);
    fprintf(f, "\n");
  }
  fclose(f);

  sprintf(fname, "%s.area", Temp); 
  f = fopen(fname, "w");
  fprintf(f, "%d\n", nTriangles);
  for (i=0; i<nTriangles; i++) {
    fprintf(f, "%d %g\n", i+1, areaTriangles[i]);
  }
  fclose(f);
}

void export_data(char* Temp, int nlhs, mxArray *plhs[])
{
  FILE *f;
  char fname[30];
  int nPoints, nTriangles, nEdges, dummy, dummy2, marker;
  double *pointsX = NULL, *pointsY = NULL, rdummy;
  int *labelsP = NULL, *labelsT = NULL, *labelsE = NULL;
  int *t_ = NULL;
  int *e = NULL;
  int i,j;
  const char *fieldname[] = { "vertices", 
			       "triangles", 
			       "edges",
		    	       "lab_vertices", 
			       "lab_triangles", 
			       "lab_edges",
                               "tri_neighbors",
                               "connectivity"};
  mxArray * p;
  
  if (nlhs < 1)
     return;
     
  if (nlhs > 1) 
    mexErrMsgTxt("Too many output arguments.");
   
  plhs[0] = mxCreateStructMatrix(1, 1, sizeof(fieldname)/sizeof(const char *), fieldname);
    
  sprintf(fname, "%s.1.node", Temp);

  f = fopen(fname, "r");
  fscanf(f, "%d %d %d %d", &nPoints, &dummy, &dummy2, &marker);
  pointsX = (double *) malloc(nPoints* sizeof(double));
  pointsY = (double *) malloc(nPoints* sizeof(double));
  if (marker) {
     labelsP = (int *) malloc(nPoints * sizeof(int));
     for (i=0; i<nPoints; i++)
       fscanf(f, "%d %lg %lg %d\n", &j, pointsX + i, pointsY + i, labelsP + i);
  }
  else {
     for (i=0; i<nPoints; i++)
       fscanf(f, "%d %lg %lg\n", &j, pointsX + i, pointsY + i);
  }
 
  p = mxCreateDoubleMatrix(nPoints, 2, mxREAL);
  if (nPoints > 0) {
      double *q = mxGetPr(p);
      for (i=0; i<nPoints; i++)
          *q++ = pointsX[i];
      for (i=0; i<nPoints; i++)
          *q++ = pointsY[i];
      }

  free(pointsX);
  free(pointsY);

  fclose(f);
  mxSetField(plhs[0], 0, fieldname[0], p);

  sprintf(fname, "%s.1.ele", Temp);

  f = fopen(fname, "r");
  fscanf(f, "%d %d %d", &nTriangles, &dummy, &marker);
  t_ = (int *) malloc(3*nTriangles* sizeof(int));
  if (marker) {
     labelsT = (int *) malloc(nTriangles * sizeof(int));
     for (i=0; i<nTriangles; i++) {
         fscanf(f, "%d %d %d %d %lg\n", 
                   &j, t_ + 3*i, t_ + 3*i+1, t_ + 3*i+2, &rdummy);
	 labelsT[i] = rdummy;
	 }
  }
  else
	for (i=0; i<nTriangles; i++)
           fscanf(f, "%d %d %d %d\n", 
                     &j, t_ + 3*i, t_ + 3*i+1, t_ + 3*i+2);
  p = mxCreateDoubleMatrix(nTriangles, 3, mxREAL);
  if (nTriangles > 0) {
      double *q = mxGetPr(p);
      for (i=0; i<nTriangles; i++)
          *q++ = t_[3*i];
      for (i=0; i<nTriangles; i++)
          *q++ = t_[3*i+1];
      for (i=0; i<nTriangles; i++)
          *q++ = t_[3*i+2];
  }


  fclose(f);
 
  mxSetField(plhs[0], 0, fieldname[1], p);

  sprintf(fname, "%s.1.edge", Temp);

  f = fopen(fname, "r");
  fscanf(f, "%d %d", &nEdges, &marker);
  e = (int *) malloc(2*nEdges* sizeof(int));
  if (marker) {
     labelsE = (int *) malloc(nEdges * sizeof(int));
     for (i=0; i<nEdges; i++)
         fscanf(f, "%d %d %d %d\n", 
                   &j, e + 2*i, e + 2*i+1, labelsE + i);
  }
  else
	for (i=0; i<nEdges; i++)
           fscanf(f, "%d %d %d\n", 
                     &j, e + 2*i, e + 2*i+1);

  p = mxCreateDoubleMatrix(nEdges, 2, mxREAL);
  if (nEdges > 0) {
      double *q = mxGetPr(p);
      for (i=0; i<nEdges; i++)
          *q++ = e[2*i];
      for (i=0; i<nEdges; i++)
          *q++ = e[2*i+1];
  }

  fclose(f);
  mxSetField(plhs[0], 0, fieldname[2], p);

  p = mxCreateDoubleMatrix(nPoints, 1, mxREAL);
  if (nPoints > 0) {
      double *q = mxGetPr(p);
      for (i=0; i<nPoints; i++)
          *q++ = labelsP ? labelsP[i] : 0;
      } 
  mxSetField(plhs[0], 0, fieldname[3], p);
  if (labelsP) free(labelsP);

  p = mxCreateDoubleMatrix(nTriangles, 1, mxREAL);
  if (nTriangles > 0) {
      double *q = mxGetPr(p);
      for (i=0; i<nTriangles; i++)
          *q++ = labelsT ? labelsT[i] : 0;
      } 
  mxSetField(plhs[0], 0, fieldname[4], p);
  if (labelsT) free(labelsT);

  p = mxCreateDoubleMatrix(nEdges, 1, mxREAL);
  if (nEdges > 0) {
      double *q = mxGetPr(p);
      for (i=0; i<nEdges; i++)
          *q++ = labelsE ? labelsE[i] : 0;
      } 
  mxSetField(plhs[0], 0, fieldname[5], p);
  if (labelsE) free(labelsE);
  
  sprintf(fname, "%s", Temp);
  unlink(fname);
  sprintf(fname, "%s.node", Temp);
  unlink(fname);
  sprintf(fname, "%s.ele", Temp);
  unlink(fname);
  sprintf(fname, "%s.area", Temp);
  unlink(fname);
  sprintf(fname, "%s.1.node", Temp);
  unlink(fname);
  sprintf(fname, "%s.1.ele", Temp);
  unlink(fname);
  sprintf(fname, "%s.1.edge", Temp);
  unlink(fname);

  {
    const char *fname[] = {
      "NodeElements",
      "NodeEdges",
      "EdgeNodes",
      "EdgeElements",
      "ElementNodes",
      "ElementEdges"
    };
    
    mxArray *p = mxCreateStructMatrix(1, 1, sizeof(fname)/sizeof(const char *), fname);
    mxArray *q, *r;
    int nP = nPoints;
    int nT = nTriangles;
    int nE = nEdges;
    
    double *d_r;
    int i, j, n;
    int *t_i;
    
    struct triangulateio i_out;
    struct triangulateio *out = &i_out;
    
    i_out.numberofpoints = nPoints;
    i_out.numberoftriangles = nTriangles;
    i_out.numberofedges = nEdges;
    i_out.trianglelist = t_;
    i_out.edgelist = e;
    
    compute_connectivities(out);
    
    q = mxCreateCellArray(1, &nP);
    
    for (i=0; i<nP; i++) {
      n = out->cNT.nlist[i];
      r = mxCreateDoubleMatrix(1, n, mxREAL);
      t_i = out->cNT.list[i];
      d_r = mxGetPr(r);
      for (j=0; j<n; j++) *d_r++ = *t_i++;
      mxSetCell(q, i, r);
    }
    
    mxSetField(p, 0, fname[0], q);
    
    q = mxCreateCellArray(1, &nP);
    
    for (i=0; i<nP; i++) {
      n = out->cNE.nlist[i];
      r = mxCreateDoubleMatrix(1, n, mxREAL);
      t_i = out->cNE.list[i];
      d_r = mxGetPr(r);
      for (j=0; j<n; j++) *d_r++ = *t_i++;
      mxSetCell(q, i, r);
    }
    
    mxSetField(p, 0, fname[1], q);
    
    q = mxCreateDoubleMatrix(nE, 2, mxREAL);
    d_r = mxGetPr(q);
    t_i = out->cEN.list;
    for (i=0; i<nE; i++, t_i+=2) *d_r++ = *t_i;
    t_i = out->cEN.list + 1;
    for (i=0; i<nE; i++, t_i+=2) *d_r++ = *t_i;
    mxSetField(p, 0, fname[2], q);
    
    q = mxCreateDoubleMatrix(nE, 2, mxREAL);
    d_r = mxGetPr(q);
    t_i = out->cET.list;
    for (i=0; i<nE; i++, t_i+=2) *d_r++ = *t_i;
    t_i = out->cEN.list + 1;
    for (i=0; i<nE; i++, t_i+=2) *d_r++ = *t_i;
    mxSetField(p, 0, fname[3], q);
    
    q = mxCreateDoubleMatrix(nT, 3, mxREAL);
    d_r = mxGetPr(q);
    t_i = out->cTN.list;
    for (i=0; i<nT; i++, t_i+=3) *d_r++ = *t_i;
    t_i = out->cTN.list + 1;
    for (i=0; i<nT; i++, t_i+=3) *d_r++ = *t_i;
    t_i = out->cTN.list + 2;
    for (i=0; i<nT; i++, t_i+=3) *d_r++ = *t_i;
    mxSetField(p, 0, fname[4], q);
    
    q = mxCreateDoubleMatrix(nT, 3, mxREAL);
    d_r = mxGetPr(q);
    t_i = out->cTE.list;
    for (i=0; i<nT; i++, t_i+=3) *d_r++ = *t_i;
    t_i = out->cTE.list + 1;
    for (i=0; i<nT; i++, t_i+=3) *d_r++ = *t_i;
    t_i = out->cTE.list + 2;
    for (i=0; i<nT; i++, t_i+=3) *d_r++ = *t_i;
    mxSetField(p, 0, fname[5], q);
    
    mxSetField(plhs[0], 0, fieldname[7], p);
  }
  free(e);
  free(t_);

}

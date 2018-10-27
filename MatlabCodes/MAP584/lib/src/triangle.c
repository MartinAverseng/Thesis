#define REAL double

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "triangle.h"
#include "tristruct.h"
#include "util.h"

#include "mex.h"

void import_data(int nrhs, const mxArray *prhs[],
		 int *nPoints, double **points, double **labelPoints,
		 int *nSegments, double **segments, double **labelSegments,
		 double *area,
		 int *nHoles, double **holes,
		 int *nRegions, double **regions, 
		 double **labelRegions, double **areaRegions,
		 char options[]);

void export_data(int nlhs, mxArray *plhs[],
                 struct triangulateio *out);
                
mxArray * export_connectivities(struct triangulateio *out);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  struct triangulateio *in, *out;
  
  int nPoints = 0, nSegments = 0;
  double *points = NULL, *segments = NULL;
  double *labelPoints = NULL, *labelSegments = NULL;
  double area = 0.0;
  double *regions = NULL, *holes = NULL, 
    *labelRegions = NULL, *areaRegions = NULL; 
  int nHoles = 0, nRegions = 0;
  char options[50];
  /*char format[50];*/

  debug_begin("triangle_debug.txt");

  import_data(nrhs, prhs, 
              &nPoints, &points, &labelPoints,
              &nSegments, &segments, &labelSegments,
              &area, 
              &nHoles, &holes, 
              &nRegions, &regions, &labelRegions, &areaRegions,
	      options);

  debug_printf("mHoles = %d Holes %p\n", nHoles, holes);
  in = create_triangulation();
  out = create_triangulation();
  
  init_triangulation(in,
		     nPoints, points, labelPoints,
		     nSegments, segments, labelSegments,
                     nHoles, holes, 
                     nRegions, regions, labelRegions, areaRegions);


  report_triangulation("in", in);
  
  triangulate(options, in, out, NULL);
  
  export_data(nlhs, plhs, out);
  report_triangulation("out", out);  

  free_triangulation(in, 1);
  free_triangulation(out, 0);

  debug_end();
}

void import_data(int nrhs, const mxArray *prhs[],
		 int *nPoints, double **points, double **labelPoints,
		 int *nSegments, double **segments, double **labelSegments,
		 double *area,
		 int *nHoles, double **holes,
		 int *nRegions, double **regions, 
		 double **labelRegions, double **areaRegions,
		 char options[50])
{
  const mxArray *geom, *p, *p2, *p3;
  char *input_buf, str[50];
  int buflen, status;

  debug_printf("nrhs = %d\n\n", nrhs);

  if ((nrhs < 1) || (nrhs > 2))
    mexErrMsgTxt("usage mesh = triangle(geom [, area])");
 
  geom = prhs[0];
  if (!mxIsStruct(geom))
    mexErrMsgTxt("first parameter is not a structure");

  p = mxGetField(geom, 0, "points");
  if (!p) 
    mexErrMsgTxt("no point array in geometry structure");

  if (mxGetN(p) != 2)
    mexErrMsgTxt("points must be in R^2");
  *nPoints = mxGetM(p);
  *points = mxGetPr(p);
  
  p = mxGetField(geom, 0, "segments");
  if (!p) 
    mexErrMsgTxt("no segment array in geometry structure");
  if (mxGetN(p) != 2)
    mexErrMsgTxt("segments must have "
                 "two and only two endpoints");
  *nSegments = mxGetM(p);
  *segments = mxGetPr(p);
    
  p = mxGetField(geom, 0, "lab_points");
  if ((p) && (mxGetM(p) > 0)) {
    if ((mxGetN(p) != 1) && (mxGetM(p) != *nPoints))
      mexErrMsgTxt("point labels must be a "
                   "vector whose size is the number of points");
    *labelPoints = mxGetPr(p);
    }

  p = mxGetField(geom, 0, "lab_segments");
  if ((p) && (mxGetM(p) > 0)) {
    if ((mxGetN(p) != 1) && (mxGetM(p) != *nSegments))
      mexErrMsgTxt("segment labels must be a "
                   "vector whose size is the number of points");
    *labelSegments = mxGetPr(p);
    }

  p = mxGetField(geom, 0, "holes");
  if (p && (mxGetM(p) > 0)) {
    if (mxGetN(p) != 2)
      mexErrMsgTxt("holes must be a 2 by n matrix (one (x,y) pair by hole)");
    *nHoles = mxGetM(p);
    *holes = mxGetPr(p);
    }

  p = mxGetField(geom, 0, "regions");
  p2 = mxGetField(geom, 0, "lab_regions");
  p3 = mxGetField(geom, 0, "area_regions");

  if (p && (mxGetM(p) > 0)) {
    if (mxGetM(p) != 2)
      mexErrMsgTxt("regions must be a 2 by n array, where n "
                   "is the number of regions (one point inside each region)");

    *nRegions = mxGetM(p);
    *regions = mxGetPr(p);

    if (p2 && (mxGetN(p2) == 1) && (mxGetM(p2) == *nRegions))
      *labelRegions = mxGetPr(p2);
    if (p3 && (mxGetN(p3) == 1) && (mxGetM(p3) == *nRegions)) {
      *areaRegions = mxGetPr(p3);
    }
  }

  input_buf = NULL;
  p = mxGetField(geom, 0, "options");
  if (p && (mxGetM(p) > 0)) {
    if(mxIsChar(p) != 1)
      mexErrMsgTxt("options must be a string.");
    if(mxGetM(p) != 1)
      mexErrMsgTxt("options must be a row vector.");
    
    buflen = (mxGetM(p) * mxGetN(p)) + 1;
    /* Allocate memory for input and output strings. */
    input_buf = mxCalloc(buflen, sizeof(char));
    /* Copy the string data from prhs[0] into a C string 
     * input_buf. If the string array contains several rows, 
     * they are copied, one column at a time, into one long 
     * string array. */
    status = mxGetString(p, input_buf, buflen);
    if(status != 0) 
      mexWarnMsgTxt("Not enough space. String is truncated.");
  }
  
  if (nrhs == 2) {
    p = prhs[1];
    if ((mxGetN(p) != 1) || (mxGetM(p) != 1))
      mexErrMsgTxt("area must be a positive scalar");
    *area = *mxGetPr(p);
  }

  if (*area > 0.0 && (*areaRegions == NULL))
    {if (input_buf != NULL){
      sprintf(str, "ApeinQa%.10f", *area);
      strcpy(options,input_buf);
      strcat(options,str);
    }
    else
      sprintf(options, "ApeinQa%.10f", *area);}
  else{
    if (input_buf != NULL){
      strcpy(options,"pAenQ");
      strcat(options,input_buf);}
    else
      strcpy(options, "pAenQ");}
  /*mexWarnMsgTxt(options);*/
}

void export_data(int nlhs, mxArray *plhs[],
                 struct triangulateio *out)
{
  double *q, *d_q;
  int i, *i_q;
  const char *fname[] = { "vertices", 
			  "triangles", 
			  "edges",
			  "lab_vertices", 
			  "lab_triangles", 
			  "lab_edges",
			  "tri_neighbors",
			  "connectivity"
  };

  if (nlhs < 1) return;
  if (!out) return;

  if (nlhs > 1) 
    mexErrMsgTxt("Too many output arguments.");


  plhs[0] = mxCreateStructMatrix(1, 1, sizeof(fname)/sizeof(const char *), fname);

  if (out->numberofpoints > 0) {
    mxArray *p = mxCreateDoubleMatrix(out->numberofpoints, 2, mxREAL);
    q = mxGetPr(p);
    
    d_q = out->pointlist;
    for (i=0; i<out->numberofpoints; i++) {
      *q++ = *d_q++; d_q++;
    }
    d_q = out->pointlist+1;
    for (i=0; i<out->numberofpoints; i++) {
      *q++ = *d_q++; d_q++;
    }
    mxSetField(plhs[0], 0, fname[0], p);
  }


  if (out->numberoftriangles > 0) {

    mxArray *p = mxCreateDoubleMatrix(out->numberoftriangles, 3, mxREAL);
    q = mxGetPr(p);
    
    i_q = out->trianglelist;
    for (i=0; i<out->numberoftriangles; i++) {
      *q++ = *i_q++; i_q += 2;
    }
    i_q = out->trianglelist+1;
    for (i=0; i<out->numberoftriangles; i++) {
      *q++ = *i_q++; i_q += 2;
    }
    i_q = out->trianglelist+2;
    for (i=0; i<out->numberoftriangles; i++) {
      *q++ = *i_q++; i_q += 2;
    }
    mxSetField(plhs[0], 0, fname[1], p);
  }
  
  if (out->numberofedges > 0) {
     
    mxArray *p = mxCreateDoubleMatrix(out->numberofedges, 2, mxREAL);
    q = mxGetPr(p);
    
    i_q = out->edgelist;
    for (i=0; i<out->numberofedges; i++) {
      *q++ = *i_q++; i_q += 1;
    }
    i_q = out->edgelist+1;
    for (i=0; i<out->numberofedges; i++) {
      *q++ = *i_q++; i_q += 1;
    }
    mxSetField(plhs[0], 0, fname[2], p);
  }

  if (out->numberofpoints > 0) {
    
    mxArray *p = mxCreateDoubleMatrix(out->numberofpoints, 1, mxREAL);
    q = mxGetPr(p);

    i_q = out->pointmarkerlist;
    for (i=0; i<out->numberofpoints; i++)
      *q++ = *i_q++;
    mxSetField(plhs[0], 0, fname[3], p);
  }

  if (out->numberoftriangles > 0) {

    mxArray *p = mxCreateDoubleMatrix(out->numberoftriangles, 1, mxREAL);
    q = mxGetPr(p);

    d_q = out->triangleattributelist;
      
    for (i=0; i<out->numberoftriangles; i++)
      *q++ = d_q ? *d_q++ : 0;

    mxSetField(plhs[0], 0, fname[4], p);
  }

  if (out->numberofedges > 0) {

    mxArray *p = mxCreateDoubleMatrix(out->numberofedges, 1, mxREAL);

    q = mxGetPr(p);
    i_q = out->edgemarkerlist;
    for (i=0; i<out->numberofedges; i++)
      *q++ = *i_q++;

    mxSetField(plhs[0], 0, fname[5], p);
  }

  if (out->neighborlist) {
    mxArray *p = mxCreateDoubleMatrix(out->numberoftriangles, 3, mxREAL);

    q = mxGetPr(p);
    i_q = out->neighborlist;
    for (i=0; i<out->numberoftriangles; i++) {
      *q++ = *i_q++; i_q += 2;
    }
    i_q = out->neighborlist+1;
    for (i=0; i<out->numberoftriangles; i++) {
      *q++ = *i_q++; i_q += 2;
    }
    i_q = out->neighborlist+2;
    for (i=0; i<out->numberoftriangles; i++) {
      *q++ = *i_q++; i_q += 2;
    } 
    mxSetField(plhs[0], 0, fname[6], p);
  }

  mxSetField(plhs[0], 0, fname[7], export_connectivities(out));
}

mxArray * export_connectivities(struct triangulateio *out)
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
  int nP = out->numberofpoints;
  int nT = out->numberoftriangles;
  int nE = out->numberofedges;

  double *d_r;
  int i, j, n;
  int *p_i;

  compute_connectivities(out);

  q = mxCreateCellArray(1, &nP);

  for (i=0; i<nP; i++) {
    n = out->cNT.nlist[i];
    r = mxCreateDoubleMatrix(1, n, mxREAL);
    p_i = out->cNT.list[i];
    d_r = mxGetPr(r);
    for (j=0; j<n; j++) *d_r++ = *p_i++;
    mxSetCell(q, i, r);
  }

  mxSetField(p, 0, fname[0], q);

  q = mxCreateCellArray(1, &nP);

  for (i=0; i<nP; i++) {
    n = out->cNE.nlist[i];
    r = mxCreateDoubleMatrix(1, n, mxREAL);
    p_i = out->cNE.list[i];
    d_r = mxGetPr(r);
    for (j=0; j<n; j++) *d_r++ = *p_i++;
    mxSetCell(q, i, r);
  }

  mxSetField(p, 0, fname[1], q);

  q = mxCreateDoubleMatrix(nE, 2, mxREAL);
  d_r = mxGetPr(q);
  p_i = out->cEN.list;
  for (i=0; i<nE; i++, p_i+=2) *d_r++ = *p_i;
  p_i = out->cEN.list + 1;
  for (i=0; i<nE; i++, p_i+=2) *d_r++ = *p_i;
  mxSetField(p, 0, fname[2], q);

  q = mxCreateDoubleMatrix(nE, 2, mxREAL);
  d_r = mxGetPr(q);
  p_i = out->cET.list;
  for (i=0; i<nE; i++, p_i+=2) *d_r++ = *p_i;
  p_i = out->cET.list + 1;
  for (i=0; i<nE; i++, p_i+=2) *d_r++ = *p_i;
  mxSetField(p, 0, fname[3], q);

  q = mxCreateDoubleMatrix(nT, 3, mxREAL);
  d_r = mxGetPr(q);
  p_i = out->cTN.list;
  for (i=0; i<nT; i++, p_i+=3) *d_r++ = *p_i;
  p_i = out->cTN.list + 1;
  for (i=0; i<nT; i++, p_i+=3) *d_r++ = *p_i;
  p_i = out->cTN.list + 2;
  for (i=0; i<nT; i++, p_i+=3) *d_r++ = *p_i;
  mxSetField(p, 0, fname[4], q);

  q = mxCreateDoubleMatrix(nT, 3, mxREAL);
  d_r = mxGetPr(q);
  p_i = out->cTE.list;
  for (i=0; i<nT; i++, p_i+=3) *d_r++ = *p_i;
  p_i = out->cTE.list + 1;
  for (i=0; i<nT; i++, p_i+=3) *d_r++ = *p_i;
  p_i = out->cTE.list + 2;
  for (i=0; i<nT; i++, p_i+=3) *d_r++ = *p_i;
  mxSetField(p, 0, fname[5], q);
 
  return p;
}

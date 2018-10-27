#include "tristruct.h"
#include "util.h"
#include <stdlib.h>
#include <stdio.h>

struct triangulateio * create_triangulation()
{
  struct triangulateio *tr = (struct triangulateio *)
    malloc(sizeof(struct triangulateio));

  tr->pointlist = NULL;
  tr->pointattributelist = NULL;
  tr->pointmarkerlist = NULL;
  tr->numberofpoints = 0;
  tr->numberofpointattributes = 0;

  tr->trianglelist = NULL;
  tr->triangleattributelist = NULL;
  tr->trianglearealist = NULL;
  tr->neighborlist = NULL;
  tr->numberoftriangles = 0;
  tr->numberofcorners = 0;
  tr->numberoftriangleattributes = 0;
  
  tr->segmentlist = NULL;
  tr->segmentmarkerlist = NULL;
  tr->numberofsegments = 0;
  
  tr->holelist = NULL;
  tr->numberofholes = 0;
  
  tr->regionlist = NULL;
  tr->numberofregions = 0;
  
  tr->edgelist = NULL;
  tr->edgemarkerlist = NULL;
  tr->normlist = NULL;
  tr->numberofedges = 0;

  tr->cTN.list = NULL;
  tr->cTE.list = NULL;

  tr->cET.list = NULL;
  tr->cEN.list = NULL;

  tr->cNE.list = NULL;
  tr->cNE.nlist = NULL;
  tr->cNT.list = NULL;
  tr->cNT.nlist = NULL;

  return tr;
}

void free_triangulation(struct triangulateio *tr, int owner)
{
  if (tr->pointlist != NULL) free(tr->pointlist);
  if (tr->pointattributelist != NULL) free(tr->pointattributelist);
  if (tr->pointmarkerlist != NULL) free(tr->pointmarkerlist);
  if (tr->trianglelist != NULL) free(tr->trianglelist);
  if (tr->triangleattributelist != NULL) free(tr->triangleattributelist);
  if (tr->trianglearealist != NULL) free(tr->trianglearealist);
  if (tr->neighborlist != NULL) free(tr->neighborlist);  
  if (tr->segmentlist != NULL) free(tr->segmentlist);
  if (tr->segmentmarkerlist != NULL) free(tr->segmentmarkerlist);
  if (owner && (tr->holelist != NULL)) free(tr->holelist);
  if (owner && (tr->regionlist != NULL)) free(tr->regionlist);
  if (tr->edgelist != NULL) free(tr->edgelist);
  if (tr->edgemarkerlist != NULL) free(tr->edgemarkerlist);
  if (tr->normlist != NULL) free(tr->normlist);
  if (tr->cNT.list != NULL) {
     free(tr->cNT.list[0]);
     free(tr->cNT.list);
     free(tr->cNT.nlist);
  }
  if (tr->cNE.list != NULL) {
     free(tr->cNE.list[0]);
     free(tr->cNE.list);
     free(tr->cNE.nlist);
  }
  if (tr->cET.list != NULL) free(tr->cET.list);
  if (tr->cTE.list != NULL) free(tr->cTE.list);
  free(tr);
}

void init_triangulation(struct triangulateio *tr, 
			int np, double *p, double *lp,
			int ns, double *s, double *ls,
                        int nh, double *h,
                        int nr, double *r, double *lr, double *ar)
{
  int i;

  tr->numberofpoints = np;
  tr->numberofsegments = ns;

  tr->pointlist = (double *) malloc(tr->numberofpoints * 2 * sizeof(double));
  for (i=0; i<tr->numberofpoints; i++)
    tr->pointlist[2*i] = *p++;
  for (i=0; i<tr->numberofpoints; i++)
    tr->pointlist[2*i+1] = *p++;

  tr->pointmarkerlist = (int *) malloc(tr->numberofpoints * sizeof(int));
  if (lp)
     for (i=0; i<tr->numberofpoints; i++)
         tr->pointmarkerlist[i] = (int) lp[i];
  else
     for (i=0; i<tr->numberofpoints; i++)
         tr->pointmarkerlist[i] = 0;

  tr->segmentlist = (int *) malloc(tr->numberofsegments * 2 * sizeof(int));
  for (i=0; i<tr->numberofsegments; i++)
    tr->segmentlist[2*i] = *s++;
  for (i=0; i<tr->numberofsegments; i++)
    tr->segmentlist[2*i+1] = *s++;

  tr->segmentmarkerlist = (int *) malloc(tr->numberofsegments * sizeof(int));
  if (ls)
     for (i=0; i<tr->numberofsegments; i++)
         tr->segmentmarkerlist[i] = (int) ls[i];
  else
     for (i=0; i<tr->numberofsegments; i++)
         tr->segmentmarkerlist[i] = 0;

  if (nr > 0) {
    tr->numberofregions = nr;
    tr->regionlist = (double *) malloc(nr * 4 * sizeof(double));
    p = tr->regionlist;
    for (i=0; i<nr; i++, p+=4) {
      p[0] = r[i];
      p[1] = r[i+nr];
      p[2] = lr ? lr[i] : i+1;
      p[3] = ar ? ar[i] : 0;
      debug_printf( "region %d : %g %g %g %g\n", i, p[0], p[1], p[2], p[3]);
    }
  }

  if (nh > 0) {
    tr->numberofholes = nh;
    tr->holelist = (double *) malloc(nh * 2 * sizeof(double));
    p = tr->holelist;
    for (i=0; i<nh; i++, p+=2) {
      p[0] = h[i];
      p[1] = h[i+nh];
      debug_printf("hole %d : %g %g\n", i, p[0], p[1]);
    }
  }
 
}

void  init_refinment      (struct triangulateio *tr,
			   int np, double *p, double *lp,
			   int ns, double *s, double *ls,
			   int nt, double *t, double *lt,
			   double *a)
{
  int i;

  tr->numberofpoints = np;
  tr->numberofsegments = ns;

  tr->pointlist = (double *) malloc(tr->numberofpoints * 2 * sizeof(double));
  for (i=0; i<tr->numberofpoints; i++)
    tr->pointlist[2*i] = *p++;
  for (i=0; i<tr->numberofpoints; i++)
    tr->pointlist[2*i+1] = *p++;

  tr->pointmarkerlist = (int *) malloc(tr->numberofpoints * sizeof(int));
  if (lp)
     for (i=0; i<tr->numberofpoints; i++)
         tr->pointmarkerlist[i] = (int) lp[i];
  else
     for (i=0; i<tr->numberofpoints; i++)
         tr->pointmarkerlist[i] = 0;

  tr->segmentlist = (int *) malloc(tr->numberofsegments * 2 * sizeof(int));
  for (i=0; i<tr->numberofsegments; i++)
    tr->segmentlist[2*i] = *s++;
  for (i=0; i<tr->numberofsegments; i++)
    tr->segmentlist[2*i+1] = *s++;

  tr->segmentmarkerlist = (int *) malloc(tr->numberofsegments * sizeof(int));
  if (ls)
     for (i=0; i<tr->numberofsegments; i++)
         tr->segmentmarkerlist[i] = (int) ls[i];
  else
     for (i=0; i<tr->numberofsegments; i++)
         tr->segmentmarkerlist[i] = 0;

  tr->trianglelist = (int *) malloc(tr->numberoftriangles * 3 * sizeof(int));
  for (i=0; i<tr->numberoftriangles; i++)
    tr->trianglelist[3*i] = *s++;
  for (i=0; i<tr->numberoftriangles; i++)
    tr->trianglelist[3*i+1] = *s++;
  for (i=0; i<tr->numberoftriangles; i++)
    tr->trianglelist[3*i+2] = *s++;

  tr->trianglearealist = (double *) 
    malloc(tr->numberoftriangles * sizeof(double));
  for (i=0; i<tr->numberoftriangles; i++)
    tr->trianglearealist[i] = a[i];

  tr->segmentmarkerlist = (int *) malloc(tr->numberofsegments * sizeof(int));
  if (ls)
     for (i=0; i<tr->numberofsegments; i++)
         tr->segmentmarkerlist[i] = (int) ls[i];
  else
     for (i=0; i<tr->numberofsegments; i++)
         tr->segmentmarkerlist[i] = 0;
}


void export_triangulation(double *p, double *t, struct triangulateio *out)
{
  double *d_p;
  int i, *i_p;

  d_p = out->pointlist;
  for (i=0; i<out->numberofpoints; i++) {
    *p++ = *d_p++; d_p++;
  }
  d_p = out->pointlist+1;
  for (i=0; i<out->numberofpoints; i++) {
    *p++ = *d_p++; d_p++;
  }

  i_p = out->trianglelist;
  for (i=0; i<out->numberoftriangles; i++) {
    *t++ = *i_p++; i_p += 2;
  }
  i_p = out->trianglelist+1;
  for (i=0; i<out->numberoftriangles; i++) {
    *t++ = *i_p++; i_p += 2;
  }
  i_p = out->trianglelist+2;
  for (i=0; i<out->numberoftriangles; i++) {
    *t++ = *i_p++; i_p += 2;
  }

}

void compute_connectivities(io)
struct triangulateio *io;
{
  int i, j, k, l, mE, npk;
  int *cPP, *p, *q, *pE, **pT, *npT, *pk;
  int **pPP;
  int *npPP, *np;

  int nP = io->numberofpoints;
  int nT = io->numberoftriangles;
  int nE = io->numberofedges;

  int *e, *t, t1, e1;

  debug_printf( "nP = %d nT = %d nE = %d\n", nP, nT, nE);

  /* Connectivite Elements -> Noeuds */

  io->cTN.list = io->trianglelist;

  /* Connectivite Aretes -> Noeuds */

  io->cEN.list = io->edgelist;

  /* Connectivite Noeuds -> Elements */

  cPP  = (int *)  malloc (sizeof(int) * 3*nT);
  npPP = (int *)  malloc (sizeof(int) * nP);
  pPP  = (int **) malloc (sizeof(int *) * nP);

  for (i=0; i<nP; i++) {
    npPP[i] = 0;
  }

  t = io->trianglelist;
  for (i=0; i<nT; i++) 
    for (j=0; j<3; j++, t++)
      npPP[*t-1]++;
 
  pPP[0] = cPP;
  for (i=0; i<nP-1; i++) {
    pPP[i+1] = pPP[i] + npPP[i];
    npPP[i] = 0;
  }
  npPP[nP-1] = 0;

  t = io->trianglelist;
  for (i=0; i<nT; i++) {  
    for (j=0; j<3; j++, t++) {
      t1 = *t - 1;
      pPP[t1][npPP[t1]] = i+1;
      npPP[t1]++;
    }
  }

  debug_printf( "Connectivite noeud -> elements\n\n");
  for (i=0; i<nP; i++) {
    debug_printf( "%6d %6d : ", i+1, npPP[i]);
    for (j=0; j<npPP[i]; j++)
      debug_printf( " %6d", pPP[i][j]);
    debug_printf( "\n");
  }
  debug_printf( "\n\n");

  io->cNT.list = pPP;
  io->cNT.nlist = npPP;

  /* Connectivite noeuds -> aretes */

  cPP  = (int *)  malloc (sizeof(int) * 2*nE);
  npPP = (int *)  malloc (sizeof(int) * nP);
  pPP  = (int **) malloc (sizeof(int *) * nP);

  for (i=0; i<nP; i++) {
    npPP[i] = 0;
  }

  e = io->edgelist;
  for (i=0; i<nE; i++) {
    npPP[e[0]-1]++; 
    npPP[e[1]-1]++;
    e+=2;
  }

  pPP[0] = cPP;
  for (i=0; i<nP-1; i++) {
    pPP[i+1] = pPP[i] + npPP[i];
    npPP[i] = 0;
  }
  npPP[nP-1] = 0;

  for (i=0; i<nP; i++) {
    npPP[i] = 0;
  }

  e = io->edgelist;
  for (i=0; i<nE; i++) {  
    for (j=0; j<2; j++, e++) {
      e1 = *e - 1;
      pPP[e1][npPP[e1]] = i+1;
      npPP[e1]++;
    }
  }

  debug_printf( "Connectivite noeud -> aretes\n\n");
  for (i=0; i<nP; i++) {
    debug_printf( "%6d %6d : ", i+1, npPP[i]);
    for (j=0; j<npPP[i]; j++)
      debug_printf( " %6d", pPP[i][j]);
    debug_printf( "\n");
  }
  debug_printf( "\n\n");
   
  io->cNE.list = pPP;
  io->cNE.nlist = npPP;

  /* Connectivite Aretes -> Elements */

  p = (int *) malloc(sizeof(int) * nE * 2);
  q = (int *) malloc(sizeof(int) * nT);

  pE = io->cEN.list;
  pT = io->cNT.list;
  npT = io->cNT.nlist;

  for (i=0; i<2*nE; i++) p[i] = -1;

  for (i=0; i<nE; i++) {
   
    for (j=0; j<nT; j++) q[j] = 0;

    for (j=0; j<2; j++) {
      k = *pE++;
      pk = pT[k-1];
      npk = npT[k-1];
      for (l=0; l<npk; l++) q[pk[l]-1]++;
    }
    mE = 0;
    for (j=0; j<nT; j++)
      if (q[j] == 2) {
	 p[2*i + mE] = j+1;
	 mE++;
      }
  }
  io->cET.list = p;
  free(q);

  debug_printf( "Connectivite arete -> elements\n\n");
  for (i=0; i<nE; i++)
    debug_printf( "%6d : %6d %6d\n", i+1, p[2*i], p[2*i+1]);
  debug_printf( "\n\n");

  /* Connectivite Elements -> Aretes */

  p = (int *) malloc(sizeof(int) * nT * 3);
  np = (int *) malloc(sizeof(int) * nT);
  q = io->cET.list;
 
  for (i=0; i<nT; i++) np[i] = 0;

  for (j=0; j<2; j++)
     for (i=0; i<nE; i++) {
       k = q[2*i+j]-1;
       if (k >= 0) {
	  p[np[k] + 3*k] = i+1;
	  np[k]++;
       }
     }

  free(np);
  io->cTE.list = p;

  debug_printf( "Connectivite element -> aretes (avant renumerotation)\n\n");
  for (i=0; i<nT; i++)
    debug_printf( "%6d : %6d %6d %6d\n", i+1, p[3*i], p[3*i+1], p[3*i+2]);
  debug_printf( "\n\n");

  {
  int e1, e2, e3, p1, p2, p3, e1p1, e2p1, e3p1, e1p2, e2p2, e3p2, ee1, ee2, ee3;
  for (i=0; i<nT; i++) {
     e1 = io->cTE.list[3*i];
     e2 = io->cTE.list[3*i+1];
     e3 = io->cTE.list[3*i+2];

     debug_printf("i = %5d e1 = %5d e2 = %5d e3 = %5d", i, e1, e2, e3);
     p1 = io->cTN.list[3*i];
     p2 = io->cTN.list[3*i+1];
     p3 = io->cTN.list[3*i+2];
     debug_printf(" p1 = %5d p2 = %5d p3 = %5d\n", p1, p2, p3);
     e1p1 = io->cEN.list[2*e1-2];   e1p2 = io->cEN.list[2*e1+1-2];
     debug_printf("         e1 = (%5d, %5d)", e1p1, e1p2);
     e2p1 = io->cEN.list[2*e2-2];   e2p2 = io->cEN.list[2*e2+1-2];
     debug_printf(" e2 = (%5d, %5d)", e2p1, e2p2);
     e3p1 = io->cEN.list[2*e3-2];   e3p2 = io->cEN.list[2*e3+1-2];
     debug_printf(" e3 = (%5d, %5d)\n", e3p1, e3p2);

     if ((p1 != e1p1) && (p1 != e1p2)) ee1 = e1;
     if ((p1 != e2p1) && (p1 != e2p2)) ee1 = e2;
     if ((p1 != e3p1) && (p1 != e3p2)) ee1 = e3;

     if ((p2 != e1p1) && (p2 != e1p2)) ee2 = e1;
     if ((p2 != e2p1) && (p2 != e2p2)) ee2 = e2;
     if ((p2 != e3p1) && (p2 != e3p2)) ee2 = e3;

     if ((p3 != e1p1) && (p3 != e1p2)) ee3 = e1;
     if ((p3 != e2p1) && (p3 != e2p2)) ee3 = e2;
     if ((p3 != e3p1) && (p3 != e3p2)) ee3 = e3;

     io->cTE.list[3*i  ] = ee1;
     io->cTE.list[3*i+1] = ee2;
     io->cTE.list[3*i+2] = ee3;
  }
  }

  debug_printf( "Connectivite element -> aretes\n\n");
  for (i=0; i<nT; i++)
    debug_printf( "%6d : %6d %6d %6d\n", i+1, p[3*i], p[3*i+1], p[3*i+2]);
  debug_printf( "\n\n");
    
}

void report_triangulation(name, io)
char *name;
struct triangulateio *io;
{
  int i, j;

  debug_printf("\nTriagulation %s : \n\n", name);
  debug_printf("io->pointlist                  = %d\n", io->pointlist);
  debug_printf("io->pointattributelist         = %p\n", io->pointattributelist);
  debug_printf("io->pointmarkerlist            = %p\n", io->pointmarkerlist);
  debug_printf("io->numberofpoints             = %d\n", io->numberofpoints);
  debug_printf("io->numberofpointattributes    = %d\n", io->numberofpointattributes);

  debug_printf("io->trianglelist               = %p\n", io->trianglelist);
  debug_printf("io->triangleattributelist      = %p\n", io->triangleattributelist);
  debug_printf("io->trianglearealist           = %p\n", io->trianglearealist);
  debug_printf("io->neighborlist               = %d\n", io->neighborlist);
  debug_printf("io->numberoftriangles          = %d\n", io->numberoftriangles);
  debug_printf("io->numberofcorners            = %d\n", io->numberofcorners);
  debug_printf("io->numberoftriangleattributes = %d\n", io->numberoftriangleattributes);
  
  debug_printf("io->segmentlist                = %p\n", io->segmentlist);
  debug_printf("io->segmentmarkerlist          = %p\n", io->segmentmarkerlist);
  debug_printf("io->numberofsegments           = %d\n", io->numberofsegments);
  
  debug_printf("io->holelist                   = %p\n", io->holelist);
  debug_printf("io->numberofholes              = %d\n", io->numberofholes);
  
  debug_printf("io->regionlist                 = %p\n", io->regionlist);
  debug_printf("io->numberofregions            = %d\n", io->numberofregions);
  
  debug_printf("io->edgelist                   = %p\n", io->edgelist);
  debug_printf("io->edgemarkerlist             = %p\n", io->edgemarkerlist);
  debug_printf("io->normlist                   = %p\n", io->normlist);
  debug_printf("io->numberofedges              = %d\n\n", io->numberofedges);

  if (io->pointlist) {

    debug_printf("\nPoints\n\n");

    for (i=0; i < io->numberofpoints; i++) {

       debug_printf("%4d  x,y :", i+1);

       for (j = 0; j < 2; j++)
           debug_printf("  %7.4g", io->pointlist[i * 2 + j]);

       if (io->numberofpointattributes > 0) {
          debug_printf("   attributes :");
 
          for (j = 0; j < io->numberofpointattributes; j++)
              debug_printf("  %10.6g",
		io->pointattributelist[i * io->numberofpointattributes + j]);
          }
       
       if (io->pointmarkerlist)
          debug_printf("   marker : %d", io->pointmarkerlist[i]);

       if (io->cNT.list) {
	  debug_printf("   triangles :");
          for (j=0; j<io->cNT.nlist[i]; j++)
	     debug_printf(" %d", io->cNT.list[i][j]);
          }

       if (io->cNE.list) {
	  debug_printf("   edges :");
          for (j=0; j<io->cNE.nlist[i]; j++)
	     debug_printf(" %d", io->cNE.list[i][j]);
          }
       debug_printf("\n");
    }
  }
  debug_printf("\n");

  if (io->trianglelist) {
     
    debug_printf("\nTriangles\n\n");
    for (i = 0; i < io->numberoftriangles; i++) {

      debug_printf("%4d  points :", i+1);

      for (j = 0; j < io->numberofcorners; j++)
        debug_printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);

      if (io->numberoftriangleattributes > 0) {
        debug_printf("   attributes :");
        for (j = 0; j < io->numberoftriangleattributes; j++)
          debug_printf("  %.6g", io->triangleattributelist
	                           [i * io->numberoftriangleattributes + j]);
      }
 
      if (io->neighborlist) {
        debug_printf("  neighbors :");
        for (j = 0; j < 3; j++)
	  debug_printf("  %4d", io->neighborlist[i * 3 + j]);
      }

      if (io->cTE.list) {
         debug_printf("  edges :");
         for (j = 0; j < 3; j++)
           debug_printf("  %4d", io->cTE.list[i * 3 + j]);
      }
      debug_printf("\n");
    }
  }

  if (io->segmentlist) {
   
    debug_printf("\nSegments\n\n");
    for (i = 0; i < io->numberofsegments; i++) {

      debug_printf("%4d  points :", i+1);
      for (j = 0; j < 2; j++)
        debug_printf(" %4d", io->segmentlist[i * 2 + j]);
    
      if (io->segmentmarkerlist)
        debug_printf("   marker : %d", io->segmentmarkerlist[i]);

      debug_printf("\n");
    }
  }

  if (io->edgelist) {
     
    debug_printf("\nEdges\n\n");
    for (i = 0; i < io->numberofedges; i++) {
  
      debug_printf("%4d  points :", i+1);
      for (j = 0; j < 2; j++)
        debug_printf("  %4d", io->edgelist[i * 2 + j]);
      if (io->normlist) {
        debug_printf("  normal :");
        for (j = 0; j < 2; j++)
	  debug_printf("  %.6g", io->normlist[i * 2 + j]);
      }
      if (io->edgemarkerlist)
        debug_printf("   marker %d", io->edgemarkerlist[i]);

      if (io->cET.list) {
        debug_printf("  elements :");
        for (j = 0; j < 2; j++)
          debug_printf("  %4d", io->cET.list[i * 2 + j]);
      }
      debug_printf("\n");
    }
  }
}

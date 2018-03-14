#ifndef _TRISTRUCT
#define _TRISTRUCT


#define REAL double
#include "triangle.h"

struct triangulateio * create_triangulation();

void free_triangulation   (struct triangulateio *tr, int owner);

void init_triangulation   (struct triangulateio *tr, 
			   int np, double *p, double *lp,
			   int ns, double *s, double *ls,
                           int nh, double *h, 
                           int nr, double *r, double *lr, double *ar);

void  init_refinment      (struct triangulateio *in,
			   int nP, double *points, double *labelPoints,
			   int nS, double *segments, double *labelSegments,
			   int nT, double *triangles, double *labelTriangles,
			   double *areaTriangles);

void compute_connectivities(struct triangulateio *in);

void report_triangulation(char*name, struct triangulateio *io);

#endif

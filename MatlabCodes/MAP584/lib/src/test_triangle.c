#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "tristruct.h"
#include "util.h"

void export_data(struct triangulateio *);
void print_array(const char *s, int n, int m, double *v);

int main()
{

  int nPoints = 12;
  double points[24] = {
    0.0,
    1.0,
    1.0,
    0.0,
    0.2,
    0.4,
    0.4,
    0.2,
    0.6,
    0.8,
    0.8,
    0.6,
    0.0,
    0.0,
    1.0,
    1.0,
    0.2,
    0.2,
    0.4,
    0.4,
    0.6,
    0.6,
    0.8,
    0.8
  };

  double * labelPoints = NULL;

  int nSegments = 12;
  double segments[24] = {
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    2,
    3,
    4,
    1,
    6,
    7,
    8,
    5,
    10,
    11,
    12,
    9
  };

  double * labelSegments = NULL;

  int nHoles = 2;
  double holes[4] = {
    0.7,
    0.3,
    0.7, 
    0.3
  };
    
  int nRegions = 0;
  double *regions = NULL;
  double *labelRegions = NULL;
  double *areaRegions = NULL;

  double area = 0.000001;

  struct triangulateio *in, *out;

  char format[40];

  in = create_triangulation();
  out = create_triangulation();
 
  init_triangulation(in,
                     nPoints, points, labelPoints,
                     nSegments, segments, labelSegments,
                     nHoles, holes,
                     nRegions, regions, labelRegions, areaRegions);
 
  print_array("Points", nPoints, 2, points);
  print_array("Segments", nSegments, 2, segments);

  if (area > 0.0 && (areaRegions == NULL))
    sprintf(format, "Apea%.15f", area);
  else
    strcpy(format, "pAe");
 
  report_triangulation("in", in);
 
  triangulate(format, in, out, NULL);
 
  report_triangulation("out", out);

  export_data(out);    
  return 0;
}

void export_data(struct triangulateio *out)
{
   FILE *f;
   int i;

   f = fopen("triangles.out", "w");
   for (i=0; i<out->numberoftriangles; i++)
     fprintf(f, "%d %d %d\n",
	     out->trianglelist[3*i],                                      
	     out->trianglelist[3*i+1],
	     out->trianglelist[3*i+2]);
   fclose(f);
   f = fopen("points.out", "w");
   for (i=0; i<out->numberofpoints; i++)
     fprintf(f, "%g %g\n",
             out->pointlist[2*i],
             out->pointlist[2*i+1]);
   fclose(f);
}

void print_array(const char *s, int n, int m, double *v)
{
  int i,j;
  printf("\n%s : \n\n", s);
  for (i=0; i<n; i++) {
    printf("%5d : ", i+1);
    for (j=0; j<m; j++)
      printf("%10.4g", *v++);
    printf("\n");
  }
  printf("\n");
}

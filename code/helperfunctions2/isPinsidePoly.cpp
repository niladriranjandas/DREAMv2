#include <iostream>
#include <math.h>
#include <mex.h>

using namespace std;

#define INF 181; // since we are using phi and psi setting extremum point as 181 suffices

#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;
#endif

struct Point
{
    double phi;
    double psi;
};


//-------------------------------------------------------------------------
// The code is taken from: http://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
//-------------------------------------------------------------------------

// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(Point p, Point q, Point r)
{
    if (q.phi <= max(p.phi, r.phi) && q.phi >= min(p.phi, r.phi) &&
            q.psi <= max(p.psi, r.psi) && q.psi >= min(p.psi, r.psi))
        return true;
    return false;
}
 
// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(Point p, Point q, Point r)
{
    int val = (q.psi - p.psi) * (r.phi - q.phi) -
              (q.phi - p.phi) * (r.psi - q.psi);
 
    if (val == 0) return 0;  // colinear
    return (val > 0)? 1: 2; // clock or counterclock wise
}
 
// The function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
bool doIntersect(Point p1, Point q1, Point p2, Point q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);
 
    // General case
    if (o1 != o2 && o3 != o4)
        return true;
 
    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;
 
    // p1, q1 and p2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;
 
    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;
 
     // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;
 
    return false; // Doesn't fall in any of the above cases
}
 
// Returns true if the point p lies inside the polygon[] with n vertices
bool isInside(Point polygon[], int n, Point p)
{
    // There must be at least 3 vertices in polygon[]
    if (n < 3)  return false;
 
    // Create a point for line segment from p to infinite
    //Point extreme = {INF, p.psi};
    Point extreme;
    extreme.phi = INF;
    extreme.psi = p.psi;
 
    // Count intersections of the above line with sides of polygon
    int count = 0, i = 0;
    do
    {
        int next = (i+1)%n;
 
        // Check if the line segment from 'p' to 'extreme' intersects
        // with the line segment from 'polygon[i]' to 'polygon[next]'
        if (doIntersect(polygon[i], polygon[next], p, extreme))
        {
            // If the point 'p' is colinear with line segment 'i-next',
            // then check if it lies on segment. If it lies, return true,
            // otherwise false
            if (orientation(polygon[i], p, polygon[next]) == 0)
               return onSegment(polygon[i], p, polygon[next]);
 
            count++;
        }
        i = next;
    } while (i != 0);
 
    // Return true if count is odd, false otherwise
    return count&1;  // Same as (count%2 == 1)
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
  //declare variables
    mxArray *a_in_m, *b_in_m, *c_out_m;
    const mwSize *dims;
    double *a, *b, *c, *d;
    int dimx, dimy, numdims;
    int i,j; 
    
    //Point region;
    
  //associate inputs
    a_in_m = mxDuplicateArray(prhs[0]);
    b_in_m = mxDuplicateArray(prhs[1]);

  //figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dimy = (int)dims[0]; dimx = (int)dims[1];

  //associate outputs
    c_out_m = plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
  //associate pointers
    a = mxGetPr(a_in_m);
    b = mxGetPr(b_in_m);
    c = mxGetPr(c_out_m);
   
//call the inside/outside polygon algorithm
   Point region[dimx],querry_pt;
   
// build the polygon out of the region data   
    for(i=0;i<dimx;i++)
    {
        for(j=0;j<dimy;j++)
        {
           // mexPrintf("element[%d][%d] = %f\n",j,i,a[i*dimy+j]);
            (j==0)? region[i].phi = a[i*dimy+j] : region[i].psi = a[i*dimy+j];
        }
    }
   
  // mexPrintf("\n (%f,%f) ",*(b),*(b+1));
   querry_pt.phi = *(b);
   querry_pt.psi = *(b+1);
   
   //for(i=0;i<dimx;i++)
   //    mexPrintf("\n %d phi: %f psi: %f",i,region[i].phi,region[i].psi);
   
   //isInside(region, dimx, querry_pt)? cout << "\n Yes \n": cout << "\n No \n";
   isInside(region, dimx, querry_pt)? *c=1: *c=0;
   
}
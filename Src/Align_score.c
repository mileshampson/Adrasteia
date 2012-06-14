// Code derived from Deacon Sweeney's. Contact sweeney.2@wright.edu

#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define MAX2(x,y)     ((x)<(y) ? (y) : (x))
#define MAX3(x,y,z)   (MAX2(x,y)<(z) ? (z) : MAX2(x,y))

typedef struct{
  char p[80];          // array of residues
  int PosCnt;           // length of chain
} Fragment ;

double gaussian(double x)
{
  return(erfc(x/sqrt(2)));
} 

double simple_align(Fragment target_fragment, Fragment template_fragment){

  int N = 0;
  if (target_fragment.PosCnt > template_fragment.PosCnt){
    N = target_fragment.PosCnt+2;
  } else {
    N = template_fragment.PosCnt+2;
  }

  // variable declarations
  double diag, down, right;
  double max, Max;
  int simple_score, i, j, tempi, tempj, x, y, topskip, bottomskip;
  int xMax, yMax;   
  char aout[N],bout[N];
  double h[N][N];
  int xTraceback[N][N], yTraceback[N][N];

  // initialize traceback array
  Max=xMax=yMax=0;
  for (i=0;i<=target_fragment.PosCnt;i++)
    for (j=0;j<=template_fragment.PosCnt;j++) {
      xTraceback[i][j]=-1;
      yTraceback[i][j]=-1;
    }

  // compute "h" local similarity array
  for (i=0;i<=target_fragment.PosCnt;i++) 
    h[i][0]=0.0;
  for (j=0;j<=template_fragment.PosCnt;j++) 
    h[0][j]=0.0;

  for (i=1;i<=target_fragment.PosCnt;i++)
    for (j=1;j<=template_fragment.PosCnt;j++) {
      if (target_fragment.p[i-1] == template_fragment.p[j-1]){
        simple_score = 1;
      } else {
        simple_score = -1;
      }
      diag    = h[i-1][j-1] + simple_score;
      down    = h[i-1][j] -2; 
      right   = h[i][j-1] -2; 
      max=MAX3(diag,down,right);
      if (max <= 0)  {
        h[i][j]=0.0;
        xTraceback[i][j]=-1;
        yTraceback[i][j]=-1;
      } else if (max == diag) {
        h[i][j]=diag;
        xTraceback[i][j]=i-1;
        yTraceback[i][j]=j-1;
      } else if (max == down) {
        h[i][j]=down;
        xTraceback[i][j]=i-1;
        yTraceback[i][j]=j;
      } else  {
        h[i][j]=right;
        xTraceback[i][j]=i;
        yTraceback[i][j]=j-1;
      }
      if (max > Max){
        Max=max;
        xMax=i;
        yMax=j;
      }
    }

  i = xMax; j = yMax;
  while ( i>0 && j>0 && h[i][j]>0){
    tempi=i;
    tempj=j;
    i=xTraceback[tempi][tempj];
    j=yTraceback[tempi][tempj];
  }

  // reset to max point to do alignment
  i=xMax; j=yMax;
  x=y=0;
  topskip = bottomskip = 1;
  while (i>0 && j>0 && h[i][j] > 0){
    if (topskip && bottomskip) {
      aout[x++]=target_fragment.p[i-1];
      bout[y++]=template_fragment.p[j-1];
    } else if (topskip) {
      aout[x++]='-';
      bout[y++]=template_fragment.p[j-1];
    } else if (bottomskip) {
      aout[x++]=target_fragment.p[i-1];
      bout[y++]='-';
    }
    topskip    = (j>yTraceback[i][j]);
    bottomskip = (i>xTraceback[i][j]);
    tempi=i;tempj=j;
    i=xTraceback[tempi][tempj];
    j=yTraceback[tempi][tempj];
  }

  return Max;
}


int score(int len, char *seg1, char *seg2)
{
  int i;
  int length = len;
  // load the two sequences into fragments
  Fragment target_fragment, template_fragment;
  
  target_fragment.PosCnt = length;
  for (i=0; i<length; i++)
    target_fragment.p[i] = seg1[i];

  template_fragment.PosCnt = length;
  for (i=0; i<length; i++)
    template_fragment.p[i] = seg2[i];

  // align the fragments
  return simple_align(target_fragment, template_fragment);
}


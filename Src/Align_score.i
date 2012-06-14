%module Align_score
%{
extern double gaussian(double x);
extern int score(int len, char *seg1, char *seg2);
%}
extern double gaussian(double x);
extern int score(int len, char *seg1, char *seg2);

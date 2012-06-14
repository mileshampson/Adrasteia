#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAXPATLEN 256
#define MAXSYM  256
/* MAXSYM * 32 (sagrep)  */
#define MEMBER_TABLE_SIZE 8192  
/* intial number of matches provided in list */
#define INITPAIRS 10	
#define INCRPAIRS 1.5
/* 2 ^ 16, i.e. combinations of 2 chars */
#define HASHTABLESIZE 65536  
#define TWOBYTES 16
/* The length of pattern above which lagrep is used instead of sagrep. 
lagrep does work for smaller pattern, so useful for testing */
//lagrep does not work on 64 bit platforms so if we are running
//on one of these use sagrep for all queries
#if defined(__amd64__) || defined(x86_64)
#define SHORT_LONG  256  
#else
#define SHORT_LONG  24  
#endif

#ifndef TRUE
#define FALSE 0
#define TRUE 1
#endif

/* typedef unsigned char boolean;  Cannot be use because not
				   currently supported by SWIG */
typedef int boolean;

/* The bundled parameters for sagrep aka agrep */
typedef struct sagrep_struct
  {
  unsigned char  SHIFT[MAXSYM ];
  unsigned char  MEMBER[MEMBER_TABLE_SIZE ];
  unsigned int Mask[MAXSYM];
  unsigned int endposition;
  int shift_1;
  signed char NErrors;
  } sagrep_struct;

/* The bundled parameters for lagrep aka a_monkey */
typedef struct lagrep_struct
  {
  unsigned int Hashmask;
  char MEMBER_1[HASHTABLESIZE];
  signed char NErrors;
  } lagrep_struct;

typedef union param_struct
{
  sagrep_struct sagrep;
  lagrep_struct lagrep;
} param_struct;


typedef struct int_pair
{
  int start, end;
} int_pair;

typedef struct int_pair_list
{
  int npairs, maxpairs;
  int_pair *pairs;
} int_pair_list;

extern int_pair_list * add_ends(int start, int end, int_pair_list *matches);
extern void printnstring(char *string, int start, int end);
extern int_pair_list *agrepy(char *pat, int patlen, char *text, int textlen,
        int gotoends, param_struct *parampt);
extern param_struct *compile(char* Pattern, int patlen, int NErrors);
extern param_struct *sagrepy_compile(char* Pattern, int patlen, signed char NErrors);
extern param_struct *lagrepy_compile(char* Pattern, int patlen, signed char NErrors);
extern int_pair_list *exec_sagrepy(char *pat, int patlen, char *text, int textlen,
                boolean gotoends, param_struct *parampt);
extern int_pair_list *exec_lagrepy(char *pat, int patlen, char *text, int textlen,
                boolean gotoends, param_struct *parampt);


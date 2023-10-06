#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <inttypes.h>
#include <math.h>
#include "sfmtrng.h"

#ifndef __LRSGLIB_H_INC__
#define __LRSGLIB_H_INC__

/* utils.c */
extern uint32_t softplus_u32(int32_t x);

extern bool strsme(const char *str1, const char *str2);
extern void strrmac(char *__strdst, const char *__strsrc, const char __chrctr);
extern void strrmuc(char *__strdst, const char *__strsrc, const char __chrctr);
extern bool strbws(char *__strsrc, const char *__strbegin);
extern size_t strtozu(const char *s);
extern uint32_t strtou32(const char *s);

extern bool __feexist(const char *fn);
extern bool __fnexist(const char *fn);
extern void __fopen(FILE **f, const char *fn, const char *md);
extern void __fread_check(size_t __frdval, size_t __frdcnt);
extern void __popen(FILE **p, const char *fn, const char *md);
extern void __fgets_row(FILE **fc, char *row);


void __challoc(void *__ptr);
void *__chMalloc(size_t n);
void *__chCalloc(size_t __nmemb, size_t __size);

extern pid_t call(char* argv[]);
extern void __wait_childs(void);

char *rand_string(char *str, size_t size);

typedef double *double_p;


extern void __make_adj_from_tmp(size_t i, size_t j, double tmp, double_p **adj);
extern void __fill_adj__(FILE **f, size_t N, double_p **adj);


#endif /* __LRSGLIB_H_INC__ */
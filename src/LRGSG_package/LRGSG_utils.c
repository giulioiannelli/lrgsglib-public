#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>
#include <inttypes.h>
#include "LRGSG_utils.h"
#include "LRGSG_customs.h"
#include "sfmtrng.h"
/**
 * @brief get the ReLU(x) (or softplus) of input unsigned 32-bit integer.
 * 
 * @param x (uint32_t) an integer
 * @return (int32_t) max(0, x).
 */
extern uint32_t softplus_u32(int32_t x)
{
    if(x > 0)
    {
        return x;
    }
    return 0;
}

/* STRING RELATED FUNCTIONS ************************************************* */
//
/**
 * @brief check if two strings are the same
 * 
 * @param str1 (const char *) the first string to compare
 * @param str2 (const char *) the second string to compare
 * @return true if the two strings are the same
 * @return false if the two strings differ for some characters
 */
extern bool strsme(const char *str1, const char *str2)
{
    if (strcmp(str1, str2))
        return false;
    else
        return true;
}
/**
 * @brief copy a number of characters from source string to destination string
 * 
 * @param nc (size_t) number of characters to copy
 * @param outs (char *) the destination string
 * @param psrc (const char *) the source string
 */
extern void strncpy_safe(size_t nc, char *outs, const char* psrc)
{
    strncpy(outs, psrc, nc);
    outs[nc - 1] = '\0';
}
/**
 * @brief create a copy of a string after a character is met in string
 * 
 * @param __strdst (char *) the destination string
 * @param __strsrc (const char *) the source string
 * @param __chrctr (const char) the caracter to lookup
 */
extern void strrmac(char *__strdst, const char *__strsrc, const char __chrctr)
{
    char *__rmstr = strrchr(__strsrc, __chrctr);
    strcpy(__strdst, __rmstr);
}
/**
 * @brief create a copy of a string until a character is met in string
 * 
 * @param __strdst (char *) the destination string
 * @param __strsrc (const char *) the source string
 * @param __chrctr (const char) the caracter to lookup
 */
extern void strrmuc(char *__strdst, const char *__strsrc, const char __chrctr)
{
    char *__rmstr = strrchr(__strsrc, __chrctr);
    strcpy(__strdst, __strsrc);
    __strdst[strlen(__strsrc)-strlen(__rmstr)] = '\0';
}
/**
 * @brief chack whether a string begin with another 
 * 
 * @param __strsrc (char *) the string where to look
 * @param __strbegin (const char *) the string to be contained
 * @return (bool) true if __strsrc begins with __strbegin
 * @return (bool) false otherwise
 */
extern bool strbws(char *__strsrc, const char *__strbegin)
{
    return (bool) (!(strncmp(__strsrc, __strbegin, strlen(__strbegin))));
}
/**
 * @brief get a string from file and check it is nor empty nor there are errors
 * with files
 * 
 * @param fc (FILE **) file from which to read
 * @param row (char *) the char pointer to store the row content
 */
extern void __fgets_row(FILE **fc, char *row)
{
    if ((fgets(row, STRL1024, *fc) == NULL))
    {
        perror(MSG_ERR_FGETS);
        exit(EXIT_FAILURE);
    }
    row[strlen(row) - 1] = '\0';
}

/**
 * @brief acquire size_t data from a string
 * 
 * @param s (const char *) the input string
 * @return (size_t) the sscanned number
 */
extern size_t strtozu(const char *s)
{
    char c;
    int scanned;
    size_t i;
    scanned = sscanf(s, "%zu%c", &i, &c);
    if (scanned == 1)
        return i;
    else if (scanned > 1)
    {
        perror("strtozu");
        return i;
    }
    if (c != '\0' || errno != 0)
    {
        perror("strtozu");
        exit(EXIT_FAILURE);
    }
    return 0;
}

/**
 * @brief acquire uint32_t data from a string
 * 
 * @param s (const char *) the input string
 * @return (uint32_t) the scanned number
 */
extern uint32_t strtou32(const char *s)
{
    char *endptr;
    uintmax_t tmp;
    tmp = strtoumax(s, &endptr, 10);
    if (*endptr != '\0')
    {
        fprintf(stderr, MSG_WRN_SCNU32 "%s\n", endptr);
        return (uint32_t) tmp;
    }
    if (tmp < UINT32_MAX)
    {
        return (uint32_t) tmp;
    }
    else
    {
        fprintf(stderr, MSG_ERR_SCNU32);
        exit(EXIT_FAILURE);
    }
    // char c;
    // int scanned;
    // uint32_t i;
    // scanned = sscanf(s, "%" SCNu32 "%c", &i, &c);
    // printf("%" PRIu32 "%c\n", i, c);
    // if (scanned == 1)
    //     return i;
    // else if (scanned > 1)
    // {
    //     perror("strtou32");
    //     return i;
    // }
    // if (c != '\0' || errno != 0)
    // {
    //     perror("strtou32");
    //     printf("%" PRIu32 "%c", i, c);
    //     exit(EXIT_FAILURE);
    // }
    // return 0;
}
/* FILES I/O FUNCTIONS ****************************************************** */
//
/**
 * @brief check that input file pointer points to an existing file
 * 
 * @param n (const char *) string containing path to file
 * @return true if file exists
 * @return false if file does not exist
 */
extern bool __feexist(const char *fn)
{
    FILE *f;
    if ((f = fopen(fn, "r")))
        fclose(f);
    else
        return false;
    return true;
}
/**
 * @brief check that input file pointer points to a non existing file
 * 
 * @param n (const char *) string containing path to file
 * @return true if file does not exist
 * @return false if file exists
 */
extern bool __fnexist(const char *fn)
{
    return (!(__feexist(fn)));
}
/**
 * @brief open a file according to an operative mode allowed by fopen
 * 
 * @param f (FILE **) FILE pointer
 * @param fn (const char *) file name string
 * @param md (const char *) opening mode
 */
extern void __fopen(FILE **f, const char *fn, const char *md)
{
    if ((*f = fopen(fn, md)) == NULL)
    {
        perror(fn);
        exit(EXIT_FAILURE);
    }
}
/**
 * @brief check an fread worked out correctly
 * 
 * @param __frdval (size_t) the number of elements read off by fread
 * @param __frdcnt (size_t) the requested number of read
 */
extern void __fread_check(size_t __frdval, size_t __frdcnt)
{
     if (__frdval != __frdcnt)
     {
        perror(MSG_ERR_FREAD);
        exit(EXIT_FAILURE);
    }
}
/**
 * @brief open a pipe according to an operative mode allowed by fopen
 * 
 * @param p (FILE **) pipe pointer
 * @param fn (const char *) file name string
 * @param md (const char *) opening mode
 */
extern void __popen(FILE **p, const char *fn, const char *md)
{
    if ((*p = popen(fn, md)) == NULL)
    {
        perror(MSG_ERR_POPEN);
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief check if allocation made with malloc, calloc and realloc has worked,
 * othewise print on stderr the errno
 * 
 * @param __ptr (void *) allocated pointer
 */
void __challoc(void *__ptr)
{
    if (__ptr == NULL)
    {
        perror("Alloc fail");
        exit(EXIT_FAILURE);
    }
}
/**
 * @brief (m)allocate array and check sucessfull allocation, return the
 * allocated pointer
 * 
 * @param n (size_t) the number of bytes to allocate
 * @return (void*) allocated pointer
 */
void *__chMalloc(size_t n)
{
    void *p = malloc(n);
    __challoc(p);
    return p;
}
/**
 * @brief (c)allocate array and check sucessfull allocation, return the
 * allocated pointer
 * 
 * @param n (size_t) the number of bytes to allocate
 * @return (void*) allocated pointer
 */
void *__chCalloc(size_t __nmemb, size_t __size)
{
    void *p = calloc(__nmemb, __size);
    __challoc(p);
    return p;
}
/**
 * @brief clone the current process, then, in child process, execve the program 
 * file to replace the child with the desired program.
 * 
 * @param argv 
 * @return pid (pid_t) 
 */
extern pid_t call(char* argv[]) {
    pid_t pid = fork();
    if (pid == 0)
    {
        char* envp[] = { NULL };
        
        execve(argv[0], argv, envp);
        perror("Error execve");
        exit(EXIT_FAILURE);
    }
    else
    {
        
        return pid;
    }
}
/**
 * @brief wait for all the children of the program to end
 * 
 */
extern void __wait_childs(void)
{
    int crps;
    int status;
    while ((crps = wait(&status)) > 0)
    {
        printf("%d: PID %d exited with status 0x%.4X\n",
               (int)getpid(), crps, status);
    }
}


/**
 * @brief generate random string
 * 
 * @param str 
 * @return size
 */
char *rand_string(char *str, size_t size)
{
    const char charset[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKILMNOPQRSTUVWXYZ";
    if (size) {
        --size;
        for (size_t n = 0; n < size; n++) {
            uint64_t key = RNG_u64() % (sizeof charset - 1);
            str[n] = charset[key];
        }
        str[size] = '\0';
    }
    return str;
}






// /* DICTIONARY IMPLEMENTATION ************************************************ */
// //
// #define HASHSIZE 101
// static struct nlist *hashtab[HASHSIZE]; /* pointer table */

// struct nlist { /* table entry: */
//     struct nlist *next; /* next entry in chain */
//     char *name; /* defined name */
//     char *defn; /* replacement text */
// };


// /* hash: form hash value for string s */
// unsigned hash(char *s)
// {
//     unsigned hashval;
//     for (hashval = 0; *s != '\0'; s++)
//       hashval = *s + 31 * hashval;
//     return hashval % HASHSIZE;
// }

// /* lookup: look for s in hashtab */
// struct nlist *lookup(char *s)
// {
//     struct nlist *np;
//     for (np = hashtab[hash(s)]; np != NULL; np = np->next)
//         if (strcmp(s, np->name) == 0)
//           return np; /* found */
//     return NULL; /* not found */
// }

// char *strdup(char *);
// /* install: put (name, defn) in hashtab */
// struct nlist *install(char *name, char *defn)
// {
//     struct nlist *np;
//     unsigned hashval;
//     if ((np = lookup(name)) == NULL) { /* not found */
//         np = (struct nlist *) malloc(sizeof(*np));
//         if (np == NULL || (np->name = strdup(name)) == NULL)
//           return NULL;
//         hashval = hash(name);
//         np->next = hashtab[hashval];
//         hashtab[hashval] = np;
//     } else /* already there */
//         free((void *) np->defn); /*free previous defn */
//     if ((np->defn = strdup(defn)) == NULL)
//        return NULL;
//     return np;
// }

// char *strdup(char *s) /* make a duplicate of s */
// {
//     char *p;
//     p = (char *) malloc(strlen(s)+1); /* +1 for ’\0’ */
//     if (p != NULL)
//        strcpy(p, s);
//     return p;
// }
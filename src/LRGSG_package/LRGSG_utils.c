#include "LRGSG_utils.h"
//

extern void print_stdout_cwd(void)
{
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) == NULL)
        perror("getcwd() error");
    else
        printf("Current working directory: %s\n", cwd);
}

void flip_spin(size_t nd, spin_tp s) {
    *(s + nd) = -*(s + nd);
}
double calc_ext_magn(size_t N, spin_tp s) {
    double m = 0.;
    for (size_t i = 0; i < N; i++)
        m += *(s + i);
    return m;
}
double calc_magn(size_t N, spin_tp s) {
    return (calc_ext_magn(N, s) / N);
}
double calc_ext_magn2(size_t N, spin_tp s) {
    double m2 = 0.;
    for (size_t i = 0; i < N; i++)
        m2 += *(s + i) * *(s + i);
    return m2;
}
double calc_clust_magn(size_t cli_l, size_tp cli, spin_tp s) {
    double clm = 0.;
    for (size_t i = 0; i < cli_l; i++)
        clm += s[cli[i]];
    return clm / cli_l;
}
double neigh_weight_magn(size_t nd, size_t n_nn, spin_tp s, size_tp *neighs,
                         double_p *edgl) {
    double sum = 0.;
    for (size_t i = 0; i < n_nn; i++)
        sum += *(*(edgl + nd) + i) * *(s + *(*(neighs + nd) + i));
    return sum / n_nn;
}
double calc_energy_full(size_t N, spin_tp s, size_tp nlen, size_tp *neighs,
                        double_p *edgl) {
    double sum = 0., tmp = 0.;
    for (size_t i = 0; i < N; i++) {
        tmp = *(s + i) * neigh_weight_magn(i, *(nlen + i), s, neighs, edgl);
        sum += tmp;
    }
    return -sum / N;
}

/** perform the sum of a floating point array 
 * @param n (size_t) the number of vector components
 * @param v (double *) the floaring point array
 * @return (double) sum(v)
 */
extern double sum_vs(size_t n, double *v)
{
    double s = 0.;
    for (size_t i = 0; i < n; i++)
        s += *(v + i);
    return s;
}

/** perform the sum of squared components of an array
 * @param n (size_t) the number of vector components
 * @param v (double *) the array
 * @return (double) sum(v)
 */
extern double sum_vs_2(size_t n, double *v)
{
    double s = 0.;
    for (size_t i = 0; i < n; i++)
        s += v[i] * v[i];
    return s;
}

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
 * @return true if the two strings are the same false if the two strings differ 
 * for some characters
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
/* Prepends t into s. Assumes s has enough space allocated
** for the combined string.
*/
void prepend(char *s, const char *t)
{
    size_t len = strlen(t);
    memmove(s + len, s, strlen(s) + 1);
    memcpy(s, t, len);
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


extern void __make_adj_from_tmp(size_t i, size_t j, double tmp, double_p **adj)
{
    *(*(*adj + j) + i) = tmp;
    *(*(*adj + i) + j) = *(*(*adj + j) + i);
}

extern void __fill_adj__(FILE **f, size_t N, double_p **adj)
{
    double tmp;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = i; j < N; j++)
        {
            __fread_check(fread(&tmp, sizeof(tmp), 1, *f), 1);
            __make_adj_from_tmp(i, j, tmp, adj);
        }
    }
}

extern void __fill_edgl_read__(FILE **f, size_t N, double_p **edgl, size_tp **neighs, size_tp *neigh_len)
{
    size_t node_i, cntr, last = 0;
    double w_ij;
    for (size_t i = 0; i < N; i++)
    {
        *(*edgl + i) = __chMalloc(1 * sizeof(**edgl));
        *(*neighs + i) = __chMalloc(1 * sizeof(**neighs));
    }
    cntr = 0;
    node_i = 0;
    for(size_t i, j; fscanf(*f, "%zu %zu %lf", &i, &j, &w_ij) != EOF;)
    {
        last = i;
        if (i != node_i)
        {
            *(*neigh_len + i-1) = cntr;
            node_i++;
            // printf("%zu, %zu, %zu\n", i-1, *(*neigh_len + i-1), cntr);
            cntr = 0;
        }
        *(*(*neighs + i) + cntr) = j;
        *(*(*edgl + i) + cntr) = w_ij;
        *(*edgl + i) = realloc(*(*edgl + i), (++cntr + 1) * sizeof(**edgl));
        *(*neighs + i) = realloc(*(*neighs + i), (cntr + 1) * sizeof(**neighs));
    }
    *(*neigh_len + last) = cntr;
}

extern void __fill_edgl_make__(FILE **f, size_t N, double_p **adj, double_p **edgl, size_tp **neighs, size_tp *neigh_len)
{
    size_t cntr;
    for (size_t i = 0; i < N; i++)
    {
        cntr = 0;
        *(*edgl + i) = __chMalloc(1 * sizeof(**edgl));
        *(*neighs + i) = __chMalloc(1 * sizeof(**neighs));
        for (size_t j = 0; j < N; j++)
        {
            if (fabs(*(*(*adj + j) + i)) > 0.)
            {
                *(*(*neighs + i) + cntr) = j;
                *(*(*edgl + i) + cntr) = *(*(*adj + j) + i);
                *(*edgl + i) = realloc(*(*edgl + i), (++cntr + 1) * sizeof(**edgl));
                *(*neighs + i) = realloc(*(*neighs + i), (cntr + 1) * sizeof(**neighs));
                fprintf(*f, "%zu %zu %lf\n", i, j, *(*(*adj + j) + i));
            }
        }
        *(*neigh_len + i) = cntr;
        // printf("%zu, %zu, %zu\n", i, *(*neigh_len + i), cntr);
    }
}

// extern void __make_adj_fromfile(size_t i, size_t j, double tmp, double_p **adj)
// {
//     *(*(*adj + i) + j) = (double) tmp;
//     *(*(*adj + j) + i) = *(*(*adj + i) + j);
// }


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
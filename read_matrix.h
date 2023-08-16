#pragma once

#include "matrix_csr.h"
#include "parse.h"

#include <limits>
#include <unistd.h>

#ifdef __cplusplus
#    define restrict __restrict
#endif

#ifdef HAVE_LIBZ
#    include <zlib.h>
#endif

#ifndef IDXTYPEWIDTH
typedef int idx_t;
#    define PRIdx       "d"
#    define IDX_T_MIN   INT_MIN
#    define IDX_T_MAX   INT_MAX
#    define parse_idx_t parse_int
#elif IDXTYPEWIDTH == 32
typedef int32_t idx_t;
#    define PRIdx       PRId32
#    define IDX_T_MIN   INT32_MIN
#    define IDX_T_MAX   INT32_MAX
#    define parse_idx_t parse_int32_t
#elif IDXTYPEWIDTH == 64
typedef int64_t idx_t;
#    define PRIdx       PRId64
#    define IDX_T_MIN   INT64_MIN
#    define IDX_T_MAX   INT64_MAX
#    define parse_idx_t parse_int64_t
#endif

enum partition {
    partition_rows,
    partition_nonzeros,
};

enum mtxobject {
    mtxmatrix,
    mtxvector,
};

enum mtxformat {
    mtxarray,
    mtxcoordinate,
};

enum mtxfield {
    mtxreal,
    mtxinteger,
    mtxpattern,
};

enum mtxsymmetry {
    mtxgeneral,
    mtxsymmetric,
};

void stream_close(enum streamtype streamtype, union stream s)
{
    if (streamtype == stream_stdio) {
        fclose(s.f);
#ifdef HAVE_LIBZ
    } else if (streamtype == stream_zlib) {
        gzclose(s.gzf);
#endif
    }
}

static int mtxfile_fread_header(enum mtxobject   *object,
                                enum mtxformat   *format,
                                enum mtxfield    *field,
                                enum mtxsymmetry *symmetry,
                                idx_t            *num_rows,
                                idx_t            *num_columns,
                                int64_t          *num_nonzeros,
                                enum streamtype   streamtype,
                                union stream      stream,
                                int64_t          *lines_read,
                                int64_t          *bytes_read)
{
    int   line_max = sysconf(_SC_LINE_MAX);
    char *linebuf  = (char *)malloc(line_max + 1);
    if (!linebuf)
        return errno;

    /* read and parse header line */
    int err = freadline(linebuf, line_max, streamtype, stream);
    if (err) {
        free(linebuf);
        return err;
    }
    char *s = linebuf;
    char *t = s;
    if (strncmp("%%MatrixMarket ", t, strlen("%%MatrixMarket ")) == 0) {
        t += strlen("%%MatrixMarket ");
    } else {
        free(linebuf);
        return EINVAL;
    }
    if (bytes_read)
        *bytes_read += t - s;
    s = t;
    if (strncmp("matrix ", t, strlen("matrix ")) == 0) {
        t += strlen("matrix ");
        *object = mtxmatrix;
    } else if (strncmp("vector ", t, strlen("vector ")) == 0) {
        t += strlen("vector ");
        *object = mtxvector;
    } else {
        free(linebuf);
        return EINVAL;
    }
    if (bytes_read)
        *bytes_read += t - s;
    s = t;
    if (strncmp("array ", t, strlen("array ")) == 0) {
        t += strlen("array ");
        *format = mtxarray;
    } else if (strncmp("coordinate ", t, strlen("coordinate ")) == 0) {
        t += strlen("coordinate ");
        *format = mtxcoordinate;
    } else {
        free(linebuf);
        return EINVAL;
    }
    if (bytes_read)
        *bytes_read += t - s;
    s = t;
    if (strncmp("real ", t, strlen("real ")) == 0) {
        t += strlen("real ");
        *field = mtxreal;
    } else if (strncmp("integer ", t, strlen("integer ")) == 0) {
        t += strlen("integer ");
        *field = mtxinteger;
    } else if (strncmp("pattern ", t, strlen("pattern ")) == 0) {
        t += strlen("pattern ");
        *field = mtxpattern;
    } else {
        free(linebuf);
        return EINVAL;
    }
    if (bytes_read)
        *bytes_read += t - s;
    s = t;
    if (strncmp("general", t, strlen("general")) == 0) {
        t += strlen("general");
        *symmetry = mtxgeneral;
    } else if (strncmp("symmetric", t, strlen("symmetric")) == 0) {
        t += strlen("symmetric");
        *symmetry = mtxsymmetric;
    } else {
        free(linebuf);
        return EINVAL;
    }
    if (bytes_read)
        *bytes_read += t - s;
    s = t;

    /* skip lines starting with '%' */
    do {
        if (lines_read)
            (*lines_read)++;
        err = freadline(linebuf, line_max, streamtype, stream);
        if (err) {
            free(linebuf);
            return err;
        }
        s = t = linebuf;
    } while (linebuf[0] == '%');

    /* parse size line */
    if (*object == mtxmatrix && *format == mtxcoordinate) {
        err = parse_idx_t(num_rows, s, &t, bytes_read);
        if (err) {
            free(linebuf);
            return err;
        }
        if (s == t || *t != ' ') {
            free(linebuf);
            return EINVAL;
        }
        if (bytes_read)
            (*bytes_read)++;
        s   = t + 1;
        err = parse_idx_t(num_columns, s, &t, bytes_read);
        if (err) {
            free(linebuf);
            return err;
        }
        if (s == t || *t != ' ') {
            free(linebuf);
            return EINVAL;
        }
        if (bytes_read)
            (*bytes_read)++;
        s   = t + 1;
        err = parse_int64_t(num_nonzeros, s, &t, bytes_read);
        if (err) {
            free(linebuf);
            return err;
        }
        if (s == t) {
            free(linebuf);
            return EINVAL;
        }
        if (lines_read)
            (*lines_read)++;
    } else if (*object == mtxvector && *format == mtxarray) {
        err = parse_idx_t(num_rows, s, &t, bytes_read);
        if (err) {
            free(linebuf);
            return err;
        }
        if (s == t) {
            free(linebuf);
            return EINVAL;
        }
        if (lines_read)
            (*lines_read)++;
    } else {
        free(linebuf);
        return EINVAL;
    }
    free(linebuf);
    return 0;
}

static int mtxfile_fread_matrix_coordinate(enum mtxfield   field,
                                           idx_t           num_rows,
                                           idx_t           num_columns,
                                           int64_t         num_nonzeros,
                                           idx_t          *rowidx,
                                           idx_t          *colidx,
                                           double         *a,
                                           enum streamtype streamtype,
                                           union stream    stream,
                                           int64_t        *lines_read,
                                           int64_t        *bytes_read)
{
    int   line_max = sysconf(_SC_LINE_MAX);
    char *linebuf  = (char *)malloc(line_max + 1);
    if (!linebuf)
        return errno;
    if (field == mtxreal || field == mtxinteger) {
        for (int64_t i = 0; i < num_nonzeros; i++) {
            int err = freadline(linebuf, line_max, streamtype, stream);
            if (err) {
                free(linebuf);
                return err;
            }
            char *s = linebuf;
            char *t = s;
            err     = parse_idx_t(&rowidx[i], s, &t, bytes_read);
            if (err) {
                free(linebuf);
                return err;
            }
            if (s == t || *t != ' ') {
                free(linebuf);
                return EINVAL;
            }
            if (bytes_read)
                (*bytes_read)++;
            s   = t + 1;
            err = parse_idx_t(&colidx[i], s, &t, bytes_read);
            if (err) {
                free(linebuf);
                return err;
            }
            if (s == t || *t != ' ') {
                free(linebuf);
                return EINVAL;
            }
            if (bytes_read)
                (*bytes_read)++;
            s   = t + 1;
            err = parse_double(&a[i], s, &t, bytes_read);
            if (err) {
                free(linebuf);
                return err;
            }
            if (s == t) {
                free(linebuf);
                return EINVAL;
            }
            if (lines_read)
                (*lines_read)++;
        }
    } else if (field == mtxinteger) {
        for (int64_t i = 0; i < num_nonzeros; i++) {
            int err = freadline(linebuf, line_max, streamtype, stream);
            if (err) {
                free(linebuf);
                return err;
            }
            char *s = linebuf;
            char *t = s;
            err     = parse_idx_t(&rowidx[i], s, &t, bytes_read);
            if (err) {
                free(linebuf);
                return err;
            }
            if (s == t || *t != ' ') {
                free(linebuf);
                return EINVAL;
            }
            if (bytes_read)
                (*bytes_read)++;
            s   = t + 1;
            err = parse_idx_t(&colidx[i], s, &t, bytes_read);
            if (err) {
                free(linebuf);
                return err;
            }
            if (s == t || *t != ' ') {
                free(linebuf);
                return EINVAL;
            }
            if (bytes_read)
                (*bytes_read)++;
            s = t + 1;
            int x;
            err = parse_int(&x, s, &t, bytes_read);
            if (err) {
                free(linebuf);
                return err;
            }
            if (s == t) {
                free(linebuf);
                return EINVAL;
            }
            a[i] = x;
            if (lines_read)
                (*lines_read)++;
        }
    } else if (field == mtxpattern) {
        for (int64_t i = 0; i < num_nonzeros; i++) {
            int err = freadline(linebuf, line_max, streamtype, stream);
            if (err) {
                free(linebuf);
                return err;
            }
            char *s = linebuf;
            char *t = s;
            err     = parse_idx_t(&rowidx[i], s, &t, bytes_read);
            if (err) {
                free(linebuf);
                return err;
            }
            if (s == t || *t != ' ') {
                free(linebuf);
                return EINVAL;
            }
            if (bytes_read)
                (*bytes_read)++;
            s   = t + 1;
            err = parse_idx_t(&colidx[i], s, &t, bytes_read);
            if (err) {
                free(linebuf);
                return err;
            }
            if (s == t) {
                free(linebuf);
                return EINVAL;
            }
            a[i] = 1;
            if (lines_read)
                (*lines_read)++;
        }
    } else {
        free(linebuf);
        return EINVAL;
    }
    free(linebuf);
    return 0;
}

static int csr_from_coo_size(enum mtxsymmetry symmetry,
                             idx_t            num_rows,
                             idx_t            num_columns,
                             int64_t          num_nonzeros,
                             const idx_t *__restrict rowidx,
                             const idx_t *__restrict colidx,
                             const double *__restrict a,
                             int64_t *__restrict rowptr,
                             int64_t *__restrict csrsize,
                             idx_t *__restrict rowsizemin,
                             idx_t *__restrict rowsizemax,
                             idx_t *__restrict diagsize,
                             bool           separate_diagonal,
                             enum partition partition)
{
#ifdef _OPENMP
#    pragma omp parallel for
#endif
    for (idx_t i = 0; i < num_rows; i++)
        rowptr[i] = 0;
    rowptr[num_rows] = 0;
    if (num_rows == num_columns && symmetry == mtxsymmetric && separate_diagonal) {
        for (int64_t k = 0; k < num_nonzeros; k++) {
            if (rowidx[k] != colidx[k]) {
                rowptr[rowidx[k]]++;
                rowptr[colidx[k]]++;
            }
        }
    } else if (num_rows == num_columns && symmetry == mtxsymmetric && !separate_diagonal) {
        for (int64_t k = 0; k < num_nonzeros; k++) {
            if (rowidx[k] != colidx[k]) {
                rowptr[rowidx[k]]++;
                rowptr[colidx[k]]++;
            } else {
                rowptr[rowidx[k]]++;
            }
        }
    } else if (num_rows == num_columns && separate_diagonal) {
        for (int64_t k = 0; k < num_nonzeros; k++) {
            if (rowidx[k] != colidx[k])
                rowptr[rowidx[k]]++;
        }
    } else {
        for (int64_t k = 0; k < num_nonzeros; k++)
            rowptr[rowidx[k]]++;
    }
    idx_t rowmin = num_rows > 0 ? rowptr[1] : 0;
    idx_t rowmax = 0;
    for (idx_t i = 1; i <= num_rows; i++) {
        rowmin = rowmin <= rowptr[i] ? rowmin : rowptr[i];
        rowmax = rowmax >= rowptr[i] ? rowmax : rowptr[i];
        rowptr[i] += rowptr[i - 1];
    }
    if (num_rows == num_columns && separate_diagonal) {
        rowmin++;
        rowmax++;
    }
    *rowsizemin = rowmin;
    *rowsizemax = rowmax;
    *csrsize    = rowptr[num_rows];
    *diagsize   = (num_rows == num_columns && separate_diagonal) ? num_rows : 0;
    return 0;
}

static int csr_from_coo(enum mtxsymmetry symmetry,
                        idx_t            num_rows,
                        idx_t            num_columns,
                        int64_t          num_nonzeros,
                        const idx_t *__restrict rowidx,
                        const idx_t *__restrict colidx,
                        const double *__restrict a,
                        int64_t *__restrict rowptr,
                        int64_t csrsize,
                        idx_t   rowsizemin,
                        idx_t   rowsizemax,
                        idx_t *__restrict csrcolidx,
                        double *__restrict csra,
                        double *__restrict csrad,
                        bool           separate_diagonal,
                        bool           sort_rows,
                        enum partition partition)
{
    if (num_rows == num_columns && symmetry == mtxsymmetric && separate_diagonal) {
        for (int64_t k = 0; k < num_nonzeros; k++) {
            if (rowidx[k] == colidx[k]) {
                csrad[rowidx[k] - 1] += a[k];
            } else {
                idx_t i = rowidx[k] - 1, j = colidx[k] - 1;
                csrcolidx[rowptr[i]] = j;
                csra[rowptr[i]]      = a[k];
                rowptr[i]++;
                csrcolidx[rowptr[j]] = i;
                csra[rowptr[j]]      = a[k];
                rowptr[j]++;
            }
        }
        for (idx_t i = num_rows; i > 0; i--)
            rowptr[i] = rowptr[i - 1];
        rowptr[0] = 0;
    } else if (num_rows == num_columns && symmetry == mtxsymmetric && !separate_diagonal) {
        for (int64_t k = 0; k < num_nonzeros; k++) {
            idx_t i = rowidx[k] - 1, j = colidx[k] - 1;
            csrcolidx[rowptr[i]] = j;
            csra[rowptr[i]]      = a[k];
            rowptr[i]++;
            if (i != j) {
                csrcolidx[rowptr[j]] = i;
                csra[rowptr[j]]      = a[k];
                rowptr[j]++;
            }
        }
        for (idx_t i = num_rows; i > 0; i--)
            rowptr[i] = rowptr[i - 1];
        rowptr[0] = 0;
    } else if (num_rows == num_columns && separate_diagonal) {
        for (int64_t k = 0; k < num_nonzeros; k++) {
            idx_t i = rowidx[k] - 1, j = colidx[k] - 1;
            if (i == j) {
                csrad[i] += a[k];
            } else {
                csrcolidx[rowptr[i]] = j;
                csra[rowptr[i]]      = a[k];
                rowptr[i]++;
            }
        }
        for (idx_t i = num_rows; i > 0; i--)
            rowptr[i] = rowptr[i - 1];
        rowptr[0] = 0;
    } else {
        /* simpler, serial version: */
        /* for (int64_t k = 0; k < num_nonzeros; k++) { */
        /*     idx_t i = rowidx[k]-1, j = colidx[k]-1; */
        /*     csrcolidx[rowptr[i]] = j; csra[rowptr[i]] = a[k]; rowptr[i]++; */
        /* } */
        /* for (idx_t i = num_rows; i > 0; i--) rowptr[i] = rowptr[i-1]; */
        /* rowptr[0] = 0; */

        int64_t *__restrict perm = (int64_t *)malloc(num_nonzeros * sizeof(int64_t));
        if (!perm) {
            return errno;
        }
#ifdef _OPENMP
#    pragma omp parallel for
#endif
        for (int64_t k = 0; k < num_nonzeros; k++)
            perm[k] = 0;
        for (int64_t k = 0; k < num_nonzeros; k++) {
            idx_t i           = rowidx[k] - 1;
            perm[rowptr[i]++] = k;
        }
        for (idx_t i = num_rows; i > 0; i--)
            rowptr[i] = rowptr[i - 1];
        rowptr[0] = 0;
#ifdef _OPENMP
#    pragma omp parallel for
#endif
        for (int64_t k = 0; k < num_nonzeros; k++) {
            csrcolidx[k] = colidx[perm[k]] - 1;
            csra[k]      = a[perm[k]];
        }
        free(perm);
    }

#if 0
    /* If requested, sort nonzeros by column within each row */
    if (sort_rows) {
        int err = rowsort(
            num_rows, num_columns,
            rowptr, rowsizemax, csrcolidx, csra);
        if (err) return err;
    }
#endif
    return 0;
}

void read_matrix(matrix_csr<double, int64_t, int> &matrix, const char *matrix_path)
{
    long pagesize = sysconf(_SC_PAGESIZE);

    FILE *f;
    if ((f = fopen(matrix_path, "r")) == NULL) {
        perror("fopen (matrix)");
        exit(EXIT_FAILURE);
    }

    enum mtxobject   object;
    enum mtxformat   format;
    enum mtxfield    field;
    enum mtxsymmetry symmetry;
    idx_t            num_rows;
    idx_t            num_columns;
    int64_t          num_nonzeros;
    int64_t          lines_read = 0;
    int64_t          bytes_read = 0;
    int              err;

    stream     stream_{f};
    streamtype streamtype_ = stream_stdio;

    err = mtxfile_fread_header(&object,
                               &format,
                               &field,
                               &symmetry,
                               &num_rows,
                               &num_columns,
                               &num_nonzeros,
                               streamtype_,
                               stream_,
                               &lines_read,
                               &bytes_read);

    size_t rowidxsize = num_nonzeros * sizeof(idx_t);
    idx_t *rowidx     = (idx_t *)aligned_alloc(pagesize, rowidxsize + pagesize - rowidxsize % pagesize);

    size_t colidxsize = num_nonzeros * sizeof(idx_t);
    idx_t *colidx     = (idx_t *)aligned_alloc(pagesize, colidxsize + pagesize - colidxsize % pagesize);

    size_t  asize = num_nonzeros * sizeof(double);
    double *a     = (double *)aligned_alloc(pagesize, asize + pagesize - asize % pagesize);

    err = mtxfile_fread_matrix_coordinate(
        field, num_rows, num_columns, num_nonzeros, rowidx, colidx, a, streamtype_, stream_, &lines_read, &bytes_read);

    stream_close(streamtype_, stream_);

    size_t   csrrowptrsize = (num_rows + 1) * sizeof(int64_t);
    int64_t *csrrowptr     = (int64_t *)aligned_alloc(pagesize, csrrowptrsize + pagesize - csrrowptrsize % pagesize);

    int64_t csrsize;
    idx_t   rowsizemin, rowsizemax;
    idx_t   diagsize;
    err = csr_from_coo_size(symmetry,
                            num_rows,
                            num_columns,
                            num_nonzeros,
                            rowidx,
                            colidx,
                            a,
                            csrrowptr,
                            &csrsize,
                            &rowsizemin,
                            &rowsizemax,
                            &diagsize,
                            false,
                            partition_rows);

    matrix.nnz  = csrsize;
    matrix.nrow = num_rows;
    matrix.ncol = num_columns;

    size_t csrcolidxsize = csrsize * sizeof(idx_t);
    idx_t *csrcolidx     = (idx_t *)aligned_alloc(pagesize, csrcolidxsize + pagesize - csrcolidxsize % pagesize);

#pragma omp parallel for
    for (idx_t i = 0; i < num_rows; i++) {
        for (int64_t k = csrrowptr[i]; k < csrrowptr[i + 1]; k++)
            csrcolidx[k] = 0;
    }

    size_t  csrasize = csrsize * sizeof(double);
    double *csra     = (double *)aligned_alloc(pagesize, csrasize + pagesize - csrasize % pagesize);

#pragma omp parallel for
    for (idx_t i = 0; i < num_rows; i++) {
        for (int64_t k = csrrowptr[i]; k < csrrowptr[i + 1]; k++)
            csra[k] = 0;
    }

    err = csr_from_coo(symmetry,
                       num_rows,
                       num_columns,
                       num_nonzeros,
                       rowidx,
                       colidx,
                       a,
                       csrrowptr,
                       csrsize,
                       rowsizemin,
                       rowsizemax,
                       csrcolidx,
                       csra,
                       nullptr,
                       false,
                       false,
                       partition_rows);

    matrix.val     = csra;
    matrix.row_ptr = csrrowptr;
    matrix.col_idx = csrcolidx;

    free(a);
    free(rowidx);
    free(colidx);
}
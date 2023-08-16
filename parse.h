// #include "enums.h"

#include <cerrno>
#include <cinttypes>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>

enum streamtype {
    stream_stdio,
#ifdef HAVE_LIBZ
    stream_zlib,
#endif
};

union stream {
    FILE *f;
#ifdef HAVE_LIBZ
    gzFile gzf;
#endif
};

static int
parse_long_long_int(const char *s, char **outendptr, int base, long long int *out_number, int64_t *bytes_read)
{
    errno = 0;
    char         *endptr;
    long long int number = strtoll(s, &endptr, base);
    if ((errno == ERANGE && (number == LLONG_MAX || number == LLONG_MIN)) || (errno != 0 && number == 0))
        return errno;
    if (outendptr)
        *outendptr = endptr;
    if (bytes_read)
        *bytes_read += endptr - s;
    *out_number = number;
    return 0;
}

int parse_int(int *x, const char *s, char **endptr, int64_t *bytes_read)
{
    long long int y;
    int           err = parse_long_long_int(s, endptr, 10, &y, bytes_read);
    if (err)
        return err;
    if (y < INT_MIN || y > INT_MAX)
        return ERANGE;
    *x = y;
    return 0;
}

int parse_int32_t(int32_t *x, const char *s, char **endptr, int64_t *bytes_read)
{
    long long int y;
    int           err = parse_long_long_int(s, endptr, 10, &y, bytes_read);
    if (err)
        return err;
    if (y < INT32_MIN || y > INT32_MAX)
        return ERANGE;
    *x = y;
    return 0;
}

int parse_int64_t(int64_t *x, const char *s, char **endptr, int64_t *bytes_read)
{
    long long int y;
    int           err = parse_long_long_int(s, endptr, 10, &y, bytes_read);
    if (err)
        return err;
    if (y < INT64_MIN || y > INT64_MAX)
        return ERANGE;
    *x = y;
    return 0;
}

int parse_double(double *x, const char *s, char **outendptr, int64_t *bytes_read)
{
    errno = 0;
    char *endptr;
    *x = strtod(s, &endptr);
    if ((errno == ERANGE && (*x == HUGE_VAL || *x == -HUGE_VAL)) || (errno != 0 && x == 0)) {
        return errno;
    }
    if (outendptr)
        *outendptr = endptr;
    if (bytes_read)
        *bytes_read += endptr - s;
    return 0;
}

/**
 * ‘freadline()’ reads a single line from a stream.
 */
static int freadline(char *linebuf, size_t line_max, enum streamtype streamtype, union stream stream)
{
    if (streamtype == stream_stdio) {
        char *s = fgets(linebuf, line_max + 1, stream.f);
        if (!s && feof(stream.f))
            return -1;
        else if (!s)
            return errno;
        int n = strlen(s);
        if (n > 0 && n == line_max && s[n - 1] != '\n')
            return EOVERFLOW;
        return 0;
#ifdef HAVE_LIBZ
    } else if (streamtype == stream_zlib) {
        char *s = gzgets(stream.gzf, linebuf, line_max + 1);
        if (!s && gzeof(stream.gzf))
            return -1;
        else if (!s)
            return errno;
        int n = strlen(s);
        if (n > 0 && n == line_max && s[n - 1] != '\n')
            return EOVERFLOW;
        return 0;
#endif
    } else {
        return EINVAL;
    }
}

#ifdef HAVE_ALIGNED_ALLOC
    long pagesize = sysconf(_SC_PAGESIZE);
#endif

    /* Set program invocation name. */
    program_invocation_name = argv[0];
    program_invocation_short_name = (
        strrchr(program_invocation_name, '/')
        ? strrchr(program_invocation_name, '/') + 1
        : program_invocation_name);

    /* 1. Parse program options. */
    struct program_options args;
    int nargs;
    err = parse_program_options(argc, argv, &args, &nargs);
    if (err) {
        fprintf(stderr, "%s: %s %s\n", program_invocation_short_name,
                strerror(err), argv[nargs]);
        return EXIT_FAILURE;
    }

#ifdef _OPENMP
    #pragma omp parallel
    {
      /*
       * This empty parallel section is used to make the OpenMP
       * runtime output its configuration now if the environment
       * variable OMP_DISPLAY_ENV is set.
       */
    }
#endif

    /* 2. Read the matrix from a Matrix Market file. */
    if (args.verbose > 0) {
        fprintf(stderr, "mtxfile_read: ");
        clock_gettime(CLOCK_MONOTONIC, &t0);
    }

    enum streamtype streamtype;
    union stream stream;
#ifdef HAVE_LIBZ
    if (!args.gzip) {
#endif
        streamtype = stream_stdio;
        if ((stream.f = fopen(args.Apath, "r")) == NULL) {
            fprintf(stderr, "%s: %s: %s\n",
                    program_invocation_short_name, args.Apath, strerror(errno));
            program_options_free(&args);
            return EXIT_FAILURE;
        }
#ifdef HAVE_LIBZ
    } else {
        streamtype = stream_zlib;
        if ((stream.gzf = gzopen(args.Apath, "r")) == NULL) {
            fprintf(stderr, "%s: %s: %s\n",
                    program_invocation_short_name, args.Apath, strerror(errno));
            program_options_free(&args);
            return EXIT_FAILURE;
        }
    }
#endif

    enum mtxobject object;
    enum mtxformat format;
    enum mtxfield field;
    enum mtxsymmetry symmetry;
    idx_t num_rows;
    idx_t num_columns;
    int64_t num_nonzeros;
    int64_t lines_read = 0;
    int64_t bytes_read = 0;
    err = mtxfile_fread_header(
        &object, &format, &field, &symmetry,
        &num_rows, &num_columns, &num_nonzeros,
        streamtype, stream, &lines_read, &bytes_read);
    if (err) {
        if (args.verbose > 0) fprintf(stderr, "\n");
        fprintf(stderr, "%s: %s:%"PRId64": %s\n",
                program_invocation_short_name,
                args.Apath, lines_read+1, strerror(err));
        stream_close(streamtype, stream);
        program_options_free(&args);
        return EXIT_FAILURE;
    }
#ifdef HAVE_ALIGNED_ALLOC
    size_t rowidxsize = num_nonzeros*sizeof(idx_t);
    idx_t * rowidx = aligned_alloc(pagesize, rowidxsize + pagesize - rowidxsize % pagesize);
#else
    idx_t * rowidx = malloc(num_nonzeros * sizeof(idx_t));
#endif
    if (!rowidx) {
        if (args.verbose > 0) fprintf(stderr, "\n");
        fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(errno));
        stream_close(streamtype, stream);
        program_options_free(&args);
        return EXIT_FAILURE;
    }
#ifdef HAVE_ALIGNED_ALLOC
    size_t colidxsize = num_nonzeros*sizeof(idx_t);
    idx_t * colidx = aligned_alloc(pagesize, colidxsize + pagesize - colidxsize % pagesize);
#else
    idx_t * colidx = malloc(num_nonzeros * sizeof(idx_t));
#endif
    if (!colidx) {
        if (args.verbose > 0) fprintf(stderr, "\n");
        fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(errno));
        free(rowidx);
        stream_close(streamtype, stream);
        program_options_free(&args);
        return EXIT_FAILURE;
    }
#ifdef HAVE_ALIGNED_ALLOC
    size_t asize = num_nonzeros*sizeof(double);
    double * a = aligned_alloc(pagesize, asize + pagesize - asize % pagesize);
#else
    double * a = malloc(num_nonzeros * sizeof(double));
#endif
    if (!a) {
        if (args.verbose > 0) fprintf(stderr, "\n");
        fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(errno));
        free(colidx); free(rowidx);
        stream_close(streamtype, stream);
        program_options_free(&args);
        return EXIT_FAILURE;
    }
    err = mtxfile_fread_matrix_coordinate(
        field, num_rows, num_columns, num_nonzeros, rowidx, colidx, a,
        streamtype, stream, &lines_read, &bytes_read);
    if (err) {
        if (args.verbose > 0) fprintf(stderr, "\n");
        fprintf(stderr, "%s: %s:%"PRId64": %s\n",
                program_invocation_short_name,
                args.Apath, lines_read+1, strerror(err));
        free(a); free(colidx); free(rowidx);
        stream_close(streamtype, stream);
        program_options_free(&args);
        return EXIT_FAILURE;
    }

    if (args.verbose > 0) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        fprintf(stderr, "%'.6f seconds (%'.1f MB/s)\n",
                timespec_duration(t0, t1),
                1.0e-6 * bytes_read / timespec_duration(t0, t1));
    }
    stream_close(streamtype, stream);

    /* 3. Convert to CSR format. */
    if (args.verbose > 0) {
        fprintf(stderr, "csr_from_coo: ");
        clock_gettime(CLOCK_MONOTONIC, &t0);
    }

#ifdef HAVE_ALIGNED_ALLOC
    size_t csrrowptrsize = (num_rows+1)*sizeof(int64_t);
    int64_t * csrrowptr = aligned_alloc(pagesize, csrrowptrsize + pagesize - csrrowptrsize % pagesize);
#else
    int64_t * csrrowptr = malloc((num_rows+1) * sizeof(int64_t));
#endif
    if (!csrrowptr) {
        if (args.verbose > 0) fprintf(stderr, "\n");
        fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(errno));
        free(a); free(colidx); free(rowidx);
        program_options_free(&args);
        return EXIT_FAILURE;
    }
    int64_t csrsize;
    idx_t rowsizemin, rowsizemax;
    idx_t diagsize;
    err = csr_from_coo_size(
        symmetry, num_rows, num_columns, num_nonzeros, rowidx, colidx, a,
        csrrowptr, &csrsize, &rowsizemin, &rowsizemax, &diagsize,
        args.separate_diagonal, args.partition);
    if (err) {
        if (args.verbose > 0) fprintf(stderr, "\n");
        fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(err));
        free(csrrowptr); free(a); free(colidx); free(rowidx);
        program_options_free(&args);
        return EXIT_FAILURE;
    }

    /* precompute per-thread partitioning of rows/columns/nonzeros */
    idx_t * startrows = NULL;
    idx_t * endrows = NULL;
    idx_t * startcolumns = NULL;
    idx_t * endcolumns = NULL;
#ifdef _OPENMP
    if (args.partition == partition_rows && args.rows_per_thread ||
        args.partition == partition_nonzeros && args.precompute_partition)
    {
        #pragma omp parallel
        #pragma omp master
        {
            int nthreads = omp_get_num_threads();
            startrows = malloc(nthreads * sizeof(idx_t));
        }
        if (!startrows) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(errno));
            free(csrrowptr); free(a); free(colidx); free(rowidx);
            program_options_free(&args);
            return EXIT_FAILURE;
        }
        #pragma omp parallel
        #pragma omp master
        {
            int nthreads = omp_get_num_threads();
            endrows = malloc(nthreads * sizeof(idx_t));
        }
        if (!endrows) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(errno));
            free(startrows);
            free(csrrowptr); free(a); free(colidx); free(rowidx);
            program_options_free(&args);
            return EXIT_FAILURE;
        }
    }
    if (args.partition == partition_rows && args.columns_per_thread)
    {
        #pragma omp parallel
        #pragma omp master
        {
            int nthreads = omp_get_num_threads();
            startcolumns = malloc(nthreads * sizeof(idx_t));
        }
        if (!startcolumns) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(errno));
            free(endrows); free(startrows);
            free(csrrowptr); free(a); free(colidx); free(rowidx);
            program_options_free(&args);
            return EXIT_FAILURE;
        }
        #pragma omp parallel
        #pragma omp master
        {
            int nthreads = omp_get_num_threads();
            endcolumns = malloc(nthreads * sizeof(idx_t));
        }
        if (!endcolumns) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(errno));
            free(startcolumns); free(endrows); free(startrows);
            free(csrrowptr); free(a); free(colidx); free(rowidx);
            program_options_free(&args);
            return EXIT_FAILURE;
        }
    }

    if (args.partition == partition_rows && args.rows_per_thread) {
        int nthreads;
        #pragma omp parallel
        #pragma omp master
        {
            nthreads = omp_get_num_threads();
            if (nthreads > 0) startrows[0] = 0;
            if (nthreads > 0) endrows[0] = args.rows_per_thread > 0 ? args.rows_per_thread[0] : 0;
            for (int p = 1; p < nthreads; p++) {
                startrows[p] = endrows[p-1];
                if (p < args.rows_per_thread_size) endrows[p] = startrows[p] + args.rows_per_thread[p];
                else endrows[p] = startrows[p];
            }
        }
        if (args.rows_per_thread_size != nthreads) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: warning: --rows-per-thread does not match the number of threads (%d)\n",
                    program_invocation_short_name, nthreads);
        }
        if (nthreads > 0 && endrows[nthreads-1] > num_rows) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: %s: the sum of --rows-per-thread (%'"PRIdx") exceeds the number of rows (%'"PRIdx")\n",
                    program_invocation_short_name, strerror(EINVAL), endrows[nthreads-1], num_rows);
            free(endrows); free(startrows);
            free(csrrowptr); free(a); free(colidx); free(rowidx);
            program_options_free(&args);
            return EXIT_FAILURE;
        } else if (nthreads > 0 && endrows[nthreads-1] < num_rows) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: warning: the sum of --rows-per-thread (%'"PRIdx") is less than the number of rows (%'"PRIdx")\n",
                    program_invocation_short_name, endrows[nthreads-1], num_rows);
        }
    } else if (args.partition == partition_nonzeros &&
               args.precompute_partition)
    {
        #pragma omp parallel
        {
            int nthreads = omp_get_num_threads();
            int p = omp_get_thread_num();
            int64_t startnz = p*(csrsize+nthreads-1)/nthreads;
            int64_t endnz = (p+1)*(csrsize+nthreads-1)/nthreads;
            if (endnz > csrsize) endnz = csrsize;
            idx_t startrow = 0;
            while (startrow < num_rows && startnz > csrrowptr[startrow+1]) startrow++;
            idx_t endrow = startrow;
            while (endrow < num_rows && endnz-1 > csrrowptr[endrow+1]) endrow++;
            startrows[p] = startrow;
            endrows[p] = endrow;
        }
    }

    if (args.partition == partition_rows && args.columns_per_thread) {
        int nthreads;
        #pragma omp parallel
        #pragma omp master
        {
            nthreads = omp_get_num_threads();
            if (nthreads > 0) startcolumns[0] = 0;
            if (nthreads > 0) endcolumns[0] = args.columns_per_thread > 0 ? args.columns_per_thread[0] : 0;
            for (int p = 1; p < nthreads; p++) {
                startcolumns[p] = endcolumns[p-1];
                if (p < args.columns_per_thread_size) endcolumns[p] = startcolumns[p] + args.columns_per_thread[p];
                else endcolumns[p] = startcolumns[p];
            }
        }
        if (args.columns_per_thread_size != nthreads) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: warning: --columns-per-thread does not match the number of threads (%d)\n",
                    program_invocation_short_name, nthreads);
        }
        if (nthreads > 0 && endcolumns[nthreads-1] > num_columns) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: %s: the sum of --columns-per-thread (%'"PRIdx") exceeds the number of columns (%'"PRIdx")\n",
                    program_invocation_short_name, strerror(EINVAL), endcolumns[nthreads-1], num_columns);
            free(endcolumns); free(startcolumns); free(endrows); free(startrows);
            free(csrrowptr); free(a); free(colidx); free(rowidx);
            program_options_free(&args);
            return EXIT_FAILURE;
        } else if (nthreads > 0 && endcolumns[nthreads-1] < num_columns) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: warning: the sum of --columns-per-thread (%'"PRIdx") is less than the number of columns (%'"PRIdx")\n",
                    program_invocation_short_name, endcolumns[nthreads-1], num_columns);
        }
    }
#endif

#ifdef HAVE_ALIGNED_ALLOC
    size_t csrcolidxsize = csrsize*sizeof(idx_t);
    idx_t * csrcolidx = aligned_alloc(pagesize, csrcolidxsize + pagesize - csrcolidxsize % pagesize);
#else
    idx_t * csrcolidx = malloc(csrsize * sizeof(idx_t));
#endif
    if (!csrcolidx) {
        if (args.verbose > 0) fprintf(stderr, "\n");
        fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(errno));
        free(endcolumns); free(startcolumns); free(endrows); free(startrows);
        free(csrrowptr); free(a); free(colidx); free(rowidx);
        program_options_free(&args);
        return EXIT_FAILURE;
    }
#ifdef _OPENMP
    if (args.partition == partition_rows && !args.rows_per_thread) {
        #pragma omp parallel for
        for (idx_t i = 0; i < num_rows; i++) {
            for (int64_t k = csrrowptr[i]; k < csrrowptr[i+1]; k++)
                csrcolidx[k] = 0;
        }
    } else if (args.partition == partition_rows) {
        #pragma omp parallel
        {
            int p = omp_get_thread_num();
            for (idx_t i = startrows[p]; i < endrows[p]; i++) {
                for (int64_t k = csrrowptr[i]; k < csrrowptr[i+1]; k++)
                    csrcolidx[k] = 0;
            }
        }
    } else if (args.partition == partition_nonzeros) {
        #pragma omp parallel for
        for (int64_t k = 0; k < csrsize; k++) csrcolidx[k] = 0;
    }
#endif
#ifdef HAVE_ALIGNED_ALLOC
    size_t csrasize = csrsize*sizeof(double);
    double * csra = aligned_alloc(pagesize, csrasize + pagesize - csrasize % pagesize);
#else
    double * csra = malloc(csrsize * sizeof(double));
#endif
    if (!csra) {
        if (args.verbose > 0) fprintf(stderr, "\n");
        fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(errno));
        free(csrcolidx);
        free(endcolumns); free(startcolumns); free(endrows); free(startrows);
        free(csrrowptr); free(a); free(colidx); free(rowidx);
        program_options_free(&args);
        return EXIT_FAILURE;
    }
#ifdef HAVE_ALIGNED_ALLOC
    size_t csradsize = diagsize*sizeof(double);
    double * csrad = aligned_alloc(pagesize, csradsize + pagesize - csradsize % pagesize);
#else
    double * csrad = malloc(diagsize * sizeof(double));
#endif
    if (!csrad) {
        if (args.verbose > 0) fprintf(stderr, "\n");
        fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(errno));
        free(csra); free(csrcolidx);
        free(endcolumns); free(startcolumns); free(endrows); free(startrows);
        free(csrrowptr); free(a); free(colidx); free(rowidx);
        program_options_free(&args);
        return EXIT_FAILURE;
    }
#ifdef _OPENMP
    if (args.partition == partition_rows && !args.rows_per_thread) {
        #pragma omp parallel for
        for (idx_t i = 0; i < num_rows; i++) {
            for (int64_t k = csrrowptr[i]; k < csrrowptr[i+1]; k++)
                csra[k] = 0;
        }
        if (diagsize > 0) {
            #pragma omp parallel for
            for (idx_t i = 0; i < num_rows; i++) csrad[i] = 0;
        }
    } else if (args.partition == partition_rows) {
        #pragma omp parallel
        {
            int p = omp_get_thread_num();
            for (idx_t i = startrows[p]; i < endrows[p]; i++) {
                for (int64_t k = csrrowptr[i]; k < csrrowptr[i+1]; k++)
                    csra[k] = 0;
            }
            if (diagsize > 0) {
                for (idx_t i = startrows[p]; i < endrows[p]; i++) csrad[i] = 0;
            }
        }
    } else if (args.partition == partition_nonzeros) {
        #pragma omp parallel for
        for (int64_t k = 0; k < csrsize; k++) csra[k] = 0;
        if (diagsize > 0) {
            #pragma omp parallel for
            for (idx_t i = 0; i < num_rows; i++) csrad[i] = 0;
        }
    }
#endif
    err = csr_from_coo(
        symmetry, num_rows, num_columns, num_nonzeros, rowidx, colidx, a,
        csrrowptr, csrsize, rowsizemin, rowsizemax, csrcolidx, csra, csrad,
        args.separate_diagonal, args.sort_rows, args.partition);
    if (err) {
        if (args.verbose > 0) fprintf(stderr, "\n");
        fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(err));
        free(csrad); free(csra); free(csrcolidx);
        free(endcolumns); free(startcolumns); free(endrows); free(startrows);
        free(csrrowptr); free(a); free(colidx); free(rowidx);
        program_options_free(&args);
        return EXIT_FAILURE;
    }
    free(a); free(colidx); free(rowidx);

    if (args.verbose > 0) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        fprintf(stderr, "%'.6f seconds, %'"PRIdx" rows, %'"PRIdx" columns, %'"PRId64" nonzeros"
                ", %'"PRIdx" to %'"PRIdx" nonzeros per row",
                timespec_duration(t0, t1), num_rows, num_columns, csrsize+diagsize, rowsizemin, rowsizemax);
#ifdef _OPENMP
        int nthreads;
        idx_t min_rows_per_thread = IDX_T_MAX;
        idx_t max_rows_per_thread = 0;
        int64_t min_nonzeros_per_thread = INT64_MAX;
        int64_t max_nonzeros_per_thread = 0;
        if (args.partition == partition_rows && !args.rows_per_thread) {
            #pragma omp parallel \
                reduction(min:min_rows_per_thread) reduction(max:max_rows_per_thread) \
                reduction(min:min_nonzeros_per_thread) reduction(max:max_nonzeros_per_thread)
            {
                nthreads = omp_get_num_threads();
                int p = omp_get_thread_num();
                min_rows_per_thread = max_rows_per_thread = num_rows/nthreads + (p < (num_rows % nthreads));
                int64_t num_nonzeros = 0;
                #pragma omp for
                for (int i = 0; i < num_rows; i++)
                    num_nonzeros += csrrowptr[i+1]-csrrowptr[i] + (diagsize > 0 ? 1 : 0);
                min_nonzeros_per_thread = num_nonzeros;
                max_nonzeros_per_thread = num_nonzeros;
            }
        } else if (args.partition == partition_rows) {
            #pragma omp parallel \
                reduction(min:min_rows_per_thread) reduction(max:max_rows_per_thread) \
                reduction(min:min_nonzeros_per_thread) reduction(max:max_nonzeros_per_thread)
            {
                nthreads = omp_get_num_threads();
                int p = omp_get_thread_num();
                idx_t startrow = startrows[p];
                idx_t endrow = endrows[p];
                min_rows_per_thread = max_rows_per_thread = endrow - startrow;
                int64_t num_nonzeros = 0;
                for (idx_t i = startrows[p]; i < endrows[p]; i++)
                    num_nonzeros += csrrowptr[i+1]-csrrowptr[i] + (diagsize > 0 ? 1 : 0);
                min_nonzeros_per_thread = num_nonzeros;
                max_nonzeros_per_thread = num_nonzeros;
            }
        } else if (args.partition == partition_nonzeros) {
            #pragma omp parallel \
                reduction(min:min_rows_per_thread) reduction(max:max_rows_per_thread) \
                reduction(min:min_nonzeros_per_thread) reduction(max:max_nonzeros_per_thread)
            {
                nthreads = omp_get_num_threads();
                int p = omp_get_thread_num();
                int64_t startnz = p*(csrsize+nthreads-1)/nthreads;
                int64_t endnz = (p+1)*(csrsize+nthreads-1)/nthreads;
                if (endnz > csrsize) endnz = csrsize;
                idx_t startrow = 0;
                if (startrows) { startrow = startrows[p]; }
                else { while (startrow < num_rows && startnz > csrrowptr[startrow+1]) startrow++; }
                idx_t endrow = startrow;
                if (endrows) { endrow = endrows[p]; }
                else { while (endrow < num_rows && endnz-1 > csrrowptr[endrow+1]) endrow++; }
                min_rows_per_thread = max_rows_per_thread = endrow - startrow;
                min_nonzeros_per_thread = max_nonzeros_per_thread = csrsize/nthreads + (p < (csrsize % nthreads));
            }
        }
        fprintf(stderr, ", %'d threads, %'"PRIdx" to %'"PRIdx" rows per thread, %'"PRId64" to %'"PRId64" nonzeros per thread",
                nthreads, min_rows_per_thread, max_rows_per_thread,
                min_nonzeros_per_thread, max_nonzeros_per_thread);
#endif
        fputc('\n', stderr);
    }

    /* 4. allocate vectors */
#ifdef HAVE_ALIGNED_ALLOC
    size_t xsize = num_columns*sizeof(double);
    double * x = aligned_alloc(pagesize, xsize + pagesize - xsize % pagesize);
#else
    double * x = malloc(num_columns * sizeof(double));
#endif
    if (!x) {
        if (args.verbose > 0) fprintf(stderr, "\n");
        fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(errno));
        free(endcolumns); free(startcolumns); free(endrows); free(startrows);
        free(csrad); free(csra); free(csrcolidx); free(csrrowptr);
        program_options_free(&args);
        return EXIT_FAILURE;
    }

#ifdef _OPENMP
    if (args.partition == partition_rows && args.columns_per_thread) {
        #pragma omp parallel
        {
            int p = omp_get_thread_num();
            for (idx_t i = startcolumns[p]; i < endcolumns[p]; i++) x[i] = 1.0;
            int nthreads = omp_get_num_threads();
            #pragma omp master
            for (idx_t i = endcolumns[nthreads-1]; i < num_columns; i++) x[i] = 1.0;
        }
    } else if (args.partition == partition_rows && args.rows_per_thread &&
               num_rows == num_columns)
    {
        #pragma omp parallel
        {
            int p = omp_get_thread_num();
            for (idx_t i = startrows[p]; i < endrows[p]; i++) x[i] = 1.0;
            int nthreads = omp_get_num_threads();
            #pragma omp master
            for (idx_t i = endrows[nthreads-1]; i < num_rows; i++) x[i] = 1.0;
        }
    } else {
        #pragma omp parallel for
        for (idx_t i = 0; i < num_columns; i++) x[i] = 1.0;
    }
#else
    for (idx_t i = 0; i < num_columns; i++) x[i] = 1.0;
#endif

    /* read x vector from a Matrix Market file */
    if (args.xpath) {
        if (args.verbose > 0) {
            fprintf(stderr, "mtxfile_read: ");
            clock_gettime(CLOCK_MONOTONIC, &t0);
        }

        enum streamtype streamtype;
        union stream stream;
#ifdef HAVE_LIBZ
        if (!args.gzip) {
#endif
            streamtype = stream_stdio;
            if ((stream.f = fopen(args.xpath, "r")) == NULL) {
                fprintf(stderr, "%s: %s: %s\n",
                        program_invocation_short_name, args.xpath, strerror(errno));
                free(x);
                free(endcolumns); free(startcolumns); free(endrows); free(startrows);
                free(csrad); free(csra); free(csrcolidx); free(csrrowptr);
                program_options_free(&args);
                return EXIT_FAILURE;
            }
#ifdef HAVE_LIBZ
        } else {
            streamtype = stream_zlib;
            if ((stream.gzf = gzopen(args.xpath, "r")) == NULL) {
                fprintf(stderr, "%s: %s: %s\n",
                        program_invocation_short_name, args.xpath, strerror(errno));
                free(x);
                free(endcolumns); free(startcolumns); free(endrows); free(startrows);
                free(csrad); free(csra); free(csrcolidx); free(csrrowptr);
                program_options_free(&args);
                return EXIT_FAILURE;
            }
        }
#endif

        enum mtxobject object;
        enum mtxformat format;
        enum mtxfield field;
        enum mtxsymmetry symmetry;
        idx_t xnum_rows;
        idx_t xnum_columns;
        int64_t xnum_nonzeros;
        int64_t lines_read = 0;
        int64_t bytes_read = 0;
        err = mtxfile_fread_header(
            &object, &format, &field, &symmetry,
            &xnum_rows, &xnum_columns, &xnum_nonzeros,
            streamtype, stream, &lines_read, &bytes_read);
        if (err) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: %s:%"PRId64": %s\n",
                    program_invocation_short_name,
                    args.xpath, lines_read+1, strerror(err));
            stream_close(streamtype, stream);
            free(x);
            free(endcolumns); free(startcolumns); free(endrows); free(startrows);
            free(csrad); free(csra); free(csrcolidx); free(csrrowptr);
            program_options_free(&args);
            return EXIT_FAILURE;
        } else if (object != mtxvector || format != mtxarray || xnum_rows != num_columns) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: %s:%"PRId64": "
                    "expected vector in array format of size %"PRIdx"\n",
                    program_invocation_short_name,
                    args.xpath, lines_read+1, num_columns);
            stream_close(streamtype, stream);
            free(x);
            free(endcolumns); free(startcolumns); free(endrows); free(startrows);
            free(csrad); free(csra); free(csrcolidx); free(csrrowptr);
            program_options_free(&args);
            return EXIT_FAILURE;
        }

        err = mtxfile_fread_vector_array(
            field, num_rows, x, streamtype, stream, &lines_read, &bytes_read);
        if (err) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: %s:%"PRId64": %s\n",
                    program_invocation_short_name,
                    args.xpath, lines_read+1, strerror(err));
            free(x);
            free(endcolumns); free(startcolumns); free(endrows); free(startrows);
            free(csrad); free(csra); free(csrcolidx); free(csrrowptr);
            stream_close(streamtype, stream);
            program_options_free(&args);
            return EXIT_FAILURE;
        }

        if (args.verbose > 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            fprintf(stderr, "%'.6f seconds (%'.1f MB/s)\n",
                    timespec_duration(t0, t1),
                    1.0e-6 * bytes_read / timespec_duration(t0, t1));
        }
        stream_close(streamtype, stream);
    }

#ifdef HAVE_ALIGNED_ALLOC
    size_t ysize = num_rows*sizeof(double);
    double * y = aligned_alloc(pagesize, ysize + pagesize - ysize % pagesize);
#else
    double * y = malloc(num_rows * sizeof(double));
#endif
    if (!y) {
        if (args.verbose > 0) fprintf(stderr, "\n");
        fprintf(stderr, "%s: %s\n", program_invocation_short_name, strerror(errno));
        free(x);
        free(endcolumns); free(startcolumns); free(endrows); free(startrows);
        free(csrad); free(csra); free(csrcolidx); free(csrrowptr);
        program_options_free(&args);
        return EXIT_FAILURE;
    }

#ifdef _OPENMP
    if (args.partition == partition_rows && !args.rows_per_thread) {
        #pragma omp parallel for
        for (idx_t i = 0; i < num_rows; i++) y[i] = 0.0;
    } else if (args.partition == partition_rows) {
        #pragma omp parallel
        {
            int p = omp_get_thread_num();
            for (idx_t i = startrows[p]; i < endrows[p]; i++) y[i] = 0.0;
            int nthreads = omp_get_num_threads();
            #pragma omp master
            for (idx_t i = endrows[nthreads-1]; i < num_rows; i++) y[i] = 0;
        }
    } else if (args.partition == partition_nonzeros) {
        #pragma omp parallel
        {
            int nthreads = omp_get_num_threads();
            int p = omp_get_thread_num();
            int64_t startnz = p*(csrsize+nthreads-1)/nthreads;
            int64_t endnz = (p+1)*(csrsize+nthreads-1)/nthreads;
            if (endnz > csrsize) endnz = csrsize;
            idx_t startrow = 0;
            if (startrows) { startrow = startrows[p]; }
            else { while (startrow < num_rows && startnz > csrrowptr[startrow+1]) startrow++; }
            idx_t endrow = startrow;
            if (endrows) { endrow = endrows[p]; }
            else { while (endrow < num_rows && endnz-1 > csrrowptr[endrow+1]) endrow++; }
            for (idx_t i = startrow; i < endrow; i++) y[i] = 0.0;
        }
    }
#else
    for (idx_t i = 0; i < num_rows; i++) y[i] = 0.0;
#endif

    /* read y vector from a Matrix Market file */
    if (args.ypath) {
        if (args.verbose > 0) {
            fprintf(stderr, "mtxfile_read: ");
            clock_gettime(CLOCK_MONOTONIC, &t0);
        }

        enum streamtype streamtype;
        union stream stream;
#ifdef HAVE_LIBZ
        if (!args.gzip) {
#endif
            streamtype = stream_stdio;
            if ((stream.f = fopen(args.ypath, "r")) == NULL) {
                fprintf(stderr, "%s: %s: %s\n",
                        program_invocation_short_name, args.ypath, strerror(errno));
                free(y); free(x);
                free(endcolumns); free(startcolumns); free(endrows); free(startrows);
                free(csrad); free(csra); free(csrcolidx); free(csrrowptr);
                program_options_free(&args);
                return EXIT_FAILURE;
            }
#ifdef HAVE_LIBZ
        } else {
            streamtype = stream_zlib;
            if ((stream.gzf = gzopen(args.ypath, "r")) == NULL) {
                fprintf(stderr, "%s: %s: %s\n",
                        program_invocation_short_name, args.ypath, strerror(errno));
                free(y); free(x);
                free(endcolumns); free(startcolumns); free(endrows); free(startrows);
                free(csrad); free(csra); free(csrcolidx); free(csrrowptr);
                program_options_free(&args);
                return EXIT_FAILURE;
            }
        }
#endif

        enum mtxobject object;
        enum mtxformat format;
        enum mtxfield field;
        enum mtxsymmetry symmetry;
        idx_t ynum_rows;
        idx_t ynum_columns;
        int64_t ynum_nonzeros;
        int64_t lines_read = 0;
        int64_t bytes_read = 0;
        err = mtxfile_fread_header(
            &object, &format, &field, &symmetry,
            &ynum_rows, &ynum_columns, &ynum_nonzeros,
            streamtype, stream, &lines_read, &bytes_read);
        if (err) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: %s:%"PRId64": %s\n",
                    program_invocation_short_name,
                    args.ypath, lines_read+1, strerror(err));
            stream_close(streamtype, stream);
            free(y); free(x);
            free(endcolumns); free(startcolumns); free(endrows); free(startrows);
            free(csrad); free(csra); free(csrcolidx); free(csrrowptr);
            program_options_free(&args);
            return EXIT_FAILURE;
        } else if (object != mtxvector || format != mtxarray || ynum_rows != num_rows) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: %s:%"PRId64": "
                    "expected vector in array format of size %'"PRIdx"\n",
                    program_invocation_short_name,
                    args.ypath, lines_read+1, num_rows);
            stream_close(streamtype, stream);
            free(y); free(x);
            free(endcolumns); free(startcolumns); free(endrows); free(startrows);
            free(csrad); free(csra); free(csrcolidx); free(csrrowptr);
            program_options_free(&args);
            return EXIT_FAILURE;
        }

        err = mtxfile_fread_vector_array(
            field, num_rows, y, streamtype, stream, &lines_read, &bytes_read);
        if (err) {
            if (args.verbose > 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s: %s:%"PRId64": %s\n",
                    program_invocation_short_name,
                    args.ypath, lines_read+1, strerror(err));
            free(y); free(x);
            free(endcolumns); free(startcolumns); free(endrows); free(startrows);
            free(csrad); free(csra); free(csrcolidx); free(csrrowptr);
            stream_close(streamtype, stream);
            program_options_free(&args);
            return EXIT_FAILURE;
        }

        if (args.verbose > 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            fprintf(stderr, "%'.6f seconds (%'.1f MB/s)\n",
                    timespec_duration(t0, t1),
                    1.0e-6 * bytes_read / timespec_duration(t0, t1));
        }
        stream_close(streamtype, stream);
    }
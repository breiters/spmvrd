#if CALCULATE_NNZ_PER_ROW_VARIANCE
double variance = 0.0;

double avg_nnzs_per_row = (double)matrix.nnz / matrix.nrow;
for (unsigned r = 0; r < matrix.nrow; ++r) {
    int nnzs = matrix.row_ptr[r + 1] - matrix.row_ptr[r];
    variance += ((double)nnzs - avg_nnzs_per_row) * ((double)nnzs - avg_nnzs_per_row);
}
variance /= matrix.nrow;
FILE *varfile = fopen("variance.csv", "a");
fprintf(varfile, "%s,%f,%f\n", matrix_path, variance, variance / avg_nnzs_per_row);
fclose(varfile);
exit(0);
#endif

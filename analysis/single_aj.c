/* Exact arithmetic acceleration of the frozen cp_b_aj recurrence.
 * No fitted hazard, coefficient, tie convention or estimator is changed.
 * Column-major hazards: jumps x 2; relative hazards: participants x 2.
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <math.h>

SEXP ff_aj_native(SEXP hazards, SEXP relative) {
    SEXP hd = getAttrib(hazards, R_DimSymbol);
    SEXP rd = getAttrib(relative, R_DimSymbol);
    if (TYPEOF(hazards) != REALSXP || TYPEOF(relative) != REALSXP ||
        length(hd) != 2 || length(rd) != 2 || INTEGER(hd)[1] != 2 ||
        INTEGER(rd)[1] != 2 || INTEGER(rd)[0] < 1)
        error("FF_AJ_NATIVE_MATRIX_DIMENSION_INVALID");
    int jumps = INTEGER(hd)[0], n = INTEGER(rd)[0];
    double *h = REAL(hazards), *rr = REAL(relative);
    for (int j = 0; j < 2 * jumps; ++j)
        if (!R_FINITE(h[j]) || h[j] < 0)
            error("FF_AJ_NATIVE_HAZARD_INVALID_NO_CLIPPING");
    for (int i = 0; i < 2 * n; ++i)
        if (!R_FINITE(rr[i]) || rr[i] <= 0)
            error("FF_AJ_NATIVE_RELATIVE_HAZARD_INVALID_NO_CLIPPING");
    SEXP out = PROTECT(allocMatrix(REALSXP, n, 3));
    double *p = REAL(out);
    for (int i = 0; i < n; ++i) {
        p[i] = 0.0; p[i + n] = 0.0; p[i + 2 * n] = 1.0;
    }
    for (int j = 0; j < jumps; ++j) {
        if (j % 128 == 0) R_CheckUserInterrupt();
        for (int i = 0; i < n; ++i) {
            double h1 = h[j] * rr[i], h2 = h[j + jumps] * rr[i + n];
            if (!R_FINITE(h1 + h2) || h1 + h2 > 1.0)
                error("B_AJ_JUMP_PROBABILITY_INVALID");
            double s = p[i + 2 * n];
            p[i] = p[i] + s * h1;
            p[i + n] = p[i + n] + s * h2;
            p[i + 2 * n] = s * (1.0 - h1 - h2);
        }
    }
    for (int i = 0; i < n; ++i)
        if (!R_FINITE(p[i] + p[i + n] + p[i + 2 * n]) ||
            p[i] < 0 || p[i + n] < 0 || p[i + 2 * n] < 0 ||
            p[i] > 1 || p[i + n] > 1 || p[i + 2 * n] > 1 ||
            fabs(p[i] + p[i + n] + p[i + 2 * n] - 1.0) >= 1e-9)
            error("B_PIECEWISE_COMPETING_PROBABILITY_FAILURE");
    UNPROTECT(1);
    return out;
}

static const R_CallMethodDef call_methods[] = {
    {"ff_aj_native", (DL_FUNC) &ff_aj_native, 2},
    {NULL, NULL, 0}
};

void R_init_ff_aj(DllInfo *dll) {
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

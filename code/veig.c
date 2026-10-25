/* Approximate lowest eigenpair of a symmetric matrix by inverse iteration (one LU, many solves).
 * Only used to choose WHAT to certify (vchol does the certifying), so plain floating point is fine.
 * Input file: n, then n*n entries (row-major, decimal). Output: lambda, then the n vector entries.
 * Build: gcc -O2 -o veig veig.c -lflint -lmpfr -lgmp
 */
#include <stdio.h>
#include <stdlib.h>
#include <flint/flint.h>
#include <flint/arb.h>
#include <flint/arb_mat.h>

int main(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "usage: veig file prec [iters]\n"); return 2; }
    slong prec = atol(argv[2]), iters = argc > 3 ? atol(argv[3]) : 20;
    FILE *f = fopen(argv[1], "r");
    static char buf[1 << 16];
    slong n;
    if (!f || !fgets(buf, sizeof buf, f) || sscanf(buf, "%ld", &n) != 1) { printf("FAIL\n"); return 1; }
    arb_mat_t A, LU, x, y;
    arb_mat_init(A, n, n); arb_mat_init(LU, n, n); arb_mat_init(x, n, 1); arb_mat_init(y, n, 1);
    for (slong i = 0; i < n * n; i++)
    {
        if (!fgets(buf, sizeof buf, f)) { printf("FAIL\n"); return 1; }
        arb_set_str(arb_mat_entry(A, i / n, i % n), buf, prec);
        arb_get_mid_arb(arb_mat_entry(A, i / n, i % n), arb_mat_entry(A, i / n, i % n));
    }
    fclose(f);
    slong *perm = flint_malloc(sizeof(slong) * n);
    if (!arb_mat_approx_lu(perm, LU, A, prec)) { printf("FAIL lu\n"); return 1; }
    arb_t nrm, t, lam;
    arb_init(nrm); arb_init(t); arb_init(lam);
    for (slong i = 0; i < n; i++) arb_one(arb_mat_entry(x, i, 0));
    for (slong it = 0; it < iters; it++)
    {
        arb_mat_approx_solve_lu_precomp(y, perm, LU, x, prec);
        arb_zero(nrm);
        for (slong i = 0; i < n; i++) arb_addmul(nrm, arb_mat_entry(y, i, 0), arb_mat_entry(y, i, 0), prec);
        arb_sqrt(nrm, nrm, prec);
        for (slong i = 0; i < n; i++)
        {
            arb_div(arb_mat_entry(x, i, 0), arb_mat_entry(y, i, 0), nrm, prec);
            arb_get_mid_arb(arb_mat_entry(x, i, 0), arb_mat_entry(x, i, 0));
        }
    }
    arb_mat_approx_mul(y, A, x, prec);  /* Rayleigh quotient x^T A x (x has unit norm) */
    arb_zero(lam);
    for (slong i = 0; i < n; i++) arb_addmul(lam, arb_mat_entry(x, i, 0), arb_mat_entry(y, i, 0), prec);
    printf("%s\n", arb_get_str(lam, 40, ARB_STR_NO_RADIUS));
    for (slong i = 0; i < n; i++) printf("%s\n", arb_get_str(arb_mat_entry(x, i, 0), 60, ARB_STR_NO_RADIUS));
    flint_free(perm);
    arb_mat_clear(A); arb_mat_clear(LU); arb_mat_clear(x); arb_mat_clear(y);
    arb_clear(nrm); arb_clear(t); arb_clear(lam);
    flint_cleanup();
    return 0;
}

/* Verified positive semidefiniteness: is A - lam*I >= 0 ?
 *
 * Input file: n, lam, then the n*n entries of A (row-major), each a ball "[mid +/- rad]" or a decimal.
 * Method (same as certify.certified_pd_py):
 *   1. R = approximate Cholesky of mid(A) - (lam+d) I, d = lam/100, on midpoints only (every entry rounded to its
 *      midpoint as it is produced; FLINT's ball Cholesky arb_mat_cho blows up at n ~ 2500).
 *   2. E = A - (lam+d) I - R R^T, in ball arithmetic (multithreaded).
 *   3. If ||E||_F < d then A - lam I = R R^T + E + d I >= 0, rigorously.
 * Prints "PD <frobenius bound>" or "FAIL <reason>".
 * Build: gcc -O2 -o vchol vchol.c -lflint -lmpfr -lgmp
 */
#include <stdio.h>
#include <stdlib.h>
#include <flint/flint.h>
#include <flint/arb.h>
#include <flint/arb_mat.h>

static int read_ball(FILE *f, arb_t x, slong prec)
{
    static char buf[1 << 16];
    if (!fgets(buf, sizeof buf, f)) return 0;
    return arb_set_str(x, buf, prec) == 0;
}

int main(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "usage: vchol file prec [threads]\n"); return 2; }
    slong prec = atol(argv[2]);
    flint_set_num_threads(argc > 3 ? atoi(argv[3]) : 1);
    FILE *f = fopen(argv[1], "r");
    if (!f) { printf("FAIL open\n"); return 1; }
    slong n; char line[256];
    if (!fgets(line, sizeof line, f) || sscanf(line, "%ld", &n) != 1) { printf("FAIL header\n"); return 1; }

    arb_t lam, d, s, fro, e;
    arb_init(lam); arb_init(d); arb_init(s); arb_init(fro); arb_init(e);
    if (!read_ball(f, lam, prec)) { printf("FAIL lam\n"); return 1; }
    arb_div_ui(d, lam, argc > 4 ? atol(argv[4]) : 100, prec);  /* d = lam/ddiv */
    arb_add(s, lam, d, prec);  /* shift = lam + d */

    arb_mat_t A, M, R, RT, P;
    arb_mat_init(A, n, n); arb_mat_init(M, n, n); arb_mat_init(R, n, n);
    arb_mat_init(RT, n, n); arb_mat_init(P, n, n);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++)
            if (!read_ball(f, arb_mat_entry(A, i, j), prec)) { printf("FAIL entry %ld %ld\n", i, j); return 1; }
    fclose(f);

    /* 1. Cholesky of the exact midpoint matrix minus the shift */
    arb_mat_get_mid(M, A);
    for (slong i = 0; i < n; i++)
    {
        arb_sub(arb_mat_entry(M, i, i), arb_mat_entry(M, i, i), s, prec);
        arb_get_mid_arb(arb_mat_entry(M, i, i), arb_mat_entry(M, i, i));
    }
    arb_t piv, tmp;
    arb_init(piv); arb_init(tmp);
    arb_mat_zero(R);
    for (slong j = 0; j < n; j++)
    {
        arb_ptr Rj = arb_mat_entry(R, j, 0);
        arb_dot(piv, arb_mat_entry(M, j, j), 1, Rj, 1, Rj, 1, j, prec);  /* M_jj - sum_k R_jk^2 */
        arb_get_mid_arb(piv, piv);
        if (!arb_is_positive(piv)) { printf("FAIL cholesky pivot %ld\n", j); return 0; }
        arb_sqrt(arb_mat_entry(R, j, j), piv, prec);
        arb_get_mid_arb(arb_mat_entry(R, j, j), arb_mat_entry(R, j, j));
        for (slong i = j + 1; i < n; i++)
        {
            arb_ptr Ri = arb_mat_entry(R, i, 0);
            arb_dot(tmp, arb_mat_entry(M, i, j), 1, Ri, 1, Rj, 1, j, prec);
            arb_div(tmp, tmp, arb_mat_entry(R, j, j), prec);
            arb_get_mid_arb(arb_mat_entry(R, i, j), tmp);
        }
    }
    arb_clear(piv); arb_clear(tmp);

    /* 2. E = A - s I - R R^T (rigorous) */
    arb_mat_transpose(RT, R);
    arb_mat_mul(P, R, RT, prec);  /* multithreaded when threads > 1 */
    arb_mat_sub(M, A, P, prec);
    for (slong i = 0; i < n; i++)
        arb_sub(arb_mat_entry(M, i, i), arb_mat_entry(M, i, i), s, prec);

    /* 3. Frobenius upper bound */
    arb_zero(fro);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++)
        {
            arb_get_abs_ubound_arf(arb_midref(e), arb_mat_entry(M, i, j), prec);
            mag_zero(arb_radref(e));
            arb_addmul(fro, e, e, prec);
        }
    arb_sqrtpos(fro, fro, prec);
    arb_get_ubound_arf(arb_midref(fro), fro, prec);
    mag_zero(arb_radref(fro));

    if (arb_lt(fro, d)) printf("PD %s\n", arb_get_str(fro, 10, 0));
    else printf("FAIL residual %s >= %s\n", arb_get_str(fro, 10, 0), arb_get_str(d, 10, 0));

    arb_mat_clear(A); arb_mat_clear(M); arb_mat_clear(R); arb_mat_clear(RT); arb_mat_clear(P);
    arb_clear(lam); arb_clear(d); arb_clear(s); arb_clear(fro); arb_clear(e);
    flint_cleanup();
    return 0;
}

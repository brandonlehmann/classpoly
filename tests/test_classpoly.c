#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "ff_poly.h"
#include "classpoly.h"
#include "cstd.h"

/*
    Test harness for classpoly serial and parallel paths.

    Runs both compute_classpoly() and compute_classpoly_parallel() on a set
    of known discriminants, then compares output coefficients against cached
    known-good reference values (from ~/finder/classpoly_cache/).

    Usage: test_classpoly <path_to_tests_data_dir>

    6-way cross-validation per discriminant:
      - serial P1 vs cache P1
      - serial P2 vs cache P2
      - parallel P1 vs cache P1
      - parallel P2 vs cache P2
      - serial P1 vs parallel P1
      - serial P2 vs parallel P2
*/

#define MAX_COEFFS  2048
#define MAX_LINE    256

struct test_vector {
	long D;
	const char *P1;
	const char *P2;
};

static struct test_vector test_vectors[] = {
	{-167995,
	 "57896044618658097711785492504343953926348427374841791341041551676248176146099",
	 "57896044618658097711785492504343953926634992332820282019728792003956564819949"},
	{-7857907,
	 "57896044618658097711785492504343953926549254372227246365156541811699034343327",
	 "57896044618658097711785492504343953926634992332820282019728792003956564819949"},
	{-4059339,
	 "57896044618658097711785492504343953926452459951711404613951846209448911595739",
	 "57896044618658097711785492504343953926634992332820282019728792003956564819949"},
	{-14141835,
	 "57896044618658097711785492504343953926543511531128800564411133339292026663159",
	 "57896044618658097711785492504343953926634992332820282019728792003956564819949"},
	{-54558307,
	 "57896044618658097711785492504343953926342357520481333358459613131403471893763",
	 "57896044618658097711785492504343953926634992332820282019728792003956564819949"},
	{-64752003,
	 "57896044618658097711785492504343953926483949935030029656422428983527059104287",
	 "57896044618658097711785492504343953926634992332820282019728792003956564819949"},
	{-92169307,
	 "57896044618658097711785492504343953926544611479355037105499359504869493802389",
	 "57896044618658097711785492504343953926634992332820282019728792003956564819949"},
};

#define NUM_TESTS (sizeof(test_vectors) / sizeof(test_vectors[0]))

/*
    Parse classpoly output file: extract coefficient strings.
    Format: lines like "12345*X^0 + \n" and "1*X^N\n" (leading term).
    Returns count of coefficients (including leading 1).
*/
static int parse_classpoly_output (const char *filename, char coeffs[][MAX_LINE], int max_coeffs)
{
	FILE *fp = fopen(filename, "r");
	char line[4096];
	int count = 0;

	if ( !fp ) { fprintf(stderr, "Cannot open %s\n", filename); return -1; }

	while ( fgets(line, sizeof(line), fp) && count < max_coeffs ) {
		char *star = strstr(line, "*X^");
		if ( !star ) continue;
		/* Extract coefficient: everything before *X^ */
		int len = (int)(star - line);
		if ( len <= 0 || len >= MAX_LINE ) continue;
		strncpy(coeffs[count], line, len);
		coeffs[count][len] = '\0';
		/* Strip trailing whitespace */
		while ( len > 0 && (coeffs[count][len-1] == ' ' || coeffs[count][len-1] == '\n') )
			coeffs[count][--len] = '\0';
		count++;
	}
	fclose(fp);
	return count;
}

/*
    Load reference coefficients from .coeffs file (one bare integer per line).
*/
static int load_reference_coeffs (const char *filename, char coeffs[][MAX_LINE], int max_coeffs)
{
	FILE *fp = fopen(filename, "r");
	char line[4096];
	int count = 0;

	if ( !fp ) { fprintf(stderr, "Cannot open reference file %s\n", filename); return -1; }

	while ( fgets(line, sizeof(line), fp) && count < max_coeffs ) {
		/* Strip trailing whitespace/newline */
		int len = strlen(line);
		while ( len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' || line[len-1] == ' ') )
			line[--len] = '\0';
		if ( len == 0 ) continue;
		if ( len >= MAX_LINE ) { fprintf(stderr, "Coefficient too long in %s\n", filename); fclose(fp); return -1; }
		strcpy(coeffs[count], line);
		count++;
	}
	fclose(fp);
	return count;
}

/*
    Compare two coefficient arrays. Returns 0 if identical, 1 if different.
*/
static int compare_coeffs (char a[][MAX_LINE], int na, char b[][MAX_LINE], int nb)
{
	if ( na != nb ) return 1;
	for ( int i = 0; i < na; i++ ) {
		if ( strcmp(a[i], b[i]) != 0 ) return 1;
	}
	return 0;
}

static int check (const char *label, char a[][MAX_LINE], int na, char b[][MAX_LINE], int nb)
{
	int fail = compare_coeffs(a, na, b, nb);
	printf("  %-26s %s\n", label, fail ? "FAIL" : "OK");
	return fail;
}

int main (int argc, char *argv[])
{
	const char *data_dir;
	int total_pass = 0, total_fail = 0;
	int i;

	if ( argc < 2 ) {
		fprintf(stderr, "Usage: %s <test_data_dir>\n", argv[0]);
		return 1;
	}
	data_dir = argv[1];

	/* Global init — same as classpoly_program.c */
	_cstd_seed = 42;
	mpz_util_init();
	dbg_level = 0;

	printf("\n=== classpoly test suite ===\n\n");

	for ( i = 0; i < (int)NUM_TESTS; i++ ) {
		struct test_vector *tv = &test_vectors[i];
		mpz_t P1, P2;
		char tmpdir[] = "/tmp/classpoly_test_XXXXXX";
		char serial_p1[512], serial_p2[512], par_p1[512], par_p2[512];
		char ref_p1_path[1024], ref_p2_path[1024];
		int failures = 0;

		/* Allocate coefficient arrays on heap (too large for stack with many tests) */
		char (*s_p1)[MAX_LINE] = malloc(MAX_COEFFS * MAX_LINE);
		char (*s_p2)[MAX_LINE] = malloc(MAX_COEFFS * MAX_LINE);
		char (*p_p1)[MAX_LINE] = malloc(MAX_COEFFS * MAX_LINE);
		char (*p_p2)[MAX_LINE] = malloc(MAX_COEFFS * MAX_LINE);
		char (*r_p1)[MAX_LINE] = malloc(MAX_COEFFS * MAX_LINE);
		char (*r_p2)[MAX_LINE] = malloc(MAX_COEFFS * MAX_LINE);
		int ns_p1, ns_p2, np_p1, np_p2, nr_p1, nr_p2;

		if ( !s_p1 || !s_p2 || !p_p1 || !p_p2 || !r_p1 || !r_p2 ) {
			fprintf(stderr, "malloc failed\n");
			return 1;
		}

		/* Create temp directory for output files */
		if ( !mkdtemp(tmpdir) ) { fprintf(stderr, "mkdtemp failed\n"); return 1; }

		snprintf(serial_p1, sizeof(serial_p1), "%s/serial_p1.txt", tmpdir);
		snprintf(serial_p2, sizeof(serial_p2), "%s/serial_p2.txt", tmpdir);
		snprintf(par_p1, sizeof(par_p1), "%s/par_p1.txt", tmpdir);
		snprintf(par_p2, sizeof(par_p2), "%s/par_p2.txt", tmpdir);

		/* Build reference file paths */
		snprintf(ref_p1_path, sizeof(ref_p1_path), "%s/D%ld_p%s.coeffs", data_dir, -tv->D, tv->P1);
		snprintf(ref_p2_path, sizeof(ref_p2_path), "%s/D%ld_p%s.coeffs", data_dir, -tv->D, tv->P2);

		/* Load reference coefficients */
		nr_p1 = load_reference_coeffs(ref_p1_path, r_p1, MAX_COEFFS);
		nr_p2 = load_reference_coeffs(ref_p2_path, r_p2, MAX_COEFFS);
		if ( nr_p1 < 0 || nr_p2 < 0 ) {
			printf("[D=%ld] SKIP — missing reference files\n\n", tv->D);
			free(s_p1); free(s_p2); free(p_p1); free(p_p2); free(r_p1); free(r_p2);
			continue;
		}

		printf("[D=%ld, h=%d]\n", tv->D, nr_p1 - 1);  /* h = coeffs - 1 (leading term is 1) */

		/* Initialize P values */
		mpz_init(P1); mpz_init(P2);
		mpz_set_str(P1, tv->P1, 10);
		mpz_set_str(P2, tv->P2, 10);

		/* Run serial for P1 */
		if ( !compute_classpoly(tv->D, 0, P1, serial_p1) ) {
			printf("  serial P1: FAILED to compute\n");
			failures += 6;
			goto next;
		}
		/* Run serial for P2 */
		if ( !compute_classpoly(tv->D, 0, P2, serial_p2) ) {
			printf("  serial P2: FAILED to compute\n");
			failures += 4;
			goto next;
		}
		/* Run parallel for both P values */
		if ( !compute_classpoly_parallel(tv->D, 0, P1, P2, par_p1, par_p2, 2) ) {
			printf("  parallel: FAILED to compute\n");
			failures += 4;
			goto next;
		}

		/* Parse output files */
		ns_p1 = parse_classpoly_output(serial_p1, s_p1, MAX_COEFFS);
		ns_p2 = parse_classpoly_output(serial_p2, s_p2, MAX_COEFFS);
		np_p1 = parse_classpoly_output(par_p1, p_p1, MAX_COEFFS);
		np_p2 = parse_classpoly_output(par_p2, p_p2, MAX_COEFFS);

		if ( ns_p1 < 0 || ns_p2 < 0 || np_p1 < 0 || np_p2 < 0 ) {
			printf("  FAILED to parse output files\n");
			failures += 6;
			goto next;
		}

		/* 6-way cross-validation */
		failures += check("serial P1 vs cache:",   s_p1, ns_p1, r_p1, nr_p1);
		failures += check("serial P2 vs cache:",   s_p2, ns_p2, r_p2, nr_p2);
		failures += check("parallel P1 vs cache:", p_p1, np_p1, r_p1, nr_p1);
		failures += check("parallel P2 vs cache:", p_p2, np_p2, r_p2, nr_p2);
		failures += check("serial vs parallel P1:", s_p1, ns_p1, p_p1, np_p1);
		failures += check("serial vs parallel P2:", s_p2, ns_p2, p_p2, np_p2);

next:
		mpz_clear(P1); mpz_clear(P2);

		total_pass += (6 - failures);
		total_fail += failures;
		printf("\n");

		/* Cleanup temp files */
		remove(serial_p1);
		remove(serial_p2);
		remove(par_p1);
		remove(par_p2);
		rmdir(tmpdir);

		free(s_p1); free(s_p2); free(p_p1); free(p_p2); free(r_p1); free(r_p2);
	}

	printf("=== Results: %d/%d passed, %d failed ===\n\n",
	       total_pass, total_pass + total_fail, total_fail);

	return total_fail ? 1 : 0;
}

/*
    Copyright 2012 Andrew V. Sutherland
    Copyright 2026 Brandon Lehmann

    This file is part of classpoly.

    classpoly is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 2 of the License.

    classpoly is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with classpoly.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <gmp.h>
#include "ff_poly.h"
#include "classpoly.h"
#include "cstd.h"

#define BIG_D		1000000
#define MAX_P_VALUES	16

struct cli_args {
	long D;
	int inv;
	int inv_set;
	mpz_t P[MAX_P_VALUES];
	int num_P;
	int nworkers;
	int verbosity;		/* -1 = auto */
};

static void print_help (void)
{
	printf("classpoly v%s — Compute Hilbert class polynomials via CRT\n\n", CLASSPOLY_VERSION_STRING);
	printf("Usage:\n");
	printf("  classpoly --D <discriminant> [--inv <invariant>] [--P <prime>]... [-j <workers>] [-v <level>]\n");
	printf("  classpoly -h | --help\n\n");

	printf("Required:\n");
	printf("  --D <value>        Imaginary quadratic discriminant (positive values auto-negated)\n\n");

	printf("Optional:\n");
	printf("  --inv <name|code>  Class invariant to use (default: auto-pick best for D)\n");
	printf("  --P <prime>        Prime modulus for reduction (repeatable; omit for computation over Z)\n");
	printf("  -j <N>             Number of parallel worker processes (requires at least one --P)\n");
	printf("  -v <level>         Verbosity level: -1=quiet, 0=normal, 1=info, 2=debug (default: auto)\n");
	printf("  -h, --help         Show this help message\n\n");

	printf("Invariant names (--inv accepts any of these, or the numeric code):\n");
	printf("  j, hilbert       Hilbert class polynomial (code 0)\n");
	printf("  f, weber         Weber f-function (code 1)\n");
	printf("  f2, f3, f4, f8   Powers of Weber f (codes 2-4, 8)\n");
	printf("  g2, gamma2       Cube root of j (code 5)\n");
	printf("  t, ramanujan     Ramanujan invariant (code 11)\n");
	printf("  t2, t6           Powers of t (codes 12-13)\n");
	printf("  w2w3e1, w3w3e1, ...  Double eta-quotients (codes 6, 9, etc.)\n");
	printf("  a<N>             Atkin invariant A_N, e.g. a71 (code 100+N)\n");
	printf("  w<N>             Single eta-quotient, e.g. w13 (code 400+N)\n");
	printf("  w<p1>w<p2>       Double eta-quotient, e.g. w2w3 (code 500+p1*p2)\n\n");

	printf("P expressions:\n");
	printf("  --P supports arithmetic expressions: 2^127-1, 2^255-19, M31 (Mersenne), etc.\n\n");

	printf("Output files:\n");
	printf("  H_D<|D|>.coeffs           when computing over Z\n");
	printf("  H_D<|D|>_p<P>.coeffs      when reducing mod P (one file per --P)\n\n");

	printf("Examples:\n");
	printf("  classpoly --D 167995\n");
	printf("  classpoly --D 167995 --inv weber --P 2^127-1\n");
	printf("  classpoly --D 167995 --P 2^127-1 --P 2^255-19\n");
	printf("  classpoly --D 167995 --inv 0 --P 2^127-1 --P 2^255-19 -j 4\n\n");

	printf("Version: classpoly %s (ff_poly %s)\n", CLASSPOLY_VERSION_STRING, FF_POLY_VERSION_STRING);
	printf("See https://arxiv.org/abs/1001.3394 for algorithm details and constraints on D.\n");
}

static int parse_args (struct cli_args *args, int argc, char *argv[])
{
	int i;

	args->D = 0;
	args->inv = -1;
	args->inv_set = 0;
	args->num_P = 0;
	args->nworkers = 0;
	args->verbosity = -1;
	for ( i = 0; i < MAX_P_VALUES; i++ ) mpz_init(args->P[i]);

	for ( i = 1; i < argc; i++ ) {
		if ( strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0 ) {
			print_help();
			exit(0);
		} else if ( strcmp(argv[i], "--D") == 0 ) {
			if ( ++i >= argc ) { err_printf("--D requires a value\n"); return 0; }
			args->D = atol(argv[i]);
			if ( args->D > 0 ) args->D = -args->D;
		} else if ( strcmp(argv[i], "--inv") == 0 ) {
			if ( ++i >= argc ) { err_printf("--inv requires a value\n"); return 0; }
			args->inv = inv_parse_name(argv[i]);
			if ( args->inv == -2 ) { err_printf("Unknown invariant: '%s'\n", argv[i]); return 0; }
			args->inv_set = 1;
		} else if ( strcmp(argv[i], "--P") == 0 ) {
			if ( ++i >= argc ) { err_printf("--P requires a value\n"); return 0; }
			if ( args->num_P >= MAX_P_VALUES ) { err_printf("Too many --P values (max %d)\n", MAX_P_VALUES); return 0; }
			mpz_eval_expr(args->P[args->num_P], argv[i]);
			args->num_P++;
		} else if ( strcmp(argv[i], "-j") == 0 ) {
			if ( ++i >= argc ) { err_printf("-j requires a value\n"); return 0; }
			args->nworkers = atoi(argv[i]);
			if ( args->nworkers < 1 ) { err_printf("-j requires a positive integer\n"); return 0; }
		} else if ( strcmp(argv[i], "-v") == 0 ) {
			if ( ++i >= argc ) { err_printf("-v requires a value\n"); return 0; }
			args->verbosity = atoi(argv[i]);
		} else {
			err_printf("Unknown option: '%s'\n", argv[i]);
			return 0;
		}
	}

	if ( !args->D ) {
		err_printf("Missing required option: --D <discriminant>\n");
		return 0;
	}

	return 1;
}

static void cli_args_clear (struct cli_args *args)
{
	int i;
	for ( i = 0; i < MAX_P_VALUES; i++ ) mpz_clear(args->P[i]);
}

/*
    Generate a deterministic output filename for a given D and P.
    If P is NULL or zero, produces "H_D<|D|>.coeffs".
    Otherwise produces "H_D<|D|>_p<P>.coeffs".
    Caller must free() the returned string.
*/
static char *make_filename (long D, mpz_t P)
{
	char *buf;
	if ( P && mpz_sgn(P) ) {
		size_t needed = mpz_sizeinbase(P, 10) + 32;
		buf = malloc(needed);
		gmp_snprintf(buf, needed, "H_D%ld_p%Zd.coeffs", -D, P);
	} else {
		buf = malloc(64);
		snprintf(buf, 64, "H_D%ld.coeffs", -D);
	}
	return buf;
}

int main (int argc, char *argv[])
{
	struct cli_args args;
	int inv, picked, ret;
	char *filename = NULL;
	char **par_filenames = NULL;
	int i;

	if ( argc < 2 ) { print_help(); return 0; }

	if ( !parse_args(&args, argc, argv) ) {
		printf("\nRun 'classpoly --help' for usage information.\n");
		cli_args_clear(&args);
		return -1;
	}

	_cstd_seed = 42;
	mpz_util_init();
	ret = 0;

	if ( !discriminant_test(args.D) ) {
		err_printf("%ld is not a valid discriminant.\n", args.D);
		ret = -1; goto cleanup;
	}

	picked = 0;
	if ( args.inv_set ) {
		inv = args.inv;
		if ( !inv_good_invariant(inv) ) { err_printf("Unknown or invalid invariant %d specified\n", inv); ret = -1; goto cleanup; }
		if ( !inv_enum(inv) ) { err_printf("Unsupported invariant %d specified\n", inv); ret = -1; goto cleanup; }
		if ( !inv_good_discriminant(args.D, inv) ) { err_printf("D=%ld is not a good discriminant for invariant %s (%d)\n", args.D, inv_string(inv), inv); ret = -1; goto cleanup; }
	} else {
		if ( args.D >= -4 ) inv = 0; else { inv = inv_pick_invariant(args.D); picked = 1; }
		printf("Using invariant %s (%d) with height factor %.1f\n", inv_string(inv), inv, inv_height_factor(inv));
	}

	dbg_level = ( args.verbosity >= 0 ) ? args.verbosity : ( -args.D > BIG_D ? 1 : 0 );

	if ( args.nworkers > 0 ) {
		/* ==================== PARALLEL MODE ==================== */
		if ( args.num_P < 1 ) {
			err_printf("Parallel mode (-j) requires at least one --P value\n");
			ret = -1; goto cleanup;
		}
		if ( args.num_P > MAX_P_VALUES ) {
			err_printf("Too many --P values for parallel mode (got %d, max %d)\n", args.num_P, MAX_P_VALUES);
			ret = -1; goto cleanup;
		}

		par_filenames = malloc(args.num_P * sizeof(char *));
		for ( i = 0; i < args.num_P; i++ )
			par_filenames[i] = make_filename(args.D, args.P[i]);

		if ( !compute_classpoly_parallel(args.D, inv, args.P, par_filenames, args.num_P, args.nworkers) ) {
			if ( !picked || !inv ) goto cleanup;
			err_printf("Unable to compute classpoly using invariant %s (%d), reverting to j-invariant\n", inv_string(inv), inv);
			if ( !compute_classpoly_parallel(args.D, 0, args.P, par_filenames, args.num_P, args.nworkers) ) {
				err_printf("Error computing Hilbert class polynomial!\n");
				ret = -1; goto cleanup;
			}
		}
	} else {
		/* ==================== SERIAL MODE ==================== */
		if ( args.num_P == 0 ) {
			mpz_t P0;
			mpz_init(P0);
			filename = make_filename(args.D, NULL);
			if ( !compute_classpoly(args.D, inv, P0, filename) ) {
				if ( !picked || !inv ) { mpz_clear(P0); goto cleanup; }
				err_printf("Unable to compute classpoly using invariant %s (%d), reverting to j-invariant\n", inv_string(inv), inv);
				if ( !compute_classpoly(args.D, 0, P0, filename) ) { err_printf("Error computing Hilbert class polynomial!\n"); mpz_clear(P0); ret = -1; goto cleanup; }
			}
			mpz_clear(P0);
		} else {
			for ( i = 0; i < args.num_P; i++ ) {
				free(filename);
				filename = make_filename(args.D, args.P[i]);
				if ( !compute_classpoly(args.D, inv, args.P[i], filename) ) {
					if ( !picked || !inv ) goto cleanup;
					err_printf("Unable to compute classpoly using invariant %s (%d), reverting to j-invariant\n", inv_string(inv), inv);
					if ( !compute_classpoly(args.D, 0, args.P[i], filename) ) { err_printf("Error computing Hilbert class polynomial!\n"); ret = -1; goto cleanup; }
				}
			}
		}
	}

cleanup:
	free(filename);
	if ( par_filenames ) {
		for ( i = 0; i < args.num_P; i++ ) free(par_filenames[i]);
		free(par_filenames);
	}
	cli_args_clear(&args);
	return ret;
}

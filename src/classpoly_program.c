/*
    Copyright 2012 Andrew V. Sutherland

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

#define BIG_D	1000000

int main (int argc, char *argv[])
{
	long D;
	int inv;
	mpz_t P;
	char *filename;
	int picked;
	int nworkers = 0;
	mpz_t P2;
	char *filename2 = NULL;
	int i;

	/* Scan for -j flag before positional arg parsing */
	for ( i = 1; i < argc; i++ ) {
		if ( strcmp(argv[i], "-j") == 0 && i + 1 < argc ) {
			nworkers = atoi(argv[i + 1]);
			if ( nworkers < 1 ) { err_printf("-j requires a positive integer\n"); return -1; }
			/* Shift remaining args down to remove -j N */
			int k;
			for ( k = i; k < argc - 2; k++ ) argv[k] = argv[k + 2];
			argc -= 2;
			break;
		}
	}

	if ( argc < 2 ) {
		printf ("classpoly D [inv P filename verbosity]\n");
		printf ("classpoly -j N D inv P1 [P2 filename1 filename2 verbosity]\n\n");
		printf ("See README.txt for the possible values of inv (and see https://arxiv.org/abs/1001.3394 for explanations and constraints on D)\n");
		printf ("The default (inv=-1) is to pick the best class invariant for the given D, use inv=0 to force j to be used (Hilbert class polynomial)\n");
		printf ("The default P=0 computes H^inv_D(X) over Z, set P to a prime to get the reduction modulo P (slightly faster and uses less memory)\n");
		printf ("-j N: use N parallel worker processes (requires P to be set)\n");
		printf ("Version %s (ff_poly version %s)\n", CLASSPOLY_VERSION_STRING, FF_POLY_VERSION_STRING);
		return 0;
	}

	_cstd_seed = 42;

	mpz_util_init();

	D = atol(argv[1]);  if ( D > 0 ) D = -D;
	if ( ! discriminant_test(D) ) { err_printf ("%ld is not a valid discriminant.\n", D); return -1; }

	picked = 0;
	if ( argc > 2 ) {
		inv = atoi(argv[2]);
		if ( inv == - 1 ) { inv = inv_pick_invariant(D); picked = 1; info_printf ("Using invariant %s (%d)\n", inv_string(inv), inv); }
		if ( ! inv_good_invariant(inv) ) { err_printf ("Unknown or invalid invariant %d specified\n", inv); return -1; }
		if ( ! inv_enum(inv) ) { err_printf ("Unsupported invariant %d specified\n", inv); return -1; }
		if ( ! inv_good_discriminant (D, inv) ) { err_printf ("D=%ld is not a good discriminant for invariant %s (%d)\n", D, inv_string(inv), inv); return -1; }
	} else {
		if ( D >= -4 ) inv = 0; else { inv = inv_pick_invariant(D); picked = 1; }
		printf ("Using invariant %s (%d) with height factor %.1f\n", inv_string(inv), inv, inv_height_factor(inv));
	}

	if ( nworkers > 0 ) {
		/*
		    Parallel mode: classpoly -j N D inv P1 [P2 filename1 filename2 verbosity]
		    After -j N is stripped, positional args are:
		      argv[1]=D  argv[2]=inv  argv[3]=P1  [argv[4]=P2  argv[5]=filename1  argv[6]=filename2  argv[7]=verbosity]
		*/
		mpz_init(P);
		mpz_init(P2);

		if ( argc > 3 ) mpz_eval_expr(P, argv[3]);
		if ( !mpz_sgn(P) ) { err_printf("Parallel mode requires P1 to be set\n"); return -1; }

		if ( argc > 4 ) mpz_eval_expr(P2, argv[4]);
		filename = ( argc > 5 ? argv[5] : NULL );
		if ( filename && strlen(filename) == 0 ) filename = NULL;
		filename2 = ( argc > 6 ? argv[6] : NULL );
		if ( filename2 && strlen(filename2) == 0 ) filename2 = NULL;
		if ( argc > 7 ) { dbg_level = atoi(argv[7]); } else { dbg_level = ( -D > BIG_D ? 1 : 0 ); }

		if ( !compute_classpoly_parallel(D, inv, P, P2, filename, filename2, nworkers) ) {
			if ( !picked || !inv ) { mpz_clear(P); mpz_clear(P2); return 0; }
			err_printf("Unable to compute classpoly using invariant %s (%d), reverting to j-invariant\n", inv_string(inv), inv);
			if ( !compute_classpoly_parallel(D, 0, P, P2, filename, filename2, nworkers) ) {
				err_printf("Error computing Hilbert class polynomial!\n");
				mpz_clear(P); mpz_clear(P2);
				return -1;
			}
		}
		mpz_clear(P);
		mpz_clear(P2);
	} else {
		/* Original serial mode */
		if ( argc > 5 ) { dbg_level = atoi (argv[5]); } else { dbg_level = ( -D > BIG_D ? 1 : 0 ); }

		mpz_init (P);
		if ( argc > 3 )  mpz_eval_expr (P, argv[3]);
		filename = ( argc > 4 ? argv[4] : 0 );
		if ( filename && strlen(filename) == 0 ) filename = 0;
		if ( ! compute_classpoly (D, inv, P, filename) ) {
			if ( ! picked || ! inv ) return 0;
			err_printf ("Unable to compute classpoly using invariant %s (%d), reverting to j-invariant\n", inv_string(inv), inv);
			if ( ! compute_classpoly (D, 0, P, filename) ) { err_printf ("Error computing Hilbert class polynomial!\n"); return-1; }
		}
		mpz_clear (P);
	}
	return 0;
}

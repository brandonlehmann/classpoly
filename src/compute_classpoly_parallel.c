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
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/mman.h>
#include <unistd.h>
#include <gmp.h>
#include "ff_poly.h"
#include "phi_poly.h"
#include "class_inv.h"
#include "classpoly.h"
#include "qform.h"
#include "pickprimes.h"
#include "classpoly_crt.h"
#include "classpoly_inv.h"
#include "crt.h"
#include "cstd.h"

#define DELTA_PERCENT		0.05
#define HEIGHT_MARGIN		256
#define PARALLEL_BATCH_SIZE	64

/* Shared state between parent and children, allocated via mmap MAP_SHARED */
struct parallel_shared {
	volatile int counter;		/* atomic work counter: next prime index to claim */
	volatile int primes_done;	/* atomic: total primes completed across all workers */
	int pcnt;			/* total prime count */
	int batch_size;			/* primes per batch */
};

/*
    Atomically claim a batch of primes.
    Returns the start index, or a value >= pcnt if no work remains.
*/
static inline int claim_batch (struct parallel_shared *sh)
{
	return __sync_fetch_and_add(&sh->counter, sh->batch_size);
}

/*
    Write one ECRT output file: header + coefficients.
    This is the same format as classpoly_crt_finish for the ECRT path.
*/
static int write_ecrt_output (ecrt_context_t ecrt, int inv, long D, mpz_t P, int ccnt, const char *filename)
{
	mpz_t X;
	FILE *fp;
	int i;

	ecrt_finalize(ecrt);
	fp = fopen(filename, "w");
	if ( !fp ) { err_printf("Error opening output file %s\n", filename); return 0; }

	gmp_fprintf(fp, "I=%d\n", inv);
	gmp_fprintf(fp, "D=%ld\n", D);
	gmp_fprintf(fp, "P=%Zd\n", P);

	mpz_init(X);
	for ( i = 0; ecrt_next_coeff(X, ecrt); i++ )
		gmp_fprintf(fp, "%Zd*X^%d + \n", X, i);
	mpz_clear(X);

	if ( i != ccnt ) {
		err_printf("ECRT postcomputation returned only %d of %d expected coefficients\n", i, ccnt);
		fclose(fp);
		return 0;
	}
	gmp_fprintf(fp, "1*X^%d\n", ccnt);
	fclose(fp);
	return 1;
}

/*
    Merge N worker accumulator buffers into the parent ecrt context.
    Each worker's Cdata/sdata is at bufs[w] offset by (k * Climbs) and (k * slimbs).
*/
static void merge_buffers (ecrt_context_t ecrt, mp_limb_t **Cbufs, mp_limb_t **Sbufs, int nworkers)
{
	int j, w;
	int k = ecrt->k;
	int Climbs = ecrt->Climbs;
	int slimbs = ecrt->slimbs;

	/* Zero parent accumulators */
	memset(ecrt->Cdata, 0, (size_t)k * Climbs * sizeof(mp_limb_t));
	memset(ecrt->sdata, 0, (size_t)k * slimbs * sizeof(mp_limb_t));

	for ( w = 0; w < nworkers; w++ ) {
		for ( j = 0; j < k; j++ ) {
			mpn_add_n(ecrt->Cdata + j * Climbs,
				  ecrt->Cdata + j * Climbs,
				  Cbufs[w] + j * Climbs,
				  Climbs);
			mpn_add_n(ecrt->sdata + j * slimbs,
				  ecrt->sdata + j * slimbs,
				  Sbufs[w] + j * slimbs,
				  slimbs);
		}
	}
}

int compute_classpoly_parallel (long D, int inv, mpz_t P1, mpz_t P2, char *filename1, char *filename2, int nworkers)
{
	time_t begin, start, end;
	char buf1[256], buf2[256];
	classpoly_t H;
	classpoly_inv_t Hinv;
	classpoly_crt_t crt;
	ecrt_context_t ecrt2;
	double bits, height_factor;
	long vfilter, maxp;
	double tbits;
	long *crt_p;
	int crt_pcnt, p1, p2;
	int H_d, ell0, N, classno;
	int have_P2;
	register int i;

	/* Shared memory regions */
	struct parallel_shared *shared;
	mp_limb_t **Cbufs1 = NULL, **Sbufs1 = NULL;
	mp_limb_t **Cbufs2 = NULL, **Sbufs2 = NULL;
	size_t Csize, Ssize;
	pid_t *pids;
	int w;

	have_P2 = (P2 && mpz_sgn(P2));
	if ( !P1 || !mpz_sgn(P1) ) { err_printf("P1 is required for parallel mode\n"); return 0; }

	if ( !filename1 ) { filename1 = buf1; sprintf(filename1, "H_%ld_P1.txt", -D); }
	if ( have_P2 && !filename2 ) { filename2 = buf2; sprintf(filename2, "H_%ld_P2.txt", -D); }

	begin = clock();

	/* ==================== PRECOMPUTATION (parent only) ==================== */

	phi_poly_setup();
	classpoly_init(H, D, inv);
	classpoly_setup_find_jinv(H);

	N = inv_level(inv);
	ell0 = ( (inv_atkin(inv) || inv_double_eta(inv)) && inv_ramified(D, inv) ) ? inv_degree(&p1, &p2, inv) : 0;
	if ( ell0 ) { info_printf("Computing square-root of classpoly for D=%ld\n", D); }
	else { info_printf("Computing classpoly for D=%ld\n", D); }

	if ( !classpoly_setup_enum_roots(H, N, ell0, inv) ) return 0;
	if ( dbg_level >= INFO_LEVEL ) classgroup_pcp_print(H->pres);
	classno = H->pres->h;
	H_d = (int) classno / (ell0 ? 2 : 1);
	height_factor = inv_enum_height_factor(inv);
	if ( !height_factor ) { fprintf(stderr, "height factor error for D=%ld\n", D); abort(); }
	bits = qform_hilbert_height_bound(D, H_d) / height_factor;
	if ( inv_single_eta(inv) ) bits = 1.1 * bits + 200;
	bits += HEIGHT_MARGIN;
	info_printf("bits = %.2f (classpoly degree %d height factor %.2f, safety margin %d bits)\n", bits, H_d, height_factor, HEIGHT_MARGIN);
	vfilter = classgroup_pcp_vfilter(H->pres);

	/* Pick primes */
	crt_pcnt = pick_primes(&crt_p, 0, 0, D, classno, bits, inv_pfilter(inv), vfilter);
	for ( i = 0, maxp = 0, tbits = 0.0; i < crt_pcnt; i++ ) {
		if ( crt_p[i] > maxp ) maxp = crt_p[i];
		tbits += log2(crt_p[i]);
	}
	info_printf("Selected %d primes, maxp=%ld, totbits = %.1f\n", crt_pcnt, maxp, tbits);

	/* CRT setup with P1 — this calls ecrt_init for P1 inside */
	classpoly_crt_start(crt, D, inv, (unsigned long *)crt_p, crt_pcnt, H_d, P1);

	/* classpoly_inv_setup uses the crt to consume a few primes for trace computation */
	if ( !classpoly_inv_setup(Hinv, H, inv, crt) ) {
		classpoly_clear(H);
		classpoly_crt_end(crt);
		return 0;
	}

	/* Initialize second ecrt for P2 if needed */
	if ( have_P2 ) {
		char prefix2[64];
		sprintf(prefix2, "H_%ld_2", -D);
		ecrt_init(ecrt2, (unsigned long *)crt_p, crt_pcnt, H_d, P2, 0, 0, prefix2);
	}

	qform_table_free();
	end = clock();

	/* Print run summary (always visible) */
	out_printf("\n=== Parallel classpoly run summary ===\n");
	out_printf("  D             = %ld\n", D);
	out_printf("  invariant     = %s (%d)\n", inv_string(inv), inv);
	out_printf("  class number  = %ld\n", (long)classno);
	out_printf("  poly degree   = %d%s\n", H_d, ell0 ? " (sqrt)" : "");
	out_printf("  height bound  = %.1f bits (margin %d)\n", bits, HEIGHT_MARGIN);
	out_printf("  CRT primes    = %d (%.1f bits total, maxp=%ld)\n", crt_pcnt, tbits, maxp);
	out_printf("  P1            = "); gmp_printf("%Zd", P1); out_printf(" (%ld bits)\n", mpz_sizeinbase(P1, 2));
	if ( have_P2 ) { out_printf("  P2            = "); gmp_printf("%Zd", P2); out_printf(" (%ld bits)\n", mpz_sizeinbase(P2, 2)); }
	out_printf("  workers       = %d\n", nworkers);
	out_printf("  batch size    = %d primes\n", PARALLEL_BATCH_SIZE);
	out_printf("  primes/worker = ~%d (avg)\n", crt_pcnt / nworkers);
	out_printf("  precomp time  = %ld msecs\n", delta_msecs(begin, end));
	out_printf("======================================\n\n");

	/* ==================== SHARED MEMORY SETUP ==================== */

	/* Atomic work counter */
	shared = mmap(NULL, sizeof(*shared), PROT_READ | PROT_WRITE,
		      MAP_SHARED | MAP_ANONYMOUS, -1, 0);
	if ( shared == MAP_FAILED ) { err_printf("mmap failed for shared counter\n"); abort(); }
	shared->counter = 0;
	shared->primes_done = 0;
	shared->pcnt = crt_pcnt;
	shared->batch_size = PARALLEL_BATCH_SIZE;

	/* Per-worker accumulator buffers */
	Csize = (size_t)H_d * crt->ecrt->Climbs * sizeof(mp_limb_t);
	Ssize = (size_t)H_d * crt->ecrt->slimbs * sizeof(mp_limb_t);

	Cbufs1 = malloc(nworkers * sizeof(mp_limb_t *));
	Sbufs1 = malloc(nworkers * sizeof(mp_limb_t *));
	if ( have_P2 ) {
		Cbufs2 = malloc(nworkers * sizeof(mp_limb_t *));
		Sbufs2 = malloc(nworkers * sizeof(mp_limb_t *));
	}

	for ( w = 0; w < nworkers; w++ ) {
		Cbufs1[w] = mmap(NULL, Csize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
		Sbufs1[w] = mmap(NULL, Ssize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
		if ( Cbufs1[w] == MAP_FAILED || Sbufs1[w] == MAP_FAILED ) { err_printf("mmap failed for worker %d P1 buffers\n", w); abort(); }
		memset(Cbufs1[w], 0, Csize);
		memset(Sbufs1[w], 0, Ssize);

		if ( have_P2 ) {
			Cbufs2[w] = mmap(NULL, Csize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
			Sbufs2[w] = mmap(NULL, Ssize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
			if ( Cbufs2[w] == MAP_FAILED || Sbufs2[w] == MAP_FAILED ) { err_printf("mmap failed for worker %d P2 buffers\n", w); abort(); }
			memset(Cbufs2[w], 0, Csize);
			memset(Sbufs2[w], 0, Ssize);
		}
	}

	/* ==================== FORK WORKERS ==================== */

	fflush(stdout);
	fflush(stderr);

	pids = malloc(nworkers * sizeof(pid_t));
	start = clock();

	for ( w = 0; w < nworkers; w++ ) {
		pids[w] = fork();
		if ( pids[w] < 0 ) { err_printf("fork() failed for worker %d\n", w); abort(); }

		if ( pids[w] == 0 ) {
			/* ============ CHILD PROCESS ============ */
			ff_t *work, J;
			int idx, batch_end;
			int my_start;
			long primes_done = 0;

			/* Point ecrt accumulators at our shared memory buffers */
			/* Free the COW-inherited Cdata/sdata first */
			free(crt->ecrt->Cdata);
			free(crt->ecrt->sdata);
			crt->ecrt->Cdata = Cbufs1[w];
			crt->ecrt->sdata = Sbufs1[w];

			if ( have_P2 ) {
				free(ecrt2->Cdata);
				free(ecrt2->sdata);
				ecrt2->Cdata = Cbufs2[w];
				ecrt2->sdata = Sbufs2[w];
			}

			work = mem_alloc(2 * H_d * sizeof(ff_t));

			while ( (my_start = claim_batch(shared)) < crt_pcnt ) {
				batch_end = my_start + shared->batch_size;
				if ( batch_end > crt_pcnt ) batch_end = crt_pcnt;

				for ( idx = my_start; idx < batch_end; idx++ ) {
					ff_setup_ui(crt_p[idx]);

					if ( !ff_classpoly_find_jinv(&J, H) ) {
						err_printf("[W%02d] ff_classpoly_find_jinv failed for p=%ld\n", w, _ff_p);
						_exit(1);
					}

					ff_classpoly_setup(H);
					ff_classpoly_inv_from_j(&J, &J, Hinv, H);

					if ( !ff_classpoly_enum_roots(H, J, inv) ) {
						err_printf("[W%02d] ff_classpoly_enum_roots failed for p=%ld\n", w, _ff_p);
						_exit(1);
					}

					if ( !ff_classpoly_inv_disambiguate(H, Hinv) ) {
						err_printf("[W%02d] ff_classpoly_inv_disambiguate failed for p=%ld, aborting\n", w, _ff_p);
						_exit(1);
					}

					ff_poly_from_roots(work, H->roots, H_d);
					for ( i = 0; i < H_d; i++ ) work[i] = _ff_get_ui(work[i]);

					ecrt_update(crt->ecrt, idx, (unsigned long *)work, H_d);
					if ( have_P2 )
						ecrt_update(ecrt2, idx, (unsigned long *)work, H_d);

					primes_done++;
					__sync_fetch_and_add(&shared->primes_done, 1);
				}
			}

			/* Don't free mmap buffers — parent reads them */
			_exit(0);
		}
	}

	/* ==================== PARENT: PROGRESS + WAIT ==================== */

	int failures = 0;
	int workers_alive = nworkers;
	{
		struct timespec ts_start, ts_now;
		clock_gettime(CLOCK_MONOTONIC, &ts_start);

		while ( workers_alive > 0 ) {
			/* Sleep 10 seconds between progress updates */
			sleep(10);

			/* Reap any finished children (non-blocking) */
			for ( w = 0; w < nworkers; w++ ) {
				if ( pids[w] == 0 ) continue; /* already reaped */
				int status;
				pid_t ret = waitpid(pids[w], &status, WNOHANG);
				if ( ret > 0 ) {
					if ( !WIFEXITED(status) || WEXITSTATUS(status) != 0 ) {
						err_printf("Worker %d failed (status %d)\n", w, status);
						failures++;
					}
					pids[w] = 0;
					workers_alive--;
				}
			}

			/* Print progress */
			clock_gettime(CLOCK_MONOTONIC, &ts_now);
			double elapsed = (ts_now.tv_sec - ts_start.tv_sec) + (ts_now.tv_nsec - ts_start.tv_nsec) / 1e9;
			int done = shared->primes_done;
			if ( done > crt_pcnt ) done = crt_pcnt;
			double pct = 100.0 * done / crt_pcnt;
			double rate = (elapsed > 0) ? done / elapsed : 0;
			int remaining = crt_pcnt - done;
			double eta = (rate > 0) ? remaining / rate : 0;
			int eta_h = (int)(eta / 3600);
			int eta_m = (int)((eta - eta_h * 3600) / 60);
			int eta_s = (int)(eta - eta_h * 3600 - eta_m * 60);
			int el_h = (int)(elapsed / 3600);
			int el_m = (int)((elapsed - el_h * 3600) / 60);
			int el_s = (int)(elapsed - el_h * 3600 - el_m * 60);
			out_printf("  [%6.2f%%] %d/%d primes  %.0f primes/s  elapsed %dh%02dm%02ds  ETA %dh%02dm%02ds  [%d workers]\n",
				   pct, done, crt_pcnt, rate, el_h, el_m, el_s, eta_h, eta_m, eta_s, workers_alive);
			fflush(stdout);
		}
	}
	free(pids);

	end = clock();
	info_printf("All workers finished in %ld msecs, %d failures\n", delta_msecs(start, end), failures);

	if ( failures ) {
		err_printf("Aborting due to worker failures\n");
		goto cleanup;
	}

	/* Merge P1 */
	info_printf("Merging %d worker buffers for P1...\n", nworkers);
	merge_buffers(crt->ecrt, Cbufs1, Sbufs1, nworkers);

	/* Write P1 output */
	start = clock();
	if ( !write_ecrt_output(crt->ecrt, inv, D, P1, H_d, filename1) ) {
		err_printf("Failed to write P1 output\n");
		failures = 1;
		goto cleanup;
	}
	end = clock();
	info_printf("P1 finalization and output in %ld msecs\n", delta_msecs(start, end));
	out_printf("Class polynomial for inv=%d, D=%ld reduced mod P1 (%ld bits) written to %s, degree %d\n",
		   inv, D, mpz_sizeinbase(P1, 2), filename1, H_d);

	/* Merge and write P2 if present */
	if ( have_P2 ) {
		info_printf("Merging %d worker buffers for P2...\n", nworkers);
		merge_buffers(ecrt2, Cbufs2, Sbufs2, nworkers);

		start = clock();
		if ( !write_ecrt_output(ecrt2, inv, D, P2, H_d, filename2) ) {
			err_printf("Failed to write P2 output\n");
			failures = 1;
			goto cleanup;
		}
		end = clock();
		info_printf("P2 finalization and output in %ld msecs\n", delta_msecs(start, end));
		out_printf("Class polynomial for inv=%d, D=%ld reduced mod P2 (%ld bits) written to %s, degree %d\n",
			   inv, D, mpz_sizeinbase(P2, 2), filename2, H_d);
	}

cleanup:
	/* Unmap shared memory */
	munmap(shared, sizeof(*shared));
	for ( w = 0; w < nworkers; w++ ) {
		munmap(Cbufs1[w], Csize);
		munmap(Sbufs1[w], Ssize);
		if ( have_P2 ) {
			munmap(Cbufs2[w], Csize);
			munmap(Sbufs2[w], Ssize);
		}
	}
	free(Cbufs1); free(Sbufs1);
	if ( have_P2 ) { free(Cbufs2); free(Sbufs2); }

	classpoly_clear(H);
	classpoly_inv_clear(Hinv);
	if ( have_P2 ) ecrt_clear(ecrt2);

	end = clock();
	out_printf("Total parallel computation time: %ld msecs\n", delta_msecs(begin, end));

	return failures ? 0 : 1;
}

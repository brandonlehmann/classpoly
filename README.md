# classpoly v1.0.3 — Parallel Fork

Computes Hilbert class polynomials H_D(X) using the Chinese Remainder Theorem, with parallelized multi-process computation via `fork()`.

Based on the algorithms described in:

> **\[Sutherland 2011\]** Andrew V. Sutherland, "Computing Hilbert class polynomials with the Chinese Remainder Theorem", Math. Comp. 80 (2011), 501-538.

> **\[Enge-Sutherland 2010\]** Andreas Enge and Andrew V. Sutherland, "Class invariants by the CRT method", ANTS IX, LNCS 6197 (2010), 142-156.

## Prerequisites

- **GMP** (version 6 or later) — <https://gmplib.org> (on Ubuntu: `libgmp-dev`)
- **Modular polynomials** — download at least `phi_j.tar` from <https://math.mit.edu/~drew/SmallModPolys.html> and extract to `$HOME/phi_files/`. For class invariants beyond j, download all from <https://math.mit.edu/~drew/phi_polys.tar>
- **64-bit OS** (Linux or macOS, x86_64 or aarch64/Apple Silicon)

**Important:**
- The `$HOME/phi_files` directory must contain the modular polynomials before running. To change this location, modify `phi_dir()` in `phi_poly.h`.
- The `$HOME/temp` directory must exist (used for intermediate CRT files when computing over Z).

## Building

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

## Usage

```bash
# Serial (original)
./classpoly D [inv P filename verbosity]

# Parallel with single P
./classpoly -j N D inv P1 0 filename1 "" [verbosity]

# Parallel with two P values
./classpoly -j N D inv P1 P2 filename1 filename2 [verbosity]
```

### Parameters

| Parameter | Description |
|-----------|-------------|
| `D` | Imaginary quadratic discriminant (negative; positive values are auto-negated) |
| `inv` | Class invariant identifier (see table below; -1 = auto-pick best, 0 = j-invariant) |
| `P` | Prime modulus (0 or omitted = compute over Z) |
| `filename` | Output file (default: `H_{-D}.txt`) |
| `verbosity` | 0 = normal, 1 or 2 = more detail, -1 = suppress |
| `-j N` | Number of parallel worker processes (use physical core count; requires P) |

### Class invariants

| inv | Invariant | Notes |
|-----|-----------|-------|
| 0 | j | Hilbert class polynomial (Sutherland 2011) |
| 1 | f (Weber function) | See Enge-Sutherland 2010, sec. 3 |
| 2 | f^2 | |
| 5 | gamma\_2 (cube-root of j) | |
| 6 | w\_{2,3} (double eta-quotient) | See Enge-Sutherland 2010, sec. 3 |
| 9 | w\_{3,3} | |
| 10 | w\_{2,5} | |
| 11 | t (Ramanujan-related) | See Enge-Sutherland 2010, sec. 4.4 |
| 12 | t^2 | |
| 14 | w\_{2,7} | |
| 15 | w\_{3,5} | |
| 21 | w\_{3,7} | |
| 23 | w\_{2,3}^2 | |
| 24 | w\_{2,5}^2 | |
| 26 | w\_{2,13} | |
| 27 | w\_{2,7}^2 | |
| 28 | w\_{3,3}^2 | |
| 100+N | A\_N (Atkin, N=3,5,7,11,13,17,19,23,29,31,41,47,59,71) | |
| 400+N | w\_N^s single-eta (N=3,5,7,13; s=24/gcd(24,N-1)) | |
| 500+p1\*p2 | w\_{p1,p2}^s double-eta | Pairs: (2,3),(2,5),(2,7),(2,13),(3,3),(3,5),(3,7),(3,13),(5,7) |

## Testing

```bash
cd build && cmake .. && make test_classpoly
./test_classpoly ../tests/data
# or: ctest --verbose
```

The test harness runs 7 discriminants through both serial and parallel paths, cross-validating 6 ways against known-good cached coefficients. 42 checks total.

## Parallel Architecture

- **fork()** gives each worker its own copy of all mutable ff\_poly globals via COW — zero changes to the finite field library
- **Atomic work queue** (`__sync_fetch_and_add`) dynamically balances load across workers — no static chunking
- **Shared memory** (`mmap MAP_SHARED`) accumulators — zero disk I/O
- **Single `ecrt_init`** in the parent before fork — workers inherit via COW, no redundant CRT precomputation
- Per-prime work (find j-invariant, enumerate roots, build polynomial) done once, then `ecrt_update` called per P value

## Benchmark Results (15 cores, cross-validated against known-good cached coefficients)

| D | h(D) | CRT primes | Serial | Parallel (15-core) | Speedup | All checks |
|---|------|------------|--------|-------------------|---------|------------|
| -167995 | 57 | 222 | 0.1s | 9ms | 13.9x | 6/6 OK |
| -7857907 | 293 | 1,103 | 3.5s | 19ms | 183.1x | 6/6 OK |
| -4059339 | 389 | 1,333 | 2.5s | 27ms | 90.9x | 6/6 OK |
| -14141835 | 589 | 1,978 | 10.6s | 28ms | 376.9x | 6/6 OK |
| -54558307 | 781 | 2,964 | 26.2s | 57ms | 460.0x | 6/6 OK |
| -137469067 | 1,198 | 4,399 | 65.6s | 76ms | 862.9x | 6/6 OK |
| -64752003 | 1,397 | 4,563 | 41.9s | 65ms | 644.4x | 6/6 OK |
| -268755307 | 1,677 | 6,161 | 2.1m | 91ms | 1,405x | 6/6 OK |
| -92169307 | 1,868 | 5,764 | 46.1s | 77ms | 599x | 6/6 OK |
| -261756723 | 2,237 | 7,910 | 2.9m | 133ms | 1,318x | 6/6 OK |
| -289977907 | 2,405 | 7,790 | 2.7m | 161ms | 989x | 6/6 OK |
| -471122787 | 2,749 | 9,973 | 6.6m | 182ms | 2,181x | 6/6 OK |
| -309744003 | 3,369 | 10,756 | 3.9m | 138ms | 1,678x | 6/6 OK |
| -947098723 | 3,531 | 12,390 | 8.1m | 183ms | 2,661x | 6/6 OK |
| -367827107 | 4,695 | 14,428 | 5.4m | 228ms | 1,434x | 6/6 OK |
| -895993939 | 4,883 | 16,030 | 8.0m | 232ms | 2,076x | 6/6 OK |
| -692013211 | 5,215 | 15,993 | 7.1m | 332ms | 1,280x | 6/6 OK |
| -500038611 | 5,425 | 16,493 | 7.6m | 249ms | 1,825x | 6/6 OK |
| -1354169499 | 6,369 | 21,250 | 17.0m | 430ms | 2,374x | 6/6 OK |
| -972154635 | 6,649 | 20,313 | 21.7m | 351ms | 3,708x | 6/6 OK |

Cross-validation: each discriminant verified 6 ways (serial vs cache P1, serial vs cache P2, parallel vs cache P1, parallel vs cache P2, serial vs parallel P1, serial vs parallel P2). 120/120 passed.


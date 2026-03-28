# classpoly v1.0.3 — Parallel Fork

Parallelized Hilbert class polynomial computation via `fork()` with atomic work queue and shared memory accumulators.

## Usage

```bash
# Serial (original)
./classpoly D [inv P filename verbosity]

# Parallel with single P
./classpoly -j N D inv P1 0 filename1 "" [verbosity]

# Parallel with two P values
./classpoly -j N D inv P1 P2 filename1 filename2 [verbosity]
```

- `-j N` — number of worker processes (use physical core count)
- `D` — negative discriminant
- `inv` — class invariant (-1 = auto, 0 = j-invariant)
- `P1`, `P2` — target prime moduli

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

## Architecture

- **fork()** gives each worker its own copy of all mutable ff_poly globals via COW — zero changes to the finite field library
- **Atomic work queue** (`__sync_fetch_and_add`) dynamically balances load across workers — no static chunking
- **Shared memory** (`mmap MAP_SHARED`) accumulators — zero disk I/O
- **Single `ecrt_init`** in the parent before fork — workers inherit via COW, no redundant CRT precomputation
- Per-prime work (find j-invariant, enumerate roots, build polynomial) done once, then `ecrt_update` called per P value

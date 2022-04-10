/*
 * Copyright 2021 Benjamin Edgington
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdlib.h> // malloc(), free(), atoi()
#include <stdio.h>  // printf()
#include <assert.h> // assert()
#include <unistd.h> // EXIT_SUCCESS/FAILURE
#include "bench_util.h"
#include "test_util.h"
#include "c_kzg_alloc.h"
#include "c_kzg.h"

// Run the benchmark for `max_seconds` and return the time per iteration in nanoseconds.
long run_bench(int scale, int max_seconds) {
    timespec_t t0, t1;
    unsigned long total_time = 0, nits = 0;

    // Fill with randomness
    uint64_t coeffs[1 << (scale - 1)];
    int poly_len = sizeof coeffs / sizeof coeffs[0];
    for (uint64_t i = 0; i < poly_len; i++) {
        coeffs[i] = rand_uint64();
    }

    // The FFT settings size
    uint64_t n = scale, n_len = (uint64_t)1 << n;

    FFTSettings fs;
    KZGSettings ks;
    FK20SingleSettings fk;
    uint64_t secrets_len = n_len + 1;
    g1_t s1[secrets_len];
    g2_t s2[secrets_len];
    poly p;
    g1_t commitment;
    g1_t *all_proofs;

    assert(n_len >= 2 * poly_len);
    assert(C_KZG_OK == new_g1_array(&all_proofs, 2 * poly_len));
    assert(C_KZG_OK == new_poly(&p, poly_len));
    for (uint64_t i = 0; i < poly_len; i++) {
        fr_from_uint64(&p.coeffs[i], coeffs[i]);
    }

    // Initialise the secrets and data structures
    generate_trusted_setup(s1, s2, &secret, secrets_len);
    assert(C_KZG_OK == new_fft_settings(&fs, n));
    assert(C_KZG_OK == new_kzg_settings(&ks, s1, s2, secrets_len, &fs));
    assert(C_KZG_OK == new_fk20_single_settings(&fk, 2 * poly_len, &ks));

    // Commit to the polynomial
    assert(C_KZG_OK == commit_to_poly(&commitment, &p, &ks));

    while (total_time < max_seconds * NANO) {
        clock_gettime(CLOCK_REALTIME, &t0);
        assert(C_KZG_OK == da_using_fk20_single(all_proofs, &p, &fk, true));
        clock_gettime(CLOCK_REALTIME, &t1);
        nits++;
        total_time += tdiff(t0, t1);
    }

    free_poly(&p);
    free(all_proofs);
    free_fft_settings(&fs);
    free_kzg_settings(&ks);
    free_fk20_single_settings(&fk);

    return total_time / nits;
}

int main(int argc, char *argv[]) {
    int nsec = 0;

    switch (argc) {
        case 1:
            nsec = NSEC;
            break;
        case 2:
            nsec = atoi(argv[1]);
            break;
        default:
            break;
    };

    if (nsec == 0) {
        printf("Usage: %s [test time in seconds > 0]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    printf("*** Benchmarking fk_single_da, %d second%s per test.\n", nsec, nsec == 1 ? "" : "s");
    for (int scale = 4; scale <= 14; scale++) {
        printf("fk_single_da/scale_%d %lu ns/op\n", scale, run_bench(scale, nsec));
    }

    return EXIT_SUCCESS;
}

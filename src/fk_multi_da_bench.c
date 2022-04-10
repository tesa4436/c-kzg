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
#include "utility.h"
#include "c_kzg.h"

// Run the benchmark for `max_seconds` and return the time per iteration in nanoseconds.
long run_bench(int scale, int max_seconds) {
    timespec_t t0, t1;
    unsigned long total_time = 0, nits = 0;

    FFTSettings fs;
    KZGSettings ks;
    FK20MultiSettings fk;
    uint64_t chunk_count, width;
    uint64_t secrets_len;
    g1_t *s1;
    g2_t *s2;
    poly p;
    uint64_t vv[] = {1, 2, 3, 4, 7, 8, 9, 10, 13, 14, 1, 15, 1, 1000, 134, 33};
    g1_t commitment;
    g1_t *all_proofs;

    int chunk_len = 16;
    int n = 1 << scale;

    assert(is_power_of_two(n));
    assert(is_power_of_two(chunk_len));
    assert(n % 16 == 0);
    assert(n >= chunk_len);

    chunk_count = n / chunk_len;
    secrets_len = 2 * n;
    width = log2_pow2(secrets_len);

    assert(C_KZG_OK == new_g1_array(&s1, secrets_len));
    assert(C_KZG_OK == new_g2_array(&s2, secrets_len));

    generate_trusted_setup(s1, s2, &secret, secrets_len);
    assert(C_KZG_OK == new_fft_settings(&fs, width));
    assert(C_KZG_OK == new_kzg_settings(&ks, s1, s2, secrets_len, &fs));
    assert(C_KZG_OK == new_fk20_multi_settings(&fk, n * 2, chunk_len, &ks));

    // Create a test polynomial of size n that's independent of chunk_len
    assert(C_KZG_OK == new_poly(&p, n));
    for (int i = 0; i < chunk_count; i++) {
        for (int j = 0; j < chunk_len; j++) {
            int p_index = i * chunk_len + j;
            int v_index = p_index % 16;
            uint64_t v = vv[v_index];
            int tmp = i * chunk_len / 16;
            if (v_index == 3) v += tmp;
            if (v_index == 5) v += tmp * tmp;
            fr_from_uint64(&p.coeffs[p_index], v);
            if (v_index == 12) fr_negate(&p.coeffs[p_index], &p.coeffs[p_index]);
            if (v_index == 14) fr_negate(&p.coeffs[p_index], &p.coeffs[p_index]);
        }
    }

    assert(C_KZG_OK == commit_to_poly(&commitment, &p, &ks));
    assert(C_KZG_OK == new_g1_array(&all_proofs, 2 * chunk_count));

    while (total_time < max_seconds * NANO) {
        clock_gettime(CLOCK_REALTIME, &t0);
        assert(C_KZG_OK == da_using_fk20_multi(all_proofs, &p, &fk));
        clock_gettime(CLOCK_REALTIME, &t1);
        nits++;
        total_time += tdiff(t0, t1);
    }

    free_poly(&p);
    free(all_proofs);
    free(s1);
    free(s2);
    free_fft_settings(&fs);
    free_kzg_settings(&ks);
    free_fk20_multi_settings(&fk);

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

    printf("*** Benchmarking fk_multi_da, %d second%s per test.\n", nsec, nsec == 1 ? "" : "s");
    for (int scale = 4; scale <= 14; scale++) {
        printf("fk_multi_da/scale_%d %lu ns/op\n", scale, run_bench(scale, nsec));
    }

    return EXIT_SUCCESS;
}

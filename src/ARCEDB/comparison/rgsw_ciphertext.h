#pragma once

#include "seal/seal.h"

namespace arcedb
{
    /**
     * 2 * l * seal::ciphertext;
    */
    using RGSWCiphertext = seal::KSwitchKeys;

    /**
    We use RNS decomposition gadget to generate RGSW ciphertext.
    g_rns[i] = (Q / q_i) * q_{special} * [(Q / q_i)^{-1} mod q_i ]
    c(x) = c_1(x) * (Q / q_1) * [(Q / q_1)^{-1} mod q_1 ] + c_2(x) * (Q / q_2) * [(Q / q_2)^{-1} mod q_2 ] + ...
    c(x) * g_rns[1] = g_rns[1] * c_1(x) * (Q / q_1) * [(Q / q_1)^{-1} mod q_1 ] + g_rns[1] *c_2(x) * (Q / q_2) * [(Q / q_2)^{-1} mod q_2 ] + ...
    c(x) * g_rns[1] mod q_1 = q_{special} * c1(x) mod q_1
    c(x) * g_rns[1] mod q_2 = 0

    An rgsw ciphertext contains 2 * L rlwe ciphertexts.
    First Row: RLWE(m*g[1]), RLWE(m*g[2]), ... RLWE(m*g[L])
    Second Row: RLWE(m*s*g[1]), RLWE(m*s*g[2]), ... RLWE(m*s*g[L])
    Note for the second rowe:
    Instead encrypt (-a*s + ms +e,a), We do (-a*s + e, a+m) to do the same thing.
    And we can get two rows of ciphers in one loop
    */
    void encrypt_rgsw_rns(std::vector<int64_t> values, RGSWCiphertext &destination, 
                        seal::Encryptor &encryptor, seal::SecretKey &secret_key, seal::SEALContext &context);

    void external_product(seal::Ciphertext &encrypted, RGSWCiphertext &rgsw_cipher, seal::Ciphertext &destination,
                        seal::SEALContext &context, seal::Evaluator &evaluator);
} // namespace arcedb

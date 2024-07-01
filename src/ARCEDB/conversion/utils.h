#pragma once

#include "seal/seal.h"
#include "lwe_ciphertext.h"

namespace arcedb
{
    void encode_coefficient(std::vector<uint64_t> &values_matrix,  seal::Plaintext &destination, seal::SEALContext &context);

    void decode_coefficient(seal::Plaintext &plain, std::vector<uint64_t> &destination, seal::SEALContext &context);

    void decrypt_and_decode_lwe(LWECiphertext &cipher, uint64_t &destination, seal::SEALContext &context, seal::Decryptor &decryptor);
} // namespace arcedb

#pragma once
#include "seal/seal.h"
#include "ARCEDB/comparison/comparable.h"
namespace arcedb
{
    void exponetial_to_base_test_vector(seal::Ciphertext &exponential_cipher, seal::Ciphertext &base_cipher, 
                            seal::SEALContext &context, seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator, seal::Decryptor &decryptor);
} // namespace arcedb

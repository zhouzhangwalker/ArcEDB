#pragma once
#include "seal/seal.h"
#include "lwe_ciphertext.h"

namespace arcedb
{
    using namespace seal;
    
    void rlwe_trace_inplace(Ciphertext &rlwe_cipher, size_t num_rep, SEALContext &context, GaloisKeys &galois_keys, Evaluator &evaluator);

    void multiply_rlwe_mono_inplace(Ciphertext &rlwe_cipher, uint64_t coefficient, bool is_positive, 
                                    SEALContext &context);
    
    void multiply_rlwe_mono(Ciphertext &rlwe_cipher, Ciphertext &destination, uint64_t coefficient, bool is_positive, 
                                    SEALContext &context);
    
    void generate_pack_galois_keys(const SEALContext &context, GaloisKeys &galois_keys, KeyGenerator &keygen);

    void pack_lwe_preprocess(LWECiphertext &lwe_cipher, SEALContext &context);

    void lwe_to_rlwe(LWECiphertext &lwe_cipher, Ciphertext &rlwe_cipher, SEALContext &context);

    void pack_lwes(std::vector<LWECiphertext> &lwe_ciphers, Ciphertext &rlwe_cipher, 
                    SEALContext &context, GaloisKeys &galois_keys, Evaluator &evaluator);
} // namespace arcedb

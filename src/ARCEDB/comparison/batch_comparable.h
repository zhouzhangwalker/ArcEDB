#pragma once

#include "seal/seal.h"
#include "rgsw_ciphertext.h"
#include "ARCEDB/utils/types.h"
#include "ARCEDB/conversion/lwe_ciphertext.h"
#include "ARCEDB/comparison/comparable.h"
// #include "ARCEDB/comparison/data.h"
#include "ARCEDB/comparison/poly_eval.h"
#include "ARCEDB/comparison/batch_bootstrap.h"

namespace arcedb
{
    
    void compare_greater_than_batch_rlwe(std::vector<ComparableCipher> &cipher1, std::vector<ComparableRGSWCipher> &cipher2, std::vector<seal::Ciphertext> &destination,
                                    seal::SEALContext &context, seal::Evaluator &evaluator, seal::Decryptor &decryptor);
    
    void compare_greater_than_batch_lwe(std::vector<ComparableCipher> &cipher1, std::vector<ComparableRGSWCipher> &cipher2, std::vector<LWECiphertext> &destination,
                                seal::SEALContext &context, seal::Evaluator &evaluator, seal::Decryptor &decryptor);
    
    void compare_greater_than_batch_optimized(std::vector<ComparableCipher> &cipher1, std::vector<ComparableRGSWCipher> &cipher2, std::vector<LWECiphertext> &destination,
                                    seal::SEALContext &context, seal::Evaluator &evaluator, seal::Decryptor &decryptor);
    
    void simd_greater_than(std::vector<ComparableCipher> &cipher1, std::vector<ComparableRGSWCipher> &cipher2, size_t modular_size, std::vector<LWECiphertext> &destination,
                            seal::SEALContext &context, seal::Evaluator &evaluator, seal::Decryptor &decryptor,
                            LinearTransformKey &eval_key, seal::SEALContext &rlwe_context, seal::BatchEncoder &rlwe_batch_encoder, 
                            seal::GaloisKeys &rlwe_galois_keys, seal::Evaluator &rlwe_evaluator, seal::RelinKeys &rlwe_relin_keys,
                            seal::Encryptor &rlwe_encryptor, seal::Decryptor &rlwe_decryptor,
                            std::vector<uint64_t> &rlwe_mux_poly, seal::KSwitchKeys &rlwe_key_switch_key, std::vector<std::vector<Plaintext>> pre_plains);
    
    void simd_greater_than_32(std::vector<ComparableCipher> &cipher1, std::vector<ComparableRGSWCipher> &cipher2, size_t modular_size, 
                            std::vector<LWECiphertext> &destination, seal::SEALContext &lwe_context, seal::Evaluator &lwe_evaluator, seal::Decryptor &lwe_decryptor,
                            LinearTransformKey &eval_key, seal::SEALContext &rlwe_context, seal::BatchEncoder &rlwe_batch_encoder, 
                            seal::GaloisKeys &rlwe_galois_keys, seal::Evaluator &rlwe_evaluator, seal::RelinKeys &rlwe_relin_keys,
                            seal::Encryptor &rlwe_encryptor, seal::Decryptor &rlwe_decryptor,
                            std::vector<uint64_t> &rlwe_mux_poly, seal::KSwitchKeys &rlwe_key_switch_key, std::vector<std::vector<Plaintext>> pre_plains);
    
    void simd_greater_than_64(std::vector<ComparableCipher> &cipher1, std::vector<ComparableRGSWCipher> &cipher2, size_t modular_size, 
                            std::vector<LWECiphertext> &destination, seal::SEALContext &lwe_context, seal::Evaluator &lwe_evaluator, seal::Decryptor &lwe_decryptor,
                            LinearTransformKey &eval_key, seal::SEALContext &rlwe_context, seal::BatchEncoder &rlwe_batch_encoder, 
                            seal::GaloisKeys &rlwe_galois_keys, seal::Evaluator &rlwe_evaluator, seal::RelinKeys &rlwe_relin_keys,
                            seal::Encryptor &rlwe_encryptor, seal::Decryptor &rlwe_decryptor,
                            std::vector<uint64_t> &rlwe_mux_poly, seal::KSwitchKeys &rlwe_key_switch_key, std::vector<std::vector<Plaintext>> pre_plains);
    
    void batch_bootstrap_mux(std::vector<LWECiphertext> &input_cipher, std::vector<LWECiphertext> &output_cipher,
                            LinearTransformKey &eval_key, seal::SEALContext &rlwe_context, seal::BatchEncoder &rlwe_batch_encoder, 
                            seal::GaloisKeys &rlwe_galois_keys, seal::Evaluator &rlwe_evaluator, seal::RelinKeys &rlwe_relin_keys,
                            seal::Encryptor &rlwe_encryptor, seal::Decryptor &rlwe_decryptor,
                            std::vector<uint64_t> &rlwe_mux_poly, seal::KSwitchKeys &rlwe_key_switch_key, std::vector<std::vector<Plaintext>> pre_plains,
                            seal::SEALContext &lwe_context);
} // namespace arcedb

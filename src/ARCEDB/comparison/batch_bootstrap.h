#pragma once

#include "seal/seal.h"
#include "ARCEDB/conversion/lwe_ciphertext.h"
#include "ARCEDB/utils/utils.h"

namespace arcedb
{
    typedef struct LinearTransformKey
    {
        seal::Ciphertext key;
        std::vector<seal::Ciphertext> rotated_keys;
    } LinearTransformKey;

    using LWESecretKey = seal::Plaintext;



    void generate_slot_to_coeff_matrix(std::vector<std::vector<uint64_t>> &matrix, seal::SEALContext &context, seal::BatchEncoder &batch_encoder);

    void read_slot_to_coeff_matrix(std::vector<std::vector<uint64_t>> &matrix, seal::BatchEncoder &batch_encoder);

    void generate_slot_to_coeff_galois_keys(seal::SEALContext&context, seal::KeyGenerator &keygen, seal::GaloisKeys &galois_keys);

    void slot_to_coeff(seal::Ciphertext &cipher, seal::Ciphertext &destination, std::vector<std::vector<uint64_t>> &UMatrix,
                    seal::SEALContext &context, seal::BatchEncoder &batch_encoder, seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator);

    void generate_slot_to_coeff_plain(seal::parms_id_type parms_id, std::vector<std::vector<uint64_t>> &UMatrix, std::vector<std::vector<seal::Plaintext>> &plains,
                                    seal::SEALContext &context, seal::BatchEncoder &batch_encoder, seal::Evaluator &evaluator);

    void slot_to_coeff_optimized(seal::Ciphertext &cipher, seal::Ciphertext &destination, std::vector<std::vector<seal::Plaintext>> &plains,
                    seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator);
    
    void pack_encode(std::vector<uint64_t> &input, seal::Plaintext &plain, seal::BatchEncoder &batch_encoder);

    void convert_secert_key(seal::SecretKey &lwe_key, seal::SEALContext &context, std::vector<int64_t> &key_values);

    void generate_linear_transform_key(LinearTransformKey &eval_key, std::vector<int64_t> &lwe_key_values, seal::SecretKey &rlwe_key, 
                                        seal::SEALContext &context, BatchEncoder &batch_encoder, Encryptor &encryptor, Evaluator &evaluator);

    void LinearTransform(seal::Ciphertext &result, std::vector<std::vector<uint64_t>> &matrix, LinearTransformKey &eval_key, 
                        seal::SEALContext &context, seal::BatchEncoder &batch_encoder, seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator);
    
    void As_plus_b(std::vector<LWECiphertext> &lwe_ciphers, seal::Ciphertext &middle_result, LinearTransformKey &eval_key, 
                seal::SEALContext &context, seal::BatchEncoder &batch_encoder, seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator);
    
    
    void generate_key_switch_key(seal::SecretKey &lwe_secret_key, seal::SEALContext &lwe_context,
                            seal::SecretKey &rlwe_secret_key, seal::SEALContext &rlwe_context, 
                            seal::SecretKey &pad_lwe_secret_key, seal::KSwitchKeys &key_switch_key);

    void modulus_switching_inplace(seal::Ciphertext &cipher, const seal::Modulus &before_modulus, const seal::Modulus &after_modulus);

    void part_extract(seal::Ciphertext &rlwe_cipher, std::vector<LWECiphertext> &lwe_ciphers, seal::SEALContext &rlwe_context, seal::SEALContext &lwe_context);

    void modulus_switching_lwe_inplace(LWECiphertext &cipher, const seal::Modulus &before_modulus, const seal::Modulus &after_modulus);
} // namespace arcedb

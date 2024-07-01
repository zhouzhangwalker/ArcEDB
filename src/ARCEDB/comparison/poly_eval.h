#pragma once

#include "seal/seal.h"
namespace arcedb
{
    std::vector<uint32_t> optimal_parameter_paterson(uint32_t degree);

    std::vector<uint32_t> ComputeDegreesPS(const uint32_t n);

    int get_significant_coeff_count_poly(uint64_t *poly, size_t coeff_count, size_t coeff_uint64_count);

    inline std::uint64_t *get_poly_coeff(std::uint64_t *poly, int coeff_index, int coeff_uint64_count)
    {
        return poly + (static_cast<std::size_t>(coeff_index) * coeff_uint64_count);
    }

    void compute_all_powers(seal::Ciphertext &encrypted, uint32_t degree, std::vector<seal::Ciphertext> &destination, 
                            seal::Evaluator &evaluator, seal::RelinKeys &relin_keys);
    
    void divide_poly_poly_coeffmod_inplace(uint64_t *numerator, uint64_t *denominator, int coeff_count, seal::Modulus &modulus, uint64_t *quotient);

    inline void divide_poly_poly_coeffmod(uint64_t *numerator, uint64_t *denominator, int coeff_count, seal::Modulus &modulus, uint64_t *quotient, uint64_t *remainder)
    {
        int coeff_uint64_count = modulus.uint64_count();
        std::memcpy(remainder, numerator, coeff_count * sizeof(uint64_t) * coeff_uint64_count);
        divide_poly_poly_coeffmod_inplace(remainder, denominator, coeff_count, modulus, quotient);
    }

    void poly_evaluation_powers(std::vector<seal::Ciphertext> &all_powers_encrypted, std::vector<uint64_t> &coeff, seal::Ciphertext &destination,
                                seal::Evaluator &evaluator, seal::Decryptor &decryptor, seal::BatchEncoder &batch_encoder);
    
    void paterson_stockmeyer(seal::Ciphertext &encrypted, std::vector<seal::Ciphertext> &baby_step, std::vector<seal::Ciphertext> &giant_step,
                        std::vector<uint64_t> &coeff, seal::Ciphertext &destination, seal::SEALContext &context,
                        seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, int k, int step,
                        seal::Decryptor &decryptor, seal::BatchEncoder &batch_encoder);
    
    void poly_evaluation_paterson_stockmeyer(seal::Ciphertext &encrypted, std::vector<uint64_t> &coeff, seal::Ciphertext &destination,
                    seal::SEALContext &context, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, seal::Encryptor &encryptor,
                    seal::Decryptor &decryptor, seal::BatchEncoder &batch_encoder);
    
    void compute_even_powers(seal::Ciphertext &encrypted, uint32_t degree, std::vector<seal::Ciphertext> &destination, 
                        seal::Evaluator &evaluator, seal::RelinKeys &relin_keys);
    
    void optimized_poly_evaluation_paterson_stockmeyer(seal::Ciphertext &encrypted, std::vector<uint64_t> &coeff, seal::Ciphertext &destination,
                    seal::SEALContext &context, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, seal::Encryptor &encryptor,
                    seal::Decryptor &decryptor, seal::BatchEncoder &batch_encoder);
    
    void compute_mux_poly(seal::Ciphertext &encrypted, std::vector<uint64_t> &coeff, seal::Ciphertext &destination,
                        seal::SEALContext &context, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, seal::Encryptor &encryptor,
                        seal::Decryptor &decryptor, seal::BatchEncoder &batch_encoder);


} // namespace arcedb

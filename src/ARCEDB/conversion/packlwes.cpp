#include "packlwes.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/numth.h"
#include "seal/util/uintarithsmallmod.h"

namespace arcedb
{
    using namespace seal;
    
    void rlwe_trace_inplace(Ciphertext &rlwe_cipher, size_t num_rep, SEALContext &context, GaloisKeys &galois_keys, Evaluator &evaluator)
    {
        size_t coeff_count = rlwe_cipher.poly_modulus_degree();
        int coeff_power = util::get_power_of_two(coeff_count);
        Ciphertext automorphism_cipher;
        for (int k = 1; k <= static_cast<int>(std::log2(coeff_count / num_rep)); k++)
        {
            uint32_t galois_elt = std::pow(2, coeff_power - k + 1) + 1;
            evaluator.apply_galois(rlwe_cipher, galois_elt, galois_keys, automorphism_cipher);
            evaluator.add_inplace(rlwe_cipher, automorphism_cipher);
        }
        
    }

    void multiply_rlwe_mono_inplace(Ciphertext &rlwe_cipher, uint64_t coefficient, bool is_positive, 
                                    SEALContext &context)
    {
        auto context_data_ptr = context.get_context_data(rlwe_cipher.parms_id());
        auto &parms = context_data_ptr->parms();
        size_t coeff_count = parms.poly_modulus_degree();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_modulus_size = coeff_modulus.size();

        uint64_t mono_coeff;
        uint64_t mono_exponent;
        bool sig_is_negative;

        if (is_positive)
        {
            mono_exponent = coefficient % coeff_count;
            sig_is_negative = (coefficient / coeff_count) & 1;
        }
        else
        {
            mono_exponent = coefficient % coeff_count;
            mono_exponent = coeff_count - mono_exponent;
            sig_is_negative = ((coefficient / coeff_count + 1) & 1);
        }

        for (size_t t = 0; t < rlwe_cipher.size(); t++)
        {
            for (size_t i = 0; i < coeff_modulus_size; i++)
            {
                uint64_t *src_ptr = rlwe_cipher.data(t) + i * coeff_count;
                uint64_t q = coeff_modulus[i].value();
                if (sig_is_negative)
                {
                    mono_coeff = coeff_modulus[i].value() - 1;
                }
                else
                {
                    mono_coeff = 1;
                }
                util::multiply_poly_scalar_coeffmod(src_ptr, coeff_count, 
                                                    mono_coeff, coeff_modulus[i], src_ptr);
                std::rotate(src_ptr, src_ptr + coeff_count - mono_exponent, src_ptr + coeff_count);
                std::transform(src_ptr, src_ptr + mono_exponent, src_ptr, 
                                [q](uint64_t v){
                                    return v == 0? 0 : q - v;
                                });
            }
        }
    }

    void multiply_rlwe_mono(Ciphertext &rlwe_cipher, Ciphertext &destination, uint64_t coefficient, bool is_positive, 
                                    SEALContext &context)
    {
        destination = rlwe_cipher;
        multiply_rlwe_mono_inplace(destination, coefficient, is_positive, context);
    }

    void generate_pack_galois_keys(const SEALContext &context, GaloisKeys &galois_keys, KeyGenerator &keygen)
    {
        size_t N = context.key_context_data()->parms().poly_modulus_degree();
        int logN = util::get_power_of_two(N);
        std::vector<uint32_t> galois_elt;
        for (uint32_t i = 1; i <= logN; i++)
        {
            galois_elt.push_back((1u << i) + 1);
        }
        galois_elt.push_back(2 * N - 1);
        keygen.create_galois_keys(galois_elt, galois_keys);
    }

    void pack_lwe_preprocess(LWECiphertext &lwe_cipher, SEALContext &context)
    {
        auto context_data_ptr = context.get_context_data(lwe_cipher.parms_id());
        auto &parms = context_data_ptr->parms();
        size_t coeff_count = parms.poly_modulus_degree();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_modulus_size = coeff_modulus.size();
        std::vector<uint64_t> invert_N(coeff_modulus_size);

        for (size_t i = 0; i < coeff_modulus_size; i++)
        {
            util::try_invert_uint_mod(coeff_count, coeff_modulus[i], invert_N[i]);
        }
        
        for (size_t j = 0; j < coeff_modulus_size; j++)
        {
            util::multiply_poly_scalar_coeffmod(lwe_cipher.dataA() + j * coeff_count, coeff_count, 
                                    invert_N[j], coeff_modulus[j], lwe_cipher.dataA() + j * coeff_count);
            *(lwe_cipher.dataB()) = util::multiply_uint_mod(invert_N[j], *(lwe_cipher.dataB()), coeff_modulus[j]); 
        }
    }

    void lwe_to_rlwe(LWECiphertext &lwe_cipher, Ciphertext &rlwe_cipher, SEALContext &context)
    {
        if (lwe_cipher.is_ntt_form() || !lwe_cipher.is_set())
        {
            throw std::invalid_argument("Invalid LWE ciphertext.");
        }
        auto parms_id = lwe_cipher.parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &parms = context_data_ptr->parms();
        auto coeff_modulus = parms.coeff_modulus();
        size_t coeff_count = lwe_cipher.poly_modulus_degree();
        size_t coeff_modulus_size = lwe_cipher.coeff_modulus_size(); 
        rlwe_cipher.resize(context, lwe_cipher.parms_id(), 2);
        uint64_t * ptr = rlwe_cipher.data(1);
        for (size_t i = 0; i < coeff_modulus_size; i++)
        {
            uint64_t q = coeff_modulus[i].value();
            rlwe_cipher[i * coeff_count] = lwe_cipher[i];
            std::copy_n(lwe_cipher.dataA() + i * coeff_count, coeff_count, ptr);
            std::transform(ptr + 1, ptr + coeff_count, ptr + 1, 
                            [q](uint64_t v)
                            {
                                return v == 0? 0 : q -v;
                            });
            std::reverse(ptr + 1, ptr + coeff_count);
            ptr += coeff_count;
        }
        
        rlwe_cipher.scale() = lwe_cipher.scale();
        rlwe_cipher.is_ntt_form() = false;

    }

    void pack_lwes(std::vector<LWECiphertext> &lwe_ciphers, Ciphertext &result, 
                    SEALContext &context, GaloisKeys &galois_keys, Evaluator &evaluator)
    {
        size_t lwes_count = lwe_ciphers.size();
        if (lwes_count & lwes_count - 1)
        {
            throw std::invalid_argument("The number of repacking LWE ciphers is not power-of-2");
        }
        
        Ciphertext even_rlwe_cipher;
        Ciphertext odd_rlwe_cipher;
        if (lwes_count == 2)
        {
            lwe_to_rlwe(lwe_ciphers[0], even_rlwe_cipher, context);
            lwe_to_rlwe(lwe_ciphers[1], odd_rlwe_cipher, context);
        }
        else
        {
            std::vector<LWECiphertext> even_lwe_ciphers, odd_lwe_ciphers;
            for (size_t i = 0; i < lwes_count; i++)
            {
                if (i % 2 == 0)
                {
                    even_lwe_ciphers.push_back(lwe_ciphers[i]);
                }
                else
                {
                    odd_lwe_ciphers.push_back(lwe_ciphers[i]);
                }
            }
            pack_lwes(even_lwe_ciphers, even_rlwe_cipher, context, galois_keys, evaluator);
            pack_lwes(odd_lwe_ciphers, odd_rlwe_cipher, context, galois_keys, evaluator);
        }
        size_t coeff_count = lwe_ciphers[0].poly_modulus_degree();
        multiply_rlwe_mono_inplace(odd_rlwe_cipher, coeff_count / lwes_count, true, context);

        evaluator.add(even_rlwe_cipher, odd_rlwe_cipher, result);
        evaluator.sub_inplace(even_rlwe_cipher, odd_rlwe_cipher);
        evaluator.apply_galois_inplace(even_rlwe_cipher, lwes_count + 1, galois_keys);
        evaluator.add_inplace(result, even_rlwe_cipher);
    }
} // namespace arcedb

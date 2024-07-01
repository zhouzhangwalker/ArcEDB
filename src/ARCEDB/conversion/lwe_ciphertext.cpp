#include "ARCEDB/conversion/lwe_ciphertext.h"
#include <algorithm>

namespace arcedb
{
    LWECiphertext &LWECiphertext::operator=(const LWECiphertext &assign)
    {
        if (this == &assign)
        {
            return *this;
        }

        parms_id_ = assign.parms_id_;
        is_ntt_form_ = assign.is_ntt_form_;
        scale_ = assign.scale_;

        resize_internal(assign.poly_modulus_degree_, assign.coeff_modulus_size_);

        std::copy(assign.data_.cbegin(), assign.data_.cend(), data_.begin());

        return *this;
    }

    void LWECiphertext::resize(const SEALContext &context, parms_id_type parms_id)
    {
        if(!context.parameters_set())
        {
            throw std::invalid_argument("encryption parameters are not set correctly");
        }

        auto context_data_ptr = context.get_context_data(parms_id);
        if (!context_data_ptr)
        {
            throw std::invalid_argument("parms_id is not valid for encryption parameters");
        }

        auto &parms = context_data_ptr->parms();
        parms_id_ = context_data_ptr->parms_id();

        resize_internal(parms.poly_modulus_degree(), parms.coeff_modulus().size());

    }

    void LWECiphertext::resize_internal(size_t poly_modulus_degree, size_t coeff_modulus_size)
    {
        size_t new_data_size = util::mul_safe(poly_modulus_degree, coeff_modulus_size);
        new_data_size += coeff_modulus_size;

        data_.resize(new_data_size);
        poly_modulus_degree_ = poly_modulus_degree;
        coeff_modulus_size_ = coeff_modulus_size;
    }

    void LWECiphertext::extract(const SEALContext &context, Ciphertext &rlwe_cipher, size_t coeff_index)
    {
        // check parameters
        if (!is_metadata_valid_for(rlwe_cipher, context))
        {
            throw std::invalid_argument("rlwe ciphertext is not valid for encryption parameters");
        }

        if (rlwe_cipher.is_ntt_form())
        {
            throw std::invalid_argument("rlwe ciphertext should in non-ntt form");
        }

        resize(context, rlwe_cipher.parms_id());
        auto context_data_ptr = context.get_context_data(rlwe_cipher.parms_id());
        auto &parms = context_data_ptr->parms();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_count = rlwe_cipher.poly_modulus_degree();
        size_t coeff_modulus_size = rlwe_cipher.coeff_modulus_size();
        ct_coeff_type *b_ptr = data_.begin();
        ct_coeff_type *a_ptr = data_.begin() + coeff_modulus_size;
        ct_coeff_type *B_ptr = rlwe_cipher.data(0);
        ct_coeff_type *A_ptr = rlwe_cipher.data(1);
        for (size_t i = 0; i < coeff_modulus_size; i++)
        {
            uint64_t q = coeff_modulus[i].value();
            // Extract b
            b_ptr[i] = B_ptr[coeff_index];

            // Extract a
            std::copy_n(A_ptr, coeff_count, a_ptr);
            std::reverse(a_ptr, a_ptr + coeff_index + 1);
            std::reverse(a_ptr + coeff_index + 1, a_ptr + coeff_count);
            std::transform(a_ptr + coeff_index + 1, a_ptr + coeff_count, 
                            a_ptr + coeff_index + 1, 
                            [q](ct_coeff_type v)
                            {
                                return v == 0 ? 0 : q - v;
                            });
            
            b_ptr += 1;
            a_ptr += coeff_count;
            B_ptr += coeff_count;
            A_ptr += coeff_count;
        }
        scale_ = rlwe_cipher.scale();
        is_ntt_form_ = false;
    }

    uint64_t lwe_decrypt_q(LWECiphertext &cipher, std::vector<int64_t> &secret_key, seal::Modulus &modulus)
    {
        std::vector<uint64_t> temp(secret_key.size());
        if (secret_key.size() != cipher.poly_modulus_degree())
        {
            std::cout << "Location: " << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")\n";
            throw std::invalid_argument("Size mismatch.");
        }
        for (size_t i = 0; i < temp.size(); i++)
        {
            temp[i] = secret_key[i] >= 0 ? secret_key[i] : modulus.value() + secret_key[i];
        }
        uint64_t value = 0;
        for (size_t i = 0; i < secret_key.size(); i++)
        {
            uint64_t tt = util::multiply_uint_mod(*(cipher.dataA() + i), temp[i], modulus);
            value = util::add_uint_mod(value, tt, modulus);
        }
        value = util::add_uint_mod(*cipher.dataB(), value, modulus);
        // uint64_t value = util::dot_product_mod(cipher.dataA(), temp.data(), secret_key.size(), modulus);
        // value += *cipher.dataB();
        // value = util::modulo_uint(&value, 1, modulus);
        return value;
    }

    uint64_t lwe_decrypt(LWECiphertext &cipher, std::vector<int64_t> &secret_key, seal::Modulus &modulus, uint64_t t)
    {
        uint64_t value = lwe_decrypt_q(cipher, secret_key, modulus);
        uint64_t bound = modulus.value() / t;
        value = std::round(static_cast<double>(value) / bound);
        value = value % t;
        return value;

    }


    // Only one modulus
    void sample_extract(Ciphertext &rlwe_cipher, LWECiphertext &lwe_cipher, size_t coeff_index, seal::Modulus &modulus)
    {
        
        size_t coeff_count = rlwe_cipher.poly_modulus_degree();
        lwe_cipher.resize_internal(coeff_count, 1);
        uint64_t *b_ptr = lwe_cipher.dataB();
        uint64_t *a_ptr = lwe_cipher.dataA();
        uint64_t *B_ptr = rlwe_cipher.data(0);
        uint64_t *A_ptr = rlwe_cipher.data(1);
        uint64_t q = modulus.value();
        // Extract b
        b_ptr[0] = B_ptr[coeff_index];

        // Extract a
        std::copy_n(A_ptr, coeff_count, a_ptr);
        std::reverse(a_ptr, a_ptr + coeff_index + 1);
        std::reverse(a_ptr + coeff_index + 1, a_ptr + coeff_count);
        std::transform(a_ptr + coeff_index + 1, a_ptr + coeff_count, 
                        a_ptr + coeff_index + 1, 
                        [q](uint64_t v)
                        {
                            return v == 0 ? 0 : q - v;
                        });
        lwe_cipher.scale() = rlwe_cipher.scale();
        lwe_cipher.is_ntt_form() = false;
    }

    // Only one modulus
    void lwe_add_inplace(LWECiphertext &cipher1, LWECiphertext &cipher2, seal::Modulus &modulus)
    {
        size_t coeff_count = cipher1.poly_modulus_degree();
        uint64_t q = modulus.value();
        for (size_t i = 0; i < coeff_count; i++)
        {
            *(cipher1.dataA() + i) +=  *(cipher2.dataA() + i);
            if (*(cipher1.dataA() + i) >= q) *(cipher1.dataA() + i) -= q;
        }
        *cipher1.dataB() +=  *cipher2.dataB();
        if (*cipher1.dataB() >= q) *cipher1.dataB() -= q;
    }
} // namespace arcedb

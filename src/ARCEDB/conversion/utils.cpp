#include "utils.h"
#include "packlwes.h"
#include "lwe_ciphertext.h"

namespace arcedb
{
    void encode_coefficient(std::vector<uint64_t> &values_matrix,  seal::Plaintext &destination, seal::SEALContext &context)
    {
        auto &context_data = *context.first_context_data();
        size_t slots = context_data.parms().poly_modulus_degree();
        size_t values_matrix_size = values_matrix.size();
        if (values_matrix_size > slots)
        {
            throw std::invalid_argument("values_matrix size is too large");
        }

        destination.resize(slots);
        destination.parms_id() = seal::parms_id_zero;
        size_t gap = slots / values_matrix_size;
        for (size_t i = 0; i < values_matrix_size; i++)
        {
            *(destination.data() + gap * i) = values_matrix[i];
        }
    }

    void decode_coefficient(seal::Plaintext &plain, std::vector<uint64_t> &destination, seal::SEALContext &context)
    {
        if (!is_valid_for(plain, context))
        {
            throw std::invalid_argument("plain is not valid for encryption parameters");
        }
        if (plain.is_ntt_form())
        {
            throw std::invalid_argument("plain cannot be in NTT form");
        }

        auto &context_data = *context.first_context_data();
        size_t slots = context_data.parms().poly_modulus_degree();
        uint64_t modulus = context_data.parms().plain_modulus().value();
        size_t plain_coeff_count = std::min(plain.coeff_count(), slots);
        // Set destination size
        destination.resize(plain_coeff_count);
        uint64_t plain_modulus_div_two = modulus >> 1;
        for (size_t i = 0; i < plain_coeff_count; i++)
        {
            uint64_t curr_value = plain[i];
            destination[i] = (curr_value > plain_modulus_div_two)
                                 ? (static_cast<int64_t>(curr_value) - static_cast<int64_t>(modulus))
                                 : static_cast<int64_t>(curr_value);
        }
    }

    void decrypt_and_decode_lwe(LWECiphertext &cipher, uint64_t &destination, seal::SEALContext &context, seal::Decryptor &decryptor)
    {
        seal::Ciphertext temp_cipher;
        lwe_to_rlwe(cipher, temp_cipher, context);
        
        seal::Plaintext plain;
        std::vector<uint64_t> values;
        decryptor.decrypt(temp_cipher, plain);
        decode_coefficient(plain, values, context);
        destination = values[0];
    }

} // namespace arcedb
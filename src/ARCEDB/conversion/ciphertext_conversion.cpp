#include "ARCEDB/conversion/ciphertext_conversion.h"
#include "ARCEDB/conversion/packlwes.h"
namespace arcedb
{
    void exponetial_to_base_test_vector(seal::Ciphertext &exponential_cipher, seal::Ciphertext &base_cipher, 
                            seal::SEALContext &context, seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator, seal::Decryptor &decryptor)
    {
        using namespace seal;
        size_t poly_modulus_degree = exponential_cipher.poly_modulus_degree();
        Ciphertext nega_cipher;
        // X^{-b}
        evaluator.apply_galois(exponential_cipher, 2 * poly_modulus_degree - 1, galois_keys, nega_cipher);
        Plaintext constant_plaintext(poly_modulus_degree);
        for (size_t i = 0; i < poly_modulus_degree; i++)
        {
            constant_plaintext[i] = i;
        }
        evaluator.multiply_plain(nega_cipher, constant_plaintext, base_cipher);
        Plaintext debug_plaintext;
        decryptor.decrypt(base_cipher, debug_plaintext);
        auto context_data_ptr = context.get_context_data(exponential_cipher.parms_id());
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        Modulus plain_modulus = parms.plain_modulus();
        uint64_t inv_N_value;
        util::try_invert_uint_mod(poly_modulus_degree, plain_modulus, inv_N_value);
        for (size_t i = 0; i < poly_modulus_degree; i++)
        {
            constant_plaintext[i] = 0;
        }
        constant_plaintext[0] = inv_N_value;
        
        rlwe_trace_inplace(base_cipher, 1, context, galois_keys, evaluator);
        evaluator.multiply_plain_inplace(base_cipher, constant_plaintext);

    }
} // namespace arcedb

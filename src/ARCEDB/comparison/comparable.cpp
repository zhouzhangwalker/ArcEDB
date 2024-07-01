#include "ARCEDB/comparison/comparable.h"
#include "ARCEDB/conversion/packlwes.h"
#include "ARCEDB/conversion/lwe_ciphertext.h"

namespace arcedb
{
    
    // Encrypt (u * X^v)
    void comparable_encrypt_bfv(uint64_t v, seal::Ciphertext &destination, uint64_t u,
                            seal::SEALContext &context, seal::Encryptor &encryptor)
    {
        using namespace seal;
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t coeff_count = parms.poly_modulus_degree();

        if (v >= coeff_count)
        {
            throw std::invalid_argument("Out of bound.");
        }
        if (parms.scheme() != scheme_type::bfv)
        {
            throw std::invalid_argument("Scheme type mismatch.");
        }
        seal::Plaintext plain(coeff_count);
        plain.set_zero();
        // u* X^{v}
        plain[v] = u;
        encryptor.encrypt_symmetric(plain, destination);
    }

    void comparable_encrypt_bfv(uint64_t v, ComparableCipher &destination, uint64_t u,
                            seal::SEALContext &context, seal::Encryptor &encryptor)
    {
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t coeff_count = parms.poly_modulus_degree();

        uint8_t cipher_size = static_cast<uint8_t>(std::ceil(bits_count(v) / log(coeff_count)));
        destination.resize(cipher_size);
        uint64_t base = coeff_count;
        for (size_t i = 0; i < cipher_size; i++)
        {
            uint64_t value = v % base;
            comparable_encrypt_bfv(value, destination[i], u, context, encryptor);
            v = (v - value) / base;
        }
        
    }

    void comparable_encrypt_bfv(std::vector<uint64_t> v, ComparableCipher &destination, uint64_t u,
                            seal::SEALContext &context, seal::Encryptor &encryptor)
    {
        destination.resize(v.size());
        for (size_t i = 0; i < v.size(); i++)
        {
            comparable_encrypt_bfv(v[i], destination[i], u, context, encryptor);
        }
        
    }

    void comparable_encrypt_bfv(uint64_t v, uint64_t precision, ComparableCipher &destination, uint64_t u, 
                                seal::SEALContext &context, seal::Encryptor &encryptor)
    {
        if (v >= (1ULL << precision))
        {
            throw std::invalid_argument("Out of range.");
        }
        using namespace seal;
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t coeff_count = parms.poly_modulus_degree();

        uint8_t cipher_size = static_cast<uint8_t>(ceil(precision/log2(coeff_count)));
        destination.resize(cipher_size);
        uint64_t base = coeff_count;
        for (size_t i = 0; i < cipher_size; i++)
        {
            uint64_t value = v % base;
            comparable_encrypt_bfv(value, destination[i], u, context, encryptor);
            v = (v - value) / base;
        }
    }


    void comparable_rgsw_encrypt(uint64_t v, RGSWCiphertext &destination, bool is_negative,
                            seal::Encryptor &encryptor, seal::SecretKey &secret_key, seal::SEALContext &context)
    {
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t coeff_count = parms.poly_modulus_degree();
        std::vector<int64_t> plain(coeff_count, 0);

        if (v >= coeff_count)
        {
            throw std::invalid_argument("Out of bound.");
        }
        if (v == 0)
        {
            plain[v] = 1;
        }
        else
        {
            if (is_negative)
            {
                plain[coeff_count - v] = -1;
            }
            else
            {
                plain[v] = 1;
            }
        }
        
        encrypt_rgsw_rns(plain, destination, encryptor, secret_key, context);
    }

    

    void comparable_rgsw_encrypt(uint64_t v, ComparableRGSWCipher &destination, bool is_negative,
                                seal::Encryptor &encryptor, seal::SecretKey &secret_key, seal::SEALContext &context)
    {
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t coeff_count = parms.poly_modulus_degree();

        uint8_t cipher_size = static_cast<uint8_t>(std::ceil(log(v) / log(coeff_count)));
        destination.resize(cipher_size);
        uint64_t base = coeff_count;
        for (size_t i = 0; i < cipher_size; i++)
        {
            uint64_t value = v % base;
            comparable_rgsw_encrypt(value, destination[i], is_negative, encryptor, secret_key, context);
            v = (v - value) / base;
        }
    }

    void comparable_rgsw_encrypt(std::vector<uint64_t> v, ComparableRGSWCipher &destination, bool is_negative,
                                seal::Encryptor &encryptor, seal::SecretKey &secret_key, seal::SEALContext &context)
    {
        destination.resize(v.size());
        for (size_t i = 0; i < v.size(); i++)
        {
            comparable_rgsw_encrypt(v[i], destination[i], is_negative, encryptor, secret_key, context);
        }
    }

    void comparable_rgsw_encrypt(uint64_t v, uint64_t precision, ComparableRGSWCipher &destination, bool is_negative,
                                seal::Encryptor &encryptor, seal::SecretKey &secret_key, seal::SEALContext &context)
    {
        if (v >= (1ULL << precision))
        {
            throw std::invalid_argument("Out of range.");
        }
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t coeff_count = parms.poly_modulus_degree();

        uint8_t cipher_size = static_cast<uint8_t>(ceil(precision/log2(coeff_count)));
        destination.resize(cipher_size);
        uint64_t base = coeff_count;
        for (size_t i = 0; i < cipher_size; i++)
        {
            uint64_t value = v % base;
            comparable_rgsw_encrypt(value, destination[i], is_negative, encryptor, secret_key, context);
            v = (v - value) / base;
        }
    }

    void comparable_decrypt(seal::Ciphertext &cipher, uint64_t &destination,
                            seal::SEALContext &context, seal::Decryptor &decryptor)
    {
        using namespace seal;
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_modulus_size = coeff_modulus.size();
        size_t coeff_count = parms.poly_modulus_degree();
        auto ntt_tables = context_data.small_ntt_tables();

        if(parms.scheme() == scheme_type::bfv)
        {
            seal::Plaintext plain;
            decryptor.decrypt(cipher, plain);
            destination = plain[0];
        }
    }

    void generate_offset_plaintext(seal::Plaintext &offset_plaintext, seal::SEALContext &context)
    {
        using namespace seal;
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t poly_modulus_degree = parms.poly_modulus_degree();
        offset_plaintext.resize(poly_modulus_degree);
        Modulus plain_modulus = parms.plain_modulus();
        uint64_t inv_2;
        util::try_invert_uint_mod(2, plain_modulus, inv_2);
        offset_plaintext[0] = inv_2;
    }

    void generate_test_vector(compare_type compare_operator, seal::Plaintext &test_plaintext, seal::SEALContext &context)
    {
        using namespace seal;
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t poly_modulus_degree = parms.poly_modulus_degree();

        test_plaintext.resize(poly_modulus_degree);
        Modulus plain_modulus = parms.plain_modulus();
        uint64_t inv_2, neg_inv_2;
        util::try_invert_uint_mod(2, plain_modulus, inv_2);
        neg_inv_2 =  plain_modulus.value() - inv_2;

        switch (compare_operator)
        {
        case compare_type::gt:
            for (size_t i = 0; i < poly_modulus_degree; i++)
            {
                test_plaintext[i] = neg_inv_2;
            }
            break;
        case compare_type::lt:
            test_plaintext[0] = neg_inv_2;
            for (size_t i = 1; i < poly_modulus_degree; i++)
            {
                test_plaintext[i] = inv_2;
            }
            break;
        case compare_type::ge:
            test_plaintext[0] = inv_2;
            for (size_t i = 1; i < poly_modulus_degree; i++)
            {
                test_plaintext[i] = neg_inv_2;
            }
            break;
        case compare_type::le:
            for (size_t i = 0; i < poly_modulus_degree; i++)
            {
                test_plaintext[i] = inv_2;
            }
            break;
        default:
            throw std::invalid_argument("Just support compare type: gt, lt, ge, le.");
        }
    }

    /**
     * Evaluate if a is greater than b
    @param[in] cipher1 The first ciphertext to compare X^{a}
    @param[in] cipher2 The second ciphertext to compare X^{-b}
    Note the precision should less than log2(poly_modulus_degree)
    */

   

    void compare_greater_than(seal::Ciphertext &cipher1, RGSWCiphertext &cipher2,
                                    seal::Ciphertext &destination, seal::SEALContext &context, seal::Evaluator &evaluator)
    {
        using namespace seal;
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t poly_modulus_degree = parms.poly_modulus_degree();
        Plaintext test_plaintext(poly_modulus_degree), offset_plaintext(poly_modulus_degree);
        Modulus plain_modulus = parms.plain_modulus();
        uint64_t inv_2, neg_inv_2;
        util::try_invert_uint_mod(2, plain_modulus, inv_2);
        neg_inv_2 =  plain_modulus.value() - inv_2;
        for (size_t i = 0; i < poly_modulus_degree; i++)
        {
            test_plaintext[i] = neg_inv_2;
        }
        offset_plaintext[0] = inv_2;
        external_product(cipher1, cipher2, destination, context, evaluator);
        evaluator.multiply_plain_inplace(destination, test_plaintext);
        evaluator.add_plain_inplace(destination, offset_plaintext);    
    }

    void compare_greater_than(ComparableCipher &cipher1, ComparableRGSWCipher &cipher2, size_t compute_size, seal::Ciphertext &destination, 
                                seal::GaloisKeys &galois_keys, seal::RelinKeys relin_keys, seal::SEALContext &context, seal::Evaluator &evaluator)
    {
        using namespace seal;
        if (compute_size == 1)
        {
            compare_greater_than(cipher1[0], cipher2[0], destination, context, evaluator);
        }
        else
        {
            auto parms_id = context.first_parms_id();
            auto context_data_ptr = context.get_context_data(parms_id);
            auto &context_data = *context_data_ptr;
            auto &parms = context_data.parms();
            size_t poly_modulus_degree = parms.poly_modulus_degree();

            Ciphertext is_equal_cipher, high_compare_result, low_compare_result;
            Plaintext test_plaintext(poly_modulus_degree), offset_plaintext(poly_modulus_degree);

            Modulus plain_modulus = parms.plain_modulus();
            uint64_t inv_2, neg_inv_2;
            util::try_invert_uint_mod(2, plain_modulus, inv_2);
            neg_inv_2 =  plain_modulus.value() - inv_2;
            for (size_t i = 0; i < poly_modulus_degree; i++)
            {
                test_plaintext[i] = neg_inv_2;
            }
            offset_plaintext[0] = inv_2;

            // Evaluate a1 == b1 ?
            external_product(cipher1[compute_size-1], cipher2[compute_size-1], is_equal_cipher, context, evaluator);

            // Evaluate a1 > b1 ?
            evaluator.multiply_plain(is_equal_cipher, test_plaintext, high_compare_result);
            evaluator.add_plain_inplace(high_compare_result, offset_plaintext); 

            // Evaluate a0 > b1 ?
            compare_greater_than(cipher1, cipher2, compute_size - 1, low_compare_result, galois_keys, relin_keys, context, evaluator);

            rlwe_trace_inplace(is_equal_cipher, 1, context, galois_keys, evaluator);
            Plaintext plain(poly_modulus_degree);
            uint64_t inv_N;
            util::try_invert_uint_mod(poly_modulus_degree, plain_modulus, inv_N);
            plain[0] = inv_N;
            evaluator.multiply_plain_inplace(is_equal_cipher, plain);
            evaluator.sub_inplace(low_compare_result, high_compare_result);
            evaluator.multiply(is_equal_cipher, low_compare_result, destination);
            evaluator.relinearize_inplace(destination, relin_keys);
            evaluator.add_inplace(destination, high_compare_result);
            
        }
    }

    // Default two modular


    

    void compare_less_than(seal::Ciphertext &cipher1, RGSWCiphertext &cipher2,
                                    seal::Ciphertext &destination, seal::SEALContext &context, seal::Evaluator &evaluator)
    {
        using namespace seal;
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t poly_modulus_degree = parms.poly_modulus_degree();
        Plaintext test_plaintext(poly_modulus_degree), offset_plaintext(poly_modulus_degree);
        if (parms.scheme() == scheme_type::bfv)
        {
            Modulus plain_modulus = parms.plain_modulus();
            uint64_t inv_2, neg_inv_2;
            util::try_invert_uint_mod(2, plain_modulus, inv_2);
            neg_inv_2 =  plain_modulus.value() - inv_2;
            test_plaintext[0] = neg_inv_2;
            for (size_t i = 1; i < poly_modulus_degree; i++)
            {
                test_plaintext[i] = inv_2;
            }
            offset_plaintext[0] = inv_2;
            external_product(cipher1, cipher2, destination, context, evaluator);
            evaluator.multiply_plain_inplace(destination, test_plaintext);
            evaluator.add_plain_inplace(destination, offset_plaintext);
            
        }
        else if (parms.scheme() == scheme_type::ckks)
        {

        }
        else
        {
            std::invalid_argument("Scheme type not support.");
        }       
    }

    void compare_greater_than_or_equal(seal::Ciphertext &cipher1, RGSWCiphertext &cipher2,
                                    seal::Ciphertext &destination, seal::SEALContext &context, seal::Evaluator &evaluator)
    {
        using namespace seal;
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t poly_modulus_degree = parms.poly_modulus_degree();
        Plaintext test_plaintext(poly_modulus_degree), offset_plaintext(poly_modulus_degree);
        if (parms.scheme() == scheme_type::bfv)
        {
            Modulus plain_modulus = parms.plain_modulus();
            uint64_t inv_2, neg_inv_2;
            util::try_invert_uint_mod(2, plain_modulus, inv_2);
            neg_inv_2 =  plain_modulus.value() - inv_2;
            test_plaintext[0] = inv_2;
            for (size_t i = 1; i < poly_modulus_degree; i++)
            {
                test_plaintext[i] = neg_inv_2;
            }
            offset_plaintext[0] = inv_2;
            external_product(cipher1, cipher2, destination, context, evaluator);
            evaluator.multiply_plain_inplace(destination, test_plaintext);
            evaluator.add_plain_inplace(destination, offset_plaintext);
            
        }
        else if (parms.scheme() == scheme_type::ckks)
        {

        }
        else
        {
            std::invalid_argument("Scheme type not support.");
        }       
    }

    void compare_less_than_or_equal(seal::Ciphertext &cipher1, RGSWCiphertext &cipher2,
                                    seal::Ciphertext &destination, seal::SEALContext &context, seal::Evaluator &evaluator)
    {
        using namespace seal;
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t poly_modulus_degree = parms.poly_modulus_degree();
        Plaintext test_plaintext(poly_modulus_degree), offset_plaintext(poly_modulus_degree);
        if (parms.scheme() == scheme_type::bfv)
        {
            Modulus plain_modulus = parms.plain_modulus();
            uint64_t inv_2, neg_inv_2;
            util::try_invert_uint_mod(2, plain_modulus, inv_2);
            neg_inv_2 =  plain_modulus.value() - inv_2;
            for (size_t i = 0; i < poly_modulus_degree; i++)
            {
                test_plaintext[i] = inv_2;
            }
            offset_plaintext[0] = inv_2;
            external_product(cipher1, cipher2, destination, context, evaluator);
            evaluator.multiply_plain_inplace(destination, test_plaintext);
            evaluator.add_plain_inplace(destination, offset_plaintext);
            
        }
        else if (parms.scheme() == scheme_type::ckks)
        {

        }
        else
        {
            std::invalid_argument("Scheme type not support.");
        }       
    }

    void compare_equal(seal::Ciphertext &cipher1, RGSWCiphertext &cipher2,
                        seal::Ciphertext &destination, seal::SEALContext &context, seal::Evaluator &evaluator)
    {
        
        external_product(cipher1, cipher2, destination, context, evaluator);
            
    }

    void compare_inequal(seal::Ciphertext &cipher1, RGSWCiphertext &cipher2,
                        seal::Ciphertext &destination, seal::SEALContext &context, seal::Evaluator &evaluator)
    {
        using namespace seal;
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t poly_modulus_degree = parms.poly_modulus_degree();
        Plaintext one_plaintext(poly_modulus_degree);
        if (parms.scheme() == scheme_type::bfv)
        {
            for (size_t i = 0; i < poly_modulus_degree; i++)
            {
                one_plaintext[i] = 1;
            }
            external_product(cipher1, cipher2, destination, context, evaluator);
            evaluator.negate_inplace(destination);
            evaluator.add_plain_inplace(destination, one_plaintext);
        }
        else if (parms.scheme() == scheme_type::ckks)
        {

        }
        else
        {
            std::invalid_argument("Scheme type not support.");
        }           
    }



    void greater_than_tfhepp(TRLWELvl1 &cipher1, TRGSWLvl1 &cipher2, TLWELvl1 &res, TFHESecretKey &sk)
    {
        TRLWELvl1 trlwelvl1, trlwe_mul;
        TFHEpp::trgswfftExternalProduct<Lvl1>(trlwelvl1, cipher1, cipher2);
        TFHEpp::Polynomial<Lvl1> test_plaintext;
        for (size_t i = 0; i < Lvl1::n; i++)
        {
            test_plaintext[i] = 1;
        }

        TFHEpp::PolyMul<Lvl1>(trlwe_mul[0], trlwelvl1[0], test_plaintext);
        TFHEpp::PolyMul<Lvl1>(trlwe_mul[1], trlwelvl1[1], test_plaintext);
        TFHEpp::SampleExtractIndex<Lvl1>(res, trlwe_mul, 0);
        for (size_t i = 0; i <= Lvl1::n; i++)
        {
            res[i] = -res[i];
        }
    }

    void greater_than_tfhepp(std::vector<TRLWELvl1> &ciphers1, std::vector<TRGSWLvl1> &ciphers2, size_t cipher_size, TLWELvl1 &res, 
                            TFHEEvalKey &ek, TFHESecretKey &sk)
    {
        if (cipher_size == 1)
        {
            greater_than_tfhepp(ciphers1[0], ciphers2[0], res, sk);
        }
        else
        {
            TLWELvl1 low_res, high_res, equal_res;
            TRLWELvl1 trlwelvl1;
            greater_than_tfhepp(ciphers1, ciphers2, cipher_size - 1, low_res, ek, sk);
            TFHEpp::trgswfftExternalProduct<Lvl1>(trlwelvl1, ciphers1[cipher_size-1], ciphers2[cipher_size-1]);
            TFHEpp::SampleExtractIndex<Lvl1>(equal_res, trlwelvl1, 0);
            greater_than_tfhepp(ciphers1[cipher_size-1], ciphers2[cipher_size-1], high_res, sk);
            for (size_t i = 0; i <= Lvl1::n; i++)
            {
                high_res[i] = high_res[i] + high_res[i];
            }

            TLWELvl1 tlwelvl1;
            uint32_t offset = Lvl1::μ >> 1;
            for (size_t i = 0; i <= Lvl1::k * Lvl1::n; i++)
            {
                tlwelvl1[i] = equal_res[i] + high_res[i] + low_res[i];
            }
            tlwelvl1[Lvl1::n] += offset;
            TLWELvl0 tlwelvl0;
            TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlwelvl1, *ek.iksklvl10);
            TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(res, tlwelvl0, *ek.bkfftlvl01, TFHEpp::μpolygen<Lvl1, Lvl1::μ>());
        }
    }

    // Just for two digit
    void greater_than_tfhepp(std::vector<TRLWELvl1> &ciphers1, std::vector<TRGSWLvl1> &ciphers2, size_t cipher_size, TLWELvl1 &res, uint32_t scale_bits,
                            TFHEEvalKey &ek, TFHESecretKey &sk)
    {
        if (cipher_size == 1)
        {
            greater_than_tfhepp(ciphers1[0], ciphers2[0], res, sk);
        }
        else
        {
            TLWELvl1 low_res, high_res, equal_res;
            TRLWELvl1 trlwelvl1;
            greater_than_tfhepp(ciphers1, ciphers2, cipher_size - 1, low_res, ek, sk);
            TFHEpp::trgswfftExternalProduct<Lvl1>(trlwelvl1, ciphers1[cipher_size-1], ciphers2[cipher_size-1]);
            TFHEpp::SampleExtractIndex<Lvl1>(equal_res, trlwelvl1, 0);
            greater_than_tfhepp(ciphers1[cipher_size-1], ciphers2[cipher_size-1], high_res, sk);
            for (size_t i = 0; i <= Lvl1::n; i++)
            {
                high_res[i] = high_res[i] + high_res[i];
            }

            TLWELvl1 tlwelvl1;
            uint32_t offset = Lvl1::μ >> 1;
            for (size_t i = 0; i <= Lvl1::k * Lvl1::n; i++)
            {
                tlwelvl1[i] = equal_res[i] + high_res[i] + low_res[i];
            }
            tlwelvl1[Lvl1::n] += offset;
            Lvl1::T c = (1ULL << (scale_bits-1));
            TLWELvl0 tlwelvl0;
            TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlwelvl1, *ek.iksklvl10);
            TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(-c));
            res[Lvl1::k * Lvl1::n] += c;
        }
    }

    //Just two digit
    void greater_than_tfhepp(std::vector<TRLWELvl1> &ciphers1, std::vector<TRGSWLvl1> &ciphers2, size_t cipher_size, TLWELvl2 &res, uint32_t scale_bits,
                            TFHEEvalKey &ek, TFHESecretKey &sk)
    {
        TLWELvl1 low_res, high_res, equal_res;
        TRLWELvl1 trlwelvl1;
        greater_than_tfhepp(ciphers1, ciphers2, cipher_size - 1, low_res, ek, sk);
        TFHEpp::trgswfftExternalProduct<Lvl1>(trlwelvl1, ciphers1[cipher_size-1], ciphers2[cipher_size-1]);
        TFHEpp::SampleExtractIndex<Lvl1>(equal_res, trlwelvl1, 0);
        greater_than_tfhepp(ciphers1[cipher_size-1], ciphers2[cipher_size-1], high_res, sk);
        for (size_t i = 0; i <= Lvl1::n; i++)
        {
            high_res[i] = high_res[i] + high_res[i];
        }

        TLWELvl1 tlwelvl1;
        uint32_t offset = Lvl1::μ >> 1;
        for (size_t i = 0; i <= Lvl1::k * Lvl1::n; i++)
        {
            tlwelvl1[i] = equal_res[i] + high_res[i] + low_res[i];
        }
        tlwelvl1[Lvl1::n] += offset;
        Lvl2::T c = (1ULL << (scale_bits-1));
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlwelvl1, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl02>(res, tlwelvl0, *ek.bkfftlvl02, μ_polygen<Lvl2>(-c));
        res[Lvl2::k * Lvl2::n] += c;
    }

    void equality_tfhepp(TRLWELvl1 &cipher1, TRGSWLvl1 &cipher2, TLWELvl1 &res, TFHESecretKey &sk)
    {
        TRLWELvl1 trlwe_mul;
        TFHEpp::trgswfftExternalProduct<Lvl1>(trlwe_mul, cipher1, cipher2);
        TFHEpp::SampleExtractIndex<Lvl1>(res, trlwe_mul, 0);
        for (size_t i = 0; i <= Lvl1::n; i++)
        {
            res[i] = 2 * res[i];
        }
        res[Lvl1::n] -= Lvl1::μ;
    }

    void equality_tfhepp(std::vector<TRLWELvl1> &ciphers1, std::vector<TRGSWLvl1> &ciphers2, size_t cipher_size, TLWELvl1 &res, 
                            TFHEEvalKey &ek, TFHESecretKey &sk)
    {
        if (cipher_size == 1)
        {
            equality_tfhepp(ciphers1[0], ciphers2[0], res, sk);
        }
        else
        {
            TLWELvl1 low_res, high_res, equal_res;
            TRLWELvl1 trlwelvl1;
            equality_tfhepp(ciphers1, ciphers2, cipher_size - 1, low_res, ek, sk);
            equality_tfhepp(ciphers1[cipher_size-1], ciphers2[cipher_size-1], high_res, sk);
            TFHEpp::HomAND(res, low_res, high_res, ek);
        }
    }

    void less_than_tfhepp(TRLWELvl1 &cipher1, TRGSWLvl1 &cipher2, TLWELvl1 &res, TFHESecretKey &sk)
    {
        TRLWELvl1 trlwelvl1, trlwe_mul;
        TFHEpp::trgswfftExternalProduct<Lvl1>(trlwelvl1, cipher1, cipher2);
        TFHEpp::Polynomial<Lvl1> test_plaintext;
        test_plaintext[0] = Lvl1::plain_modulus - 1;
        for (size_t i = 1; i < Lvl1::n; i++)
        {
            test_plaintext[i] = 1;
        }

        TFHEpp::PolyMul<Lvl1>(trlwe_mul[0], trlwelvl1[0], test_plaintext);
        TFHEpp::PolyMul<Lvl1>(trlwe_mul[1], trlwelvl1[1], test_plaintext);
        TFHEpp::SampleExtractIndex<Lvl1>(res, trlwe_mul, 0);
    }

    void less_than_tfhepp(std::vector<TRLWELvl1> &ciphers1, std::vector<TRGSWLvl1> &ciphers2, size_t cipher_size, TLWELvl1 &res, 
                            TFHEEvalKey &ek, TFHESecretKey &sk)
    {
        if (cipher_size == 1)
        {
            less_than_tfhepp(ciphers1[0], ciphers2[0], res, sk);
        }
        else
        {
            TLWELvl1 low_res, high_res, equal_res;
            TRLWELvl1 trlwelvl1;
            less_than_tfhepp(ciphers1, ciphers2, cipher_size - 1, low_res, ek, sk);
            TFHEpp::trgswfftExternalProduct<Lvl1>(trlwelvl1, ciphers1[cipher_size-1], ciphers2[cipher_size-1]);
            TFHEpp::SampleExtractIndex<Lvl1>(equal_res, trlwelvl1, 0);
            less_than_tfhepp(ciphers1[cipher_size-1], ciphers2[cipher_size-1], high_res, sk);
            for (size_t i = 0; i <= Lvl1::n; i++)
            {
                high_res[i] = high_res[i] + high_res[i];
            }

            TLWELvl1 tlwelvl1;
            uint32_t offset = Lvl1::μ >> 1;
            for (size_t i = 0; i <= Lvl1::k * Lvl1::n; i++)
            {
                tlwelvl1[i] = equal_res[i] + high_res[i] + low_res[i];
            }
            tlwelvl1[Lvl1::n] += offset;
            TLWELvl0 tlwelvl0;
            TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlwelvl1, *ek.iksklvl10);
            TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(res, tlwelvl0, *ek.bkfftlvl01, TFHEpp::μpolygen<Lvl1, Lvl1::μ>());
        }
    }

    void lift_and_and(TLWELvl1 &cipher1, TLWELvl1 &cipher2, TLWELvl1 &res, uint32_t scale_bits, TFHEpp::EvalKey &ek, TFHEpp::SecretKey &sk)
    {
        using namespace TFHEpp;
        TLWELvl1 temp;
        for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
            temp[i] = cipher1[i] + cipher2[i];
        temp[Lvl1::k * Lvl1::n] -= Lvl1::μ;
        Lvl1::T c = (1ULL << (scale_bits-1));
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, temp, ek.getiksk<Lvl10>());
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(res, tlwelvl0, ek.getbkfft<Lvl01>(), μ_polygen<Lvl1>(-c));
        res[Lvl1::k * Lvl1::n] += c;
    }

    void lift_and_and(TLWELvl1 &cipher1, TLWELvl1 &cipher2, TLWELvl2 &res, uint32_t scale_bits, TFHEpp::EvalKey &ek, TFHEpp::SecretKey &sk)
    {
        using namespace TFHEpp;
        TLWELvl1 temp;
        for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
            temp[i] = cipher1[i] + cipher2[i];
        temp[Lvl1::k * Lvl1::n] -= Lvl1::μ;
        Lvl2::T c = (1ULL << (scale_bits-1));
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, temp, ek.getiksk<Lvl10>());
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl02>(res, tlwelvl0, ek.getbkfft<Lvl02>(), μ_polygen<Lvl2>(-c));
        res[Lvl2::k * Lvl2::n] += c;
    }
    
} // namespace arcedb

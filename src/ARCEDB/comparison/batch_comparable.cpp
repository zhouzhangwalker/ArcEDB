
#include "batch_comparable.h"

namespace arcedb
{

    // Only consider the two modular
    void compare_greater_than_batch_rlwe(std::vector<ComparableCipher> &cipher1, std::vector<ComparableRGSWCipher> &cipher2, std::vector<seal::Ciphertext> &destination,
                                    seal::SEALContext &context, seal::Evaluator &evaluator, seal::Decryptor &decryptor)
    {
        using namespace seal;

        size_t cipher_size = cipher1.size();
        destination.resize(cipher_size);
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t poly_modulus_degree = parms.poly_modulus_degree();
        Modulus plain_modulus = parms.plain_modulus();
        Plaintext debug_plain;
        Plaintext offset_plaintext(poly_modulus_degree),test_plaintext(poly_modulus_degree);

        for (size_t i = 0; i < poly_modulus_degree; i++)
        {
            test_plaintext[i] = plain_modulus.value() - 1;
        }

        seal::Ciphertext compare_result_low, is_equal_ciphertext, compare_result_high, rlwe_result;
        for (size_t i = 0; i < cipher_size; i++)
        {
            // get compare_result_low = a0 > b0 (0 or q/4)
            external_product(cipher1[i][0], cipher2[i][0], compare_result_low, context, evaluator);
            evaluator.multiply_plain_inplace(compare_result_low, test_plaintext); // may be optimized
            offset_plaintext[0] = 1;
            evaluator.add_plain_inplace(compare_result_low, offset_plaintext);
            decryptor.decrypt(compare_result_low, debug_plain);

            // get is_equal_ciphertext = a1 == b1 (0 or q/8)
            external_product(cipher1[i][1], cipher2[i][1], is_equal_ciphertext, context, evaluator);

            // get compare_result_high = a1 > b1 (0 or q/2)
            evaluator.multiply_plain(is_equal_ciphertext, test_plaintext, compare_result_high); // may be optimized
            offset_plaintext[0] = 1;
            evaluator.add_plain_inplace(compare_result_high, offset_plaintext);
            evaluator.add_inplace(compare_result_high, compare_result_high);


            // Add compare_result_low + compare_result_high + is_equal_ciphertext * 2
            evaluator.add(is_equal_ciphertext, is_equal_ciphertext, destination[i]);
            evaluator.add_inplace(destination[i], compare_result_high);
            evaluator.add_inplace(destination[i], compare_result_low);
            evaluator.mod_switch_to_next_inplace(destination[i]);
        }
    }

    void compare_greater_than_batch_lwe(std::vector<ComparableCipher> &cipher1, std::vector<ComparableRGSWCipher> &cipher2, std::vector<LWECiphertext> &destination,
                                    seal::SEALContext &context, seal::Evaluator &evaluator, seal::Decryptor &decryptor)
    {
        using namespace seal;

        size_t cipher_size = cipher1.size();
        destination.resize(cipher_size);
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t poly_modulus_degree = parms.poly_modulus_degree();
        Modulus plain_modulus = parms.plain_modulus();
        Plaintext offset_plaintext(poly_modulus_degree),test_plaintext(poly_modulus_degree);

        for (size_t i = 0; i < poly_modulus_degree; i++)
        {
            test_plaintext[i] = plain_modulus.value() - 1;
        }

        seal::Ciphertext compare_result_low, is_equal_ciphertext, compare_result_high, rlwe_result;

        for (size_t i = 0; i < cipher_size; i++)
        {
            // get compare_result_low = a0 > b0 (0 or q/4)
            external_product(cipher1[i][0], cipher2[i][0], compare_result_low, context, evaluator);
            evaluator.multiply_plain_inplace(compare_result_low, test_plaintext); // may be optimized
            offset_plaintext[0] = 1;
            evaluator.add_plain_inplace(compare_result_low, offset_plaintext);

            // get is_equal_ciphertext = a1 == b1 (0 or q/8)
            external_product(cipher1[i][1], cipher2[i][1], is_equal_ciphertext, context, evaluator);

            // get compare_result_high = a1 > b1 (0 or q/2)
            evaluator.multiply_plain(is_equal_ciphertext, test_plaintext, compare_result_high); // may be optimized
            offset_plaintext[0] = 1;
            evaluator.add_plain_inplace(compare_result_high, offset_plaintext);
            evaluator.add_inplace(compare_result_high, compare_result_high);


            // Add compare_result_low + compare_result_high + is_equal_ciphertext * 2
            evaluator.add(is_equal_ciphertext, is_equal_ciphertext, rlwe_result);
            evaluator.add_inplace(rlwe_result, compare_result_high);
            evaluator.add_inplace(rlwe_result, compare_result_low);
            evaluator.mod_switch_to_next_inplace(rlwe_result);
            destination[i].extract(context, rlwe_result, 0);
        }

    }

    void simd_greater_than(std::vector<ComparableCipher> &cipher1, std::vector<ComparableRGSWCipher> &cipher2, size_t modular_size, 
                            std::vector<LWECiphertext> &destination, seal::SEALContext &lwe_context, seal::Evaluator &lwe_evaluator, seal::Decryptor &lwe_decryptor,
                            LinearTransformKey &eval_key, seal::SEALContext &rlwe_context, seal::BatchEncoder &rlwe_batch_encoder, 
                            seal::GaloisKeys &rlwe_galois_keys, seal::Evaluator &rlwe_evaluator, seal::RelinKeys &rlwe_relin_keys,
                            seal::Encryptor &rlwe_encryptor, seal::Decryptor &rlwe_decryptor,
                            std::vector<uint64_t> &rlwe_mux_poly, seal::KSwitchKeys &rlwe_key_switch_key, std::vector<std::vector<Plaintext>> pre_plains)
    {
        std::vector<LWECiphertext> lwe_result;
        compare_greater_than_batch_lwe(cipher1, cipher2, lwe_result, lwe_context, lwe_evaluator, lwe_decryptor);
        batch_bootstrap_mux(lwe_result, destination, eval_key, rlwe_context, rlwe_batch_encoder, rlwe_galois_keys, rlwe_evaluator, rlwe_relin_keys, 
                            rlwe_encryptor, rlwe_decryptor, rlwe_mux_poly, rlwe_key_switch_key, pre_plains, lwe_context);
    }

    void simd_greater_than_32(std::vector<ComparableCipher> &cipher1, std::vector<ComparableRGSWCipher> &cipher2, size_t modular_size, 
                            std::vector<LWECiphertext> &destination, seal::SEALContext &lwe_context, seal::Evaluator &lwe_evaluator, seal::Decryptor &lwe_decryptor,
                            LinearTransformKey &eval_key, seal::SEALContext &rlwe_context, seal::BatchEncoder &rlwe_batch_encoder, 
                            seal::GaloisKeys &rlwe_galois_keys, seal::Evaluator &rlwe_evaluator, seal::RelinKeys &rlwe_relin_keys,
                            seal::Encryptor &rlwe_encryptor, seal::Decryptor &rlwe_decryptor,
                            std::vector<uint64_t> &rlwe_mux_poly, seal::KSwitchKeys &rlwe_key_switch_key, std::vector<std::vector<Plaintext>> pre_plains)
    {
        std::vector<LWECiphertext> low_lwe_result;

        // low compare result
        compare_greater_than_batch_lwe(cipher1, cipher2, low_lwe_result, lwe_context, lwe_evaluator, lwe_decryptor);
        batch_bootstrap_mux(low_lwe_result, low_lwe_result, eval_key, rlwe_context, rlwe_batch_encoder, rlwe_galois_keys, rlwe_evaluator, rlwe_relin_keys, 
                            rlwe_encryptor, rlwe_decryptor, rlwe_mux_poly, rlwe_key_switch_key, pre_plains, lwe_context);

        size_t cipher_size = cipher1.size();
        destination.resize(cipher_size);
        auto parms_id = lwe_context.first_parms_id();
        auto context_data_ptr = lwe_context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t poly_modulus_degree = parms.poly_modulus_degree();
        Modulus plain_modulus = parms.plain_modulus();
        Plaintext offset_plaintext(poly_modulus_degree),test_plaintext(poly_modulus_degree);

        for (size_t i = 0; i < poly_modulus_degree; i++)
        {
            test_plaintext[i] = plain_modulus.value() - 1;
        }

        seal::Ciphertext is_equal_ciphertext, compare_result_high, rlwe_result;

        for (size_t i = 0; i < cipher_size; i++)
        {
            // get is_equal_ciphertext = a2 == b2 (0 or q/8)
            external_product(cipher1[i][2], cipher2[i][2], is_equal_ciphertext, lwe_context, lwe_evaluator);

            // get compare_result_high = a2 > b2 (0 or q/2)
            lwe_evaluator.multiply_plain(is_equal_ciphertext, test_plaintext, compare_result_high); // may be optimized
            offset_plaintext[0] = 1;
            lwe_evaluator.add_plain_inplace(compare_result_high, offset_plaintext);
            lwe_evaluator.add_inplace(compare_result_high, compare_result_high);


            // Add compare_result_low + compare_result_high + is_equal_ciphertext * 2
            lwe_evaluator.add(is_equal_ciphertext, is_equal_ciphertext, rlwe_result);
            lwe_evaluator.add_inplace(rlwe_result, compare_result_high);
            lwe_evaluator.mod_switch_to_next_inplace(rlwe_result);
            destination[i].extract(lwe_context, rlwe_result, 0);
            lwe_add_inplace(destination[i], low_lwe_result[i], plain_modulus);
        }
        
        batch_bootstrap_mux(destination, destination, eval_key, rlwe_context, rlwe_batch_encoder, rlwe_galois_keys, rlwe_evaluator, rlwe_relin_keys, 
                            rlwe_encryptor, rlwe_decryptor, rlwe_mux_poly, rlwe_key_switch_key, pre_plains, lwe_context);
        
    }
    void simd_greater_than_64(std::vector<ComparableCipher> &cipher1, std::vector<ComparableRGSWCipher> &cipher2, size_t modular_size, 

                            std::vector<LWECiphertext> &destination, seal::SEALContext &lwe_context, seal::Evaluator &lwe_evaluator, seal::Decryptor &lwe_decryptor,
                            LinearTransformKey &eval_key, seal::SEALContext &rlwe_context, seal::BatchEncoder &rlwe_batch_encoder, 
                            seal::GaloisKeys &rlwe_galois_keys, seal::Evaluator &rlwe_evaluator, seal::RelinKeys &rlwe_relin_keys,
                            seal::Encryptor &rlwe_encryptor, seal::Decryptor &rlwe_decryptor,
                            std::vector<uint64_t> &rlwe_mux_poly, seal::KSwitchKeys &rlwe_key_switch_key, std::vector<std::vector<Plaintext>> pre_plains)
    {
        std::vector<LWECiphertext> low_lwe_result;

        // low compare result
        compare_greater_than_batch_lwe(cipher1, cipher2, low_lwe_result, lwe_context, lwe_evaluator, lwe_decryptor);
        batch_bootstrap_mux(low_lwe_result, low_lwe_result, eval_key, rlwe_context, rlwe_batch_encoder, rlwe_galois_keys, rlwe_evaluator, rlwe_relin_keys, 
                            rlwe_encryptor, rlwe_decryptor, rlwe_mux_poly, rlwe_key_switch_key, pre_plains, lwe_context);

        size_t cipher_size = cipher1.size();
        destination.resize(cipher_size);
        auto parms_id = lwe_context.first_parms_id();
        auto context_data_ptr = lwe_context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t poly_modulus_degree = parms.poly_modulus_degree();
        Modulus plain_modulus = parms.plain_modulus();
        Plaintext offset_plaintext(poly_modulus_degree),test_plaintext(poly_modulus_degree);

        for (size_t i = 0; i < poly_modulus_degree; i++)
        {
            test_plaintext[i] = plain_modulus.value() - 1;
        }

        seal::Ciphertext is_equal_ciphertext, compare_result_high, rlwe_result;

        for (size_t i = 0; i < cipher_size; i++)
        {
            // get is_equal_ciphertext = a2 == b2 (0 or q/8)
            external_product(cipher1[i][2], cipher2[i][2], is_equal_ciphertext, lwe_context, lwe_evaluator);

            // get compare_result_high = a2 > b2 (0 or q/2)
            lwe_evaluator.multiply_plain(is_equal_ciphertext, test_plaintext, compare_result_high); // may be optimized
            offset_plaintext[0] = 1;
            lwe_evaluator.add_plain_inplace(compare_result_high, offset_plaintext);
            lwe_evaluator.add_inplace(compare_result_high, compare_result_high);


            // Add compare_result_low + compare_result_high + is_equal_ciphertext * 2
            lwe_evaluator.add(is_equal_ciphertext, is_equal_ciphertext, rlwe_result);
            lwe_evaluator.add_inplace(rlwe_result, compare_result_high);
            lwe_evaluator.mod_switch_to_next_inplace(rlwe_result);
            destination[i].extract(lwe_context, rlwe_result, 0);
            lwe_add_inplace(destination[i], low_lwe_result[i], plain_modulus);
        }
        
        batch_bootstrap_mux(destination, low_lwe_result, eval_key, rlwe_context, rlwe_batch_encoder, rlwe_galois_keys, rlwe_evaluator, rlwe_relin_keys, 
                            rlwe_encryptor, rlwe_decryptor, rlwe_mux_poly, rlwe_key_switch_key, pre_plains, lwe_context);
        
        for (size_t i = 0; i < cipher_size; i++)
        {
            // get is_equal_ciphertext = a3 == b3 (0 or q/8)
            external_product(cipher1[i][3], cipher2[i][3], is_equal_ciphertext, lwe_context, lwe_evaluator);

            // get compare_result_high = a3 > b3 (0 or q/2)
            lwe_evaluator.multiply_plain(is_equal_ciphertext, test_plaintext, compare_result_high); // may be optimized
            offset_plaintext[0] = 1;
            lwe_evaluator.add_plain_inplace(compare_result_high, offset_plaintext);
            lwe_evaluator.add_inplace(compare_result_high, compare_result_high);


            // Add compare_result_low + compare_result_high + is_equal_ciphertext * 2
            lwe_evaluator.add(is_equal_ciphertext, is_equal_ciphertext, rlwe_result);
            lwe_evaluator.add_inplace(rlwe_result, compare_result_high);
            lwe_evaluator.mod_switch_to_next_inplace(rlwe_result);
            destination[i].extract(lwe_context, rlwe_result, 0);
            lwe_add_inplace(destination[i], low_lwe_result[i], plain_modulus);
        }

        batch_bootstrap_mux(destination, low_lwe_result, eval_key, rlwe_context, rlwe_batch_encoder, rlwe_galois_keys, rlwe_evaluator, rlwe_relin_keys, 
                            rlwe_encryptor, rlwe_decryptor, rlwe_mux_poly, rlwe_key_switch_key, pre_plains, lwe_context);
        
        for (size_t i = 0; i < cipher_size; i++)
        {
            // get is_equal_ciphertext = a4 == b4 (0 or q/8)
            external_product(cipher1[i][4], cipher2[i][4], is_equal_ciphertext, lwe_context, lwe_evaluator);

            // get compare_result_high = a4 > b4 (0 or q/2)
            lwe_evaluator.multiply_plain(is_equal_ciphertext, test_plaintext, compare_result_high); // may be optimized
            offset_plaintext[0] = 1;
            lwe_evaluator.add_plain_inplace(compare_result_high, offset_plaintext);
            lwe_evaluator.add_inplace(compare_result_high, compare_result_high);


            // Add compare_result_low + compare_result_high + is_equal_ciphertext * 2
            lwe_evaluator.add(is_equal_ciphertext, is_equal_ciphertext, rlwe_result);
            lwe_evaluator.add_inplace(rlwe_result, compare_result_high);
            lwe_evaluator.mod_switch_to_next_inplace(rlwe_result);
            destination[i].extract(lwe_context, rlwe_result, 0);
            lwe_add_inplace(destination[i], low_lwe_result[i], plain_modulus);
        }

        batch_bootstrap_mux(destination, low_lwe_result, eval_key, rlwe_context, rlwe_batch_encoder, rlwe_galois_keys, rlwe_evaluator, rlwe_relin_keys, 
                            rlwe_encryptor, rlwe_decryptor, rlwe_mux_poly, rlwe_key_switch_key, pre_plains, lwe_context);
        
        for (size_t i = 0; i < cipher_size; i++)
        {
            // get is_equal_ciphertext = a5 == b5 (0 or q/8)
            external_product(cipher1[i][5], cipher2[i][5], is_equal_ciphertext, lwe_context, lwe_evaluator);

            // get compare_result_high = a5 > b5 (0 or q/2)
            lwe_evaluator.multiply_plain(is_equal_ciphertext, test_plaintext, compare_result_high); // may be optimized
            offset_plaintext[0] = 1;
            lwe_evaluator.add_plain_inplace(compare_result_high, offset_plaintext);
            lwe_evaluator.add_inplace(compare_result_high, compare_result_high);


            // Add compare_result_low + compare_result_high + is_equal_ciphertext * 2
            lwe_evaluator.add(is_equal_ciphertext, is_equal_ciphertext, rlwe_result);
            lwe_evaluator.add_inplace(rlwe_result, compare_result_high);
            lwe_evaluator.mod_switch_to_next_inplace(rlwe_result);
            destination[i].extract(lwe_context, rlwe_result, 0);
            lwe_add_inplace(destination[i], low_lwe_result[i], plain_modulus);
        }

        batch_bootstrap_mux(destination, destination, eval_key, rlwe_context, rlwe_batch_encoder, rlwe_galois_keys, rlwe_evaluator, rlwe_relin_keys, 
                            rlwe_encryptor, rlwe_decryptor, rlwe_mux_poly, rlwe_key_switch_key, pre_plains, lwe_context);
        
    }

    void batch_bootstrap_mux(std::vector<LWECiphertext> &input_cipher, std::vector<LWECiphertext> &output_cipher,
                            LinearTransformKey &eval_key, seal::SEALContext &rlwe_context, seal::BatchEncoder &rlwe_batch_encoder, 
                            seal::GaloisKeys &rlwe_galois_keys, seal::Evaluator &rlwe_evaluator, seal::RelinKeys &rlwe_relin_keys,
                            seal::Encryptor &rlwe_encryptor, seal::Decryptor &rlwe_decryptor, 
                            std::vector<uint64_t> &rlwe_mux_poly, seal::KSwitchKeys &rlwe_key_switch_key, std::vector<std::vector<Plaintext>> pre_plains,
                            seal::SEALContext &lwe_context)
    {
        using namespace seal;
        size_t cipher_size = input_cipher.size();
        auto rlwe_parms_id = rlwe_context.first_parms_id();
        auto rlwe_context_data_ptr = rlwe_context.get_context_data(rlwe_parms_id);
        auto &rlwe_context_data = *rlwe_context_data_ptr;
        auto &rlwe_parms = rlwe_context_data.parms();
        size_t rlwe_poly_modulus_degree = rlwe_parms.poly_modulus_degree();

        auto lwe_parms_id = lwe_context.first_parms_id();
        auto lwe_context_data_ptr = lwe_context.get_context_data(lwe_parms_id);
        auto &lwe_context_data = *lwe_context_data_ptr;
        auto &lwe_parms = lwe_context_data.parms();
        auto lwe_coeff_modulus = lwe_parms.coeff_modulus();
        size_t lwe_poly_modulus_degree = lwe_parms.poly_modulus_degree();
        std::chrono::system_clock::time_point start, end;
        double compare_time;
        // As + b 
        Ciphertext asb_cipher;
        try
        {
            // start = std::chrono::system_clock::now();
            As_plus_b(input_cipher, asb_cipher, eval_key, rlwe_context, rlwe_batch_encoder, rlwe_galois_keys, rlwe_evaluator);
            // end = std::chrono::system_clock::now();
            // compare_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            // std::cout << "As_plus_b time: " << compare_time << std::endl;
        }
        catch(const std::exception& e)
        {
            std::cout << "Location: " << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")\n";
            std::cerr << e.what() << '\n';
        }
        
        


        // reset
        Plaintext offset_plain(rlwe_poly_modulus_degree);
        offset_plain[0] = lwe_coeff_modulus[0].value() / 8;
        rlwe_evaluator.sub_plain_inplace(asb_cipher, offset_plain);

        Ciphertext poly_cipher;
        try
        {
            // start = std::chrono::system_clock::now();
            compute_mux_poly(asb_cipher, rlwe_mux_poly, poly_cipher, rlwe_context, rlwe_evaluator, rlwe_relin_keys, rlwe_encryptor, rlwe_decryptor, rlwe_batch_encoder);
            // end = std::chrono::system_clock::now();
            // compare_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            // std::cout << "mux_poly time: " << compare_time << std::endl;
            
        }
        catch(const std::exception& e)
        {
            std::cout << "Location: " << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")\n";
            std::cerr << e.what() << '\n';
        }
        

        // Slot to Coeff
        Ciphertext coeff_cipher;
        

        try
        {
            // start = std::chrono::system_clock::now();
            slot_to_coeff_optimized(poly_cipher, coeff_cipher, pre_plains, rlwe_galois_keys, rlwe_evaluator);
            // end = std::chrono::system_clock::now();
            // compare_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            // std::cout << "Slot-to-coeff time: " << compare_time << std::endl;
            
        }
        catch(const std::exception& e)
        {
            std::cout << "Location: " << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")\n";
            std::cerr << e.what() << '\n';
        }
        

        // Key switch
        Ciphertext cipher_copy = coeff_cipher;
        util::set_zero_poly(rlwe_poly_modulus_degree, coeff_cipher.coeff_modulus_size(), coeff_cipher.data(1));
        util::RNSIter ct_a(cipher_copy.data(1), rlwe_poly_modulus_degree);
        try
        {
            // start = std::chrono::system_clock::now();
            rlwe_evaluator.switch_key_inplace(coeff_cipher, ct_a, static_cast<const KSwitchKeys &>(rlwe_key_switch_key), 0);
            // end = std::chrono::system_clock::now();
            // compare_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            // std::cout << "switch_key_inplace time: " << compare_time << std::endl;
            
        }
        catch(const std::exception& e)
        {
            std::cout << "Location: " << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")\n";
            std::cerr << e.what() << '\n';
        }
        

        Ciphertext key_switch_cipher = coeff_cipher;
        
        rlwe_evaluator.mod_switch_to_inplace(key_switch_cipher, rlwe_context.last_parms_id());

        Ciphertext modulus_switching_cipher = key_switch_cipher;
        seal::Modulus last_modulus = rlwe_parms.coeff_modulus()[0];
        try
        {
            // start = std::chrono::system_clock::now();
            modulus_switching_inplace(modulus_switching_cipher, last_modulus, lwe_coeff_modulus[0]);
            // end = std::chrono::system_clock::now();
            // compare_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            // std::cout << "modulus_switching_inplace time: " << compare_time << std::endl;
            
        }
        catch(const std::exception& e)
        {
            std::cout << "Location: " << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")\n";
            std::cerr << e.what() << '\n';
        }
        

        output_cipher.resize(cipher_size);
        LWECiphertext temp;
        try
        {
            // start = std::chrono::system_clock::now();
            for (size_t i = 0; i < cipher_size; i++)
            {
                sample_extract(modulus_switching_cipher, temp, i, lwe_coeff_modulus[0]);
                output_cipher[i].resize(lwe_context);
                output_cipher[i].parms_id() = lwe_context.last_parms_id();
                *(output_cipher[i].dataB()) = *(temp.dataB());
                std::copy_n(temp.dataA(), lwe_poly_modulus_degree, output_cipher[i].dataA());
            }
            // end = std::chrono::system_clock::now();
            // compare_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            // std::cout << "sample_extract time: " << compare_time << std::endl;
            
        }
        catch(const std::exception& e)
        {
            std::cout << "Location: " << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")\n";
            std::cerr << e.what() << '\n';
        }
        
        
    }



    // may be less noise
    void compare_greater_than_batch_optimized(std::vector<ComparableCipher> &cipher1, std::vector<ComparableRGSWCipher> &cipher2, std::vector<LWECiphertext> &destination,
                                    seal::SEALContext &context, seal::Evaluator &evaluator, seal::Decryptor &decryptor)
    {
        using namespace seal;

        size_t cipher_size = cipher1.size();
        destination.resize(cipher_size);
        auto parms_id = context.first_parms_id();
        auto context_data_ptr = context.get_context_data(parms_id);
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t poly_modulus_degree = parms.poly_modulus_degree();
        Modulus plain_modulus = parms.plain_modulus();
        Plaintext offset_plaintext(poly_modulus_degree),test_plaintext(poly_modulus_degree);

        for (size_t i = 0; i < poly_modulus_degree; i++)
        {
            test_plaintext[i] = plain_modulus.value() - 1;
        }

        seal::Ciphertext compare_result1, is_equal_ciphertext, compare_result2, rlwe_result;
        for (size_t i = 0; i < cipher_size; i++)
        {
            // get compare_result1 = a0 > b0 (0 or q/6)
            external_product(cipher1[i][0], cipher2[i][0], compare_result1, context, evaluator);
            evaluator.multiply_plain_inplace(compare_result1, test_plaintext); // may be optimized
            offset_plaintext[0] = 1;
            evaluator.add_plain_inplace(compare_result1, offset_plaintext);

            // get is_equal_ciphertext = a1 == b1 (0 or q/12)
            external_product(cipher1[i][1], cipher2[i][1], is_equal_ciphertext, context, evaluator);

            // get compare_result2 = a1 > b1 (0 or q/6)
            evaluator.multiply_plain(is_equal_ciphertext, test_plaintext, compare_result2); // may be optimized
            offset_plaintext[0] = 1;
            evaluator.add_plain_inplace(compare_result2, offset_plaintext);

            // Compute 6 * (a1 == b1) + (a1 > b1 - a0 > b0)
            offset_plaintext[0] = 6;
            evaluator.multiply_plain(is_equal_ciphertext, offset_plaintext, rlwe_result);
            evaluator.add_inplace(rlwe_result, compare_result2);
            evaluator.sub_inplace(rlwe_result, compare_result1);
            evaluator.mod_switch_to_next_inplace(rlwe_result);
            destination[i].extract(context, rlwe_result, 0);
        }

    }

} // namespace arcedb

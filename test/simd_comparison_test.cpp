#include "ARCEDB/comparison/batch_bootstrap.h"
#include "seal/util/uintarithsmallmod.h"
#include "ARCEDB/comparison/comparable.h"
#include "ARCEDB/comparison/batch_comparable.h"
#include "ARCEDB/comparison/data.h"
#include "ARCEDB/utils/multi_thread.h"
#include "ARCEDB/conversion/packlwes.h"
#include "ARCEDB/comparison/poly_eval.h"
#include <chrono>
#include <random>
#include <cassert>
using namespace seal;
using namespace arcedb;


void compare_batch_16_test()
{
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "16-bit precision comparison operator SIMDArbHCMP" << std::endl;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    EncryptionParameters lwe_parms(scheme_type::bfv);
    size_t lwe_poly_modulus_degree = 4096;

    // LWE (12, 4096, {17, 32, 60})
    lwe_parms.set_poly_modulus_degree(lwe_poly_modulus_degree);
    std::vector<Modulus> lwe_coeff_modulus = CoeffModulus::Create(lwe_poly_modulus_degree, {32,60});
    lwe_coeff_modulus.insert(lwe_coeff_modulus.begin(), Modulus(65537));
    lwe_parms.set_coeff_modulus(lwe_coeff_modulus);
    lwe_parms.set_plain_modulus(8);
    SEALContext lwe_context(lwe_parms);
    KeyGenerator lwe_keygen(lwe_context);
    SecretKey lwe_secret_key = lwe_keygen.secret_key();

    Encryptor lwe_encryptor(lwe_context, lwe_secret_key);
    Evaluator lwe_evaluator(lwe_context);
    Decryptor lwe_decryptor(lwe_context, lwe_secret_key);

    std::uniform_int_distribution<uint64_t> lwe_message_1(0, lwe_poly_modulus_degree - 1);
    std::uniform_int_distribution<uint64_t> lwe_message_2(1, lwe_poly_modulus_degree - 1);
    
    uint64_t num = 32768;
    std::vector<ComparableCipher> cipher1(num);
    std::vector<ComparableRGSWCipher> cipher2(num);
    std::vector<uint64_t> linear_result_values(num);
    std::vector<uint64_t> linear_compare_values(num);
    std::vector<uint64_t> compare_result_values(num);
    #pragma omp parallel for num_threads(48)
    for (size_t i = 0; i < num; i++)
    {

        uint64_t value11 = lwe_message_1(engine); 
        uint64_t value12 = lwe_message_2(engine);
        // uint64_t value12 = 0;
        uint64_t value21 = lwe_message_1(engine); 
        uint64_t value22 = lwe_message_2(engine);
        linear_result_values[i] = (value12 == value22) + (value12 > value22) * 2 + (value11 > value21);
        linear_result_values[i] = linear_result_values[i] % 4 ;
        linear_compare_values[i] = linear_result_values[i] <=1 ? 0 : 1;
        compare_result_values[i] = value11 + value12 * lwe_poly_modulus_degree > value21 + value22 * lwe_poly_modulus_degree ? 1 : 0;
        comparable_encrypt_bfv(value11 + value12 * lwe_poly_modulus_degree, cipher1[i], 1, lwe_context, lwe_encryptor);
        comparable_rgsw_encrypt(value21 + value22 * lwe_poly_modulus_degree, cipher2[i], true, lwe_encryptor, lwe_secret_key, lwe_context);
    }
    
    // Batch Bootstrap

    size_t rlwe_poly_modulus_degree = 32768;
    EncryptionParameters rlwe_parms(scheme_type::bfv);
    
    rlwe_parms.set_poly_modulus_degree(rlwe_poly_modulus_degree);
    rlwe_parms.set_coeff_modulus(CoeffModulus::Create(rlwe_poly_modulus_degree, {60,60,60,60,60,60,60,60,60,60,60,60}));
    rlwe_parms.set_plain_modulus(lwe_coeff_modulus[0].value());
    SEALContext rlwe_context(rlwe_parms);
    KeyGenerator rlwe_keygen(rlwe_context);
    SecretKey rlwe_secret_key = rlwe_keygen.secret_key();
    GaloisKeys rlwe_galois_keys;
    RelinKeys rlwe_relin_keys;
    rlwe_keygen.create_galois_keys(rlwe_galois_keys);
    rlwe_keygen.create_relin_keys(rlwe_relin_keys);
    Encryptor rlwe_encryptor(rlwe_context, rlwe_secret_key);
    Evaluator rlwe_evaluator(rlwe_context);
    Decryptor rlwe_decryptor(rlwe_context, rlwe_secret_key);
    BatchEncoder rlwe_batch_encoder(rlwe_context);

    // Linear Transform key generate
    LinearTransformKey eval_key;
    std::vector<int64_t> lwe_key_values;
    convert_secert_key(lwe_secret_key, lwe_context, lwe_key_values);
    generate_linear_transform_key(eval_key, lwe_key_values, rlwe_secret_key, rlwe_context, rlwe_batch_encoder, rlwe_encryptor, rlwe_evaluator);

    std::vector<std::vector<uint64_t>> matrix;
    generate_slot_to_coeff_matrix(matrix, rlwe_context, rlwe_batch_encoder);
    // read_slot_to_coeff_matrix(matrix, rlwe_batch_encoder);
    std::vector<std::vector<Plaintext>> pre_plains;

    auto context_data = rlwe_context.first_context_data();
    while (context_data)
    {
        if (context_data->parms().coeff_modulus().size() == 3)
        {
            break;
        }
        context_data = context_data->next_context_data();
    }
    generate_slot_to_coeff_plain(context_data->parms_id(), matrix, pre_plains, rlwe_context, rlwe_batch_encoder, rlwe_evaluator);
    matrix.clear();


    SecretKey pad_lwe_secret_key;
    KSwitchKeys rlwe_keys_switch_key;
    try
    {
        generate_key_switch_key(lwe_secret_key, lwe_context, rlwe_secret_key, rlwe_context, pad_lwe_secret_key, rlwe_keys_switch_key);
        std::cout << "generate_key_switch_key passed" <<std::endl;
    }
    catch(const std::exception& e)
    {
        std::cout << "Location: " << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")\n";
        std::cerr << e.what() << '\n';
    }

    
    
    std::vector<LWECiphertext> lwe_compare_results;
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    simd_greater_than(cipher1, cipher2, 2, lwe_compare_results, lwe_context, lwe_evaluator, lwe_decryptor, eval_key, 
                        rlwe_context, rlwe_batch_encoder, rlwe_galois_keys, rlwe_evaluator, rlwe_relin_keys, rlwe_encryptor, rlwe_decryptor, 
                        mux_poly, rlwe_keys_switch_key, pre_plains);
    end = std::chrono::system_clock::now();
    double compare_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Comparing " << num << " ciphertexts takes " << compare_time << " ms"<< std::endl;
    std::cout << "The amortized compare time is: " << compare_time / rlwe_poly_modulus_degree << std::endl;
    

    // // Debug
    std::vector<int64_t> lwe_secret_values_debug(lwe_poly_modulus_degree);
    seal::Modulus lwe_q = lwe_parms.coeff_modulus()[0];
    std::vector<int64_t> part_extract_values_debug;
    convert_secert_key(lwe_secret_key, lwe_context, lwe_secret_values_debug);
    for (size_t i = 0; i < num; i++)
    {
        uint64_t debug_value = lwe_decrypt(lwe_compare_results[i], lwe_secret_values_debug, lwe_q, 4);
        // if (debug_value != compare_result_values[i])
        // {
        //     std::cout << "Error " << debug_value << " " << compare_result_values[i] << std::endl;
        // }
    }
    std::cout << "Comparison passed" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
}


int main()
{
    compare_batch_16_test();
    return 0;
}
#include "ARCEDB/comparison/batch_bootstrap.h"
#include "ARCEDB/comparison/poly_eval.h"
#include <chrono>
#include <iostream>
#include <string>
#include <fstream>

namespace arcedb
{

    void generate_slot_to_coeff_matrix(std::vector<std::vector<uint64_t>> &matrix, seal::SEALContext &context, seal::BatchEncoder &batch_encoder)
    {
        size_t coeff_count = batch_encoder.slot_count();
        matrix.resize(coeff_count);
        for (size_t i = 0; i < coeff_count; i++)
        {
            matrix[i].resize(coeff_count);
        }

        
        
        for (size_t i = 0; i < coeff_count; i++)
        {
            std::vector<uint64_t> values(coeff_count), decode_values(coeff_count);
            seal::Plaintext plain(coeff_count);
            plain[i] = 1;
            batch_encoder.decode(plain, values);
            for (size_t j = 0; j < coeff_count; j++)
            {
                matrix[j][i] = values[j];
            }
            
        }

        std::string filename = "../data/UMatrix_" + std::to_string(coeff_count) + ".txt";
        std::ofstream outFile(filename);
        if (outFile.is_open()) 
        {
            for (const auto& row : matrix)
            {
                for (uint64_t value : row) 
                {
                    outFile << value << ' ';
                }
                outFile << '\n';
            }
            outFile.close();
            std::cout << "Data has been written to the file.\n";
        } 
        else 
        {
            throw std::invalid_argument("Unable to open the file for writing.");
        }
        
        
    }

    void read_slot_to_coeff_matrix(std::vector<std::vector<uint64_t>> &matrix, seal::BatchEncoder &batch_encoder)
    {
        size_t coeff_count = batch_encoder.slot_count();
        matrix.resize(coeff_count);
        std::string filename = "../data/UMatrix_" + std::to_string(coeff_count) + ".txt";
        // Read the vector from the text file
        std::ifstream inFile(filename);
        size_t index = 0;
        if (inFile.is_open()) 
        {
            std::cout << "Yes" << std::endl;
            uint64_t value;
            for (size_t i = 0; i < coeff_count; i++)
            {
                matrix[i].resize(coeff_count);
                for (size_t j = 0; j < coeff_count; j++)
                {
                    if (inFile.peek() == '\n') 
                    {
                        inFile.ignore();
                        break;
                    }
                    inFile >> value;
                    matrix[i][j] = value;

                }
                
            }
        }
        inFile.close();

    }

    // The specical case of linear_transform, and the bsgs seems doesn't work, as bfv rotate by N/2
    // But the rotation times is equal in complexity when square matrix.
    void slot_to_coeff(seal::Ciphertext &cipher, seal::Ciphertext &destination, std::vector<std::vector<uint64_t>> &UMatrix,
                        seal::SEALContext &context, seal::BatchEncoder &batch_encoder, seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator)
    {
        using namespace seal;
        size_t coeff_count = UMatrix.size();
        size_t coeff_count_div2 = coeff_count / 2;
        size_t sqrt_ct_size = std::sqrt(coeff_count_div2);
        std::vector<Ciphertext> rotate_ciphers(sqrt_ct_size);
        rotate_ciphers[0] = cipher;
        std::chrono::system_clock::time_point start, end;
        double slot_to_coeff_time = 0;
        start = std::chrono::system_clock::now();
        for (size_t i = 1; i < sqrt_ct_size; i++)
        {
            evaluator.rotate_rows(rotate_ciphers[i-1], 1, galois_keys, rotate_ciphers[i]);
        }
        end = std::chrono::system_clock::now();
        slot_to_coeff_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << slot_to_coeff_time << " ms" << std::endl;
        std::vector<uint64_t> diag(coeff_count);
        Plaintext plain(coeff_count);
        Ciphertext destination1, destination2;
        for (size_t i = 0; i < sqrt_ct_size; i++)
        {
            Ciphertext temp1, sum1, temp2, sum2;
            for (size_t j = 0; j < sqrt_ct_size; j++)
            {
                // Get diagonal
                size_t index = i * sqrt_ct_size + j;
                for (size_t r = 0; r < coeff_count_div2; r++)
                {
                    diag[r] = UMatrix[r % coeff_count_div2][(r + index) % coeff_count_div2];
                }
                for (size_t r = 0; r < coeff_count_div2; r++)
                {
                    diag[r + coeff_count_div2] = UMatrix[r % coeff_count_div2 + coeff_count_div2][(r + index) % coeff_count_div2 + coeff_count_div2];
                }

                std::rotate(diag.rbegin(), diag.rbegin() + i * sqrt_ct_size, diag.rbegin() + coeff_count_div2);
                std::rotate(diag.rbegin() + coeff_count_div2, diag.rbegin() + coeff_count_div2 + i * sqrt_ct_size, diag.rend());
                batch_encoder.encode(diag, plain);
                if (j == 0)
                {
                    evaluator.multiply_plain(rotate_ciphers[j], plain, sum1);
                }
                else
                {
                    evaluator.multiply_plain(rotate_ciphers[j], plain, temp1);
                    evaluator.add_inplace(sum1, temp1);
                }
            }
            for (size_t j = 0; j < sqrt_ct_size; j++)
            {
                // Get diagonal
                size_t index = i * sqrt_ct_size + j;
                for (size_t r = 0; r < coeff_count_div2; r++)
                {
                    diag[r] = UMatrix[r % coeff_count_div2 + coeff_count_div2][(r + index) % coeff_count_div2];
                }
                for (size_t r = 0; r < coeff_count_div2; r++)
                {
                    diag[r + coeff_count_div2] = UMatrix[r % coeff_count_div2][(r + index) % coeff_count_div2 + coeff_count_div2];
                }

                std::rotate(diag.rbegin(), diag.rbegin() + i * sqrt_ct_size, diag.rbegin() + coeff_count_div2);
                std::rotate(diag.rbegin() + coeff_count_div2, diag.rbegin() + coeff_count_div2 + i * sqrt_ct_size, diag.rend());
                batch_encoder.encode(diag, plain);
                if (j == 0)
                {
                    evaluator.multiply_plain(rotate_ciphers[j], plain, sum2);
                }
                else
                {
                    evaluator.multiply_plain(rotate_ciphers[j], plain, temp2);
                    evaluator.add_inplace(sum2, temp2);
                }
            }
            if (i == 0)
            {
                destination1 = sum1;
                destination2 = sum2;
            }
            else
            {
                evaluator.rotate_rows_inplace(sum1, i * sqrt_ct_size, galois_keys);
                evaluator.add_inplace(destination1, sum1);
                evaluator.rotate_rows_inplace(sum2, i * sqrt_ct_size, galois_keys);
                evaluator.add_inplace(destination2, sum2);
            }
        }
        evaluator.rotate_columns_inplace(destination2, galois_keys);
        evaluator.add(destination1, destination2, destination);   
    }

    void generate_slot_to_coeff_plain(seal::parms_id_type parms_id, std::vector<std::vector<uint64_t>> &UMatrix, std::vector<std::vector<seal::Plaintext>> &plains,
                                    seal::SEALContext &context, seal::BatchEncoder &batch_encoder, seal::Evaluator &evaluator)
    {
        using namespace seal;
        size_t coeff_count = UMatrix.size();
        size_t coeff_count_div2 = coeff_count / 2;
        size_t sqrt_ct_size = std::sqrt(coeff_count_div2);
        plains.resize(2);
        plains[0].resize(coeff_count_div2);
        plains[1].resize(coeff_count_div2);
        std::vector<uint64_t> diag(coeff_count);
        for (size_t i = 0; i < sqrt_ct_size; i++)
        {
            for (size_t j = 0; j < sqrt_ct_size; j++)
            {
                // Get diagonal
                size_t index = i * sqrt_ct_size + j;
                for (size_t r = 0; r < coeff_count_div2; r++)
                {
                    diag[r] = UMatrix[r % coeff_count_div2][(r + index) % coeff_count_div2];
                }
                for (size_t r = 0; r < coeff_count_div2; r++)
                {
                    diag[r + coeff_count_div2] = UMatrix[r % coeff_count_div2 + coeff_count_div2][(r + index) % coeff_count_div2 + coeff_count_div2];
                }

                std::rotate(diag.rbegin(), diag.rbegin() + i * sqrt_ct_size, diag.rbegin() + coeff_count_div2);
                std::rotate(diag.rbegin() + coeff_count_div2, diag.rbegin() + coeff_count_div2 + i * sqrt_ct_size, diag.rend());
                batch_encoder.encode(diag, plains[0][index]);
                evaluator.transform_to_ntt_inplace(plains[0][index], parms_id);
            }
            for (size_t j = 0; j < sqrt_ct_size; j++)
            {
                // Get diagonal
                size_t index = i * sqrt_ct_size + j;
                for (size_t r = 0; r < coeff_count_div2; r++)
                {
                    diag[r] = UMatrix[r % coeff_count_div2 + coeff_count_div2][(r + index) % coeff_count_div2];
                }
                for (size_t r = 0; r < coeff_count_div2; r++)
                {
                    diag[r + coeff_count_div2] = UMatrix[r % coeff_count_div2][(r + index) % coeff_count_div2 + coeff_count_div2];
                }

                std::rotate(diag.rbegin(), diag.rbegin() + i * sqrt_ct_size, diag.rbegin() + coeff_count_div2);
                std::rotate(diag.rbegin() + coeff_count_div2, diag.rbegin() + coeff_count_div2 + i * sqrt_ct_size, diag.rend());
                batch_encoder.encode(diag, plains[1][index]);
                evaluator.transform_to_ntt_inplace(plains[1][index], parms_id);
            }
        }


    }

    void slot_to_coeff_optimized(seal::Ciphertext &cipher, seal::Ciphertext &destination, std::vector<std::vector<seal::Plaintext>> &plains,
                        seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator)
    {
        using namespace seal;
        size_t coeff_count = cipher.poly_modulus_degree();
        size_t coeff_count_div2 = coeff_count / 2;
        size_t sqrt_ct_size = std::sqrt(coeff_count_div2);
        std::vector<Ciphertext> rotate_ciphers(sqrt_ct_size);
        rotate_ciphers[0] = cipher;
        std::chrono::system_clock::time_point start, end;
        double slot_to_coeff_time = 0;
        for (size_t i = 1; i < sqrt_ct_size; i++)
        {
            evaluator.rotate_rows(rotate_ciphers[i-1], 1, galois_keys, rotate_ciphers[i]);
            
        }
        for (size_t i = 0; i < sqrt_ct_size; i++)
        {
            evaluator.transform_to_ntt_inplace(rotate_ciphers[i]);
        }
        
        std::vector<uint64_t> diag(coeff_count);
        Plaintext plain(coeff_count);
        Ciphertext destination1, destination2;
        for (size_t i = 0; i < sqrt_ct_size; i++)
        {
            Ciphertext temp1, sum1, temp2, sum2;
            for (size_t j = 0; j < sqrt_ct_size; j++)
            {
                size_t index = i * sqrt_ct_size + j;
                if (j == 0)
                {
                    evaluator.mod_switch_to_inplace(plains[0][index], rotate_ciphers[j].parms_id());
                    evaluator.multiply_plain(rotate_ciphers[j], plains[0][index], sum1);
                }
                else
                {
                    evaluator.mod_switch_to_inplace(plains[0][index], rotate_ciphers[j].parms_id());
                    evaluator.multiply_plain(rotate_ciphers[j], plains[0][index], temp1);
                    evaluator.add_inplace(sum1, temp1);
                }
            }
            for (size_t j = 0; j < sqrt_ct_size; j++)
            {
                size_t index = i * sqrt_ct_size + j;
                if (j == 0)
                {
                    evaluator.mod_switch_to_inplace(plains[1][index], rotate_ciphers[j].parms_id());
                    evaluator.multiply_plain(rotate_ciphers[j], plains[1][index], sum2);
                }
                else
                {
                    evaluator.mod_switch_to_inplace(plains[1][index], rotate_ciphers[j].parms_id());
                    evaluator.multiply_plain(rotate_ciphers[j], plains[1][index], temp2);
                    evaluator.add_inplace(sum2, temp2);
                }
            }
            evaluator.transform_from_ntt_inplace(sum1);
            evaluator.transform_from_ntt_inplace(sum2);
            if (i == 0)
            {
                destination1 = sum1;
                destination2 = sum2;
            }
            else
            {
                evaluator.rotate_rows_inplace(sum1, i * sqrt_ct_size, galois_keys);
                evaluator.add_inplace(destination1, sum1);
                evaluator.rotate_rows_inplace(sum2, i * sqrt_ct_size, galois_keys);
                evaluator.add_inplace(destination2, sum2);
            }
        }
        evaluator.rotate_columns_inplace(destination2, galois_keys);
        evaluator.add(destination1, destination2, destination);  
    }

    void generate_slot_to_coeff_galois_keys(seal::SEALContext&context, seal::KeyGenerator &keygen, seal::GaloisKeys &galois_keys)
    {
        size_t N = context.key_context_data()->parms().poly_modulus_degree();
        size_t n_div_2 = N / 2;
        size_t sqrt_size = std::sqrt(n_div_2);
        std::vector<int> step;
        for (size_t i = 0; i < sqrt_size; i++)
        {
            step.push_back(i);
            step.push_back(i * sqrt_size);
        }
        keygen.create_galois_keys(step, galois_keys);
        
    }

    void pack_encode(std::vector<uint64_t> &input, seal::Plaintext &plain, seal::BatchEncoder &batch_encoder)
    {
        size_t slot_count = batch_encoder.slot_count();
        size_t input_size = input.size();
        if (input_size <= slot_count)
        {
            int step_size = slot_count / input_size;
            std::vector<uint64_t> plain_input(slot_count, 0);
            for (size_t i = 0; i < slot_count; i++)
            {
                plain_input[i] = input[i % input_size];
            }
            batch_encoder.encode(plain_input, plain);
        }
        else
        {
            throw std::invalid_argument("Out of size.");
        }
    }

    void convert_secert_key(seal::SecretKey &lwe_key, seal::SEALContext &context, std::vector<int64_t> &key_values)
    {
        auto &context_data = *context.get_context_data(context.key_parms_id());
        size_t poly_modulus_degree = context_data.parms().poly_modulus_degree();
        auto key_ntt_tables = iter(context_data.small_ntt_tables());
        std::vector<uint64_t> non_ntt(poly_modulus_degree);
        std::copy_n(lwe_key.data().data(), poly_modulus_degree, non_ntt.data());
        seal::util::inverse_ntt_negacyclic_harvey_lazy(non_ntt.data(), key_ntt_tables[0]);
        key_values.resize(poly_modulus_degree);
        uint64_t q = key_ntt_tables[0].modulus().value();
        uint64_t half_q = q >> 1;
        for (size_t i = 0; i < poly_modulus_degree; i++)
        {
            key_values[i] = non_ntt[i] >= half_q ? non_ntt[i] - q : non_ntt[i];
        }
        
    }

    // Is not for general, just for t == q
    void generate_linear_transform_key(LinearTransformKey &eval_key, std::vector<int64_t> &lwe_key_values, seal::SecretKey &rlwe_key, 
                                        seal::SEALContext &context, BatchEncoder &batch_encoder, Encryptor &encryptor, Evaluator &evaluator)
    {
        using namespace seal;
        eval_key.key.resize(context, 2);

        size_t n = lwe_key_values.size();

        std::vector<uint64_t> slots(n, 0);
        // The secret key is ntt form
        auto &context_data = *context.get_context_data(context.key_parms_id());
        auto plain_modulus = context_data.plain_ntt_tables()->modulus();

        for (size_t i = 0; i < n; i++)
        {
            slots[i] = lwe_key_values[i] >= 0 ? lwe_key_values[i] : lwe_key_values[i] + plain_modulus.value();
        }
        
        Plaintext plain_sk;
        pack_encode(slots, plain_sk, batch_encoder);
        encryptor.encrypt_symmetric(plain_sk, eval_key.key);


        size_t g = CeilSqrt(n);
        eval_key.rotated_keys.resize(g);
        eval_key.rotated_keys[0] = eval_key.key;
        for (size_t j = 1; j < g; ++j)
        {
            std::rotate(slots.begin(), slots.begin() + 1, slots.end());
            pack_encode(slots, plain_sk, batch_encoder);
            encryptor.encrypt_symmetric(plain_sk, eval_key.rotated_keys[j]);
        }

        for (size_t i = 0; i < g; i++)
        {
            evaluator.transform_to_ntt_inplace(eval_key.rotated_keys[i]);
        }
        
        
        seal::util::seal_memzero(slots.data(), slots.size() * sizeof(slots[0]));

    }

    // It seems the columns of the matrix should less than N / 2
    void LinearTransform(seal::Ciphertext &result, std::vector<std::vector<uint64_t>> &matrix, LinearTransformKey &eval_key, 
                        seal::SEALContext &context, seal::BatchEncoder &batch_encoder, seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator)
    {
        using namespace seal;
        using namespace seal::util;
        
        MemoryPoolHandle pool = MemoryManager::GetPool(mm_prof_opt::mm_force_thread_local, true);
        size_t rows = matrix.size();
        size_t columns = matrix.front().size();
        size_t slot_counts = batch_encoder.slot_count();
        if (columns > slot_counts)
        {
            throw std::invalid_argument("Convert LWE ciphers out of size.");
        }

        // BSGS Parameters
        size_t max_len = std::max(rows, columns);
        size_t min_len = std::min(rows, columns);
        size_t g_tilde = CeilSqrt(min_len);
        size_t b_tilde = CeilDiv(min_len, g_tilde);

        // Baby-Step
        if (eval_key.rotated_keys.size() < g_tilde)
        {
            std::cout << "There migh be some improvement." << std::endl;
            eval_key.rotated_keys.resize(g_tilde);
            eval_key.rotated_keys[0] = eval_key.key;
            for (size_t i = 1; i < g_tilde; i++)
            {
                evaluator.rotate_rows(eval_key.rotated_keys[i - 1], 1, galois_keys, eval_key.rotated_keys[i]);
            }
        }
        auto &context_data = *context.get_context_data(eval_key.key.parms_id());
        auto &parms = context_data.parms();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_count = parms.poly_modulus_degree();
        size_t coeff_modulus_size = coeff_modulus.size();

        // Giant-Step
        std::vector<uint64_t> diag(max_len);
        seal::Plaintext plain;
        DynArray<uint64_t> poly_lazy1(2 * coeff_count * coeff_modulus_size);
        uint64_t* poly_lazy1_ptr = poly_lazy1.begin();
        DynArray<uint64_t> poly_lazy2(2 * coeff_count * coeff_modulus_size);
        uint64_t* poly_lazy2_ptr = poly_lazy2.begin();
        PolyIter accumulator_iter1(poly_lazy1_ptr, 2, coeff_count);
        PolyIter accumulator_iter2(poly_lazy2_ptr, 2, coeff_count);

        for (size_t b = 0; b < b_tilde && g_tilde * b < min_len; b++)
        {
            seal::Ciphertext temp, sum;
            size_t lazy_reduction_summand_bound = size_t(SEAL_MULTIPLY_ACCUMULATE_USER_MOD_MAX);
            size_t lazy_reduction_counter = lazy_reduction_summand_bound;
            for (size_t g = 0; g < g_tilde && b * g_tilde + g < min_len; g++)
            {
                size_t lazy_reduction_counter = lazy_reduction_summand_bound;
                // Get diagonal
                size_t j = b * g_tilde + g;
                for (size_t r = 0; r < max_len; r++)
                {
                    diag[r] = matrix[r % rows][(r + j) % columns];
                }
                std::rotate(diag.rbegin(), diag.rbegin() + b * g_tilde, diag.rbegin() + max_len / 2);
                std::rotate(diag.rbegin() + max_len / 2, diag.rbegin() + max_len / 2 + b * g_tilde, diag.rend());
                plain.release();
                pack_encode(diag, plain, batch_encoder);
                evaluator.transform_to_ntt_inplace(plain, eval_key.rotated_keys[g].parms_id());
                unsigned long long qword[2]{ 0, 0 };
                for (size_t i = 0; i < coeff_modulus_size; i++)
                {
                    for (size_t j = 0; j < coeff_count; j++)
                    {
                        
                        if (g == 0)
                        {
                            multiply_uint64(*(eval_key.rotated_keys[g].data(0) + i * coeff_count + j), *(plain.data() + i * coeff_count + j), qword);
                            *(poly_lazy1_ptr + i * 2 * coeff_count + 2 * j) = qword[0];
                            *(poly_lazy1_ptr + i * 2 * coeff_count + 2 * j + 1) = qword[1];

                             multiply_uint64(*(eval_key.rotated_keys[g].data(1) + i * coeff_count + j), *(plain.data() + i * coeff_count + j), qword);
                            *(poly_lazy2_ptr + i * 2 * coeff_count + 2 * j) = qword[0];
                            *(poly_lazy2_ptr + i * 2 * coeff_count + 2 * j + 1) = qword[1];
                        }
                        else
                        {
                            multiply_uint64(*(eval_key.rotated_keys[g].data(0) + i * coeff_count + j), *(plain.data() + i * coeff_count + j), qword);
                            util::add_uint128(qword, poly_lazy1_ptr + i * 2 * coeff_count + 2 * j, qword);
                            *(poly_lazy1_ptr + i * 2 * coeff_count + 2 * j) = qword[0];
                            *(poly_lazy1_ptr + i * 2 * coeff_count + 2 * j + 1) = qword[1];

                            multiply_uint64(*(eval_key.rotated_keys[g].data(1) + i * coeff_count + j), *(plain.data() + i * coeff_count + j), qword);
                            util::add_uint128(qword, poly_lazy2_ptr + i * 2 * coeff_count + 2 * j, qword);
                            *(poly_lazy2_ptr + i * 2 * coeff_count + 2 * j) = qword[0];
                            *(poly_lazy2_ptr + i * 2 * coeff_count + 2 * j + 1) = qword[1];
                        }
                    }
                }
                // if (g == 0)
                // {
                    
                //     evaluator.multiply_plain(eval_key.rotated_keys[g], plain, sum);
                // }
                // else
                // {
                //     evaluator.multiply_plain(eval_key.rotated_keys[g], plain, temp);
                //     evaluator.add_inplace(sum, temp);
                // }
            }
            sum.resize(context, eval_key.key.parms_id(), 2);
            sum.is_ntt_form() = true;
            for (size_t i = 0; i < coeff_modulus_size; i++)
            {
                for (size_t j = 0; j < coeff_count; j++)
                {
                    *(sum.data(0) + i * coeff_count + j) =  barrett_reduce_128(poly_lazy1_ptr + i * 2 * coeff_count + 2 * j, coeff_modulus[i]);
                    *(sum.data(1) + i * coeff_count + j) =  barrett_reduce_128(poly_lazy2_ptr + i * 2 * coeff_count + 2 * j, coeff_modulus[i]); 
                }
            }
            evaluator.transform_from_ntt_inplace(sum);
            if (b == 0)
            {
                result = sum;
            }
            else
            {
                evaluator.rotate_rows_inplace(sum, b * g_tilde, galois_keys);
                evaluator.add_inplace(result, sum);
            }
        }

        if (rows < columns)
        {
            size_t gama = std::log2(columns / rows);
            for (size_t j = 0; j < gama; j++)
            {
                seal::Ciphertext temp = result;
                evaluator.rotate_rows_inplace(temp, (1U << j) * rows, galois_keys);
                evaluator.add_inplace(result, temp);
            }
            
        }
    }

    void As_plus_b(std::vector<LWECiphertext> &lwe_ciphers, seal::Ciphertext &middle_result, LinearTransformKey &eval_key, 
                    seal::SEALContext &context, seal::BatchEncoder &batch_encoder, seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator)
    {
        size_t lwe_ciphers_size = lwe_ciphers.size();
        size_t lwe_n = lwe_ciphers[0].poly_modulus_degree();
        std::vector<std::vector<uint64_t>> A(lwe_ciphers_size);
        std::vector<uint64_t> b(lwe_ciphers_size);
        for (size_t i = 0; i < lwe_ciphers_size; i++)
        {
            A[i].resize(lwe_n);
            for (size_t j = 0; j < lwe_n; j++)
            {
                A[i][j] = *(lwe_ciphers[i].dataA() + j);
            }
            b[i] = *lwe_ciphers[i].dataB();
        }
        LinearTransform(middle_result, A, eval_key, context, batch_encoder, galois_keys, evaluator);
        Plaintext b_plain;
        batch_encoder.encode(b, b_plain);
        evaluator.add_plain_inplace(middle_result, b_plain);
    }


    void generate_key_switch_key(seal::SecretKey &lwe_secret_key, seal::SEALContext &lwe_context,
                                seal::SecretKey &rlwe_secret_key, seal::SEALContext &rlwe_context, 
                                seal::SecretKey &pad_lwe_secret_key, seal::KSwitchKeys &key_switch_key)
    {   
        using namespace seal;
        std::vector<int64_t> key_values;
        convert_secert_key(lwe_secret_key, lwe_context, key_values);
        auto &rlwe_context_data = *rlwe_context.key_context_data();
        auto &rlwe_parms = rlwe_context_data.parms();
        auto &rlwe_coeff_modulus = rlwe_parms.coeff_modulus();
        auto ntt_tables = rlwe_context_data.small_ntt_tables();
        size_t rlwe_coeff_count = rlwe_parms.poly_modulus_degree();
        size_t rlwe_coeff_modulus_size = rlwe_coeff_modulus.size();
        std::vector<int64_t> padd_key_values(rlwe_coeff_count, 0);
        
        std::copy_n(key_values.data(), key_values.size(), padd_key_values.data());
        
        pad_lwe_secret_key.data().resize(util::mul_safe(rlwe_coeff_count, rlwe_coeff_modulus_size));

        uint64_t *ptr = pad_lwe_secret_key.data().data();
        for (size_t i = 0; i < rlwe_coeff_modulus_size; i++)
        {
            for (size_t j = 0; j < rlwe_coeff_count; j++)
            {
                *(ptr + i * rlwe_coeff_count + j) = padd_key_values[j] >= 0 ? padd_key_values[j] : rlwe_coeff_modulus[i].value() + padd_key_values[j];
            }
            util::ntt_negacyclic_harvey(ptr + i * rlwe_coeff_count, ntt_tables[i]);
        }
        
        pad_lwe_secret_key.parms_id() = rlwe_context_data.parms_id();

        seal::util::ConstPolyIter secret_key_before(rlwe_secret_key.data().data(), rlwe_coeff_count, rlwe_coeff_modulus_size);
        KeyGenerator lwe_key_generator(rlwe_context, pad_lwe_secret_key);
        lwe_key_generator.generate_kswitch_keys(secret_key_before, 1, key_switch_key, false);
        key_switch_key.parms_id() = rlwe_context.key_parms_id();
    }

    void key_switch_inplace(seal::Ciphertext &cipher, seal::KSwitchKeys &key_switch_keys, seal::Evaluator &evaluator)
    {

        Ciphertext cipher_copy = cipher;
        util::set_zero_poly(cipher.poly_modulus_degree(), cipher.coeff_modulus_size(), cipher.data(1));
        util::RNSIter ct_a(cipher_copy.data(1), cipher.poly_modulus_degree());
        evaluator.switch_key_inplace(cipher, ct_a, static_cast<const KSwitchKeys &>(key_switch_keys), 0);
        cipher_copy.release();
    }

    // Only consider the one modulus
    // Coefficients of cipher ct multiplied by coeff_modulus q2, divided by coeff_modulus q1,
    // and rounded to the nearest integer (rounded up in case of a tie). Equivalent to
    // floor((q2 * m + floor((q1+1) / 2)) / q1).
    void modulus_switching_inplace(seal::Ciphertext &cipher, const seal::Modulus &before_modulus, const seal::Modulus &after_modulus)
    {
        using namespace seal;
        uint64_t modulus_upper_half = (before_modulus.value() + 1) / 2;
        uint64_t *ptr = cipher.data();
        uint64_t q1 = before_modulus.value();
        uint64_t q2 = after_modulus.value();
        for (size_t i = 0; i < cipher.size(); i++)
        {
            for (size_t j = 0; j < cipher.poly_modulus_degree(); j++)
            {
                unsigned long long prod[2]{ 0, 0 };
                uint64_t numerator[2]{ 0, 0 };
                util::multiply_uint64(*(ptr + i * cipher.poly_modulus_degree() + j), q2, prod);
                unsigned char carry = util::add_uint64(*prod, modulus_upper_half, numerator);
                numerator[1] = static_cast<uint64_t>(prod[1]) + static_cast<uint64_t>(carry);

                uint64_t fix[2] = { 0, 0 };
                util::divide_uint128_inplace(numerator, q1, fix);
                *(ptr +  i * cipher.poly_modulus_degree() + j) = fix[0] >= q2 ? fix[0] - q2 : fix[0];
            }
        }
        
        
    }

    void part_extract(seal::Ciphertext &rlwe_cipher, std::vector<LWECiphertext> &lwe_ciphers, seal::SEALContext &rlwe_context, seal::SEALContext &lwe_context)
    {
        using namespace seal;
        auto parms_id = rlwe_cipher.parms_id();
        auto &rlwe_context_data = *rlwe_context.get_context_data(parms_id);
        auto &rlwe_parms = rlwe_context_data.parms();
        size_t rlwe_coeff_count = rlwe_parms.poly_modulus_degree();
        size_t lwe_coeff_count = lwe_context.last_context_data()->parms().poly_modulus_degree();
        lwe_ciphers.resize(rlwe_coeff_count);
        LWECiphertext temp;
        for (size_t i = 0; i < rlwe_coeff_count; i++)
        {
            temp.extract(rlwe_context, rlwe_cipher, i);
            lwe_ciphers[i].resize(lwe_context);
            lwe_ciphers[i].parms_id() = lwe_context.last_parms_id();
            *(lwe_ciphers[i].dataB()) = *(temp.dataB());
            std::copy_n(temp.dataA(), lwe_coeff_count, lwe_ciphers[i].dataA());
        }
    }

    void modulus_switching_lwe_inplace(LWECiphertext &cipher, const seal::Modulus &before_modulus, const seal::Modulus &after_modulus)
    {
        using namespace seal;
        if (cipher.coeff_modulus_size() > 1)
        {
            std::cout << "Location: " << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")\n";
            std::invalid_argument("Invalid LWE cipher, only support one modulus.");
        }
        uint64_t modulus_upper_half = (before_modulus.value() + 1) / 2;
        uint64_t *ptr = cipher.data();
        uint64_t q1 = before_modulus.value();
        uint64_t q2 = after_modulus.value();
        for (size_t i = 0; i < cipher.poly_modulus_degree() + 1; i++)
        {
            unsigned long long prod[2]{ 0, 0 };
            uint64_t numerator[2]{ 0, 0 };
            util::multiply_uint64(*(ptr + i), q2, prod);
            unsigned char carry = util::add_uint64(*prod, modulus_upper_half, numerator);
            numerator[1] = static_cast<uint64_t>(prod[1]) + static_cast<uint64_t>(carry);

            uint64_t fix[2] = { 0, 0 };
            util::divide_uint128_inplace(numerator, q1, fix);
            *(ptr +  i) = fix[0] >= q2 ? fix[0] - q2 : fix[0];
        }
    }
} // namespace arcedb

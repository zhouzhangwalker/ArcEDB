#include "ARCEDB/comparison/rgsw_ciphertext.h"
#include "seal/util/defines.h"
#include "seal/util/rlwe.h"
#include "seal/util/uintarithsmallmod.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/uintcore.h"

namespace arcedb
{
    // Default destination is ntt-form
    void encrypt_rgsw_rns(std::vector<int64_t> values, RGSWCiphertext &destination, 
                        seal::Encryptor &encryptor, seal::SecretKey &secret_key, seal::SEALContext &context)
    {   
        using namespace seal;
        using namespace seal::util;
        MemoryPoolHandle pool = MemoryManager::GetPool(mm_prof_opt::mm_force_thread_local, true);
        auto &context_data = *context.key_context_data();
        auto key_ntt_tables = iter(context_data.small_ntt_tables());
        auto &parms = context_data.parms();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_count = parms.poly_modulus_degree();
        size_t coeff_modulus_size = coeff_modulus.size();
        size_t decomp_mod_count = context.first_context_data()->parms().coeff_modulus().size();
        destination.data().resize(2);
        destination.data()[0].resize(decomp_mod_count);
        destination.data()[1].resize(decomp_mod_count);
        SEAL_ALLOCATE_GET_RNS_ITER(temp, coeff_count, decomp_mod_count, pool);
        for (size_t i = 0; i < decomp_mod_count; i++)
        {
            for (size_t j = 0; j < coeff_count; j++)
            {
                temp[i][j] = values[j] >= 0 ? values[j] : coeff_modulus[i].value() + values[j];
            }
            ntt_negacyclic_harvey_lazy(temp[i], key_ntt_tables[i]); 
        }
        for (size_t i = 0; i < decomp_mod_count; i++)
        {
            util::encrypt_zero_symmetric(secret_key, context, context_data.parms_id(), true, false, destination.data(0)[i].data());
            util::encrypt_zero_symmetric(secret_key, context, context_data.parms_id(), true, false, destination.data(1)[i].data());
            uint64_t factor = util::barrett_reduce_64(coeff_modulus.back().value(), coeff_modulus[i]);
            util::multiply_poly_scalar_coeffmod(temp[i], coeff_count, factor, coeff_modulus[i], temp[i]);
            uint64_t *cipher0_ptr = destination.data(0)[i].data().data(0) + i * coeff_count;
            uint64_t *cipher1_ptr = destination.data(1)[i].data().data(1) + i * coeff_count;
            util::add_poly_coeffmod(cipher0_ptr, temp[i], coeff_count, coeff_modulus[i], cipher0_ptr);
            util::add_poly_coeffmod(cipher1_ptr, temp[i], coeff_count, coeff_modulus[i], cipher1_ptr);
        }
        destination.parms_id() = context_data.parms_id();   
    }
    /**
     * Based on the switch_key_inplace in seal, but has some constraint
    */
    void external_product(seal::Ciphertext &encrypted, RGSWCiphertext &rgsw_cipher, seal::Ciphertext &destination,
                        seal::SEALContext &context, seal::Evaluator &evaluator)
    {
        using namespace seal;
        using namespace seal::util;
        destination.resize(context, encrypted.parms_id(), encrypted.size());
        
        destination.is_ntt_form() = encrypted.is_ntt_form();
        destination.scale() = encrypted.scale();
        MemoryPoolHandle pool = MemoryManager::GetPool(mm_prof_opt::mm_force_thread_local, true);
        auto parms_id = encrypted.parms_id();
        auto &context_data = *context.get_context_data(parms_id);
        auto &parms = context_data.parms();
        
        auto &key_context_data = *context.key_context_data();
        auto &key_parms = key_context_data.parms();
        auto scheme = parms.scheme();
        size_t rgsw_size = rgsw_cipher.size();
        size_t cipher_size = encrypted.size();

        // Verify parameters.
        if (!is_metadata_valid_for(encrypted, context) || !is_buffer_valid(encrypted))
        {
            throw std::invalid_argument("encrypted is not valid for encryption parameters");
        }
        if (!context.using_keyswitching())
        {
            throw std::logic_error("keyswitching is not supported by the context");
        }

        // Don't validate all of kswitch_keys but just check the parms_id.
        if (rgsw_cipher.parms_id() != context.key_parms_id())
        {
            std::cout << encrypted.size();
            std::cout << rgsw_cipher.data()[0].size();
            throw std::invalid_argument("RGSW parameter mismatch");
        }

        if (rgsw_size != 2)
        {
            std::invalid_argument("Invalid rgsw cipher.");
        }

        if (cipher_size != 2)
        {
            throw std::invalid_argument("Invalid cipher size, relinearize first.");
        }

        // Extract encryption parameters.
        size_t coeff_count = parms.poly_modulus_degree();
        size_t decomp_modulus_size = parms.coeff_modulus().size();
        util::set_zero_poly(coeff_count, decomp_modulus_size * 2, destination.data());
        auto &key_modulus = key_parms.coeff_modulus();
        size_t key_modulus_size = key_modulus.size();
        size_t rns_modulus_size = decomp_modulus_size + 1;
        auto key_ntt_tables = iter(key_context_data.small_ntt_tables());
        auto modswitch_factors = key_context_data.rns_tool()->inv_q_last_mod_q();
        bool is_ntt_form = encrypted.is_ntt_form();

        // Size check
        if (!product_fits_in(coeff_count, rns_modulus_size, size_t(2)))
        {
            throw std::logic_error("invalid parameters");
        }

        // Do two key switching and add
        

        // Prepare input
        auto &key_vector = rgsw_cipher.data();
        size_t key_component_count = key_vector[0][0].data().size();


        PolyIter target_iter(encrypted.data(0), coeff_count, decomp_modulus_size);

        // Create a copy of target_iter
        SEAL_ALLOCATE_GET_POLY_ITER(t_target, cipher_size, coeff_count, decomp_modulus_size, pool);
        set_uint(target_iter, cipher_size * decomp_modulus_size * coeff_count, t_target);

        // In CKKS or BGV, t_target is in NTT form; switch back to normal form
        if (is_ntt_form)
        {
            inverse_ntt_negacyclic_harvey(t_target, cipher_size, key_ntt_tables);
        }

        // Temporary result
        auto t_poly_prod(allocate_zero_poly_array(cipher_size, coeff_count, rns_modulus_size, pool));
        SEAL_ITERATE(iter(size_t(0)), rns_modulus_size, [&](auto I) {
            size_t key_index = (I == decomp_modulus_size ? key_modulus_size - 1 : I);

            // Product of two numbers is up to 60 + 60 = 120 bits, so we can sum up to 256 of them without reduction.
            size_t lazy_reduction_summand_bound = size_t(SEAL_MULTIPLY_ACCUMULATE_USER_MOD_MAX);
            size_t lazy_reduction_counter = lazy_reduction_summand_bound;

            // Allocate memory for a lazy accumulator (128-bit coefficients)
            auto t_poly_lazy(allocate_zero_poly_array(key_component_count, coeff_count, 2, pool));

            // Semantic misuse of PolyIter; this is really pointing to the data for a single RNS factor
            PolyIter accumulator_iter(t_poly_lazy.get(), 2, coeff_count);
            SEAL_ITERATE(iter(size_t(0)), cipher_size, [&](auto T){
                // Multiply with keys and perform lazy reduction on product's coefficients
                SEAL_ITERATE(iter(size_t(0)), decomp_modulus_size, [&](auto J) {
                    SEAL_ALLOCATE_GET_COEFF_ITER(t_ntt, coeff_count, pool);
                    ConstCoeffIter t_operand;

                    // RNS-NTT form exists in input
                    if (is_ntt_form && (I == J))
                    {
                        t_operand = target_iter[T][J];
                    }
                    // Perform RNS-NTT conversion
                    else
                    {
                        // No need to perform RNS conversion (modular reduction)
                        if (key_modulus[J] <= key_modulus[key_index])
                        {
                            set_uint(t_target[T][J], coeff_count, t_ntt);
                        }
                        // Perform RNS conversion (modular reduction)
                        else
                        {
                            modulo_poly_coeffs(t_target[T][J], coeff_count, key_modulus[key_index], t_ntt);
                        }
                        // NTT conversion lazy outputs in [0, 4q)
                        ntt_negacyclic_harvey_lazy(t_ntt, key_ntt_tables[key_index]);
                        t_operand = t_ntt;
                    }

                    // Multiply with keys and modular accumulate products in a lazy fashion
                    SEAL_ITERATE(iter(key_vector[T][J].data(), accumulator_iter), key_component_count, [&](auto K) {
                        if (!lazy_reduction_counter)
                        {
                            SEAL_ITERATE(iter(t_operand, get<0>(K)[key_index], get<1>(K)), coeff_count, [&](auto L) {
                                unsigned long long qword[2]{ 0, 0 };
                                multiply_uint64(get<0>(L), get<1>(L), qword);

                                // Accumulate product of t_operand and t_key_acc to t_poly_lazy and reduce
                                add_uint128(qword, get<2>(L).ptr(), qword);
                                get<2>(L)[0] = barrett_reduce_128(qword, key_modulus[key_index]);
                                get<2>(L)[1] = 0;
                            });
                        }
                        else
                        {
                            // Same as above but no reduction
                            SEAL_ITERATE(iter(t_operand, get<0>(K)[key_index], get<1>(K)), coeff_count, [&](auto L) {
                                unsigned long long qword[2]{ 0, 0 };
                                multiply_uint64(get<0>(L), get<1>(L), qword);
                                add_uint128(qword, get<2>(L).ptr(), qword);
                                get<2>(L)[0] = qword[0];
                                get<2>(L)[1] = qword[1];
                            });
                        }
                    });

                    if (!--lazy_reduction_counter)
                    {
                        lazy_reduction_counter = lazy_reduction_summand_bound;
                    }
                });
            });


            // PolyIter pointing to the destination t_poly_prod, shifted to the appropriate modulus
            PolyIter t_poly_prod_iter(t_poly_prod.get() + (I * coeff_count), coeff_count, rns_modulus_size);

            // Final modular reduction
            SEAL_ITERATE(iter(accumulator_iter, t_poly_prod_iter), key_component_count, [&](auto K) {
                if (lazy_reduction_counter == lazy_reduction_summand_bound)
                {
                    SEAL_ITERATE(iter(get<0>(K), *get<1>(K)), coeff_count, [&](auto L) {
                        get<1>(L) = static_cast<uint64_t>(*get<0>(L));
                    });
                }
                else
                {
                    // Same as above except need to still do reduction
                    SEAL_ITERATE(iter(get<0>(K), *get<1>(K)), coeff_count, [&](auto L) {
                        get<1>(L) = barrett_reduce_128(get<0>(L).ptr(), key_modulus[key_index]);
                    });
                }
            });
        });
        
        
        // Accumulated products are now stored in t_poly_prod

        // Perform modulus switching with scaling
        PolyIter t_poly_prod_iter(t_poly_prod.get(), coeff_count, rns_modulus_size);
        SEAL_ITERATE(iter(destination, t_poly_prod_iter), key_component_count, [&](auto I) {
            if (scheme == scheme_type::bgv)
            {
                const Modulus &plain_modulus = parms.plain_modulus();
                // qk is the special prime
                uint64_t qk = key_modulus[key_modulus_size - 1].value();
                uint64_t qk_inv_qp = context.key_context_data()->rns_tool()->inv_q_last_mod_t();

                // Lazy reduction; this needs to be then reduced mod qi
                CoeffIter t_last(get<1>(I)[decomp_modulus_size]);
                inverse_ntt_negacyclic_harvey(t_last, key_ntt_tables[key_modulus_size - 1]);

                SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(k, coeff_count, pool);
                modulo_poly_coeffs(t_last, coeff_count, plain_modulus, k);
                negate_poly_coeffmod(k, coeff_count, plain_modulus, k);
                if (qk_inv_qp != 1)
                {
                    multiply_poly_scalar_coeffmod(k, coeff_count, qk_inv_qp, plain_modulus, k);
                }

                SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(delta, coeff_count, pool);
                SEAL_ALLOCATE_ZERO_GET_COEFF_ITER(c_mod_qi, coeff_count, pool);
                SEAL_ITERATE(iter(I, key_modulus, modswitch_factors, key_ntt_tables), decomp_modulus_size, [&](auto J) {
                    // delta = k mod q_i
                    modulo_poly_coeffs(k, coeff_count, get<1>(J), delta);
                    // delta = k * q_k mod q_i
                    multiply_poly_scalar_coeffmod(delta, coeff_count, qk, get<1>(J), delta);

                    // c mod q_i
                    modulo_poly_coeffs(t_last, coeff_count, get<1>(J), c_mod_qi);
                    // delta = c + k * q_k mod q_i
                    // c_{i} = c_{i} - delta mod q_i
                    SEAL_ITERATE(iter(delta, c_mod_qi), coeff_count, [&](auto K) {
                        get<0>(K) = add_uint_mod(get<0>(K), get<1>(K), get<1>(J));
                    });
                    ntt_negacyclic_harvey(delta, get<3>(J));
                    SEAL_ITERATE(iter(delta, get<0, 1>(J)), coeff_count, [&](auto K) {
                        get<1>(K) = sub_uint_mod(get<1>(K), get<0>(K), get<1>(J));
                    });

                    multiply_poly_scalar_coeffmod(get<0, 1>(J), coeff_count, get<2>(J), get<1>(J), get<0, 1>(J));

                    add_poly_coeffmod(get<0, 1>(J), get<0, 0>(J), coeff_count, get<1>(J), get<0, 0>(J));
                });
            }
            else
            {
                // Lazy reduction; this needs to be then reduced mod qi
                CoeffIter t_last(get<1>(I)[decomp_modulus_size]);
                inverse_ntt_negacyclic_harvey_lazy(t_last, key_ntt_tables[key_modulus_size - 1]);

                // Add (p-1)/2 to change from flooring to rounding.
                uint64_t qk = key_modulus[key_modulus_size - 1].value();
                uint64_t qk_half = qk >> 1;
                SEAL_ITERATE(t_last, coeff_count, [&](auto &J) {
                    J = barrett_reduce_64(J + qk_half, key_modulus[key_modulus_size - 1]);
                });

                SEAL_ITERATE(iter(I, key_modulus, key_ntt_tables, modswitch_factors), decomp_modulus_size, [&](auto J) {
                    SEAL_ALLOCATE_GET_COEFF_ITER(t_ntt, coeff_count, pool);

                    // (ct mod 4qk) mod qi
                    uint64_t qi = get<1>(J).value();
                    if (qk > qi)
                    {
                        // This cannot be spared. NTT only tolerates input that is less than 4*modulus (i.e. qk <=4*qi).
                        modulo_poly_coeffs(t_last, coeff_count, get<1>(J), t_ntt);
                    }
                    else
                    {
                        set_uint(t_last, coeff_count, t_ntt);
                    }

                    // Lazy substraction, results in [0, 2*qi), since fix is in [0, qi].
                    uint64_t fix = qi - barrett_reduce_64(qk_half, get<1>(J));
                    SEAL_ITERATE(t_ntt, coeff_count, [fix](auto &K) { K += fix; });

                    uint64_t qi_lazy = qi << 1; // some multiples of qi
                    if (scheme == scheme_type::ckks)
                    {
                        // This ntt_negacyclic_harvey_lazy results in [0, 4*qi).
                        ntt_negacyclic_harvey_lazy(t_ntt, get<2>(J));
#if SEAL_USER_MOD_BIT_COUNT_MAX > 60
                        // Reduce from [0, 4qi) to [0, 2qi)
                        SEAL_ITERATE(
                            t_ntt, coeff_count, [&](auto &K) { K -= SEAL_COND_SELECT(K >= qi_lazy, qi_lazy, 0); });
#else
                        // Since SEAL uses at most 60bit moduli, 8*qi < 2^63.
                        qi_lazy = qi << 2;
#endif
                    }
                    else if (scheme == scheme_type::bfv)
                    {
                        inverse_ntt_negacyclic_harvey_lazy(get<0, 1>(J), get<2>(J));
                    }

                    // ((ct mod qi) - (ct mod qk)) mod qi with output in [0, 2 * qi_lazy)
                    SEAL_ITERATE(
                        iter(get<0, 1>(J), t_ntt), coeff_count, [&](auto K) { get<0>(K) += qi_lazy - get<1>(K); });

                    // qk^(-1) * ((ct mod qi) - (ct mod qk)) mod qi
                    multiply_poly_scalar_coeffmod(get<0, 1>(J), coeff_count, get<3>(J), get<1>(J), get<0, 1>(J));
                    add_poly_coeffmod(get<0, 1>(J), get<0, 0>(J), coeff_count, get<1>(J), get<0, 0>(J));
                });
            }
        });

    }
} // namespace arcedb

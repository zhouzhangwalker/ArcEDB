#include "poly_eval.h"
#include "seal/util/uintarith.h"
#include "seal/util/uintarithsmallmod.h"
#include "seal/util/uintcore.h"
#include "seal/util/pointer.h"
#include "seal/util/uintarithsmallmod.h"
#include <chrono>
#include <iostream>
#include <string>
#include <fstream>

namespace arcedb
{
    // Compute optimal k,m such that n < k * (2^ m - 1), k is close to sqrt(n / 2)
	// the depth = ceil(log2(k))+m is minimized. the number of homomorphic multiplications = k - 1 +2 * (m - 1) + 2 ^ (m-1) -1 is minimized.
    std::vector<uint32_t> optimal_parameter_paterson(uint32_t degree)
    {
        uint32_t minimum = (1 << 30);
        uint32_t min_k = 0;
        uint32_t min_m;
        uint32_t upper_bound = 2 * static_cast<uint32_t>(floor(sqrt(static_cast<double>(degree))));

        for (uint32_t k = 1; k < upper_bound; k++)
        {
            uint32_t m = 1;
            while (1)
            {
                uint32_t degree_prime = k * ((1 << m) - 1);
                if (degree_prime > degree)
                {
                    uint32_t number_of_mult = (k - 1) + 2 * (m - 1) + ((1 << (m - 1)) - 1);
                    if (number_of_mult < minimum)
                    {
                        minimum = number_of_mult;
                        min_k = k;
                        min_m = m;
                    }
                    break;
                }
                else
                {
                    m++;
                }
            }
        }
        return std::vector<uint32_t>{min_k, min_m};
    }

    std::vector<uint32_t> ComputeDegreesPS(const uint32_t n) 
    {
    // heuristic for larger degrees
        std::vector<uint32_t> klist;
        std::vector<uint32_t> mlist;
        std::vector<uint32_t> multlist;

        for (uint32_t k = 1; k <= n; k++) {
            for (uint32_t m = 1; m <= std::ceil(log2(n / k) + 1) + 1; m++) {
                if (int32_t(n) - int32_t(k * ((1 << m) - 1)) < 0) {
                    if (std::abs(std::floor(log2(k)) - std::floor(log2(sqrt(n / 2)))) <= 1) {
                        klist.push_back(k);
                        mlist.push_back(m);
                        multlist.push_back(k + 2 * m + (1 << (m - 1)) - 4);
                    }
                }
            }
        }
        uint32_t minIndex = std::min_element(multlist.begin(), multlist.end()) - multlist.begin();

        return std::vector<uint32_t>{klist[minIndex], mlist[minIndex]};
    }

    int get_significant_coeff_count_poly(uint64_t *poly, size_t coeff_count, size_t coeff_uint64_count)
    {
        poly += static_cast<std::size_t>(coeff_count - 1) * coeff_uint64_count;
        for (int i = coeff_count; i > 0; i--)
        {
            if (!seal::util::is_zero_uint(poly, coeff_uint64_count))
            {
                return i;
            }
            poly -= coeff_uint64_count;
        }
        return 0;
    }

    // That is to compute X, X^2, X^3, ..., X^k
    void compute_all_powers(seal::Ciphertext &encrypted, uint32_t degree, std::vector<seal::Ciphertext> &destination, 
                            seal::Evaluator &evaluator, seal::RelinKeys &relin_keys)
    {
        using namespace seal;
        destination.resize(degree);
        destination[0] = encrypted;
        Ciphertext temp;
        std::vector<uint64_t> noise_acc(degree, 0);
        noise_acc[0] = 0;
        for (size_t i = 2; i <= degree; i++)
        {
            
            if (!(i & (i - 1)))
            {
                // Compute all X^{2^j}
                evaluator.square(destination[i / 2 - 1], destination[i -1]);
                evaluator.relinearize_inplace(destination[i -1], relin_keys);
                noise_acc[i - 1] = noise_acc[i / 2 - 1] + 1;
                if (noise_acc[i - 1] % 2 == 0)
                {
                    evaluator.mod_switch_to_next_inplace(destination[i-1]);
                }
            }
            else
            {
                // Compute all X*{i} when i is not power of two
                // Question: is there any difference in X^{15} = X{7} * X{8} and X{15} = X{3} * X{12}
                int64_t powerOf2 = 1 << (int64_t)std::floor(std::log2(i));
                int64_t rem      = i % powerOf2; 
                evaluator.mod_switch_to(destination[rem - 1], destination[powerOf2-1].parms_id(), temp);
                evaluator.multiply(destination[powerOf2-1], temp, destination[i-1]);
                evaluator.relinearize_inplace(destination[i-1], relin_keys);
                noise_acc[i - 1] = std::max(noise_acc[powerOf2 - 1], noise_acc[rem - 1]) + 1;
                if (noise_acc[i - 1] % 2 == 0)
                {
                    evaluator.mod_switch_to_next_inplace(destination[i-1]);
                }
            }
        }
        
    }

        // Here, each coefficient is one uint64_t
    // It seems new seal only support coeff_uint64_count = 1?
    void divide_poly_poly_coeffmod_inplace(uint64_t *numerator, uint64_t *denominator, int coeff_count, seal::Modulus &modulus, uint64_t *quotient)
    {
        using namespace seal;
        // Clear quotient
        int coeff_uint64_count = modulus.uint64_count();
        std::memset(quotient, 0, coeff_count * sizeof(uint64_t) * coeff_uint64_count);
        // std::fill(quotient, quotient + coeff_count, 0);

        // Determine most significant coefficients of numerator and denominator.
        int numerator_coeffs = get_significant_coeff_count_poly(numerator, coeff_count, coeff_uint64_count);
        int denominator_coeffs = get_significant_coeff_count_poly(denominator, coeff_count, coeff_uint64_count);

        if (numerator_coeffs < denominator_coeffs)
        {
            return;
        }
        
        int intermediate_uint64_count = coeff_uint64_count * 2;

        // Create scalar to store value that makes denominator monic.
        uint64_t monic_denominator_scalar;

        // Create temporary scalars used during calculation of quotient.
        // Both are purposely twice as wide to store intermediate product prior to modulo operation.
        uint64_t temp_quotient;
        uint64_t subtrahend;

        // Determine scalar necessary to make denominator monic.
        const uint64_t *leading_denominator_coeff = get_poly_coeff(denominator, denominator_coeffs - 1, coeff_uint64_count);
        uint64_t invert_temp;
        if (!util::try_invert_uint_mod(*leading_denominator_coeff, modulus, monic_denominator_scalar))
        {
            throw std::invalid_argument("coeff_modulus is not coprime with leading denominator coefficient");
        }

        // Perform coefficient-wise division algorithm
        while (numerator_coeffs >= denominator_coeffs)
        {
            // Determine leading numerator coefficient.
            const uint64_t *leading_numerator_coeff = get_poly_coeff(numerator, numerator_coeffs - 1, coeff_uint64_count);

            // If leading numerator coefficient is not zero, then need to make zero by subtraction.
            if (!util::is_zero_uint(leading_numerator_coeff, coeff_uint64_count))
            {
                // Determine shift necesarry to bring significant coefficients in alignment.
                int denominator_shift = numerator_coeffs - denominator_coeffs;

                // Determine quotient's coefficient, which is scalar that makes denominator's leading coefficient one
                // multiplied by leading coefficient of denominator (which when subtracted will zero out the topmost
                // denominator coefficient).
                uint64_t *quotient_coeff = get_poly_coeff(quotient, denominator_shift, coeff_uint64_count);
                temp_quotient = util::multiply_uint_mod(monic_denominator_scalar, *leading_numerator_coeff, modulus);
                *quotient_coeff = temp_quotient;

                // Subtract numerator and quotient*denominator (shifted by denominator_shift).
                for (int denominator_coeff_index = 0; denominator_coeff_index < denominator_coeffs; ++denominator_coeff_index)
                {
                    // Multiply denominator's coefficient by quotient.
                    const uint64_t *denominator_coeff = get_poly_coeff(denominator, denominator_coeff_index, coeff_uint64_count);
                    subtrahend = util::multiply_uint_mod(temp_quotient, *denominator_coeff, modulus);

                    // Subtract numerator with resulting product, appropriately shifted by denominator shift.
                    uint64_t *numerator_coeff = get_poly_coeff(numerator, denominator_coeff_index + denominator_shift, coeff_uint64_count);
                    *numerator_coeff = util::sub_uint_mod(*numerator_coeff, subtrahend, modulus);
                }
            }

            // Top numerator coefficient must now be zero, so adjust coefficient count.
            numerator_coeffs--;
        }
    }
    // {x,x^2, ..., x^k} {c_0+c_1*x+c_2*x^2+...+c_k*x^k}
    void poly_evaluation_powers(std::vector<seal::Ciphertext> &all_powers_encrypted, std::vector<uint64_t> &coeff, seal::Ciphertext &destination,
                                seal::Evaluator &evaluator, seal::Decryptor &decryptor, seal::BatchEncoder &batch_encoder)
    {
        using namespace seal;
        size_t degree = coeff.size() - 1;
        // Check whether we have enough all powers
        if (all_powers_encrypted.size() < degree)
        {
            throw std::invalid_argument("there are not enough powers_encrypted");
        }

        auto parms_id = all_powers_encrypted[all_powers_encrypted.size() - 1].parms_id();
        size_t coeff_count = all_powers_encrypted[all_powers_encrypted.size() - 1].poly_modulus_degree();

        if (degree < 2)
        {
            Plaintext plain(coeff_count); 
            plain[0] = coeff[1];
            evaluator.multiply_plain(all_powers_encrypted[0], plain, destination);
            plain[0] = coeff[0];
            evaluator.add_plain_inplace(destination, plain);
            return ;
        }
        int index = all_powers_encrypted.size() - 2;
        
        size_t coeff_modulus_size = all_powers_encrypted[index].coeff_modulus_size();
        size_t encrypted_size = all_powers_encrypted[index].size();
        // Initialize destination

        destination = all_powers_encrypted[index];
        util::set_zero_uint(encrypted_size * coeff_count * coeff_modulus_size, destination.data());
        Plaintext plain(coeff_count);
        Ciphertext temp;
        Plaintext debug_plain;
        std::vector<uint64_t> result_matrix;
        // Add c_0 
        plain[0] = coeff[0];
        evaluator.add_plain_inplace(destination, plain);

        // decryptor.decrypt(destination, debug_plain);
        // batch_encoder.decode(debug_plain, result_matrix);
        // Add to c^{k-1} x^{k-1}
        for (size_t i = 1; i < degree; i++)
        {
            if (coeff[i] != 0)
            {
                plain[0] = coeff[i];
                evaluator.multiply_plain(all_powers_encrypted[i - 1], plain, temp);
                evaluator.mod_switch_to_inplace(temp, destination.parms_id());
                evaluator.add_inplace(destination, temp);
            }
        }
        
        // Add c^{k} x ^{k}
        plain[0] = coeff[degree];
        evaluator.multiply_plain(all_powers_encrypted[degree - 1], plain, temp);
        evaluator.mod_switch_to_inplace(destination, all_powers_encrypted[degree-1].parms_id());
        evaluator.add_inplace(destination, temp);
    }

    /***
    Recursive polynomial evaution methdod,
    First, f(x) = ((X^{k2^{step-1}}) + a(x)) * q(x) + X^(k2^{step-1}-1) + b(x)
    Second, Compute q(x) and X^(k2^{step-1}-1) + b(x) using recursive function call
    Third, the step = 1, just do dot product and return
    The code is based on the implementation from
    https://github.com/openfheorg/openfhe-development/blob/94fd76a1d965cfde13f2a540d78ce64146fc2700/src/pke/lib/scheme/ckksrns/ckksrns-advancedshe.cpp#L369
    */
    void paterson_stockmeyer(seal::Ciphertext &encrypted, std::vector<seal::Ciphertext> &baby_step, std::vector<seal::Ciphertext> &giant_step,
                            std::vector<uint64_t> &coeff, seal::Ciphertext &destination, seal::SEALContext &context,
                            seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, int k, int step,
                            seal::Decryptor &decryptor, seal::BatchEncoder &batch_encoder)
    {
        
        using namespace seal;
        Plaintext debug_plain;
        std::vector<uint64_t> debug_result;
        uint64_t debug_one_result;
        // X^{(2^step -1) * k}
        int current_degree = ((1 << step) - 1) * k;
        if (coeff.size() != current_degree + 1)
        {
            throw std::invalid_argument("Incorrect degree for paterson_stockmeyer");
        }

        if (step == 1)
        {
            poly_evaluation_powers(baby_step, coeff, destination, evaluator, decryptor, batch_encoder);
            return ;
        }
        auto parms_id = encrypted.parms_id();
        auto &context_data = *context.get_context_data(parms_id);
        auto &parms = context_data.parms();
        auto plain_modulus = parms.plain_modulus();
        uint64_t plain_p = parms.plain_modulus().value();

        // Set q(x) and r(x) to be a f(x) = X^{2^{step-1} * k} * q(x) + r(x)
        // deg(q(x)) <= deg(r(x))
        uint64_t temp_degree = (1ULL << (step-1)) * k - 1; // r(x)
        uint64_t next_degree = ((1ULL << (step-1)) - 1) * k; // q(x)
        std::vector<uint64_t> qx(next_degree + 1);
        std::vector<uint64_t> rx(temp_degree + 1); 
        for (size_t i = 0; i < coeff.size(); i++)
        {
            if (i < temp_degree + 1)
            {
                rx[i] = coeff[i];
            }
            else
            {
                qx[i - temp_degree - 1] = coeff[i];
            }
        }



        // r(x) = r(x) - X^{(2^{step-1}-1) * k}
        if (rx[next_degree] > 0)
        {
            rx[next_degree] -= 1;
        }
        else
        {
            rx[next_degree] += plain_p - 1;
        }

        // reset qx format to get the same dimension of rx
        std::vector<uint64_t> poly_qx(rx.size(), 0);
        util::set_uint(qx.data(), qx.size(), poly_qx.data());

        std::vector<uint64_t> poly_cx(rx.size(), 0);
        std::vector<uint64_t> poly_sx(rx.size(), 0);

        // r(x) = c(x) * q(x) + s(x)
        std::chrono::system_clock::time_point start, end;
        divide_poly_poly_coeffmod(rx.data(), poly_qx.data(), rx.size(), plain_modulus, poly_cx.data(), poly_sx.data());
        end = std::chrono::system_clock::now();
        

        int cx_size = get_significant_coeff_count_poly(poly_cx.data(), poly_cx.size(), plain_modulus.uint64_count());
        std::vector<uint64_t> cx(cx_size, 0);
        util::set_uint(poly_cx.data(), cx_size, cx.data());

        // Evaluate cx using baby step
        Ciphertext encrypted_cx;
        
        start = std::chrono::system_clock::now();
        poly_evaluation_powers(baby_step, cx, encrypted_cx, evaluator, decryptor, batch_encoder);
        end = std::chrono::system_clock::now();
        

        // Compute cx + X^{2^{step-1} * k}
        Ciphertext encrypted_cx_add_power_of_two(encrypted_cx);
        evaluator.mod_switch_to_inplace(encrypted_cx_add_power_of_two, giant_step[step - 2].parms_id());
        evaluator.add_inplace(encrypted_cx_add_power_of_two, giant_step[step - 2]);
        // decryptor.decrypt(encrypted_cx_add_power_of_two, debug_plain);
        // batch_encoder.decode(debug_plain, debug_result);
        // Compute sx = sx + X^{(2^{step-1} - 1) * k}
        std::vector<uint64_t> sx(next_degree + 1, 0);
        int sx_size = get_significant_coeff_count_poly(poly_sx.data(), poly_cx.size(), plain_modulus.uint64_count());
        util::set_uint(poly_sx.data(), sx_size, sx.data());
        sx[next_degree] = 1;

        // Evaluate q(x) recursively
        Ciphertext encrypted_qx;
        paterson_stockmeyer(encrypted, baby_step, giant_step, qx, encrypted_qx, context, evaluator, relin_keys, k, step - 1, decryptor, batch_encoder);

        Ciphertext encrypted_sx;
        paterson_stockmeyer(encrypted, baby_step, giant_step, sx, encrypted_sx, context, evaluator, relin_keys, k, step - 1, decryptor, batch_encoder);


        Ciphertext encrypted_part;
        if (encrypted_qx.coeff_modulus_size() >= encrypted_cx_add_power_of_two.coeff_modulus_size())
        {
            evaluator.mod_switch_to_inplace(encrypted_qx, encrypted_cx_add_power_of_two.parms_id());
            
        }
        else
        {
            evaluator.mod_switch_to_inplace(encrypted_cx_add_power_of_two, encrypted_qx.parms_id());
            std::cout << "Wow" << std::endl;
        }
        evaluator.multiply(encrypted_cx_add_power_of_two, encrypted_qx, encrypted_part);
        evaluator.relinearize_inplace(encrypted_part, relin_keys);

        //Compute f(x) = [((X^{k2^{step-1}}) + c(x)) * q(x)] + X^(k2^{step-1}-1) + s(x)
        destination = encrypted_part;
        if (encrypted_sx.coeff_modulus_size() >= destination.coeff_modulus_size())
        {
            evaluator.mod_switch_to_inplace(encrypted_sx, destination.parms_id());

        }
        else
        {
            evaluator.mod_switch_to_inplace(destination, encrypted_sx.parms_id());
            std::cout << "Wow" << std::endl;
        }
        evaluator.add_inplace(destination, encrypted_sx);



    }

    void poly_evaluation_paterson_stockmeyer(seal::Ciphertext &encrypted, std::vector<uint64_t> &coeff, seal::Ciphertext &destination,
                        seal::SEALContext &context, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, seal::Encryptor &encryptor,
                        seal::Decryptor &decryptor, seal::BatchEncoder &batch_encoder)
    {
        using namespace seal;
        uint32_t degree = coeff.size() - 1;
        Plaintext debug_plain;
        std::vector<uint64_t> debug_result;
        // Compute parameter k, m

        std::vector<uint32_t> degs = optimal_parameter_paterson(degree);
        uint32_t k = degs[0];
        uint32_t m = degs[1];
        int degree_prime = ((1 << m) - 1) * k;
        // If degree is low, just compute all powers and product
        uint32_t number_of_mult = (k - 1) + 2 * (m - 1) + (1 << (m - 1)) - 1;
        if (degree < number_of_mult || degree < 3)
        {
            std::vector<seal::Ciphertext> power_of_all(degree);
            // Compute all powers X, X^2, ..., X^{degree-1}
            compute_all_powers(encrypted, degree, power_of_all, evaluator, relin_keys);
            // Polynomial evaluation
            poly_evaluation_powers(power_of_all, coeff, destination, evaluator, decryptor, batch_encoder);

        }

        // Compute baby step and giant step cipher
        // baby step (x, x^2, ..., x^k), giant step (x^{2k}, ...,x^{2^{m-1}k})
        
        std::vector<Ciphertext> baby_step(k); // k+1?k?
        // Compute encryption of baby step (x, x^2, ..., x^k)
        std::chrono::system_clock::time_point start, end;
        start = std::chrono::system_clock::now();
        compute_all_powers(encrypted, k, baby_step, evaluator, relin_keys);
        end = std::chrono::system_clock::now();
        double compute_all_powers_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        
        start = std::chrono::system_clock::now();
        uint64_t noise_acc_xk = std::ceil(std::log2(k));
        std::vector<uint64_t> noise_acc(m-1, 0);
        std::vector<Ciphertext> giant_step(m -1);
        // Compute encryption of giant step (x^{2k}, ...,x^{2^{m-1}k})
        evaluator.square(baby_step[k-1], giant_step[0]);
        evaluator.relinearize_inplace(giant_step[0], relin_keys);
        noise_acc[0] = noise_acc_xk + 1;
        if (noise_acc[0] % 2 == 0)
        {
            evaluator.mod_switch_to_next_inplace(giant_step[0]);
        }
        for (size_t i = 0; i < m -2; i++)
        {
            evaluator.square(giant_step[i], giant_step[i+1]);
            evaluator.relinearize_inplace(giant_step[i+1], relin_keys);
            noise_acc[i + 1] = noise_acc[i] + 1;
            if (noise_acc[i + 1] % 2 == 0)
            {
                evaluator.mod_switch_to_next_inplace(giant_step[i + 1]);
            }
        }

        end = std::chrono::system_clock::now();
        double giant_step_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::vector<uint64_t> f_prime(degree_prime + 1);
        for (size_t i = 0; i < degree_prime + 1; i++)
        {
            if (i < coeff.size())
            {
                f_prime[i] = coeff[i];
            }
            else if (i < degree_prime)
            {
                f_prime[i] = 0;
            }
            else
            {
                f_prime[i] = 1;
            }
        }  

        // Evaluate f_prime
        Ciphertext encrypted_f_prime;
        start = std::chrono::system_clock::now();
        paterson_stockmeyer(encrypted, baby_step, giant_step, f_prime, encrypted_f_prime, context, evaluator, relin_keys, k, m, decryptor, batch_encoder);
        end = std::chrono::system_clock::now();
        double paterson_stockmeyer_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        // std::cout << paterson_stockmeyer_time << " ms" << std::endl;
        //Compute encryption of X^degree_prime = X^k * X^2k ... X^{2^{m-1}k}
        Ciphertext power_of_degree_prime;
        evaluator.mod_switch_to_inplace(baby_step[k-1], giant_step[0].parms_id());
        evaluator.multiply(baby_step[k-1], giant_step[0], power_of_degree_prime); // k ? k - 1?
        evaluator.relinearize_inplace(power_of_degree_prime, relin_keys);
        for (size_t i = 0; i < m -2; i++)
        {
            evaluator.mod_switch_to_inplace(power_of_degree_prime, giant_step[i + 1].parms_id());
            evaluator.multiply_inplace(power_of_degree_prime, giant_step[i + 1]);
            evaluator.relinearize_inplace(power_of_degree_prime, relin_keys);
        }
        if (noise_acc[m - 2] % 2 == 1)
        {
            evaluator.mod_switch_to_next_inplace(power_of_degree_prime);
        }
        //Substract X^degree_prime from temp
        evaluator.mod_switch_to_inplace(encrypted_f_prime, power_of_degree_prime.parms_id());
        evaluator.sub(encrypted_f_prime, power_of_degree_prime, destination);
        
    }

    // just compute even powers {x^2, x^4, ..., x^degree} 
    // equal to {a^1, a^2, ..., a^{degree/2}}, a = x^2
    void compute_even_powers(seal::Ciphertext &encrypted, uint32_t degree, std::vector<seal::Ciphertext> &destination, 
                            seal::Evaluator &evaluator, seal::RelinKeys &relin_keys)
    {
        using namespace seal;
        destination.resize(degree / 2);
        Ciphertext temp;
        std::vector<uint64_t> noise_acc(degree / 2, 0);
        evaluator.square(encrypted, destination[0]);
        evaluator.relinearize_inplace(destination[0], relin_keys);
        noise_acc[0] = 1;
        for (size_t i = 2; i <= degree / 2; i++)
        {
            
            if (!(i & (i - 1)))
            {
                // Compute all X^{2^j}
                evaluator.square(destination[i / 2 - 1], destination[i -1]);
                evaluator.relinearize_inplace(destination[i -1], relin_keys);
                noise_acc[i - 1] = noise_acc[i / 2 - 1] + 1;
                if (noise_acc[i - 1] % 2 == 0)
                {
                    evaluator.mod_switch_to_next_inplace(destination[i-1]);
                }
            }
            else
            {
                // Compute all X*{i} when i is not power of two
                // Question: is there any difference in X^{15} = X{7} * X{8} and X{15} = X{3} * X{12}
                int64_t powerOf2 = 1 << (int64_t)std::floor(std::log2(i));
                int64_t rem      = i % powerOf2; 
                evaluator.mod_switch_to(destination[rem - 1], destination[powerOf2-1].parms_id(), temp);
                evaluator.multiply(destination[powerOf2-1], temp, destination[i-1]);
                evaluator.relinearize_inplace(destination[i-1], relin_keys);
                noise_acc[i - 1] = std::max(noise_acc[powerOf2 - 1], noise_acc[rem - 1]) + 1;
                if (noise_acc[i - 1] % 2 == 0)
                {
                    evaluator.mod_switch_to_next_inplace(destination[i-1]);
                }
            }
        }
        
    }
    // For special case?
    void optimized_poly_evaluation_paterson_stockmeyer(seal::Ciphertext &encrypted, std::vector<uint64_t> &coeff, seal::Ciphertext &destination,
                        seal::SEALContext &context, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, seal::Encryptor &encryptor,
                        seal::Decryptor &decryptor, seal::BatchEncoder &batch_encoder)
    {
        using namespace seal;
        uint32_t degree = coeff.size() - 1;
        Plaintext debug_plain;
        std::vector<uint64_t> debug_result;
        // Compute parameter k, m, suppose k and m is even
        uint32_t k = 256;
        uint32_t m = 256;

        std::vector<Ciphertext> baby_step, giant_step; // k+1?k?
        // Compute encryption of baby step (x^2, x^4, ..., x^k-2)
        std::chrono::system_clock::time_point start, end;
        start = std::chrono::system_clock::now();
        compute_all_powers(encrypted, k, baby_step, evaluator, relin_keys);
        end = std::chrono::system_clock::now();
        double compute_all_powers_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << compute_all_powers_time << " ms" << std::endl;
        start = std::chrono::system_clock::now();
        compute_all_powers(baby_step[baby_step.size() - 1], m - 1, giant_step, evaluator, relin_keys);
        end = std::chrono::system_clock::now();
        compute_all_powers_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << compute_all_powers_time << " ms" << std::endl;
        // for (size_t i = 0; i < baby_step.size(); i++)
        // {
        //     std::cout << decryptor.invariant_noise_budget(baby_step[i]) << " ";
        // }
        // std::cout<< std::endl;
        size_t coeff_cout = encrypted.poly_modulus_degree();
        Plaintext temp(coeff_cout);
        Ciphertext sum0, temp_relin;
        for (size_t i = 0; i < m; i++)
        {
            Ciphertext level_sum;
            bool flag = false;
            for (size_t j = 0; j < k; j++)
            {
                if (coeff[i * k + j] != 0)
                {
                    temp[0] = coeff[i * k + j];
                    if (!flag)
                    {
                        evaluator.multiply_plain(baby_step[j-1], temp, level_sum);
                        flag = false;
                    }
                    else
                    {
                        Ciphertext cipher_temp;
                        evaluator.multiply_plain(baby_step[j - 1], temp, cipher_temp);
                        evaluator.mod_switch_to_inplace(level_sum, cipher_temp.parms_id());
                        evaluator.add_inplace(level_sum, cipher_temp);
                    }
                }
            }
            if (i == 0)
            {
                destination = level_sum;
            }
            else if (i == 1)
            {
                evaluator.mod_switch_to_inplace(level_sum, giant_step[i - 1].parms_id());
                evaluator.multiply(level_sum, giant_step[i-1], temp_relin);
            }
            else
            {
                evaluator.mod_switch_to_inplace(level_sum, giant_step[i - 1].parms_id());
                evaluator.multiply_inplace(level_sum, giant_step[i-1]);
                evaluator.mod_switch_to_inplace(temp_relin, level_sum.parms_id());
                evaluator.add_inplace(temp_relin, level_sum);
            }
                
                
        }
        evaluator.relinearize_inplace(temp_relin, relin_keys);
        evaluator.mod_switch_to_inplace(destination, temp_relin.parms_id());
        evaluator.add_inplace(destination, temp_relin);
        
        
        
    }

    void compute_mux_poly(seal::Ciphertext &encrypted, std::vector<uint64_t> &coeff, seal::Ciphertext &destination,
                        seal::SEALContext &context, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, seal::Encryptor &encryptor,
                        seal::Decryptor &decryptor, seal::BatchEncoder &batch_encoder)
    {
        using namespace seal;
        uint32_t degree = coeff.size() - 1;
        Plaintext debug_plain;
        std::vector<uint64_t> debug_result;
        // Compute parameter k, m, suppose k and m is even
        uint32_t k = 256;
        uint32_t m = 256;

        std::vector<Ciphertext> baby_step, giant_step; //
        std::chrono::system_clock::time_point start, end;
        start = std::chrono::system_clock::now();
        // {x^2, x^4 ,...,x^256}
        compute_even_powers(encrypted, k, baby_step, evaluator, relin_keys);
        
        end = std::chrono::system_clock::now();
        double compute_all_powers_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        // std::cout << compute_all_powers_time << " ms" << std::endl;
        start = std::chrono::system_clock::now();
        // {x^256, x^512 ,...,x^255 * 256}
        compute_all_powers(baby_step[baby_step.size() - 1], m - 1, giant_step, evaluator, relin_keys);
        decryptor.decrypt(giant_step[1], debug_plain);
        batch_encoder.decode(debug_plain, debug_result);
        // std::cout << debug_result[0] << " " << debug_result[1] << " " << debug_result[2] << " " << debug_result[3] << std::endl;
        end = std::chrono::system_clock::now();
        compute_all_powers_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        // std::cout << compute_all_powers_time << " ms" << std::endl;
        // for (size_t i = 0; i < baby_step.size(); i++)
        // {
        //     std::cout << decryptor.invariant_noise_budget(baby_step[i]) << " ";
        // }
        // std::cout<< std::endl;
        size_t coeff_cout = encrypted.poly_modulus_degree();
        Plaintext temp(coeff_cout);
        Ciphertext sum0, temp_relin;
        for (size_t i = 0; i < m; i++)
        {
            Ciphertext level_sum;
            bool flag = false;
            for (size_t j = 1; j <= k; j++)
            {
                if (coeff[i * k + j] != 0)
                {
                    temp[0] = coeff[i * k + j];
                    if (!flag)
                    {
                        evaluator.multiply_plain(baby_step[j/2-1], temp, level_sum);
                        flag = true;
                    }
                    else
                    {
                        Ciphertext cipher_temp;
                        evaluator.multiply_plain(baby_step[j/2 - 1], temp, cipher_temp);
                        evaluator.mod_switch_to_inplace(level_sum, cipher_temp.parms_id());
                        evaluator.add_inplace(level_sum, cipher_temp);
                    }
                }
            }
            if (i == 0)
            {
                destination = level_sum;
            }
            else if (i == 1)
            {
                evaluator.mod_switch_to_inplace(level_sum, giant_step[i - 1].parms_id());
                evaluator.multiply(level_sum, giant_step[i-1], temp_relin);
            }
            else
            {
                evaluator.mod_switch_to_inplace(level_sum, giant_step[i - 1].parms_id());
                evaluator.multiply_inplace(level_sum, giant_step[i-1]);
                evaluator.mod_switch_to_inplace(temp_relin, level_sum.parms_id());
                evaluator.add_inplace(temp_relin, level_sum);
            }
                
                
        }
        evaluator.relinearize_inplace(temp_relin, relin_keys);
        evaluator.mod_switch_to_inplace(destination, temp_relin.parms_id());
        evaluator.add_inplace(destination, temp_relin);
    }

    

} // namespace arcedb

#include "ARCEDB/comparison/rgsw_ciphertext.h"
#include "ARCEDB/comparison/comparable.h"
#include "ARCEDB/conversion/packlwes.h"
#include "ARCEDB/comparison/batch_bootstrap.h"
#include "ARCEDB/utils/serialize.h"
#include "ARCEDB/conversion/repack.h"
#include <iomanip>
#include <random>
#include <chrono>
#include <omp.h>
#include <unistd.h>

using namespace arcedb;
using namespace seal;

/***
 * TPC-H Query 1
    select
        l_returnflag,
        l_linestatus,
        sum(l_quantity) as sum_qty,
        sum(l_extendedprice) as sum_base_price,
        sum(l_extendedprice * (1 - l_discount)) as sum_disc_price,
        sum(l_extendedprice * (1 - l_discount) * (1 + l_tax)) as sum_charge,
        avg(l_quantity) as avg_qty,
        avg(l_extendedprice) as avg_price,
        avg(l_discount) as avg_disc,
        count(*) as count_order
    from
        lineitem
    where
        l_shipdate <= date '1998-12-01' - interval ':1' day (3)
    group by
        l_returnflag,
        l_linestatus
    order by
        l_returnflag,
        l_linestatus;


*/

void relational_query1(size_t num)
{
    std::cout << "Relational SQL Query1 Test: "<< std::endl;
    std::cout << "--------------------------------------------------------"<< std::endl;
    std::cout << "Records: " << num << std::endl;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = Lvl1;
    TFHESecretKey sk;
    TFHEEvalKey ek;
    using bkP = Lvl01;
    using iksP = Lvl10;
    std::uniform_int_distribution<uint32_t> shipdate_message(10592-100, 10592+100);
    std::uniform_int_distribution<uint64_t> quantity_message(0, 8);
    std::uniform_real_distribution<double> extendedprice_message(0., 8.);
    std::uniform_real_distribution<double> discount_message(0., 1.);
    std::uniform_real_distribution<double> tax_message(0., 1.);
    std::uniform_int_distribution<uint32_t> bianry_message(0, 1);
    ek.emplacebkfft<bkP>(sk);
    ek.emplaceiksk<iksP>(sk);
    ek.emplacebkfft<Lvl02>(sk);

    // Filtering
    std::vector<uint64_t> ship_date(num);
    std::vector<uint64_t> returnflag(num);
    std::vector<uint64_t> linestatus(num);
    std::vector<ComparableLvl1> ship_data_ciphers(num);
    std::vector<TRLWELvl1> returnflag_ciphers(num);
    std::vector<TRLWELvl1> linestatus_ciphers(num);

    std::vector<TRGSWLvl1> shipdate_predicate(2);
    TRGSWLvl1 returnflag_predicate_Y, returnflag_predicate_N, linestatus_predicate_Y, linestatus_predicate_N;
    exponent_encrypt_rgsw<P>(10592, 20, shipdate_predicate, sk, true);
    exponent_encrypt_rgsw<P>(1, returnflag_predicate_Y, sk, true);
    exponent_encrypt_rgsw<P>(0, returnflag_predicate_N, sk, true);
    exponent_encrypt_rgsw<P>(1, linestatus_predicate_Y, sk, true);
    exponent_encrypt_rgsw<P>(0, linestatus_predicate_N, sk, true);


    std::vector<double> quantity(num), extendedprice(num), discount(num), tax(num);
    seal::Ciphertext quantity_cipher, extendedprice_cipher, discount_cipher, tax_cipher;
    for (size_t i = 0; i < num; i++)
    {
        quantity[i] = quantity_message(engine);
        extendedprice[i] = extendedprice_message(engine);
        discount[i] = discount_message(engine);
        tax[i] = tax_message(engine);
    }

    // Start sql evaluation
    std::vector<TLWELvl1> filter_res_YY(num), filter_res_YN(num), filter_res_NY(num), filter_res_NN(num);
    std::vector<TLWELvl2> aggregation_res(num);
    TLWELvl2 count_res;
    
    for (size_t i = 0; i < num; i++)
    {
        // Generate data
        ship_date[i] = shipdate_message(engine);
        returnflag[i] = bianry_message(engine);
        linestatus[i] = bianry_message(engine);
        
        // Encrypt
        exponent_encrypt<P>(ship_date[i], 20, ship_data_ciphers[i], sk);
        exponent_encrypt<P>(returnflag[i], returnflag_ciphers[i], sk);
        exponent_encrypt<P>(linestatus[i], linestatus_ciphers[i], sk);
    }

    std::chrono::system_clock::time_point start, end;
    double filtering_time = 0, aggregation_time;
    start = std::chrono::system_clock::now();

    #pragma omp parallel for num_threads(48)
    for (size_t i = 0; i < num; i++)
    {
        
        TLWELvl1 pre_res_YY, pre_res_YN, pre_res_NY, pre_res_NN;

        // returnflag = Y, linestatus = Y
        less_than_tfhepp(ship_data_ciphers[i], shipdate_predicate, 2, filter_res_YY[i], ek, sk);
        equality_tfhepp(returnflag_ciphers[i], returnflag_predicate_Y, pre_res_YY, sk);
        TFHEpp::HomAND(filter_res_YY[i], pre_res_YY, filter_res_YY[i], ek);
        equality_tfhepp(linestatus_ciphers[i], linestatus_predicate_Y, pre_res_YY, sk);
        lift_and_and(filter_res_YY[i], pre_res_YY, filter_res_YY[i], 29, ek, sk);

        // returnflag = Y, linestatus = N
        less_than_tfhepp(ship_data_ciphers[i], shipdate_predicate, 2, filter_res_YN[i], ek, sk);
        equality_tfhepp(returnflag_ciphers[i], returnflag_predicate_Y, pre_res_YN, sk);
        TFHEpp::HomAND(filter_res_YN[i], pre_res_YN, filter_res_YN[i], ek);
        equality_tfhepp(linestatus_ciphers[i], linestatus_predicate_N, pre_res_YN, sk);
        lift_and_and(filter_res_YN[i], pre_res_YN, filter_res_YN[i], 29, ek, sk);

        // returnflag = N, linestatus = Y
        less_than_tfhepp(ship_data_ciphers[i], shipdate_predicate, 2, filter_res_NY[i], ek, sk);
        equality_tfhepp(returnflag_ciphers[i], returnflag_predicate_N, pre_res_NY, sk);
        TFHEpp::HomAND(filter_res_NY[i], pre_res_NY, filter_res_NY[i], ek);
        equality_tfhepp(linestatus_ciphers[i], linestatus_predicate_Y, pre_res_NY, sk);
        lift_and_and(filter_res_NY[i], pre_res_NY, filter_res_NY[i], 29, ek, sk);

        // returnflag = N, linestatus = N
        less_than_tfhepp(ship_data_ciphers[i], shipdate_predicate, 2, filter_res_NN[i], ek, sk);
        equality_tfhepp(returnflag_ciphers[i], returnflag_predicate_N, pre_res_NN, sk);
        TFHEpp::HomAND(filter_res_NN[i], pre_res_NN, filter_res_NN[i], ek);
        equality_tfhepp(linestatus_ciphers[i], linestatus_predicate_N, pre_res_NN, sk);
        lift_and_and(filter_res_NN[i], pre_res_NN, filter_res_NN[i], 29, ek, sk);
        
    }

    end = std::chrono::system_clock::now();

    filtering_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::vector<uint64_t> plain_filter_res_YY(num, 0), plain_filter_res_YN(num, 0), plain_filter_res_NY(num, 0), plain_filter_res_NN(num, 0);
    uint64_t plain_agg_res_YY = 0, plain_agg_res_YN = 0, plain_agg_res_NY = 0,plain_agg_res_NN = 0;
    for (size_t i = 0; i < num; i++)
    {
        if (ship_date[i] < 10592)
        {
            if (returnflag[i] == 1)
            {
                if (linestatus[i] == 1)
                {
                    plain_filter_res_YY[i] = 1;
                    plain_agg_res_YY += quantity[i];
                }
                else
                {
                    plain_filter_res_YN[i] = 1;
                    plain_agg_res_YN += quantity[i];
                }
            }
            else
            {
                if (linestatus[i] == 1)
                {
                    plain_filter_res_NY[i] = 1;
                    plain_agg_res_NY += quantity[i];
                }
                else
                {
                    plain_filter_res_NN[i] = 1;
                    plain_agg_res_NN += quantity[i];
                }
            }
        }

    }

    std::cout << "Filtering finish" << std::endl;

    std::cout << "Aggregation :" << std::endl;
    uint64_t scale_bits = 29;
    uint64_t modq_bits = 32;
    uint64_t modulus_bits = 45;
    uint64_t repack_scale_bits = modulus_bits + scale_bits - modq_bits;
    uint64_t slots_count = filter_res_YY.size();
    std::cout << "Generating Parameters..." << std::endl;
    seal::EncryptionParameters parms(seal::scheme_type::ckks);
    size_t poly_modulus_degree = 65536;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(seal::CoeffModulus::Create(poly_modulus_degree, {59, 42, 42, 42, 42, 42, 42, 42, 42, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 59}));
    double scale = std::pow(2.0, scale_bits);

    //context instance
    seal::SEALContext context(parms, true, seal::sec_level_type::none);

    //key generation
    seal::KeyGenerator keygen(context);
    seal::SecretKey seal_secret_key = keygen.secret_key();
    seal::RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    seal::GaloisKeys galois_keys;
    keygen.create_galois_keys(galois_keys);
    

    //utils
    seal::Encryptor encryptor(context, seal_secret_key);
    seal::Evaluator evaluator(context);
    seal::Decryptor decryptor(context, seal_secret_key);

    //encoder
    seal::CKKSEncoder ckks_encoder(context);

    

    //generate evaluation key
    std::cout << "Generating Conversion Key..." << std::endl;
    LTPreKey pre_key;
    LWEsToRLWEKeyGen(pre_key, std::pow(2., modulus_bits), seal_secret_key, sk, P::n, ckks_encoder, encryptor, context);


    // conversion
    std::cout << "Starting Conversion..." << std::endl;
    seal::Ciphertext resultYY, resultYN, resultNY, resultNN;
    start = std::chrono::system_clock::now();

    LWEsToRLWE(resultYY, filter_res_YY, pre_key, scale, std::pow(2., modq_bits), std::pow(2., modulus_bits - modq_bits), ckks_encoder, galois_keys, relin_keys, evaluator, context);
    HomRound(resultYY, resultYY.scale(), ckks_encoder, relin_keys, evaluator, decryptor, context);

    LWEsToRLWE(resultYN, filter_res_YN, pre_key, scale, std::pow(2., modq_bits), std::pow(2., modulus_bits - modq_bits), ckks_encoder, galois_keys, relin_keys, evaluator, context);
    HomRound(resultYN, resultYN.scale(), ckks_encoder, relin_keys, evaluator, decryptor, context);
    
    LWEsToRLWE(resultNY, filter_res_NY, pre_key, scale, std::pow(2., modq_bits), std::pow(2., modulus_bits - modq_bits), ckks_encoder, galois_keys, relin_keys, evaluator, context);
    HomRound(resultNY, resultNY.scale(), ckks_encoder, relin_keys, evaluator, decryptor, context);

    LWEsToRLWE(resultNN, filter_res_NN, pre_key, scale, std::pow(2., modq_bits), std::pow(2., modulus_bits - modq_bits), ckks_encoder, galois_keys, relin_keys, evaluator, context);
    HomRound(resultNN, resultNN.scale(), ckks_encoder, relin_keys, evaluator, decryptor, context);
    end = std::chrono::system_clock::now();
    aggregation_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    seal::Plaintext plainYY, plainYN, plainNY, plainNN;
    std::vector<double> computedYY(slots_count), computedYN(slots_count), computedNY(slots_count), computedNN(slots_count);
    decryptor.decrypt(resultYY, plainYY);
    seal::pack_decode(computedYY, plainYY, ckks_encoder);

    decryptor.decrypt(resultYN, plainYN);
    seal::pack_decode(computedYN, plainYN, ckks_encoder);

    decryptor.decrypt(resultNY, plainNY);
    seal::pack_decode(computedNY, plainNY, ckks_encoder);

    decryptor.decrypt(resultNN, plainNN);
    seal::pack_decode(computedNN, plainNN, ckks_encoder);

    double errYY = 0., errYN = 0., errNY = 0., errNN = 0.;
    
    for (size_t i = 0; i < slots_count; ++i)
    {
        errYY += std::abs(computedYY[i] - plain_filter_res_YY[i]);

        errYN += std::abs(computedYN[i] - plain_filter_res_YN[i]);

        errNY += std::abs(computedNY[i] - plain_filter_res_NY[i]);

        errNN += std::abs(computedNN[i] - plain_filter_res_NN[i]);
    }

    printf("Repack YY average error = %f ~ 2^%.1f\n", errYY / slots_count, std::log2(errYY / slots_count));
    printf("Repack YN average error = %f ~ 2^%.1f\n", errYN / slots_count, std::log2(errYN / slots_count));
    printf("Repack NY average error = %f ~ 2^%.1f\n", errNY / slots_count, std::log2(errNY / slots_count));
    printf("Repack NN average error = %f ~ 2^%.1f\n", errNN / slots_count, std::log2(errNN / slots_count));


    /*
        sum(l_quantity) as sum_qty,
        sum(l_extendedprice) as sum_base_price,
        sum(l_extendedprice * (1 - l_discount)) as sum_disc_price,
        sum(l_extendedprice * (1 - l_discount) * (1 + l_tax)) as sum_charge,
        avg(l_quantity) as avg_qty,
        avg(l_extendedprice) as avg_price,
        avg(l_discount) as avg_disc,
        count(*) as count_order
    */
    
    
    Plaintext plain;
    double qd = parms.coeff_modulus()[resultYY.coeff_modulus_size() - 1].value();
    scale = std::pow(2., 42);
    seal::pack_encode(quantity, qd, plain, ckks_encoder);
    encryptor.encrypt_symmetric(plain, quantity_cipher);

    seal::pack_encode(extendedprice, scale, plain, ckks_encoder);
    encryptor.encrypt_symmetric(plain, extendedprice_cipher);

    seal::pack_encode(discount, scale, plain, ckks_encoder);
    encryptor.encrypt_symmetric(plain, discount_cipher);

    seal::pack_encode(tax, scale, plain, ckks_encoder);
    encryptor.encrypt_symmetric(plain, tax_cipher);

    std::cout << "Aggregating quanlity, Taking SUM(quantlity) as an example.." << std::endl;

    Ciphertext sum_qty_cipher_YY, sum_qty_cipher_YN, sum_qty_cipher_NY, sum_qty_cipher_NN;
    start = std::chrono::system_clock::now();
    seal::multiply_and_relinearize(resultYY, quantity_cipher, sum_qty_cipher_YY, evaluator, relin_keys);
    evaluator.rescale_to_next_inplace(sum_qty_cipher_YY);

    seal::multiply_and_relinearize(resultYN, quantity_cipher, sum_qty_cipher_YN, evaluator, relin_keys);
    evaluator.rescale_to_next_inplace(sum_qty_cipher_YN);

    seal::multiply_and_relinearize(resultNY, quantity_cipher, sum_qty_cipher_NY, evaluator, relin_keys);
    evaluator.rescale_to_next_inplace(sum_qty_cipher_NY);

    seal::multiply_and_relinearize(resultNN, quantity_cipher, sum_qty_cipher_NN, evaluator, relin_keys);
    evaluator.rescale_to_next_inplace(sum_qty_cipher_NN);
    int logrow = log2(num);
    
    seal::Ciphertext temp;
    size_t step;
    for (size_t i = 0; i < logrow; i++)
    {
        temp = sum_qty_cipher_YY;
        step = 1 << (logrow - i - 1);
        evaluator.rotate_vector_inplace(temp, step, galois_keys);
        evaluator.add_inplace(sum_qty_cipher_YY, temp);

        temp = sum_qty_cipher_YN;
        step = 1 << (logrow - i - 1);
        evaluator.rotate_vector_inplace(temp, step, galois_keys);
        evaluator.add_inplace(sum_qty_cipher_YN, temp);

        temp = sum_qty_cipher_NY;
        step = 1 << (logrow - i - 1);
        evaluator.rotate_vector_inplace(temp, step, galois_keys);
        evaluator.add_inplace(sum_qty_cipher_NY, temp);

        temp = sum_qty_cipher_NN;
        step = 1 << (logrow - i - 1);
        evaluator.rotate_vector_inplace(temp, step, galois_keys);
        evaluator.add_inplace(sum_qty_cipher_NN, temp);
    }
    end = std::chrono::system_clock::now();
    aggregation_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::vector<double> agg_resultYY(slots_count), agg_resultYN(slots_count), agg_resultNY(slots_count), agg_resultNN(slots_count);
    decryptor.decrypt(sum_qty_cipher_YY, plain);
    seal::pack_decode(agg_resultYY, plain, ckks_encoder);

    decryptor.decrypt(sum_qty_cipher_YN, plain);
    seal::pack_decode(agg_resultYN, plain, ckks_encoder);

    decryptor.decrypt(sum_qty_cipher_NY, plain);
    seal::pack_decode(agg_resultNY, plain, ckks_encoder);

    decryptor.decrypt(sum_qty_cipher_NN, plain);
    seal::pack_decode(agg_resultNN, plain, ckks_encoder);

    std::cout << "--------------------------------------------------------"<< std::endl;
    std::cout << "Query Evaluation Time: " << filtering_time + aggregation_time << " ms" << std::endl;

    std::cout << "Encrypted query result: " << std::endl;
    std::cout << std::setw(16) <<"returnfalg" << "|" << std::setw(16) << "linestatus" << "|" << std::setw(16) << "sum_qty" << std::endl;
    std::cout << std::setw(16) <<"Y" << "|" << std::setw(16) << "Y" << "|" << std::setw(16) << std::round(agg_resultYY[0]) << std::endl;
    std::cout << std::setw(16) <<"Y" << "|" << std::setw(16) << "N" << "|" << std::setw(16) << std::round(agg_resultYN[0]) << std::endl;
    std::cout << std::setw(16) <<"N" << "|" << std::setw(16) << "Y" << "|" << std::setw(16) << std::round(agg_resultNY[0]) << std::endl;
    std::cout << std::setw(16) <<"N" << "|" << std::setw(16) << "N" << "|" << std::setw(16) << std::round(agg_resultNN[0]) << std::endl;

    std::cout << "Plain query result: " << std::endl;
    std::cout << std::setw(16) <<"returnfalg" << "|" << std::setw(16) << "linestatus" << "|" << std::setw(16) << "sum_qty" << std::endl;
    std::cout << std::setw(16) <<"Y" << "|" << std::setw(16) << "Y" << "|" << std::setw(16) << plain_agg_res_YY << std::endl;
    std::cout << std::setw(16) <<"Y" << "|" << std::setw(16) << "N" << "|" << std::setw(16) << plain_agg_res_YN << std::endl;
    std::cout << std::setw(16) <<"N" << "|" << std::setw(16) << "Y" << "|" << std::setw(16) << plain_agg_res_NY << std::endl;
    std::cout << std::setw(16) <<"N" << "|" << std::setw(16) << "N" << "|" << std::setw(16) << plain_agg_res_NN << std::endl;\

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    
    
    
}

/***
 * TPC-H Query 6
 * select
        sum(l_extendedprice * l_discount) as revenue
    from
	    lineitem
    where
        l_shipdate >= date ':1'
        and l_shipdate < date ':1' + interval '1' year
        and l_discount between :2 - 0.01 and :2 + 0.01
        and l_quantity < :3;
    
    consider data \in [10592~10957]
*/

void relational_query6(size_t num)
{
    std::cout << "Relational SQL Query6 Test: "<< std::endl;
    std::cout << "--------------------------------------------------------"<< std::endl;
    std::cout << "Records: " << num << std::endl;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = Lvl1;
    TFHESecretKey sk;
    TFHEEvalKey ek;
    using bkP = Lvl01;
    using iksP = Lvl10;
    std::uniform_int_distribution<uint32_t> shipdate_message(10000, 20000);
    std::uniform_int_distribution<uint32_t> discount_message(19000, 21000);
    std::uniform_int_distribution<uint32_t> quantity_message(20000, 40000);
    std::uniform_int_distribution<uint64_t> revenue_message(0, 100);
    ek.emplacebkfft<bkP>(sk);
    ek.emplaceiksk<iksP>(sk);
    ek.emplacebkfft<Lvl02>(sk);

    // Filtering
    std::vector<uint64_t> ship_date(num);
    std::vector<uint64_t> discount(num), quantity(num);
    std::vector<ComparableLvl1> shipdate_ciphers(num), discount_ciphers(num), quantity_ciphers(num);


    std::vector<TRGSWLvl1> predicate1_cipher(2), predicate2_cipher(2), predicate3_cipher(2), predicate4_cipher(2), predicate5_cipher(2);
    uint64_t predicate1_value = 10592, predicate2_value = 10957;
    uint64_t predicate3_value = 19900, predicate4_value = 20100, predicate5_value = 30000;
    exponent_encrypt_rgsw<P>(predicate1_value, 16, predicate1_cipher, sk, true);
    exponent_encrypt_rgsw<P>(predicate2_value, 16, predicate2_cipher, sk, true);
    exponent_encrypt_rgsw<P>(predicate3_value, 16, predicate3_cipher, sk, true);
    exponent_encrypt_rgsw<P>(predicate4_value, 16, predicate4_cipher, sk, true);
    exponent_encrypt_rgsw<P>(predicate5_value, 16, predicate5_cipher, sk, true);


    // Start sql evaluation
    std::vector<TLWELvl1> filter_res(num);
    std::vector<TLWELvl2> aggregation_res(num);
    TLWELvl2 count_res;

    std::vector<double> revenue(num);

    for (size_t i = 0; i < num; i++)
    {
        revenue[i] = revenue_message(engine);
    }
    
    for (size_t i = 0; i < num; i++)
    {
        // Generate data
        ship_date[i] = shipdate_message(engine);
        discount[i] = discount_message(engine);
        quantity[i] = quantity_message(engine);
        exponent_encrypt<P>(ship_date[i], 16, shipdate_ciphers[i], sk);
        exponent_encrypt<P>(discount[i], 16, discount_ciphers[i], sk);
        exponent_encrypt<P>(quantity[i], 16, quantity_ciphers[i], sk);
    }

    std::chrono::system_clock::time_point start, end;
    double filtering_time = 0, aggregation_time;
    start = std::chrono::system_clock::now();

    #pragma omp parallel for num_threads(48)
    for (size_t i = 0; i < num; i++)
    {
        
        TLWELvl1 pre_res;
        greater_than_tfhepp(shipdate_ciphers[i], predicate1_cipher, shipdate_ciphers[i].size(), filter_res[i], ek, sk);
        less_than_tfhepp(shipdate_ciphers[i], predicate2_cipher, shipdate_ciphers[i].size(), pre_res, ek, sk);
        TFHEpp::HomAND(filter_res[i], pre_res, filter_res[i], ek);
        greater_than_tfhepp(discount_ciphers[i], predicate3_cipher, discount_ciphers[i].size(), pre_res, ek, sk);
        TFHEpp::HomAND(filter_res[i], pre_res, filter_res[i], ek);
        less_than_tfhepp(discount_ciphers[i], predicate4_cipher, discount_ciphers[i].size(), pre_res, ek, sk);
        TFHEpp::HomAND(filter_res[i], pre_res, filter_res[i], ek);
        less_than_tfhepp(quantity_ciphers[i], predicate5_cipher, quantity_ciphers[i].size(), pre_res, ek, sk);
        lift_and_and(filter_res[i], pre_res, filter_res[i], 29, ek, sk);
        
    }
    end = std::chrono::system_clock::now();

    filtering_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::vector<uint64_t> plain_filter_res(num);
    uint64_t plain_agg_res = 0;
    for (size_t i = 0; i < num; i++)
    {
        if (ship_date[i] > predicate1_value && ship_date[i] < predicate2_value && discount[i] > predicate3_value && discount[i] < predicate4_value && quantity[i] < predicate5_value)
        {
            plain_filter_res[i] = 1.;
            plain_agg_res += revenue[i];
        }
        else
        {
            plain_filter_res[i] = 0.;
        }

    }

    std::cout << "Filtering finish" << std::endl;

    std::cout << "Aggregation :" << std::endl;
    uint64_t scale_bits = 29;
    uint64_t modq_bits = 32;
    uint64_t modulus_bits = 45;
    uint64_t repack_scale_bits = modulus_bits + scale_bits - modq_bits;
    uint64_t slots_count = filter_res.size();
    std::cout << "Generating Parameters..." << std::endl;
    seal::EncryptionParameters parms(seal::scheme_type::ckks);
    size_t poly_modulus_degree = 65536;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(seal::CoeffModulus::Create(poly_modulus_degree, {59, 42, 42, 42, 42, 42, 42, 42, 42, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 59}));
    double scale = std::pow(2.0, scale_bits);

    //context instance
    seal::SEALContext context(parms, true, seal::sec_level_type::none);

    //key generation
    seal::KeyGenerator keygen(context);
    seal::SecretKey seal_secret_key = keygen.secret_key();
    seal::RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    seal::GaloisKeys galois_keys;
    keygen.create_galois_keys(galois_keys);
    

    //utils
    seal::Encryptor encryptor(context, seal_secret_key);
    seal::Evaluator evaluator(context);
    seal::Decryptor decryptor(context, seal_secret_key);

    //encoder
    seal::CKKSEncoder ckks_encoder(context);

    

    //generate evaluation key
    std::cout << "Generating Conversion Key..." << std::endl;
    LTPreKey pre_key;
    LWEsToRLWEKeyGen(pre_key, std::pow(2., modulus_bits), seal_secret_key, sk, P::n, ckks_encoder, encryptor, context);


    // conversion
    std::cout << "Starting Conversion..." << std::endl;
    seal::Ciphertext result;
    start = std::chrono::system_clock::now();
    LWEsToRLWE(result, filter_res, pre_key, scale, std::pow(2., modq_bits), std::pow(2., modulus_bits - modq_bits), ckks_encoder, galois_keys, relin_keys, evaluator, context);
    HomRound(result, result.scale(), ckks_encoder, relin_keys, evaluator, decryptor, context);
    end = std::chrono::system_clock::now();
    aggregation_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    seal::Plaintext plain;
    std::vector<double> computed(slots_count);
    decryptor.decrypt(result, plain);
    seal::pack_decode(computed, plain, ckks_encoder);

    double err = 0.;
    
    for (size_t i = 0; i < slots_count; ++i)
    {
        err += std::abs(computed[i] - plain_filter_res[i]);
    }

    printf("Repack average error = %f ~ 2^%.1f\n", err / slots_count, std::log2(err / slots_count));


    // Filter result * data
    seal::Ciphertext revenue_cipher;
    double qd = parms.coeff_modulus()[result.coeff_modulus_size() - 1].value();
    seal::pack_encode(revenue, qd, plain, ckks_encoder);
    encryptor.encrypt_symmetric(plain, revenue_cipher);

    std::cout << "Aggregating price and discount .." << std::endl;
    start = std::chrono::system_clock::now();
    seal::multiply_and_relinearize(result, revenue_cipher, result, evaluator, relin_keys);
    evaluator.rescale_to_next_inplace(result);
    std::cout << "Remian modulus: " << result.coeff_modulus_size() << std::endl;
    int logrow = log2(num);
    
    seal::Ciphertext temp;
    for (size_t i = 0; i < logrow; i++)
    {
        temp = result;
        size_t step = 1 << (logrow - i - 1);
        evaluator.rotate_vector_inplace(temp, step, galois_keys);
        evaluator.add_inplace(result, temp);
    }
    end = std::chrono::system_clock::now();
    aggregation_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::vector<double> agg_result(slots_count);
    decryptor.decrypt(result, plain);
    seal::pack_decode(agg_result, plain, ckks_encoder);

    std::cout << "Query Evaluation Time: " << filtering_time + aggregation_time << " ms" << std::endl;

    std::cout << "Encrypted query result: " << std::endl;
    std::cout << std::setw(12) <<"revenue" << std::endl;
    std::cout << std::setw(12) << std::round(agg_result[0]) << std::endl;
    std::cout << "Plain query result: " << std::endl;
    std::cout << std::setw(12) <<"revenue" << std::endl;
    std::cout << std::setw(12) << plain_agg_res << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

}


/*
    TPC-H Query 12
    select
        l_shipmode,
        sum(case
            when o_orderpriority = '1-URGENT'
                or o_orderpriority = '2-HIGH'
                then 1
            else 0
        end) as high_line_count,
        sum(case
            when o_orderpriority <> '1-URGENT'
                and o_orderpriority <> '2-HIGH'
                then 1
            else 0
        end) as low_line_count
    from
        orders,
        lineitem
    where
        o_orderkey = l_orderkey
        and l_shipmode in (':1', ':2')
        and l_commitdate < l_receiptdate
        and l_shipdate < l_commitdate
        and l_receiptdate >= date ':3'
        and l_receiptdate < date ':3' + interval '1' year
    group by
        l_shipmode
    order by
        l_shipmode;
    Consider the joined table
*/
void relational_query12(size_t num)
{
    std::cout << "Relational SQL Query12 Test: "<< std::endl;
    std::cout << "--------------------------------------------------------"<< std::endl;
    std::cout << "Records: " << num << std::endl;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = Lvl1;
    TFHESecretKey sk;
    TFHEEvalKey ek;
    using bkP = Lvl01;
    using iksP = Lvl10;
    std::uniform_int_distribution<uint32_t> shipmode_message(1, 10);
    std::uniform_int_distribution<uint32_t> shipdate_message(10000, 20000);
    std::uniform_int_distribution<uint32_t> receiptdate_message(10000, 20000);
    std::uniform_int_distribution<uint32_t> commitdate_message(10000, 20000);
    // orderpriority \in ('1-URGENT', '2-HIGH', '3-MEDIUM', '4-NOT SPECIFIED', '5-LOW')
    std::uniform_int_distribution<uint64_t> orderpriority_message(1, 5);
    ek.emplacebkfft<bkP>(sk);
    ek.emplaceiksk<iksP>(sk);
    ek.emplacebkfft<Lvl02>(sk);

    // Filtering
    std::vector<uint64_t> shipdate(num), commitdate(num), receiptdate(num);
    std::vector<uint64_t> orderpriority(num), shipmode(num);
    std::vector<ComparableLvl1> shipdate_ciphers(num), commitdate_ciphers(num), receiptdate_ciphers(num);
    std::vector<TRLWELvl1> shipmode_ciphers(num), orderpriority_ciphers(num);
    std::vector<ComparbleRGSWLvl1> receiptdate_rgsw_ciphers(num), commitdate_rgsw_ciphers(num);

    TRGSWLvl1 predicate_mail_cipher, predicate_ship_cipher; // 'MAIL', 'SHIP'
    TRGSWLvl1 predicate_urgent_cipher, predicate_high_cipher; // '1-URGENT', '2-HIGH'
    std::vector<TRGSWLvl1> predicate_date_cipher1(2), predicate_date_cipher2(2);
    uint64_t predicate_mail = 1, predicate_ship= 2, predicate_urgent = 1, predicate_high = 2;
    uint64_t predicate_date1 = 13000, predicate_date2 = 18000;
    exponent_encrypt_rgsw<P>(predicate_mail, predicate_mail_cipher, sk, true);
    exponent_encrypt_rgsw<P>(predicate_ship, predicate_ship_cipher, sk, true);
    exponent_encrypt_rgsw<P>(predicate_urgent, predicate_urgent_cipher, sk, true);
    exponent_encrypt_rgsw<P>(predicate_high, predicate_high_cipher, sk, true);
    exponent_encrypt_rgsw<P>(predicate_date1, 16, predicate_date_cipher1, sk, true);
    exponent_encrypt_rgsw<P>(predicate_date2, 16, predicate_date_cipher2, sk, true);


    // Start sql evaluation
    std::vector<TLWELvl1> filter_res_mail(num), filter_res_ship(num), order_res(num);
    std::vector<TLWELvl1> res_mail_order(num), res_ship_order(num);
    std::vector<TLWELvl2> count_mail(num), count_ship(num), count_mail_order(num), count_ship_order(num);
    TLWELvl2 agg_mail, agg_ship, agg_mail_order, agg_ship_order;
    

    for (size_t i = 0; i < num; i++)
    {
        // Generate data
        shipdate[i] = shipdate_message(engine);
        commitdate[i] = commitdate_message(engine);
        receiptdate[i] = receiptdate_message(engine);
        shipmode[i] = shipmode_message(engine);
        orderpriority[i] = orderpriority_message(engine);
        exponent_encrypt<P>(shipdate[i], 16, shipdate_ciphers[i], sk);
        exponent_encrypt<P>(commitdate[i], 16, commitdate_ciphers[i], sk);
        exponent_encrypt<P>(receiptdate[i], 16, receiptdate_ciphers[i], sk);
        exponent_encrypt<P>(shipmode[i], shipmode_ciphers[i], sk);
        exponent_encrypt<P>(orderpriority[i], orderpriority_ciphers[i], sk);
        exponent_encrypt_rgsw<P>(receiptdate[i], 16, receiptdate_rgsw_ciphers[i], sk, true);
        exponent_encrypt_rgsw<P>(commitdate[i], 16, commitdate_rgsw_ciphers[i], sk, true);
    }

    std::chrono::system_clock::time_point start, end;
    double filtering_time = 0, aggregation_time;
    start = std::chrono::system_clock::now();

    #pragma omp parallel for num_threads(48)
    for (size_t i = 0; i < num; i++)
    {
        
        TLWELvl1 pre_res;
        less_than_tfhepp(commitdate_ciphers[i], receiptdate_rgsw_ciphers[i], commitdate_ciphers[i].size(), filter_res_mail[i], ek, sk);
        less_than_tfhepp(shipdate_ciphers[i], commitdate_rgsw_ciphers[i], shipdate_ciphers[i].size(), pre_res, ek, sk);
        TFHEpp::HomAND(filter_res_mail[i], pre_res, filter_res_mail[i], ek);
        greater_than_tfhepp(receiptdate_ciphers[i], predicate_date_cipher1, receiptdate_ciphers[i].size(), pre_res, ek, sk);
        TFHEpp::HomAND(filter_res_mail[i], pre_res, filter_res_mail[i], ek);
        less_than_tfhepp(receiptdate_ciphers[i], predicate_date_cipher2, receiptdate_ciphers[i].size(), pre_res, ek, sk);
        TFHEpp::HomAND(filter_res_mail[i], pre_res, filter_res_mail[i], ek);
        filter_res_ship[i] = filter_res_mail[i];
        equality_tfhepp(shipmode_ciphers[i], predicate_mail_cipher, pre_res,sk);
        TFHEpp::HomAND(filter_res_mail[i], pre_res, filter_res_mail[i], ek);

        equality_tfhepp(shipmode_ciphers[i], predicate_ship_cipher, pre_res,sk);
        TFHEpp::HomAND(filter_res_ship[i], pre_res, filter_res_ship[i], ek);

        equality_tfhepp(orderpriority_ciphers[i], predicate_urgent_cipher, order_res[i], sk);
        equality_tfhepp(orderpriority_ciphers[i], predicate_high_cipher, pre_res, sk);
        TFHEpp::HomOR(order_res[i], pre_res, order_res[i], ek);

        lift_and_and(filter_res_mail[i], order_res[i], count_mail_order[i], 48, ek, sk);
        lift_and_and(filter_res_ship[i], order_res[i], count_ship_order[i], 48, ek, sk);
        lift_and_and(filter_res_mail[i], filter_res_mail[i], count_mail[i], 48, ek, sk);
        lift_and_and(filter_res_ship[i], filter_res_ship[i], count_ship[i], 48, ek, sk);
        
    }

    for (size_t i = 0; i < num; i++)
    {
        if (i == 0)
        {
            agg_mail = count_mail[0];
            agg_ship = count_ship[0];
            agg_mail_order = count_mail_order[0];
            agg_ship_order = count_ship_order[0];
        }
        else
        {
            for (size_t j = 0; j <= Lvl2::n; j++)
            {
                agg_mail[j] += count_mail[i][j];
                agg_ship[j] += count_ship[i][j];
                agg_mail_order[j] += count_mail_order[i][j];
                agg_ship_order[j] += count_ship_order[i][j];
            }
            
        }
    }
    
    end = std::chrono::system_clock::now();

    filtering_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::vector<uint64_t> plain_filter_res_mail(num, 0), plain_filter_res_ship(num, 0), plain_filter_order(num, 0);
    std::vector<uint64_t> plain_res_mail_order(num, 0), plain_res_ship_order(num, 0);
    uint64_t agg_mail_res = 0, agg_mail_order_res = 0, agg_ship_res = 0, agg_ship_order_res = 0;
    bool ress;
    for (size_t i = 0; i < num; i++)
    {
        if (commitdate[i] < receiptdate[i] && shipdate[i] < commitdate[i] && receiptdate[i] > predicate_date1 && receiptdate[i] < predicate_date2)
        {
            ress = true;
        }
        else
        {
            ress = false;
        }

        if (orderpriority[i] == 1 || orderpriority[i] == 2)
        {
            plain_filter_order[i] = 1;
        }

        if (ress && shipmode[i] == 1)
        {
            plain_filter_res_mail[i] = 1;
            agg_mail_res += 1;
            if (plain_filter_order[i] == 1)
            {
                plain_res_mail_order[i] = 1;
                agg_mail_order_res += 1;
            }
        }
        if (ress && shipmode[i] == 2)
        {
            plain_filter_res_ship[i] = 1;
            agg_ship_res += 1;
            if (plain_filter_order[i] == 1)
            {
                plain_res_ship_order[i] = 1;
                agg_ship_order_res += 1;
            }
        }
        

    }

    std::cout << "Filtering finish" << std::endl;
    uint64_t query_res_mail = tlweSymInt32Decrypt<Lvl2>(agg_mail, std::pow(2.,48), sk.key.get<Lvl2>());
    uint64_t query_res_ship = tlweSymInt32Decrypt<Lvl2>(agg_ship, std::pow(2.,48), sk.key.get<Lvl2>());
    uint64_t query_res_mail_order = tlweSymInt32Decrypt<Lvl2>(agg_mail_order, std::pow(2.,48), sk.key.get<Lvl2>());
    uint64_t query_res_ship_order = tlweSymInt32Decrypt<Lvl2>(agg_ship_order, std::pow(2.,48), sk.key.get<Lvl2>());

    std::cout << "Query Evaluation Time: " << filtering_time << " ms" << std::endl;
    std::cout << "Encrypted result: " << std::endl;
    std::cout << std::setw(12) <<"shipmode" << "|" << std::setw(16) << "high_line_count" << "|" << std::setw(16) << "low_line_count" << std::endl;
    std::cout << std::setw(12) <<"MAIL" << "|" << std::setw(16) << query_res_mail_order << "|" << std::setw(16) << query_res_mail - query_res_mail_order << std::endl;
    std::cout << std::setw(12) <<"SHIP" << "|" << std::setw(16) << query_res_ship_order << "|" << std::setw(16) << query_res_ship - query_res_ship_order << std::endl;

    std::cout << "Plain result: " << std::endl;
    std::cout << std::setw(12) <<"shipmode" << "|" << std::setw(16) << "high_line_count" << "|" << std::setw(16) << "low_line_count" << std::endl;
    std::cout << std::setw(12) <<"MAIL" << "|" << std::setw(16) << agg_mail_order_res << "|" << std::setw(16) << agg_mail_res - agg_mail_order_res << std::endl;
    std::cout << std::setw(12) <<"SHIP" << "|" << std::setw(16) << agg_ship_order_res << "|" << std::setw(16) << agg_ship_res - agg_ship_order_res << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

}

/*
    select
        100.00 * sum(case
            when p_type like 'PROMO%'
                then l_extendedprice * (1 - l_discount)
            else 0
        end) / sum(l_extendedprice * (1 - l_discount)) as promo_revenue
    from
        lineitem,
        part
    where
        l_partkey = p_partkey
        and l_shipdate >= date ':1'
        and l_shipdate < date ':1' + interval '1' month;
    Consider the joined table
*/
void relational_query14(size_t num)
{
        std::cout << "Relational SQL Query14 Test: "<< std::endl;
    std::cout << "--------------------------------------------------------"<< std::endl;
    std::cout << "Records: " << num << std::endl;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = Lvl1;
    TFHESecretKey sk;
    TFHEEvalKey ek;
    using bkP = Lvl01;
    using iksP = Lvl10;
    std::uniform_int_distribution<uint32_t> shipdate_message(10000, 20000);
    std::uniform_int_distribution<uint32_t> revenue_message(0, 100);
    std::uniform_int_distribution<uint32_t> ptype_message(0, 100);
    ek.emplacebkfft<bkP>(sk);
    ek.emplaceiksk<iksP>(sk);
    ek.emplacebkfft<Lvl02>(sk);

    // Filtering
    std::vector<uint64_t> ship_date(num);
    std::vector<uint64_t> ptype(num);
    std::vector<ComparableLvl1> shipdate_ciphers(num);
    std::vector<TRLWELvl1> ptype_ciphers(num);


    std::vector<TRGSWLvl1> predicate1_cipher(2), predicate2_cipher(2);
    TRGSWLvl1 predicate3_cipher, predicate4_cipher;
    uint64_t predicate1_value = 10592, predicate2_value = 10957;
    uint64_t predicate3_value = 30, predicate4_value = 70;
    exponent_encrypt_rgsw<P>(predicate1_value, 16, predicate1_cipher, sk, true);
    exponent_encrypt_rgsw<P>(predicate2_value, 16, predicate2_cipher, sk, true);
    exponent_encrypt_rgsw<P>(predicate3_value, predicate3_cipher, sk, true);
    exponent_encrypt_rgsw<P>(predicate4_value, predicate4_cipher, sk, true);


    // Start sql evaluation
    std::vector<TLWELvl1> filter_res(num), filter_case_res(num);
    std::vector<TLWELvl2> aggregation_res(num);
    TLWELvl2 count_res;

    std::vector<double> revenue(num);

    for (size_t i = 0; i < num; i++)
    {
        revenue[i] = revenue_message(engine);
    }
    
    #pragma omp parallel for num_threads(48)
    for (size_t i = 0; i < num; i++)
    {
        // Generate data
        ship_date[i] = shipdate_message(engine);
        ptype[i] = ptype_message(engine);
        exponent_encrypt<P>(ship_date[i], 16, shipdate_ciphers[i], sk);
        exponent_encrypt<P>(ptype[i], ptype_ciphers[i], sk);
    }

    std::chrono::system_clock::time_point start, end;
    double filtering_time = 0, aggregation_time;
    start = std::chrono::system_clock::now();

    #pragma omp parallel for num_threads(48)
    for (size_t i = 0; i < num; i++)
    {
        
        TLWELvl1 pre_res;
        greater_than_tfhepp(shipdate_ciphers[i], predicate1_cipher, shipdate_ciphers[i].size(), filter_res[i], ek, sk);
        less_than_tfhepp(shipdate_ciphers[i], predicate2_cipher, shipdate_ciphers[i].size(), pre_res, ek, sk);
        TFHEpp::HomAND(filter_res[i], pre_res, filter_res[i], ek);
        greater_than_tfhepp(ptype_ciphers[i], predicate3_cipher, pre_res, sk);
        TFHEpp::HomAND(filter_case_res[i], pre_res, filter_res[i], ek);
        less_than_tfhepp(ptype_ciphers[i], predicate4_cipher, pre_res, sk);
        lift_and_and(filter_case_res[i], pre_res, filter_case_res[i], 29, ek, sk);
        lift_and_and(filter_res[i], filter_res[i], filter_res[i], 29, ek, sk);
        
    }
    end = std::chrono::system_clock::now();

    filtering_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::vector<uint64_t> plain_filter_res(num), plain_filter_case_res(num);
    uint64_t plain_agg_res = 0, plain_agg_case_res = 0;
    for (size_t i = 0; i < num; i++)
    {
        if (ship_date[i] > predicate1_value && ship_date[i] < predicate2_value)
        {
            plain_filter_res[i] = 1;
            plain_agg_res += revenue[i];
            if (ptype[i] > predicate3_value && ptype[i] < predicate4_value)
            {
                plain_filter_case_res[i] = 1;
                plain_agg_case_res += revenue[i];
            }
            else
            {
                plain_filter_case_res[i] = 0;
            }
            
        }
        else
        {
            plain_filter_res[i] = 0;
            plain_filter_case_res[i] = 0;
        }

    }

    std::cout << "Filtering finish" << std::endl;

    std::cout << "Aggregation :" << std::endl;
    uint64_t scale_bits = 29;
    uint64_t modq_bits = 32;
    uint64_t modulus_bits = 45;
    uint64_t repack_scale_bits = modulus_bits + scale_bits - modq_bits;
    uint64_t slots_count = filter_res.size();
    std::cout << "Generating Parameters..." << std::endl;
    seal::EncryptionParameters parms(seal::scheme_type::ckks);
    size_t poly_modulus_degree = 65536;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(seal::CoeffModulus::Create(poly_modulus_degree, {59, 42, 42, 42, 42, 42, 42, 42, 42, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 59}));
    double scale = std::pow(2.0, scale_bits);

    //context instance
    seal::SEALContext context(parms, true, seal::sec_level_type::none);

    //key generation
    seal::KeyGenerator keygen(context);
    seal::SecretKey seal_secret_key = keygen.secret_key();
    seal::RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    seal::GaloisKeys galois_keys;
    keygen.create_galois_keys(galois_keys);
    

    //utils
    seal::Encryptor encryptor(context, seal_secret_key);
    seal::Evaluator evaluator(context);
    seal::Decryptor decryptor(context, seal_secret_key);

    //encoder
    seal::CKKSEncoder ckks_encoder(context);

    

    //generate evaluation key
    std::cout << "Generating Conversion Key..." << std::endl;
    LTPreKey pre_key;
    LWEsToRLWEKeyGen(pre_key, std::pow(2., modulus_bits), seal_secret_key, sk, P::n, ckks_encoder, encryptor, context);


    // conversion
    std::cout << "Starting Conversion..." << std::endl;
    seal::Ciphertext result, result_case;
    start = std::chrono::system_clock::now();
    LWEsToRLWE(result, filter_res, pre_key, scale, std::pow(2., modq_bits), std::pow(2., modulus_bits - modq_bits), ckks_encoder, galois_keys, relin_keys, evaluator, context);
    HomRound(result, result.scale(), ckks_encoder, relin_keys, evaluator, decryptor, context);

    LWEsToRLWE(result_case, filter_case_res, pre_key, scale, std::pow(2., modq_bits), std::pow(2., modulus_bits - modq_bits), ckks_encoder, galois_keys, relin_keys, evaluator, context);
    HomRound(result_case, result_case.scale(), ckks_encoder, relin_keys, evaluator, decryptor, context);
    end = std::chrono::system_clock::now();
    aggregation_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    seal::Plaintext plain;
    std::vector<double> computed(slots_count), computed_case(slots_count);
    decryptor.decrypt(result, plain);
    seal::pack_decode(computed, plain, ckks_encoder);

    decryptor.decrypt(result_case, plain);
    seal::pack_decode(computed_case, plain, ckks_encoder);

    double err1 = 0., err2 = 0.;
    
    for (size_t i = 0; i < slots_count; ++i)
    {
        err1 += std::abs(computed[i] - plain_filter_res[i]);
        err2 += std::abs(computed_case[i] - plain_filter_case_res[i]);
    }

    printf("Repack average error = %f ~ 2^%.1f\n", err1 / slots_count, std::log2(err1 / slots_count));
    printf("Repack average error = %f ~ 2^%.1f\n", err2 / slots_count, std::log2(err2 / slots_count));


    // Filter result * data
    seal::Ciphertext revenue_cipher;
    double qd = parms.coeff_modulus()[result.coeff_modulus_size() - 1].value();
    seal::pack_encode(revenue, qd, plain, ckks_encoder);
    encryptor.encrypt_symmetric(plain, revenue_cipher);

    std::cout << "Aggregating price and discount .." << std::endl;
    start = std::chrono::system_clock::now();
    seal::multiply_and_relinearize(result, revenue_cipher, result, evaluator, relin_keys);
    seal::multiply_and_relinearize(result_case, revenue_cipher, result_case, evaluator, relin_keys);
    evaluator.rescale_to_next_inplace(result);
    evaluator.rescale_to_next_inplace(result_case);
    std::cout << "Remian modulus: " << result.coeff_modulus_size() << std::endl;
    int logrow = log2(num);
    
    seal::Ciphertext temp;
    size_t step;
    for (size_t i = 0; i < logrow; i++)
    {
        temp = result;
        step = 1 << (logrow - i - 1);
        evaluator.rotate_vector_inplace(temp, step, galois_keys);
        evaluator.add_inplace(result, temp);

        temp = result_case;
        evaluator.rotate_vector_inplace(temp, step, galois_keys);
        evaluator.add_inplace(result_case, temp);
    }
    end = std::chrono::system_clock::now();
    aggregation_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::vector<double> agg_result(slots_count), agg_case_result(slots_count);
    decryptor.decrypt(result, plain);
    seal::pack_decode(agg_result, plain, ckks_encoder);

    decryptor.decrypt(result_case, plain);
    seal::pack_decode(agg_case_result, plain, ckks_encoder);

    std::cout << "Query Evaluation Time: " << filtering_time + aggregation_time << " ms" << std::endl;

    std::cout << "Encrypted query result: " << std::endl;
    std::cout << std::setw(12) <<"promo_revenue" << std::endl;
    std::cout << std::setw(12) << agg_case_result[0] / agg_result[0] << std::endl;
    std::cout << "Plain query result: " << std::endl;
    std::cout << std::setw(12) <<"promo_revenue" << std::endl;
    std::cout << std::setw(12) << (plain_agg_case_res + 0.) / plain_agg_res << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    
}



int main()
{
    size_t num = 32768;
    relational_query1(num);
    relational_query6(num);
    relational_query12(num);
    relational_query14(num);
}



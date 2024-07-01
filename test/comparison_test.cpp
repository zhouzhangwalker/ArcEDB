#include "ARCEDB/comparison/rgsw_ciphertext.h"
#include "ARCEDB/comparison/comparable.h"
#include "ARCEDB/conversion/packlwes.h"
#include "ARCEDB/comparison/batch_bootstrap.h"
#include "ARCEDB/utils/serialize.h"
#include <random>
#include <chrono>
#include <unistd.h>

using namespace arcedb;
using namespace seal;

/**
    @param[in] num_test The test number.
    @param[in] p The precision
*/
void comparison_test_relational(size_t num_test, size_t p)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = Lvl1;
    TFHESecretKey sk;
    TFHEEvalKey ek;
    using bkP = Lvl01;
    using iksP = Lvl10;
    std::uniform_int_distribution<typename P::T> message(0, (1ULL << p) - 1);
    std::uniform_int_distribution<typename P::T> message64(0, 0ULL - 1);
    ek.emplacebkfft<bkP>(sk);
    ek.emplaceiksk<iksP>(sk);
    uint32_t gres, p1, p2, dgres;
    
    ComparableLvl1 cipher1;
    ComparbleRGSWLvl1 cipher2;
    std::chrono::system_clock::time_point start, end;
    // Relational Comparison
    double relational_comparison_time = 0;
    double relational_error_time = 0;
    std::cout << "-------------------------------------------------------------" << std::endl;
    std::cout << "Test " << p << "-bits relational comparison operator ArbHCMP:" << std::endl;
    for (size_t i = 0; i < num_test; i++)
    {
        if (p == 64)
        {
            p1 = message64(engine);
            p2 = message64(engine);
        }
        else
        {
            p1 = message(engine);
            p2 = message(engine);
        }
        
        if (p1 < p2) gres = 1;
        else gres = 0;
        exponent_encrypt<P>(p1, p, cipher1, sk);
        exponent_encrypt_rgsw<P>(p2, p, cipher2, sk, true);
        TLWELvl1 cres;
        start = std::chrono::system_clock::now();
        // greater_than_tfhepp(cipher1, cipher2, cipher1.size(), cres, ek, sk);
        less_than_tfhepp(cipher1, cipher2, cipher1.size(), cres, ek, sk);
        end = std::chrono::system_clock::now();
        dgres = TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.get<Lvl1>());
        if (dgres != gres) relational_error_time += 1;
        if (p <= 10) relational_comparison_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        else relational_comparison_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
    }
    std::cout << "Error Time" << " : " << relational_error_time << std::endl;
    if (p <= 10) std::cout << "Comparison Time" << " : " << relational_comparison_time / num_test / 1000 << "ms" << std::endl;
    else std::cout << "Comparison Time" << " : " << relational_comparison_time / num_test << "ms" << std::endl;
    std::cout << "-------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;




}

/**
    @param[in] num_test The test number.
    @param[in] p The precision
*/
void comparison_test_equality(size_t num_test, size_t p)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = Lvl1;
    TFHESecretKey sk;
    TFHEEvalKey ek;
    using bkP = Lvl01;
    using iksP = Lvl10;
    std::uniform_int_distribution<typename P::T> message(0, (1ULL << p - 1));;
    ek.emplacebkfft<bkP>(sk);
    ek.emplaceiksk<iksP>(sk);
    uint32_t gres, p1, p2, dgres;
    
    ComparableLvl1 cipher1;
    ComparbleRGSWLvl1 cipher2;
    std::chrono::system_clock::time_point start, end;
    // Equality Comparison
    double equality_comparison_time = 0;
    double equality_error_time = 0;
    std::cout << "Test " << p << "-bits equality comparison:" << std::endl;
    for (size_t i = 0; i < num_test; i++)
    {
        p1 = message(engine);
        if (i % 2 == 0) p2 = p1;
        else p2 = message(engine);
        if (p1 == p2) gres = 1;
        else gres = 0;
        exponent_encrypt<P>(p1, p, cipher1, sk);
        exponent_encrypt_rgsw<P>(p2, p, cipher2, sk, true);
        TLWELvl1 cres;
        start = std::chrono::system_clock::now();
        equality_tfhepp(cipher1, cipher2, cipher1.size(), cres, ek, sk);
        end = std::chrono::system_clock::now();
        dgres = TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.get<Lvl1>());
        if (dgres != gres) equality_error_time += 1;
        if (p <= 10) equality_comparison_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        else equality_comparison_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
    }
    std::cout << "Error Time" << " : " << equality_error_time << std::endl;
    if (p <= 10) std::cout << "Comparison Time" << " : " << equality_comparison_time / num_test / 1000 << "ms" << std::endl;
    else std::cout << "Comparison Time" << " : " << equality_comparison_time / num_test << "ms" << std::endl;
    
}

int main()
{
    size_t test_num = 100;
    comparison_test_relational(test_num,4);
    comparison_test_relational(test_num,8);
    comparison_test_relational(test_num,16);
    comparison_test_relational(test_num,32);
    comparison_test_relational(test_num,64);
    comparison_test_equality(test_num,4);
    comparison_test_equality(test_num,8);
    comparison_test_equality(test_num,16);
    comparison_test_equality(test_num,32);
    comparison_test_equality(test_num,64);
    
    
    
    
}
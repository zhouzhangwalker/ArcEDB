#pragma once

#include "seal/seal.h"
#include "rgsw_ciphertext.h"
#include "ARCEDB/utils/types.h"
#include "ARCEDB/conversion/lwe_ciphertext.h"

namespace arcedb
{
    enum class compare_type : std::uint8_t
    {
        // equality
        eq = 0x0,

        // greater than
        gt = 0x1,

        // less than
        lt = 0x2,

        // greater than or equal
        ge = 0x3,

        // less than or equal
        le = 0x4, 

        // inequality
        ie = 0x5
    };

    inline size_t bits_count(uint64_t v)
    {
        if (v == 0)
        {
            return 1;
        }
        else
        {
            return log2(v) + 1;
        }
        
    }

    using ComparableLvl1 = std::vector<TRLWELvl1>;
    using ComparbleRGSWLvl1 = std::vector<TRGSWLvl1>;

    using ComparableCipher = std::vector<seal::Ciphertext>;
    using ComparableRGSWCipher = std::vector<RGSWCiphertext>;
    void comparable_encrypt_bfv(uint64_t v, seal::Ciphertext &destination, uint64_t u,
                            seal::SEALContext &context, seal::Encryptor &encryptor);
    
    void comparable_encrypt_bfv(uint64_t v, ComparableCipher &destination, uint64_t u,
                            seal::SEALContext &context, seal::Encryptor &encryptor);
            
    void comparable_encrypt_bfv(std::vector<uint64_t> v, ComparableCipher &destination, uint64_t u,
                            seal::SEALContext &context, seal::Encryptor &encryptor);
    
    void comparable_encrypt_bfv(uint64_t v, uint64_t precision, ComparableCipher &destination, uint64_t u, 
                                seal::SEALContext &context, seal::Encryptor &encryptor);

    void comparable_rgsw_encrypt(uint64_t v, RGSWCiphertext &destination, bool is_negative, 
                            seal::Encryptor &encryptor, seal::SecretKey &secret_key, seal::SEALContext &context);

    void comparable_rgsw_encrypt(uint64_t v, ComparableRGSWCipher &destination, bool is_negative,
                                seal::Encryptor &encryptor, seal::SecretKey &secret_key, seal::SEALContext &context);
    
    void comparable_rgsw_encrypt(std::vector<uint64_t> v, ComparableRGSWCipher &destination, bool is_negative,
                                seal::Encryptor &encryptor, seal::SecretKey &secret_key, seal::SEALContext &context);

    void comparable_rgsw_encrypt(uint64_t v, uint64_t precision, ComparableRGSWCipher &destination, bool is_negative,
                                seal::Encryptor &encryptor, seal::SecretKey &secret_key, seal::SEALContext &context);

    void comparable_decrypt(seal::Ciphertext &cipher, uint64_t &destination,
                            seal::SEALContext &context, seal::Decryptor &decryptor);
    
    void compare_greater_than(seal::Ciphertext &cipher1, RGSWCiphertext &cipher2,
                                    seal::Ciphertext &destination, seal::SEALContext &context, seal::Evaluator &evaluator);
    
    void compare_greater_than(ComparableCipher &cipher1, ComparableRGSWCipher &cipher2, size_t compute_size, seal::Ciphertext &destination, 
                                seal::GaloisKeys &galois_keys, seal::RelinKeys relin_keys, seal::SEALContext &context, seal::Evaluator &evaluator);
    
    
    
    // The ciphertext is +- q/16
    void cmux_lut(TLWELvl1 &sig, TLWELvl1 &cipher1, TLWELvl1 &cipher2, TLWELvl1 &res, uint32_t scale_bits, const TFHEEvalKey &ek);

    void cmux_lut(TLWELvl2 &sig, TLWELvl2 &cipher1, TLWELvl2 &cipher2, TLWELvl1 &res, uint32_t scale_bits, const TFHEEvalKey &ek);

    // mux = three input gate
    template <class P>
    TFHEpp::Polynomial<P> muxpolygen(uint32_t scale_bits)
    {
        TFHEpp::Polynomial<P> poly;
        uint32_t basic_value = 1;
        basic_value = basic_value << scale_bits;
        uint32_t padding_bits = P :: nbit - 3;
        std::array<int, 8> value_list ={0, 0, 1, 1, 1, 0, 1, 0};
        for (uint32_t i = 0; i < P::n; i++)
        {
            uint32_t index = i >> padding_bits;
            poly[i] = value_list[index] ? basic_value : -basic_value;
        }
        return poly;
    }

    template <class P>
    TFHEpp::TLWE<P> tlweSymInt32Encrypt(const typename P::T p, const double α, const double scale, const TFHEpp::Key<P> &key)
    {
        using namespace TFHEpp;
        std::uniform_int_distribution<typename P::T> Torusdist(0, std::numeric_limits<typename P::T>::max());
        TLWE<P> res = {};
        res[P::k * P::n] =
            ModularGaussian<P>(static_cast<typename P::T>(p * scale), α);
        // res[P::k * P::n] = static_cast<typename P::T>(p * scale);
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++) {
                res[k * P::n + i] = Torusdist(generator);
                // res[k * P::n + i] = 1;
                res[P::k * P::n] += res[k * P::n + i] * key[k * P::n + i];
            }
        return res;
    }

    template <class P>
    typename P::T tlweSymInt32Decrypt(const TFHEpp::TLWE<P> &c, const double scale, const TFHEpp::Key<P> &key)
    {
        typename P::T phase = c[P::k * P::n];
        typename P::T plain_modulus = (1ULL << (std::numeric_limits<typename P::T>::digits -1)) / scale;
        plain_modulus *= 2;
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++)
                phase -= c[k * P::n + i] * key[k * P::n + i];
        typename P::T res = 
        static_cast<typename P::T>(std::round(phase / scale)) % plain_modulus;
        return res;
    }

    template <class P>
    TFHEpp::TRLWE<P> trlweSymInt32Encrypt(const std::array<typename P::T, P::n> &p, const double α, const double scale, const TFHEpp::Key<P> &key)
    {
        using namespace TFHEpp;
        TRLWE<P> c = trlweSymEncryptZero<P>(α, key);
        for (int i = 0; i < P::n; i++)
            c[P::k][i] += static_cast<typename P::T>(scale* p[i]);
        return c;
    }

    template <class P>
    TFHEpp::Polynomial<P> trlweSymInt32Decrypt(const TFHEpp::TRLWE<P> &c, double scale, const TFHEpp::Key<P> &key)
    {
        using namespace TFHEpp;
        Polynomial<P> phase = c[P::k];
        typename P::T plain_modulus = (1ULL << (std::numeric_limits<typename P::T>::digits -1)) / scale;
        plain_modulus *= 2;
        for (int k = 0; k < P::k; k++) {
            Polynomial<P> mulres;
            std::array<typename P::T, P::n> partkey;
            for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
            PolyMul<P>(mulres, c[k], partkey);
            for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];
        }

        Polynomial<P> p;
        for (int i = 0; i < P::n; i++)
            p[i] = static_cast<typename P::T>(std::round(phase[i] / scale)) % plain_modulus;
        return p;
    }

    template <class P>
    void exponent_encrypt(uint64_t value, TFHEpp::TRLWE<P> &cipher, TFHESecretKey &sk)
    {
        if (value >= P :: n)
        {
            throw std::invalid_argument("Out of range.");
        }
        std::array<typename P::T, P::n> plain {};
        plain[value] = P::μ ;

        cipher = TFHEpp::trlweSymEncrypt<P>(plain, P::α, sk.key.get<P>());
    }

    template <class P>
    void exponent_encrypt(uint64_t value, uint64_t scale,TFHEpp::TRLWE<P> &cipher, TFHESecretKey &sk)
    {
        if (value >= P :: n)
        {
            throw std::invalid_argument("Out of range.");
        }
        std::array<typename P::T, P::n> plain {};
        plain[value] = scale ;

        cipher = TFHEpp::trlweSymEncrypt<P>(plain, P::α, sk.key.get<P>());
    }

    template<class P>
    void exponent_encrypt(uint64_t value, uint32_t precision, std::vector<TFHEpp::TRLWE<P>> &ciphers, TFHESecretKey &sk)
    {
        if (precision !=64 && value >= (1ULL << precision))
        {
            throw std::invalid_argument("Out of range.");
        }
        uint8_t cipher_size = static_cast<uint8_t>(ceil(precision/log2(P::n)));
        ciphers.resize(cipher_size);

        
        for (size_t i = 0; i < cipher_size; i++)
        {
            uint64_t base_value = value % P::n;
            exponent_encrypt<P>(base_value, ciphers[i], sk);
            value = (value - base_value) / P::n;
        }
    }

    template <class P>
    void exponent_encrypt(uint64_t value, std::vector<TFHEpp::TRLWE<P>> &ciphers, TFHESecretKey &sk)
    {
        uint8_t cipher_size = static_cast<uint8_t>(ceil(bits_count(value) / log2(P::n)));
        ciphers.resize(cipher_size);
        
        for (size_t i = 0; i < cipher_size; i++)
        {
            uint64_t base_value = value % P::n;
            exponent_encrypt<P>(base_value, ciphers[i], sk);
            value = (value - base_value) / P::n;
        }
    }

    template <class P>
    void exponent_encrypt_rgsw(uint64_t value, TFHEpp::TRGSWFFT<P> &cipher, TFHESecretKey &sk, bool is_negative)
    {
        if (value >= P :: n)
        {
            throw std::invalid_argument("Out of range.");
        }
        TFHEpp::Polynomial<P> plainpoly = {};
        constexpr std::array<typename P::T, P::l> h = TFHEpp::hgen<P>();
        TFHEpp::TRGSW<P> trgsw;
        for (TFHEpp::TRLWE<P> &trlwe : trgsw) trlwe = TFHEpp::trlweSymEncryptZero<P>(P::α, sk.key.get<P>());
        if (value != 0)
        {
            for (int i = 0; i < P::l; i++) 
            {
                for (int k = 0; k < P::k + 1; k++) 
                {
                    if (is_negative)
                    {
                        trgsw[i + k * P::l][k][P::n - value] -= h[i];
                    }
                    else
                    {
                        trgsw[i + k * P::l][k][value] += h[i];
                    }
                }
            }
        }
        else
        {
            for (int i = 0; i < P::l; i++) 
            {
                for (int k = 0; k < P::k + 1; k++) 
                {
                    if (is_negative)
                    {
                        trgsw[i + k * P::l][k][value] += h[i];
                    }
                    else
                    {
                        trgsw[i + k * P::l][k][value] += h[i];
                    }
                }
            }
        }
        
        cipher = TFHEpp::ApplyFFT2trgsw<P>(trgsw);
    }

    template <class P>
    void exponent_encrypt_rgsw(uint64_t value, std::vector<TFHEpp::TRGSWFFT<P>> &ciphers, TFHESecretKey &sk, bool is_negative)
    {
        uint8_t cipher_size = static_cast<uint8_t>(ceil(bits_count(value) / log2(P::n)));
        ciphers.resize(cipher_size);
        for (size_t i = 0; i < cipher_size; i++)
        {
            uint64_t base_value = value % P::n;
            exponent_encrypt_rgsw<P>(base_value, ciphers[i], sk, is_negative);
            value = (value - base_value) / P::n;
        }
    }

    template <class P>
    void exponent_encrypt_rgsw(uint64_t value, uint32_t precision,
                                std::vector<TFHEpp::TRGSWFFT<P>> &ciphers, TFHESecretKey &sk, bool is_negative)
    {
        
        if (precision !=64 && value >= (1ULL << precision))
        {
            throw std::invalid_argument("Out of range.");
        }
        uint8_t cipher_size = static_cast<uint8_t>(ceil(precision/log2(P::n)));
        ciphers.resize(cipher_size);
        ciphers.resize(cipher_size);
        for (size_t i = 0; i < cipher_size; i++)
        {
            uint64_t base_value = value % P::n;
            exponent_encrypt_rgsw<P>(base_value, ciphers[i], sk, is_negative);
            value = (value - base_value) / P::n;
        }
    }


    void greater_than_tfhepp(TRLWELvl1 &cipher1, TRGSWLvl1 &cipher2, TLWELvl1 &res, TFHESecretKey &sk);

    void greater_than_tfhepp(TRLWELvl2 &cipher1, TRGSWLvl2 &cipher2, TLWELvl2 &res, TFHESecretKey &sk, TFHEEvalKey &ek);
    
    void greater_than_mux(TLWELvl1 &low_res, TLWELvl1 &high_res, TLWELvl1 &equal_res, TLWELvl1 &res, TFHEEvalKey &ek);

    void greater_than_tfhepp(std::vector<TRLWELvl1> &ciphers1, std::vector<TRGSWLvl1> &ciphers2, size_t cipher_size, TLWELvl1 &res, 
                            TFHEEvalKey &ek, TFHESecretKey &sk);
    
    void greater_than_tfhepp(std::vector<TRLWELvl1> &ciphers1, std::vector<TRGSWLvl1> &ciphers2, size_t cipher_size, TLWELvl1 &res, uint32_t scale_bits,
                        TFHEEvalKey &ek, TFHESecretKey &sk);
    
    void greater_than_tfhepp(std::vector<TRLWELvl1> &ciphers1, std::vector<TRGSWLvl1> &ciphers2, size_t cipher_size, TLWELvl2 &res, uint32_t scale_bits,
                            TFHEEvalKey &ek, TFHESecretKey &sk);
    
    template <class P>
    constexpr TFHEpp::Polynomial<P> μ_polygen(typename P::T μ)
    {
        TFHEpp::Polynomial<P> poly;
        for (typename P::T &p : poly) p = -μ;
        return poly;
    }

    void equality_tfhepp(TRLWELvl1 &cipher1, TRGSWLvl1 &cipher2, TLWELvl1 &res, TFHESecretKey &sk);

    void equality_tfhepp(std::vector<TRLWELvl1> &ciphers1, std::vector<TRGSWLvl1> &ciphers2, size_t cipher_size, TLWELvl1 &res, 
                            TFHEEvalKey &ek, TFHESecretKey &sk);
    
    void less_than_tfhepp(TRLWELvl1 &cipher1, TRGSWLvl1 &cipher2, TLWELvl1 &res, TFHESecretKey &sk);
    
    void less_than_tfhepp(std::vector<TRLWELvl1> &ciphers1, std::vector<TRGSWLvl1> &ciphers2, size_t cipher_size, TLWELvl1 &res, 
                            TFHEEvalKey &ek, TFHESecretKey &sk);

    void lift_and_and(TLWELvl1 &cipher1, TLWELvl1 &cipher2, TLWELvl1 &res, uint32_t scale_bits, TFHEpp::EvalKey &ek, TFHEpp::SecretKey &sk);
    
    void lift_and_and(TLWELvl1 &cipher1, TLWELvl1 &cipher2, TLWELvl2 &res, uint32_t scale_bits, TFHEpp::EvalKey &ek, TFHEpp::SecretKey &sk);
} // namespace arcedb
 
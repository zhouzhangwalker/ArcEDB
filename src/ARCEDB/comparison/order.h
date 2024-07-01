#pragma once

#include "seal/seal.h"
#include "ARCEDB/conversion/lwe_ciphertext.h"
#include "ARCEDB/comparison/comparable.h"

namespace arcedb
{

    template <class P>
    void syn_column(std::vector<TFHEpp::TRLWE<P>> &ciphers, std::vector<TFHEpp::TRGSWFFT<P>> &indexs, std::vector<TFHEpp::TLWE<P>> &res,
                    TFHEpp::SecretKey &sk)
    {
        TFHEpp::TRLWE<P> sum, temp;
        if (ciphers.size() != indexs.size())
        {
            std::cout << "Location: " << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")\n";
            throw std::invalid_argument("Size mismatch.");
        }
        for (size_t i = 0; i < ciphers.size(); i++)
        {
            TFHEpp::trgswfftExternalProduct<P>(temp, ciphers[i], indexs[i]);
            // TFHEpp::Polynomial<P> debug_res = trlweSymInt32Decrypt<P>(temp, std::pow(2.0, 24), sk.key.get<P>());
            if (i == 0)
            {
                sum = temp;
            }
            else
            {
                for (size_t t = 0; t <= P::k; t++)
                {
                    for (size_t j = 0; j < P::n; j++)
                    {
                        sum[t][j] = sum[t][j] + temp[t][j];
                    }
                }      
            }
        }
        res.resize(ciphers.size());
        for (size_t i = 0; i < ciphers.size(); i++)
        {
            TFHEpp::SampleExtractIndex<P>(res[i], sum, i);
        }
    }

    template <class P>
    void generate_mux_cipher(TFHEpp::TRLWE<P> &coeff_cipher, size_t num, std::vector<TFHEpp::TRLWE<P>> &ciphers)
    {
        ciphers.resize(num);
        ciphers[0] = coeff_cipher;
        for (size_t i = 1; i < num; i++)
        {
            for (size_t t = 0; t <= P::k; t++)
            {
                for (size_t j = 0; j < P::n; j++)
                {
                    if (j == 0) 
                    {
                        ciphers[i][t][j] = -ciphers[i-1][t][P::n - 1];
                    }
                    else
                    {
                        ciphers[i][t][j] = ciphers[i-1][t][j-1];
                    }
                }
            }            
        }
        
    }

    void bit_decompose(TLWELvl1 &ciphers, std::vector<TLWELvl1> &result, uint32_t scale_bits, TFHEEvalKey &ek, TFHESecretKey &sk);

    void bit_decompose(TLWELvl2 &ciphers, std::vector<TLWELvl2> &result, uint32_t scale_bits, TFHEEvalKey &ek, TFHESecretKey &sk);

    void order_by(std::vector<ComparableLvl1> ciphers, std::vector<ComparbleRGSWLvl1> rgsw_ciphers, uint32_t scale_bits,
                std::vector<TLWELvl1> &result, TFHEEvalKey &ek, TFHESecretKey &sk);
    
    void order_by(std::vector<ComparableLvl1> ciphers, std::vector<ComparbleRGSWLvl1> rgsw_ciphers, uint32_t scale_bits,
                std::vector<TLWELvl2> &result, TFHEEvalKey &ek, TFHESecretKey &sk);
    
    void combine_exponent(std::vector<TLWELvl1> &bit_deomps, TRLWELvl2 &exponent_cipher, std::vector<TRLWELvl2> &mux_ciphers,
                            uint32_t scale_bits, TFHEEvalKey &ek, TFHESecretKey &sk);

    void combine_exponent(std::vector<TLWELvl2> &bit_deomps, TRLWELvl2 &exponent_cipher, std::vector<TRLWELvl2> &mux_ciphers,
                        uint32_t scale_bits, TFHEEvalKey &ek, TFHESecretKey &sk);

    void combine_exponent(std::vector<TLWELvl1> &bit_deomps, std::vector<TRLWELvl2> &exponent_cipher, std::vector<std::vector<TRLWELvl2>> &mux_ciphers,
                            uint32_t scale_bits, TFHEEvalKey &ek, TFHESecretKey &sk);
    
    void combine_exponent(std::vector<TLWELvl2> &bit_deomps, std::vector<TRLWELvl2> &exponent_cipher, std::vector<std::vector<TRLWELvl2>> &mux_ciphers,
                            uint32_t scale_bits, TFHEEvalKey &ek, TFHESecretKey &sk);
} // namespace arcedb

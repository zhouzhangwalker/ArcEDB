#include "order.h"
#include "ARCEDB/conversion/packlwes.h"
#include "ARCEDB/comparison/comparable.h"
#include <chrono>

namespace arcedb
{
    

    void combine_exponent(std::vector<TLWELvl1> &bit_deomps, TRLWELvl2 &exponent_cipher, std::vector<TRLWELvl2> &mux_ciphers,
                            uint32_t scale_bits, TFHEEvalKey &ek, TFHESecretKey &sk)
    {
        using iksP = Lvl10;
        using bkP = Lvl02;
        int logN = bit_deomps.size();
        std::vector<TRGSWLvl2> bootedTGSWs(bit_deomps.size());
        std::chrono::system_clock::time_point start, end;
        // start = std::chrono::system_clock::now();
        for (size_t i = 0; i < logN; i++)
        {
            TFHEpp::CircuitBootstrappingFFT<iksP, bkP, Lvl22>(bootedTGSWs[i], bit_deomps[i], ek);
        }
        // end = std::chrono::system_clock::now();
        // double conversion_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        // std::cout << "Conversion Time: " << conversion_time << " ms" <<std::endl;
        TFHEpp::Polynomial<Lvl2> plainpoly = {0ULL};
        uint32_t mux_time = 0;
        for (size_t i = 0; i < logN; i++)
        {
            for (size_t j = 0; j < (1 << (logN - i - 1)); j++)
            {
                TRLWELvl2 temp;
                size_t index1 = (1 << (i + 1)) * j;
                size_t index2 = (1 << (i + 1)) * j + (1 << i);
                TFHEpp::CMUXFFT<Lvl2>(temp, bootedTGSWs[i], mux_ciphers[index2], mux_ciphers[index1]);
                // mux_time += 1;
                mux_ciphers[index1] = temp;
            }
            
        }
        exponent_cipher = mux_ciphers[0];
        // std::cout << "mux time: " << mux_time << std::endl;
    }

    void combine_exponent(std::vector<TLWELvl2> &bit_deomps, TRLWELvl2 &exponent_cipher, std::vector<TRLWELvl2> &mux_ciphers,
                            uint32_t scale_bits, TFHEEvalKey &ek, TFHESecretKey &sk)
    {
        using iksP = Lvl20;
        using bkP = Lvl02;
        int logN = bit_deomps.size();
        std::vector<TRGSWLvl2> bootedTGSWs(bit_deomps.size());
        std::chrono::system_clock::time_point start, end;
        // start = std::chrono::system_clock::now();
        for (size_t i = 0; i < logN; i++)
        {
            TFHEpp::CircuitBootstrappingFFT<iksP, bkP, Lvl22>(bootedTGSWs[i], bit_deomps[i], ek);
        }
        // end = std::chrono::system_clock::now();
        // double conversion_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        // std::cout << "Conversion Time: " << conversion_time << " ms" <<std::endl;
        TFHEpp::Polynomial<Lvl2> plainpoly = {0ULL};
        uint32_t mux_time = 0;
        for (size_t i = 0; i < logN; i++)
        {
            for (size_t j = 0; j < (1 << (logN - i - 1)); j++)
            {
                TRLWELvl2 temp;
                size_t index1 = (1 << (i + 1)) * j;
                size_t index2 = (1 << (i + 1)) * j + (1 << i);
                TFHEpp::CMUXFFT<Lvl2>(temp, bootedTGSWs[i], mux_ciphers[index2], mux_ciphers[index1]);
                // mux_time += 1;
                mux_ciphers[index1] = temp;
            }
            
        }
        exponent_cipher = mux_ciphers[0];
        // std::cout << "mux time: " << mux_time << std::endl;
    }

    void combine_exponent(std::vector<TLWELvl1> &bit_deomps, std::vector<TRLWELvl2> &exponent_cipher, std::vector<std::vector<TRLWELvl2>> &mux_ciphers,
                            uint32_t scale_bits, TFHEEvalKey &ek, TFHESecretKey &sk)
    {
        using iksP = Lvl10;
        using bkP = Lvl02;
        int logN = bit_deomps.size();
        std::vector<TRGSWLvl2> bootedTGSWs(bit_deomps.size());
        std::chrono::system_clock::time_point start, end;
        // start = std::chrono::system_clock::now();
        for (size_t i = 0; i < logN; i++)
        {
            TFHEpp::CircuitBootstrappingFFT<iksP, bkP, Lvl22>(bootedTGSWs[i], bit_deomps[i], ek);
        }
        // end = std::chrono::system_clock::now();
        // double conversion_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        // std::cout << "Conversion Time: " << conversion_time << " ms" <<std::endl;
        TFHEpp::Polynomial<Lvl2> plainpoly = {0ULL};
        uint32_t mux_time = 0;
        for (size_t t = 0; t < exponent_cipher.size(); t++)
        {
            for (size_t i = 0; i < logN; i++)
            {
                for (size_t j = 0; j < (1 << (logN - i - 1)); j++)
                {
                    TRLWELvl2 temp;
                    size_t index1 = (1 << (i + 1)) * j;
                    size_t index2 = (1 << (i + 1)) * j + (1 << i);
                    TFHEpp::CMUXFFT<Lvl2>(temp, bootedTGSWs[i], mux_ciphers[t][index2], mux_ciphers[t][index1]);
                    // mux_time += 1;
                    mux_ciphers[t][index1] = temp;
                }
                
            }
            exponent_cipher[t] = mux_ciphers[t][0];
        }
        
        // std::cout << "mux time: " << mux_time << std::endl;
    }

    void combine_exponent(std::vector<TLWELvl2> &bit_deomps, std::vector<TRLWELvl2> &exponent_cipher, std::vector<std::vector<TRLWELvl2>> &mux_ciphers,
                            uint32_t scale_bits, TFHEEvalKey &ek, TFHESecretKey &sk)
    {
        using iksP = Lvl20;
        using bkP = Lvl02;
        int logN = bit_deomps.size();
        std::vector<TRGSWLvl2> bootedTGSWs(bit_deomps.size());
        std::chrono::system_clock::time_point start, end;
        // start = std::chrono::system_clock::now();
        for (size_t i = 0; i < logN; i++)
        {
            TFHEpp::CircuitBootstrappingFFT<iksP, bkP, Lvl22>(bootedTGSWs[i], bit_deomps[i], ek);
        }
        // end = std::chrono::system_clock::now();
        // double conversion_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        // std::cout << "Conversion Time: " << conversion_time << " ms" <<std::endl;
        TFHEpp::Polynomial<Lvl2> plainpoly = {0ULL};
        uint32_t mux_time = 0;
        for (size_t t = 0; t < exponent_cipher.size(); t++)
        {
            for (size_t i = 0; i < logN; i++)
            {
                for (size_t j = 0; j < (1 << (logN - i - 1)); j++)
                {
                    TRLWELvl2 temp;
                    size_t index1 = (1 << (i + 1)) * j;
                    size_t index2 = (1 << (i + 1)) * j + (1 << i);
                    TFHEpp::CMUXFFT<Lvl2>(temp, bootedTGSWs[i], mux_ciphers[t][index2], mux_ciphers[t][index1]);
                    // mux_time += 1;
                    mux_ciphers[t][index1] = temp;
                }
                
            }
            exponent_cipher[t] = mux_ciphers[t][0];
        }
        
        // std::cout << "mux time: " << mux_time << std::endl;
    }


    void bit_decompose(TLWELvl1 &ciphers, std::vector<TLWELvl1> &result, uint32_t scale_bits, TFHEEvalKey &ek, TFHESecretKey &sk)
    {
        size_t num_bits = std::numeric_limits<uint32_t>::digits - scale_bits;
        result.resize(num_bits);
        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            result[0][i] = ciphers[i] << (num_bits - 1);
        }
        TLWELvl1 temp = ciphers;
        for (size_t i = 1; i < num_bits; i++)
        {
            TLWELvl1 tlwelvl1 = result[i - 1];
            tlwelvl1[Lvl1::n] += (Lvl1::μ << 1);
            Lvl1::T c = (1ULL << (scale_bits -2 + i));
            TLWELvl0 tlwelvl0;
            TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlwelvl1, *ek.iksklvl10);
            TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(result[i], tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(c));
            result[i][Lvl1::k * Lvl1::n] += c;
            for (size_t t = 0; t <= Lvl1::n; t++)
            {
                temp[t] -= result[i][t];
            }
            uint32_t debug_res = tlweSymInt32Decrypt<Lvl1>(temp, std::pow(2., scale_bits), sk.key.get<Lvl1>());
            for (size_t t = 0; t <= Lvl1 :: n; t++)
            {
                result[i][t] = temp[t] << (num_bits - i - 1);
            }
            
        }
        for (size_t i = 0; i < num_bits; i++)
        {
            result[i][Lvl1::n] -= (Lvl1::μ << 1);
        }
        
        
    }

    void order_by(std::vector<ComparableLvl1> ciphers, std::vector<ComparbleRGSWLvl1> rgsw_ciphers, uint32_t scale_bits,
                std::vector<TLWELvl1> &result, TFHEEvalKey &ek, TFHESecretKey &sk)
    {
        size_t num = ciphers.size();
        if ((num & (num - 1)))
        {
            std::cout << "Location: " << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")\n";
            throw std::invalid_argument("The size is not power of two.");
        }

        std::vector<std::vector<TLWELvl1>> compare_results(num);
        for (size_t i = 0; i < num; i++)
        {
            compare_results[i].resize(num);
            for (size_t j = i + 1; j < num; j++)
            {
                greater_than_tfhepp(ciphers[i], rgsw_ciphers[j], ciphers[i].size(),compare_results[i][j], scale_bits, ek, sk);
            }
            for (size_t j = 0; j < i; j++)
            {
                for (size_t t = 0; t <= Lvl1::n; t++)
                {
                    compare_results[i][j][t] = -compare_results[j][i][t];
                }
                compare_results[i][j][Lvl1::n] += (1ULL << scale_bits);
            }
        }

        // DEBUg
        std::vector<std::vector<uint32_t>> debug_res(num);
        for (size_t i = 0; i < num; i++)
        {
            debug_res[i].resize(num);
            for (size_t j = 0; j < num; j++)
            {
                debug_res[i][j] = tlweSymInt32Decrypt<Lvl1>(compare_results[i][j], std::pow(2, scale_bits), sk.key.get<Lvl1>());
            }
        }
        

        std::vector<TLWELvl1> index_ciphers(num);
        for (size_t i = 0; i < num; i++)
        {
            index_ciphers[i] = compare_results[i][0];
            for (size_t j = 1; j < num; j++)
            {
                for (size_t t = 0; t <= Lvl1::n; t++)
                {
                    index_ciphers[i][t] += compare_results[i][j][t];
                }
                
            }
            result[i] = index_ciphers[i];
            
        }
    }

    void order_by(std::vector<ComparableLvl1> ciphers, std::vector<ComparbleRGSWLvl1> rgsw_ciphers, uint32_t scale_bits,
                std::vector<TLWELvl2> &result, TFHEEvalKey &ek, TFHESecretKey &sk)
    {
        size_t num = ciphers.size();
        if ((num & (num - 1)))
        {
            std::cout << "Location: " << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")\n";
            throw std::invalid_argument("The size is not power of two.");
        }

        std::vector<std::vector<TLWELvl2>> compare_results(num);
        for (size_t i = 0; i < num; i++)
        {
            compare_results[i].resize(num);
            for (size_t j = i + 1; j < num; j++)
            {
                greater_than_tfhepp(ciphers[i], rgsw_ciphers[j], ciphers[i].size(),compare_results[i][j], scale_bits, ek, sk);
            }
            for (size_t j = 0; j < i; j++)
            {
                for (size_t t = 0; t <= Lvl2::n; t++)
                {
                    compare_results[i][j][t] = -compare_results[j][i][t];
                }
                compare_results[i][j][Lvl2::n] += (1ULL << scale_bits);
            }
        }

        // DEBUg
        std::vector<std::vector<uint32_t>> debug_res(num);
        for (size_t i = 0; i < num; i++)
        {
            debug_res[i].resize(num);
            for (size_t j = 0; j < num; j++)
            {
                debug_res[i][j] = tlweSymInt32Decrypt<Lvl2>(compare_results[i][j], std::pow(2, scale_bits), sk.key.get<Lvl2>());
            }
        }
        

        std::vector<TLWELvl2> index_ciphers(num);
        for (size_t i = 0; i < num; i++)
        {
            index_ciphers[i] = compare_results[i][0];
            for (size_t j = 1; j < num; j++)
            {
                for (size_t t = 0; t <= Lvl2::n; t++)
                {
                    index_ciphers[i][t] += compare_results[i][j][t];
                }
                
            }
            result[i] = index_ciphers[i];
            
        } 
    }

    void bit_decompose(TLWELvl2 &ciphers, std::vector<TLWELvl2> &result, uint32_t scale_bits, TFHEEvalKey &ek, TFHESecretKey &sk)
    {
        size_t num_bits = std::numeric_limits<uint64_t>::digits - scale_bits;
        result.resize(num_bits);
        using P = Lvl2;
        for (size_t i = 0; i <= P :: n; i++)
        {
            result[0][i] = ciphers[i] << (num_bits - 1);
        }
        TLWELvl2 temp = ciphers;
        for (size_t i = 1; i < num_bits; i++)
        {
            TLWELvl2 tlwelvl2 = result[i - 1];
            tlwelvl2[P::n] += (P::μ << 1);
            P::T c = (1ULL << (scale_bits -2 + i));
            TLWELvl0 tlwelvl0;
            TFHEpp::IdentityKeySwitch<Lvl20>(tlwelvl0, tlwelvl2, *ek.iksklvl20);
            TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl02>(result[i], tlwelvl0, *ek.bkfftlvl02, μ_polygen<Lvl2>(c));
            result[i][P::k * P::n] += c;
            for (size_t t = 0; t <= P::n; t++)
            {
                temp[t] -= result[i][t];
            }
            uint32_t debug_res = tlweSymInt32Decrypt<P>(temp, std::pow(2., scale_bits), sk.key.get<P>());
            for (size_t t = 0; t <= P :: n; t++)
            {
                result[i][t] = temp[t] << (num_bits - i - 1);
            }
            
        }
        for (size_t i = 0; i < num_bits; i++)
        {
            result[i][P::n] -= (P::μ << 1);
        }
        
        
    }
} // namespace arcedb

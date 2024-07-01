#include "ARCEDB/utils/serialize.h"
#include <fstream>
#include <iostream>

namespace arcedb
{
    void generate_secret_key()
    {
        TFHEpp::SecretKey sk;
        {
            std::ofstream ofs{"./../evalkey/secret.key", std::ios::binary};
            cereal::PortableBinaryOutputArchive ar(ofs);
            sk.serialize(ar);
        }
    }

    void load_secret_key(TFHEpp::SecretKey &sk)
    {
        {
            std::ifstream ifs{"./../evalkey/secret.key", std::ios::binary};
            cereal::PortableBinaryInputArchive ar(ifs);
            sk.serialize(ar);
        }
    }

    void generate_filtering_key()
    {
        TFHEpp::SecretKey sk;
        {
            std::ifstream ifs{"./../evalkey/secret.key", std::ios::binary};
            cereal::PortableBinaryInputArchive ar(ifs);
            sk.serialize(ar);
        }

        TFHEpp::EvalKey ek;
        ek.emplacebkfft<TFHEpp::lvl01param>(sk);
        ek.emplaceiksk<TFHEpp::lvl10param>(sk);
        {
            std::ofstream ofs{"./../evalkey/eval.key", std::ios::binary};
            cereal::PortableBinaryOutputArchive ar(ofs);
            ek.serialize(ar);
        }
    }

    void load_filtering_key(TFHEpp::EvalKey &ek)
    {
        {
            std::ifstream ifs{"./../evalkey/eval.key", std::ios::binary};
            cereal::PortableBinaryInputArchive ar(ifs);
            ek.serialize(ar);
        }
    }

    void generate_time_series_key()
    {
        TFHEpp::SecretKey sk;
        {
            std::ifstream ifs{"./../evalkey/secret.key", std::ios::binary};
            cereal::PortableBinaryInputArchive ar(ifs);
            sk.serialize(ar);
        }

        TFHEpp::EvalKey ek;
        ek.emplacebkfft<TFHEpp::lvl01param>(sk);
        ek.emplacebkfft<TFHEpp::lvl02param>(sk);
        ek.emplaceiksk<TFHEpp::lvl10param>(sk);
        {
            std::ofstream ofs{"./../evalkey/eval_time.key", std::ios::binary};
            cereal::PortableBinaryOutputArchive ar(ofs);
            ek.serialize(ar);
        }
    }

    void load_time_series_key(TFHEpp::EvalKey &ek)
    {
        {
            std::ifstream ifs{"./../evalkey/eval_time.key", std::ios::binary};
            cereal::PortableBinaryInputArchive ar(ifs);
            ek.serialize(ar);
        }
    }

} // namespace HEDB
#pragma once
#include <tfhe++.hpp>


namespace arcedb
{
    void generate_secret_key();

    void load_secret_key(TFHEpp::SecretKey &sk);

    void generate_filtering_key();

    void load_filtering_key(TFHEpp::EvalKey &ek);

    void generate_time_series_key();

    void load_time_series_key(TFHEpp::EvalKey &ek);


} // namespace HEDB

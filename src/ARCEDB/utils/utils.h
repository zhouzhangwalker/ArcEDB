#pragma once

#include "multi_thread.h"
#include "serialize.h"
#include "types.h"
#include "ThreadPool.h"


namespace arcedb
{
    template <typename T>
    static inline T CeilSqrt(T val) 
    {
        return static_cast<T>(std::ceil(std::sqrt(1. * val)));
    }

    template <typename T>
    static inline T CeilDiv(T a, T b) 
    {
        return (a + b - 1) / b;
    }
} // namespace arcedb

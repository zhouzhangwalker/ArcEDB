#pragma once

#include "seal/seal.h"
using namespace seal;

namespace arcedb
{
    // SEAL:: (-a*s+\delta*m+e, a)
    // TFHEPP (a, a*s + \delta * m + e)
    class LWECiphertext{
        
    public:
        using ct_coeff_type = std::uint64_t;

        LWECiphertext(MemoryPoolHandle pool = MemoryManager::GetPool()): data_(std::move(pool))
        {}

        explicit LWECiphertext(const SEALContext &context, Ciphertext &rlwe_cipher, size_t coeff_index, MemoryPoolHandle pool = MemoryManager::GetPool())
            : data_(std::move(pool))
        {
            extract(context, rlwe_cipher, coeff_index);
        }

        LWECiphertext(const LWECiphertext &copy) = default;

        LWECiphertext(LWECiphertext &&source) = default;

        LWECiphertext(const LWECiphertext &copy, MemoryPoolHandle pool) : LWECiphertext(std::move(pool))
        {
            *this = copy;
        }

        void resize(const SEALContext &context, parms_id_type parm_id);

        inline void resize(const SEALContext &context)
        {
            auto parms_id = context.first_parms_id();
            resize(context, parms_id);
        }

        inline void release() noexcept
        {
            parms_id_ = parms_id_zero;
            is_ntt_form_ = false;
            poly_modulus_degree_ = 0;
            coeff_modulus_size_ = 0;
            scale_ = 1.0;
            data_.release();
        }

        void extract(const SEALContext &context, Ciphertext &rlwe_cipher, size_t coeff_index);


        LWECiphertext &operator=(const LWECiphertext &assign);

        LWECiphertext &operator=(LWECiphertext &&assign) = default;


        SEAL_NODISCARD inline ct_coeff_type *dataA()
        {
            return data_.begin() + coeff_modulus_size_;
        }

        SEAL_NODISCARD inline const ct_coeff_type *dataA() const
        {
            return data_.cbegin() + coeff_modulus_size_;
        }

        SEAL_NODISCARD inline ct_coeff_type *dataB()
        {
            return data_.begin();
        }

        SEAL_NODISCARD inline const ct_coeff_type *dataB() const
        {
            return data_.cbegin();
        }

        SEAL_NODISCARD inline ct_coeff_type *data() noexcept
        {
            return data_.begin();
        }

        SEAL_NODISCARD inline const ct_coeff_type *data() const noexcept
        {
            return data_.cbegin();
        }

        SEAL_NODISCARD inline ct_coeff_type &operator[](std::size_t coeff_index)
        {
            return data_.at(coeff_index);
        }

        SEAL_NODISCARD inline const ct_coeff_type &operator[](std::size_t coeff_index) const
        {
            return data_.at(coeff_index);
        }

        SEAL_NODISCARD inline std::size_t poly_modulus_degree() const noexcept
        {
            return poly_modulus_degree_;
        }

        SEAL_NODISCARD inline std::size_t coeff_modulus_size() const noexcept
        {
            return coeff_modulus_size_;
        }

        SEAL_NODISCARD inline bool is_ntt_form() const noexcept
        {
            return is_ntt_form_;
        }

        SEAL_NODISCARD inline bool &is_ntt_form() noexcept
        {
            return is_ntt_form_;
        }

        SEAL_NODISCARD inline parms_id_type &parms_id() noexcept
        {
            return parms_id_;
        }

        SEAL_NODISCARD inline const parms_id_type &parms_id() const noexcept
        {
            return parms_id_;
        }

        SEAL_NODISCARD inline double &scale() noexcept
        {
            return scale_;
        }

        SEAL_NODISCARD inline const double &scale() const noexcept
        {
            return scale_;
        }

        SEAL_NODISCARD inline MemoryPoolHandle pool() const noexcept
        {
            return data_.pool();
        }

        SEAL_NODISCARD inline bool is_set() noexcept
        {
            return poly_modulus_degree_ > 0;
        }

        void resize_internal(size_t poly_modulus_degree, size_t coeff_modulus_size);

        
    
    private:

        parms_id_type parms_id_ = parms_id_zero;

        bool is_ntt_form_ = false;

        std::size_t poly_modulus_degree_ = 0;

        std::size_t coeff_modulus_size_ = 0;

        double scale_ = 1.0;

        DynArray<ct_coeff_type> data_;

    };

    uint64_t lwe_decrypt_q(LWECiphertext &cipher, std::vector<int64_t> &secret_key, seal::Modulus &modulus);

    uint64_t lwe_decrypt(LWECiphertext &cipher, std::vector<int64_t> &secret_key, seal::Modulus &modulus, uint64_t t);

    void sample_extract(Ciphertext &rlwe_cipher, LWECiphertext &lwe_cipher, size_t coeff_index, seal::Modulus &modulus);

    void lwe_add_inplace(LWECiphertext &cipher1, LWECiphertext &cipher2, seal::Modulus &modulus);

} // namespace arcedb

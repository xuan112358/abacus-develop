#include "module_hsolver/kernels/math_kernel_op.h"

#include <iomanip>
#include <iostream>

namespace hsolver
{

template <typename T>
struct line_minimize_with_block_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(T* grad_out,
                    T* hgrad_out,
                    T* psi_out,
                    T* hpsi_out,
                    const int& n_basis,
                    const int& n_basis_max,
                    const int& n_band)
    {
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Real epsilo_0 = 0.0, epsilo_1 = 0.0, epsilo_2 = 0.0;
            Real theta = 0.0, cos_theta = 0.0, sin_theta = 0.0;
            auto A = reinterpret_cast<const Real*>(grad_out + band_idx * n_basis_max);
            Real norm = BlasConnector::dot(2 * n_basis, A, 1, A, 1);
            norm = 1.0 / sqrt(norm);
            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                grad_out[item] *= norm;
                hgrad_out[item] *= norm;
                epsilo_0 += std::real(hpsi_out[item] * std::conj(psi_out[item]));
                epsilo_1 += std::real(grad_out[item] * std::conj(hpsi_out[item]));
                epsilo_2 += std::real(grad_out[item] * std::conj(hgrad_out[item]));
            }
            theta = 0.5 * std::abs(std::atan(2 * epsilo_1 / (epsilo_0 - epsilo_2)));
            cos_theta = std::cos(theta);
            sin_theta = std::sin(theta);
            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                psi_out[item] = psi_out[item] * cos_theta + grad_out[item] * sin_theta;
                hpsi_out[item] = hpsi_out[item] * cos_theta + hgrad_out[item] * sin_theta;
            }
        }
    }
};

template <typename T>
struct calc_grad_with_block_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const Real* prec_in,
                    Real* err_out,
                    Real* beta_out,
                    T* psi_out,
                    T* hpsi_out,
                    T* grad_out,
                    T* grad_old_out,
                    const int& n_basis,
                    const int& n_basis_max,
                    const int& n_band)
    {
        for (int band_idx = 0; band_idx < n_band; band_idx++)
        {
            Real err = 0.0;
            Real beta = 0.0;
            Real epsilo = 0.0;
            Real grad_2 = {0.0};
            T grad_1 = {0.0, 0.0};
            auto A = reinterpret_cast<const Real*>(psi_out + band_idx * n_basis_max);
            Real norm = BlasConnector::dot(2 * n_basis, A, 1, A, 1);
            norm = 1.0 / sqrt(norm);
            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                psi_out[item] *= norm;
                hpsi_out[item] *= norm;
                epsilo += std::real(hpsi_out[item] * std::conj(psi_out[item]));
            }
            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                grad_1 = hpsi_out[item] - epsilo * psi_out[item];
                grad_2 = std::norm(grad_1);
                err += grad_2;
                beta += grad_2 / prec_in[basis_idx]; /// Mark here as we should div the prec?
            }
            for (int basis_idx = 0; basis_idx < n_basis; basis_idx++)
            {
                auto item = band_idx * n_basis_max + basis_idx;
                grad_1 = hpsi_out[item] - epsilo * psi_out[item];
                grad_out[item] = -grad_1 / prec_in[basis_idx] + beta / beta_out[band_idx] * grad_old_out[item];
            }
            beta_out[band_idx] = beta;
            err_out[band_idx] = sqrt(err);
        }
    }
};

template <typename FPTYPE>
struct dot_real_op<FPTYPE, base_device::DEVICE_CPU>
{
    FPTYPE operator()(const base_device::DEVICE_CPU* d,
                      const int& dim,
                      const FPTYPE* psi_L,
                      const FPTYPE* psi_R,
                      const bool reduce)
    {
        FPTYPE result = BlasConnector::dot(dim, psi_L, 1, psi_R, 1);
        if (reduce)
        {
            Parallel_Reduce::reduce_pool(result);
        }
        return result;
    }
};

// CPU specialization of actual computation.
template <typename FPTYPE>
struct dot_real_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>
{
    FPTYPE operator()(const base_device::DEVICE_CPU* d,
                      const int& dim,
                      const std::complex<FPTYPE>* psi_L,
                      const std::complex<FPTYPE>* psi_R,
                      const bool reduce)
    {
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // qianrui modify 2021-3-14
        // Note that  ddot_(2*dim,a,1,b,1) = REAL( zdotc_(dim,a,1,b,1) )
        const FPTYPE* pL = reinterpret_cast<const FPTYPE*>(psi_L);
        const FPTYPE* pR = reinterpret_cast<const FPTYPE*>(psi_R);
        FPTYPE result = BlasConnector::dot(2 * dim, pL, 1, pR, 1);
        if (reduce)
        {
            Parallel_Reduce::reduce_pool(result);
        }
        return result;
    }
};

template <typename T>
struct vector_div_constant_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const base_device::DEVICE_CPU* d, const int dim, T* result, const T* vector, const Real constant)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(Real))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector[i] / constant;
        }
    }
};

template <typename T>
struct vector_mul_vector_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const base_device::DEVICE_CPU* d, const int& dim, T* result, const T* vector1, const Real* vector2)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(Real))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] * vector2[i];
        }
    }
};

template <typename T>
struct vector_div_vector_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const base_device::DEVICE_CPU* d, const int& dim, T* result, const T* vector1, const Real* vector2)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(Real))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] / vector2[i];
        }
    }
};

template <typename T>
struct constantvector_addORsub_constantVector_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const base_device::DEVICE_CPU* d,
                    const int& dim,
                    T* result,
                    const T* vector1,
                    const Real constant1,
                    const T* vector2,
                    const Real constant2)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 8192 / sizeof(T))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] * constant1 + vector2[i] * constant2;
        }
    }
};

template <typename FPTYPE>
struct scal_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* /*ctx*/,
                    const int& N,
                    const std::complex<FPTYPE>* alpha,
                    std::complex<FPTYPE>* X,
                    const int& incx)
    {
        BlasConnector::scal(N, *alpha, X, incx);
    }
};

template <typename T>
struct gemv_op<T, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* d,
                    const char& trans,
                    const int& m,
                    const int& n,
                    const T* alpha,
                    const T* A,
                    const int& lda,
                    const T* X,
                    const int& incx,
                    const T* beta,
                    T* Y,
                    const int& incy)
    {
        BlasConnector::gemv(trans, m, n, *alpha, A, lda, X, incx, *beta, Y, incy);
    }
};

template <typename T>
struct axpy_op<T, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* /*ctx*/,
                    const int& dim,
                    const T* alpha,
                    const T* X,
                    const int& incX,
                    T* Y,
                    const int& incY)
    {
        BlasConnector::axpy(dim, *alpha, X, incX, Y, incY);
    }
};

template <typename T>
struct gemm_op<T, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* /*ctx*/,
                    const char& transa,
                    const char& transb,
                    const int& m,
                    const int& n,
                    const int& k,
                    const T* alpha,
                    const T* a,
                    const int& lda,
                    const T* b,
                    const int& ldb,
                    const T* beta,
                    T* c,
                    const int& ldc)
    {
        BlasConnector::gemm(transb, transa, n, m, k, *alpha, b, ldb, a, lda, *beta, c, ldc);
    }
};

#ifdef __DSP
template <typename T>
struct gemm_op_mt<T, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* /*ctx*/,
                    const char& transa,
                    const char& transb,
                    const int& m,
                    const int& n,
                    const int& k,
                    const T* alpha,
                    const T* a,
                    const int& lda,
                    const T* b,
                    const int& ldb,
                    const T* beta,
                    T* c,
                    const int& ldc)
    {
        BlasConnector::gemm(transb, transa, n, m, k, *alpha, b, ldb, a, lda, *beta, c, ldc, base_device::AbacusDevice_t::DspDevice);
    }
};
#endif

template <typename T>
struct matrixTranspose_op<T, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* d,
                    const int& row,
                    const int& col,
                    const T* input_matrix,
                    T* output_matrix)
    {
        T* temp = nullptr;
        base_device::memory::resize_memory_op<T, base_device::DEVICE_CPU>()(d, temp, row * col, "MTransOp");
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 8192 / sizeof(T))
#endif
        for (int j = 0; j < col; j++)
        {
            for (int i = 0; i < row; i++)
            {
                temp[j * row + i] = input_matrix[i * col + j];
            }
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 8192 / sizeof(T))
#endif
        for (int i = 0; i < row * col; i++)
        {
            output_matrix[i] = temp[i];
        }
        base_device::memory::delete_memory_op<T, base_device::DEVICE_CPU>()(d, temp);
    }
};

template <typename T>
struct matrixSetToAnother<T, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* d, const int& n, const T* A, const int& LDA, T* B, const int& LDB)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 8192 / sizeof(T))
#endif
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < LDA; j++)
            {
                B[i * LDB + j] = A[i * LDA + j];
            }
        }
    }
};

// Explicitly instantiate functors for the types of functor registered.
template struct scal_op<float, base_device::DEVICE_CPU>;
template struct axpy_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct gemv_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct gemv_op<float, base_device::DEVICE_CPU>;
template struct gemm_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct gemm_op<float, base_device::DEVICE_CPU>;
template struct dot_real_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct vector_div_constant_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct vector_mul_vector_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct vector_div_vector_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct constantvector_addORsub_constantVector_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct matrixTranspose_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct matrixSetToAnother<std::complex<float>, base_device::DEVICE_CPU>;
template struct calc_grad_with_block_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct line_minimize_with_block_op<std::complex<float>, base_device::DEVICE_CPU>;

template struct scal_op<double, base_device::DEVICE_CPU>;
template struct axpy_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct gemv_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct gemv_op<double, base_device::DEVICE_CPU>;
template struct gemm_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct gemm_op<double, base_device::DEVICE_CPU>;
template struct dot_real_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct vector_div_constant_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct vector_mul_vector_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct vector_div_vector_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct constantvector_addORsub_constantVector_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct matrixTranspose_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct matrixSetToAnother<std::complex<double>, base_device::DEVICE_CPU>;
template struct calc_grad_with_block_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct line_minimize_with_block_op<std::complex<double>, base_device::DEVICE_CPU>;

#ifdef __LCAO
template struct axpy_op<double, base_device::DEVICE_CPU>;
template struct dot_real_op<double, base_device::DEVICE_CPU>;
template struct vector_mul_vector_op<double, base_device::DEVICE_CPU>;
template struct vector_div_constant_op<double, base_device::DEVICE_CPU>;
template struct vector_div_vector_op<double, base_device::DEVICE_CPU>;
template struct matrixTranspose_op<double, base_device::DEVICE_CPU>;
template struct matrixSetToAnother<double, base_device::DEVICE_CPU>;
template struct constantvector_addORsub_constantVector_op<double, base_device::DEVICE_CPU>;
#endif
#ifdef __DSP
template struct gemm_op_mt<std::complex<float>, base_device::DEVICE_CPU>;
template struct gemm_op_mt<std::complex<double>, base_device::DEVICE_CPU>;
#endif
} // namespace hsolver
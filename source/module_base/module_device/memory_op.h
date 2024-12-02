#ifndef MODULE_DEVICE_MEMORY_H_
#define MODULE_DEVICE_MEMORY_H_

#include "types.h"

#include <complex>
#include <cstddef>

namespace base_device
{

namespace memory
{

template <typename FPTYPE, typename Device>
struct resize_memory_op
{
    /// @brief Allocate memory for a given pointer. Note this op will free the pointer first.
    ///
    /// Input Parameters
    /// \param dev : the type of computing device
    /// \param size : array size
    /// \param record_string : label for memory record
    ///
    /// Output Parameters
    /// \param arr : allocated array
    void operator()(const Device* dev, FPTYPE*& arr, const size_t size, const char* record_in = nullptr);
};

template <typename FPTYPE, typename Device>
struct set_memory_op
{
    /// @brief memset for multi-device
    ///
    /// Input Parameters
    /// \param dev : the type of computing device
    /// \param var : the specified constant value
    /// \param size : array size
    ///
    /// Output Parameters
    /// \param arr : output array initialized by the input value
    void operator()(const Device* dev, FPTYPE* arr, const int var, const size_t size);
};

template <typename FPTYPE, typename Device_out, typename Device_in>
struct synchronize_memory_op
{
    /// @brief memcpy for multi-device
    ///
    /// Input Parameters
    /// \param dev_out : the type of computing device of arr_out
    /// \param dev_in : the type of computing device of arr_in
    /// \param arr_in : input array
    /// \param size : array size
    ///
    /// Output Parameters
    /// \param arr_out : output array initialized by the input array
    void operator()(const Device_out* dev_out,
                    const Device_in* dev_in,
                    FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size);
};

template <typename FPTYPE_out, typename FPTYPE_in, typename Device_out, typename Device_in>
struct cast_memory_op
{
    /// @brief memcpy for multi-device
    ///
    /// Input Parameters
    /// \param dev_out : the type of computing device of arr_out
    /// \param dev_in : the type of computing device of arr_in
    /// \param arr_in : input array
    /// \param size : array size
    ///
    /// Output Parameters
    /// \param arr_out : output array initialized by the input array
    void operator()(const Device_out* dev_out,
                    const Device_in* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size);
};

template <typename FPTYPE, typename Device>
struct delete_memory_op
{
    /// @brief free memory for multi-device
    ///
    /// Input Parameters
    /// \param dev : the type of computing device
    /// \param arr : the input array
    void operator()(const Device* dev, FPTYPE* arr);
};


#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
// Partially specialize operator for base_device::GpuDevice.
template <typename FPTYPE>
struct resize_memory_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* dev,
                    FPTYPE*& arr,
                    const size_t size,
                    const char* record_in = nullptr);
};

template <typename FPTYPE>
struct set_memory_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* dev, FPTYPE* arr, const int var, const size_t size);
};

template <typename FPTYPE>
struct synchronize_memory_op<FPTYPE, base_device::DEVICE_CPU, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_CPU* dev_out,
                    const base_device::DEVICE_GPU* dev_in,
                    FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size);
};
template <typename FPTYPE>
struct synchronize_memory_op<FPTYPE, base_device::DEVICE_GPU, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_GPU* dev_out,
                    const base_device::DEVICE_CPU* dev_in,
                    FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size);
};
template <typename FPTYPE>
struct synchronize_memory_op<FPTYPE, base_device::DEVICE_GPU, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* dev_out,
                    const base_device::DEVICE_GPU* dev_in,
                    FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size);
};

template <typename FPTYPE>
struct delete_memory_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* dev, FPTYPE* arr);
};
#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM

#ifdef __DSP

template <typename FPTYPE, typename Device>
struct resize_memory_op_mt
{
    /// @brief Allocate memory for a given pointer. Note this op will free the pointer first.
    ///
    /// Input Parameters
    /// \param dev : the type of computing device
    /// \param size : array size
    /// \param record_string : label for memory record
    ///
    /// Output Parameters
    /// \param arr : allocated array
    void operator()(const Device* dev, FPTYPE*& arr, const size_t size, const char* record_in = nullptr);
};

template <typename FPTYPE, typename Device>
struct delete_memory_op_mt
{
    /// @brief free memory for multi-device
    ///
    /// Input Parameters
    /// \param dev : the type of computing device
    /// \param arr : the input array
    void operator()(const Device* dev, FPTYPE* arr);
};

#endif // __DSP

} // end of namespace memory
} // end of namespace base_device

using resmem_sh_op = base_device::memory::resize_memory_op<float, base_device::DEVICE_CPU>;
using resmem_dh_op = base_device::memory::resize_memory_op<double, base_device::DEVICE_CPU>;
using resmem_ch_op = base_device::memory::resize_memory_op<std::complex<float>, base_device::DEVICE_CPU>;
using resmem_zh_op = base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_CPU>;

using resmem_sd_op = base_device::memory::resize_memory_op<float, base_device::DEVICE_GPU>;
using resmem_dd_op = base_device::memory::resize_memory_op<double, base_device::DEVICE_GPU>;
using resmem_cd_op = base_device::memory::resize_memory_op<std::complex<float>, base_device::DEVICE_GPU>;
using resmem_zd_op = base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>;

using setmem_sh_op = base_device::memory::set_memory_op<float, base_device::DEVICE_CPU>;
using setmem_dh_op = base_device::memory::set_memory_op<double, base_device::DEVICE_CPU>;
using setmem_ch_op = base_device::memory::set_memory_op<std::complex<float>, base_device::DEVICE_CPU>;
using setmem_zh_op = base_device::memory::set_memory_op<std::complex<double>, base_device::DEVICE_CPU>;

using setmem_sd_op = base_device::memory::set_memory_op<float, base_device::DEVICE_GPU>;
using setmem_dd_op = base_device::memory::set_memory_op<double, base_device::DEVICE_GPU>;
using setmem_cd_op = base_device::memory::set_memory_op<std::complex<float>, base_device::DEVICE_GPU>;
using setmem_zd_op = base_device::memory::set_memory_op<std::complex<double>, base_device::DEVICE_GPU>;

using delmem_sh_op = base_device::memory::delete_memory_op<float, base_device::DEVICE_CPU>;
using delmem_dh_op = base_device::memory::delete_memory_op<double, base_device::DEVICE_CPU>;
using delmem_ch_op = base_device::memory::delete_memory_op<std::complex<float>, base_device::DEVICE_CPU>;
using delmem_zh_op = base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_CPU>;

using delmem_sd_op = base_device::memory::delete_memory_op<float, base_device::DEVICE_GPU>;
using delmem_dd_op = base_device::memory::delete_memory_op<double, base_device::DEVICE_GPU>;
using delmem_cd_op = base_device::memory::delete_memory_op<std::complex<float>, base_device::DEVICE_GPU>;
using delmem_zd_op = base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>;

using syncmem_s2s_h2h_op
    = base_device::memory::synchronize_memory_op<float, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
using syncmem_s2s_h2d_op
    = base_device::memory::synchronize_memory_op<float, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
using syncmem_s2s_d2h_op
    = base_device::memory::synchronize_memory_op<float, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;
using syncmem_d2d_h2h_op
    = base_device::memory::synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
using syncmem_d2d_h2d_op
    = base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
using syncmem_d2d_d2h_op
    = base_device::memory::synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;

using syncmem_c2c_h2h_op
    = base_device::memory::synchronize_memory_op<std::complex<float>, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
using syncmem_c2c_h2d_op
    = base_device::memory::synchronize_memory_op<std::complex<float>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
using syncmem_c2c_d2h_op
    = base_device::memory::synchronize_memory_op<std::complex<float>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;
using syncmem_z2z_h2h_op
    = base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
using syncmem_z2z_h2d_op
    = base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
using syncmem_z2z_d2h_op
    = base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;

using castmem_s2d_h2h_op
    = base_device::memory::cast_memory_op<double, float, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
using castmem_s2d_h2d_op
    = base_device::memory::cast_memory_op<double, float, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
using castmem_s2d_d2h_op
    = base_device::memory::cast_memory_op<double, float, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;
using castmem_d2s_h2h_op
    = base_device::memory::cast_memory_op<float, double, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
using castmem_d2s_h2d_op
    = base_device::memory::cast_memory_op<float, double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
using castmem_d2s_d2h_op
    = base_device::memory::cast_memory_op<float, double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;

using castmem_c2z_h2h_op = base_device::memory::
    cast_memory_op<std::complex<double>, std::complex<float>, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
using castmem_c2z_h2d_op = base_device::memory::
    cast_memory_op<std::complex<double>, std::complex<float>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
using castmem_c2z_d2h_op = base_device::memory::
    cast_memory_op<std::complex<double>, std::complex<float>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;
using castmem_z2c_h2h_op = base_device::memory::
    cast_memory_op<std::complex<float>, std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
using castmem_z2c_h2d_op = base_device::memory::
    cast_memory_op<std::complex<float>, std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
using castmem_z2c_d2h_op = base_device::memory::
    cast_memory_op<std::complex<float>, std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;

static base_device::DEVICE_CPU* cpu_ctx = {};
static base_device::DEVICE_GPU* gpu_ctx = {};
#endif // MODULE_DEVICE_MEMORY_H_
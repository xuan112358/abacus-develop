#ifndef HAMILTPW_H
#define HAMILTPW_H

#include "module_base/macros.h"
#include "module_cell/klist.h"
#include "module_elecstate/potentials/potential_new.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "module_hsolver/kernels/math_kernel_op.h"

namespace hamilt
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class HamiltPW : public Hamilt<T, Device>
{
  private:
    // Note GetTypeReal<T>::type will 
    // return T if T is real type(float, double), 
    // otherwise return the real type of T(complex<float>, complex<double>)
    using Real = typename GetTypeReal<T>::type;
  public:
    HamiltPW(elecstate::Potential* pot_in, ModulePW::PW_Basis_K* wfc_basis, K_Vectors* p_kv, pseudopot_cell_vnl* nlpp);
    template<typename T_in, typename Device_in = Device>
    explicit HamiltPW(const HamiltPW<T_in, Device_in>* hamilt);
    ~HamiltPW();

    // for target K point, update consequence of hPsi() and matrix()
    void updateHk(const int ik) override;

    void sPsi(const T* psi_in, // psi
              T* spsi,         // spsi
              const int nrow,  // dimension of spsi: nbands * nrow
              const int npw,   // number of plane waves
              const int nbands // number of bands
    ) const override;

protected:
    // used in sPhi, which are calculated in hPsi or sPhi
    const pseudopot_cell_vnl* ppcell = nullptr;
    mutable T* vkb = nullptr;
    Real* qq_nt = nullptr;
    T* qq_so = nullptr;

    Device* ctx = {};
    using gemv_op = hsolver::gemv_op<T, Device>;
    using gemm_op = hsolver::gemm_op<T, Device>;
    using setmem_complex_op = base_device::memory::set_memory_op<T, Device>;
    using resmem_complex_op = base_device::memory::resize_memory_op<T, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<T, Device>;
    using syncmem_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
};

} // namespace hamilt

#endif
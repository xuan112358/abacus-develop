#ifndef STRESS_FUNC_H
#define STRESS_FUNC_H

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/parallel_reduce.h"
#include "module_base/realarray.h"
#include "module_base/vector3.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/stress_op.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_psi/psi.h"

//-------------------------------------------------------------------
// mohan reconstruction note: 2021-02-07
// the stress code needs reconstructions (by Daye Zheng)
// 1) add explanations for each function, each variable, for
// main procedures, make the code readable
// 2) divide the stress class into several files, each file
// deals with only one part of the stress, it is convenient for
// the next-step reconstruction, for example, we want to make
// pw as an external variable instead of a global variable
// 3) for PW and LCAO, keeps only one copy of the code, for example
// the ewald term needs only one copy, it will reduce the
// time to maintain both copies of stress codes.
// 4) remain openning interfaces for others to contribute stress
// codes, for example, molecular dynamics will have an ionic stress
// term, +U? exx? may have other stress terms.
// 5) delete useless comments and tests, if you have a useless code,
// please explicitly explain why you want to keep the test
// 6) format should be beautiful! code should be readable like a
// note (let readers be comfortable)
//-------------------------------------------------------------------

//----------------------------------------------------------------
// compute the stress terms in terms of the plane wave basis set
// the stress terms include:
// 1) the stress from the electron kinetic energy
// 2) the stress from the local pseudopotentials
// 3) the stress from the non-local pseudopotentials
// 4) the stress from the Hartree term
// 5) the stress from the non-linear core correction (if any)
// 6) the strees from the exchange-correlation functional term
// 7) the stress from the ewald term (ion-ion intraction under
//		periodic boundary conditions).
// 8) the stress from ionic contributions (for molecular dynamics)
//----------------------------------------------------------------

template <typename FPTYPE, typename Device = base_device::DEVICE_CPU>
class Stress_Func
{
  public:
    Stress_Func(){};
    ~Stress_Func(){};

    // stress functions
    //  1) the stress from the electron kinetic energy
    void stress_kin(ModuleBase::matrix& sigma,
                    const ModuleBase::matrix& wg,
                    ModuleSymmetry::Symmetry* p_symm,
                    K_Vectors* p_kv,
                    ModulePW::PW_Basis_K* wfc_basis,
                    const UnitCell& ucell_in,
                    const psi::Psi<complex<FPTYPE>, Device>* psi_in = nullptr); // electron kinetic part in PW basis

    // 2) the stress from the Hartree term
    void stress_har(const UnitCell& ucell,
                    ModuleBase::matrix& sigma,
                    ModulePW::PW_Basis* rho_basis,
                    const bool is_pw,
                    const Charge* const chr); // hartree part in PW or LCAO basis

    // 3) the stress from the ewald term (ion-ion intraction under
    //		periodic boundary conditions).
    void stress_ewa(const UnitCell& ucell,
                    ModuleBase::matrix& sigma,
                    ModulePW::PW_Basis* rho_basis,
                    const bool is_pw); // ewald part in PW or LCAO basis

    // 4) the stress from the local pseudopotentials
    void stress_loc(const UnitCell& ucell,
                    ModuleBase::matrix& sigma,
                    ModulePW::PW_Basis* rho_basis,
                    const ModuleBase::matrix& vloc,
                    const Structure_Factor* p_sf,
                    const bool is_pw,
                    const Charge* const chr); // local pseudopotential part in PW or LCAO

    void dvloc_of_g(const int& msh,
                    const FPTYPE* rab,
                    const FPTYPE* r,
                    const FPTYPE* vloc_at,
                    const FPTYPE& zp,
                    FPTYPE* dvloc,
                    ModulePW::PW_Basis* rho_basis,
                    const UnitCell& ucell_in); // used in local pseudopotential stress

    /**
     * @brief compute the derivative of the coulomb potential in reciprocal space
     *        D V(g^2) / D g^2 = 4pi e^2/omegai /G^4
     *
     */
    void dvloc_coulomb(const UnitCell& ucell,
                       const FPTYPE& zp,
                       FPTYPE* dvloc,
                       ModulePW::PW_Basis* rho_basis); // used in local pseudopotential stress

    // 5) the stress from the non-linear core correction (if any)
    void stress_cc(ModuleBase::matrix& sigma,
                   ModulePW::PW_Basis* rho_basis,
                   UnitCell& ucell,
                   const Structure_Factor* p_sf,
                   const bool is_pw,
                   const bool *numeric,
                   const Charge* const chr); // nonlinear core correction stress in PW or LCAO basis

    void deriv_drhoc(const bool& numeric,
                     const double& omega,
                     const double& tpiba2,
                     const int mesh,
                     const FPTYPE* r,
                     const FPTYPE* rab,
                     const FPTYPE* rhoc,
                     FPTYPE* drhocg,
                     ModulePW::PW_Basis* rho_basis,
                     int type); // used in nonlinear core correction stress
                     
    // 6) the stress from the exchange-correlation functional term
    void stress_gga(const UnitCell& ucell,
                    ModuleBase::matrix& sigma,
                    ModulePW::PW_Basis* rho_basis,
                    const Charge* const chr); // gga part in both PW and LCAO basis
    void stress_mgga(const UnitCell& ucell,
                     ModuleBase::matrix& sigma,
                     const ModuleBase::matrix& wg,
                     const ModuleBase::matrix& v_ofk,
                     const Charge* const chr,
                     K_Vectors* p_kv,
                     ModulePW::PW_Basis_K* wfc_basis,
                     const psi::Psi<complex<FPTYPE>, Device>* psi_in); // gga part in PW basis

    // 7) the stress from the non-local pseudopotentials
    /**
     * @brief This routine computes the atomic force of non-local pseudopotential
     *    Stress^{NL}_{ij} = -1/\Omega \sum_{n,k}f_{nk}\sum_I \sum_{lm,l'm'}D_{l,l'}^{I} [
     *               \sum_G \langle c_{nk}(\mathbf{G+K})|\beta_{lm}^I(\mathbf{G+K})\rangle *
     *               \sum_{G'}\langle \partial \beta_{lm}^I(\mathbf{G+K})/\partial \varepsilon_{ij}
     * |c_{nk}(\mathbf{G+K})\rangle ] there would be three parts in the above equation: (1) sum over becp and dbecp with
     * D_{l,l'}^{I} ----- first line in the above equation (2) calculate becp = <psi | beta> ----- second line in the
     * above equation (3) calculate dbecp = <psi | dbeta> ----- third line in the above equation
     */
    void stress_nl(ModuleBase::matrix& sigma,
                   const ModuleBase::matrix& wg,
                   const ModuleBase::matrix& ekb,
                   Structure_Factor* p_sf,
                   K_Vectors* p_kv,
                   ModuleSymmetry::Symmetry* p_symm,
                   ModulePW::PW_Basis_K* wfc_basis,
                   const psi::Psi<complex<FPTYPE>, Device>* psi_in,
                   const pseudopot_cell_vnl& nlpp_in,
                   const UnitCell& ucell_in); // nonlocal part in PW basis
    // 8) the stress from the DFT+U and DeltaSpin calculations
    /**
     * @brief This routine computes the stress contribution from the DFT+U and DeltaSpin calculations
     *    Stress^{NL}_{ij} = -1/\Omega \sum_{n,k}f_{nk}\sum_I \sum_{lm,l'm'}(V^U_{lmm'\sigma\sigma'} +
     * f(\lambda,\sigma\sigma')) [ \sum_G \langle c_{nk}(\mathbf{G+K})|\alpha_{lm}^I(\mathbf{G+K})\rangle *
     *               \sum_{G'}\langle \partial \alpha_{lm}^I(\mathbf{G+K})/\partial \varepsilon_{ij}
     * |c_{nk}(\mathbf{G+K})\rangle ] there would be three parts in the above equation: (1) sum over becp and dbecp with
     * f(U+\lambda, \sigma\sigma', lmm')^{I} ----- first line in the above equation (2) calculate becp = <psi | alpha>
     * ----- second line in the above equation (3) calculate dbecp = <psi | dalpha> ----- third line in the above
     * equation
     */
    void stress_onsite(ModuleBase::matrix& sigma,
                       const ModuleBase::matrix& wg,
                       const ModulePW::PW_Basis_K* wfc_basis,
                       const UnitCell& ucell_in,
                       const psi::Psi<complex<FPTYPE>, Device>* psi_in,
                       ModuleSymmetry::Symmetry* p_symm); // nonlocal part in PW basis

    void get_dvnl1(ModuleBase::ComplexMatrix& vkb,
                   const int ik,
                   const int ipol,
                   Structure_Factor* p_sf,
                   ModulePW::PW_Basis_K* wfc_basis); // used in nonlocal part in PW basis
    void get_dvnl2(ModuleBase::ComplexMatrix& vkb,
                   const int ik,
                   Structure_Factor* p_sf,
                   ModulePW::PW_Basis_K* wfc_basis); // used in nonlocal part in PW basis

    FPTYPE Polynomial_Interpolation_nl(const ModuleBase::realArray& table,
                                       const int& dim1,
                                       const int& dim2,
                                       const int& dim3,
                                       const FPTYPE& table_interval,
                                       const FPTYPE& x);

    /**
     * @brief Compute the derivatives of the radial Fourier transform of the Q functions
     *
     * This routine computes the derivatives of the Fourier transform of
     * the Q function needed in stress assuming that the radial fourier
     * transform is already computed and stored in table qrad.
     * The formula implemented here is:
     *
     *   dq(g,i,j) = sum_lm (-i)^l ap(lm,i,j) *
     *              ( yr_lm(g^) dqrad(g,l,i,j) + dyr_lm(g^) qrad(g,l,i,j))
     *
     * @param ih [in] the first index of Q
     * @param jh [in] the second index of Q
     * @param itype [in] the atomic type
     * @param ipol [in] the polarization of the derivative
     * @param ng [in] the number of G vectors
     * @param g [in] the G vectors
     * @param qnorm [in] the norm of q+g vectors
     * @param tpiba [in] 2pi/a factor, multiplies G vectors
     * @param ylmk0 [in] the real spherical harmonics
     * @param dylmk0 [in] derivetives of spherical harmonics
     * @param dqg [out] the Fourier transform of interest
     */
    void dqvan2(const pseudopot_cell_vnl& ppcell_in,
                const int ih,
                const int jh,
                const int itype,
                const int ipol,
                const int ng,
                const ModuleBase::Vector3<FPTYPE>* g,
                const FPTYPE* qnorm,
                const FPTYPE& tpiba,
                const ModuleBase::matrix& ylmk0,
                const ModuleBase::matrix& dylmk0,
                std::complex<FPTYPE>* dqg);

  protected:
    Device* ctx = {};
    base_device::DEVICE_CPU* cpu_ctx = {};
    base_device::AbacusDevice_t device = {};
  private:
    using gemm_op = hsolver::gemm_op<std::complex<FPTYPE>, Device>;
    using cal_stress_nl_op = hamilt::cal_stress_nl_op<FPTYPE, Device>;
    using cal_dbecp_noevc_nl_op = hamilt::cal_dbecp_noevc_nl_op<FPTYPE, Device>;

    using resmem_complex_op = base_device::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using resmem_complex_h_op = base_device::memory::resize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>;
    using setmem_complex_op = base_device::memory::set_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_h_op = base_device::memory::delete_memory_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>;
    using syncmem_complex_h2d_op
        = base_device::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, base_device::DEVICE_CPU>;
    using syncmem_complex_d2h_op
        = base_device::memory::synchronize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_CPU, Device>;

    using resmem_var_op = base_device::memory::resize_memory_op<FPTYPE, Device>;
    using resmem_var_h_op = base_device::memory::resize_memory_op<FPTYPE, base_device::DEVICE_CPU>;
    using setmem_var_op = base_device::memory::set_memory_op<FPTYPE, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<FPTYPE, Device>;
    using delmem_var_h_op = base_device::memory::delete_memory_op<FPTYPE, base_device::DEVICE_CPU>;
    using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<FPTYPE, Device, base_device::DEVICE_CPU>;
    using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<FPTYPE, base_device::DEVICE_CPU, Device>;

    using resmem_int_op = base_device::memory::resize_memory_op<int, Device>;
    using delmem_int_op = base_device::memory::delete_memory_op<int, Device>;
    using syncmem_int_h2d_op = base_device::memory::synchronize_memory_op<int, Device, base_device::DEVICE_CPU>;

    using cal_vq_op = hamilt::cal_vq_op<FPTYPE, Device>;
    using cal_vq_deri_op = hamilt::cal_vq_deri_op<FPTYPE, Device>;

    using cal_vkb_op = hamilt::cal_vkb_op<FPTYPE, Device>;
    using cal_vkb_deri_op = hamilt::cal_vkb_deri_op<FPTYPE, Device>;

  protected:
    pseudopot_cell_vnl* nlpp = nullptr;
    const UnitCell* ucell = nullptr;
};

#endif

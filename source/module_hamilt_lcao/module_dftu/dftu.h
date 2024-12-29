//==========================================================
// Author: Xin Qu
// DATE : 2019-12-10
//==========================================================
#ifndef DFTU_H
#define DFTU_H

#include "module_cell/klist.h"
#include "module_cell/unitcell.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#ifdef __LCAO
#include "module_elecstate/module_charge/charge_mixing.h"
#include "module_hamilt_general/hamilt.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/force_stress_arrays.h" // mohan add 2024-06-15
#endif

#include <string>
#include <vector>

//==========================================================
// CLASS :
// NAME : DTFU (DFT+U)
//==========================================================
namespace ModuleDFTU
{

class DFTU
{

  public:
    DFTU(); // constructor
    ~DFTU(); // deconstructor

    //=============================================================
    // In dftu.cpp
    // Initialization & Calculating energy
    //=============================================================
  public:
    // allocate relevant data strcutures
    void init(UnitCell& cell, // unitcell class
              const Parallel_Orbitals* pv,
              const int nks
#ifdef __LCAO
              , const LCAO_Orbitals* orb = nullptr
#endif
              );
    
    static DFTU* get_instance();

    // calculate the energy correction
    void cal_energy_correction(const UnitCell& ucell, const int istep);
    double get_energy(){return EU;}
    void uramping_update(); // update U by uramping
    bool u_converged(); // check if U is converged

    std::vector<double> U = {}; // U (Hubbard parameter U)
    std::vector<double> U0; // U0 (target Hubbard parameter U0)
    std::vector<int> orbital_corr = {}; //
    double uramping; // increase U by uramping, default is -1.0
    int omc; // occupation matrix control
    int mixing_dftu; //whether to mix locale

    double EU; //+U energy
  private:
    const Parallel_Orbitals* paraV = nullptr;
    int cal_type = 3; // 1:dftu_tpye=1, dc=1; 2:dftu_type=1, dc=2; 3:dftu_tpye=2, dc=1; 4:dftu_tpye=2, dc=2;

    // FIXME: the following variable does not have static lifetime;
    // while the present class is used via a global variable. This has
    // potential to cause dangling pointer issues.
#ifdef __LCAO
    const LCAO_Orbitals* ptr_orb_ = nullptr;
    std::vector<double> orb_cutoff_;
#endif
    
    // transform between iwt index and it, ia, L, N and m index
    std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>
        iatlnmipol2iwt; // iatlnm2iwt[iat][l][n][m][ipol]

#ifdef __LCAO
    //=============================================================
    // In dftu_hamilt.cpp
    // For calculating contribution to Hamiltonian matrices
    //=============================================================
  public:
    void cal_eff_pot_mat_complex(const int ik, std::complex<double>* eff_pot, const std::vector<int>& isk, const std::complex<double>* sk);
    void cal_eff_pot_mat_real(const int ik, double* eff_pot, const std::vector<int>& isk, const double* sk);
    void cal_eff_pot_mat_R_double(const int ispin, double* SR, double* HR);
    void cal_eff_pot_mat_R_complex_double(const int ispin, std::complex<double>* SR, std::complex<double>* HR);
#endif

    //=============================================================
    // In dftu_occup.cpp
    // For calculating occupation matrix and saving to locale
    // and other operations of locale: copy,zero out,mix
    //=============================================================
  public:
    /// interface for PW base
    /// calculate the local occupation number matrix for PW based wave functions
    void cal_occ_pw(const int iter, const void* psi_in, const ModuleBase::matrix& wg_in, const UnitCell& cell, const double& mixing_beta);
    /// calculate the local DFT+U effective potential matrix for PW base.
    void cal_VU_pot_pw(const int spin);
    /// get effective potential matrix for PW base
    const std::complex<double>* get_eff_pot_pw(const int iat) const { return &(eff_pot_pw[this->eff_pot_pw_index[iat]]); }
    int get_size_eff_pot_pw() const { return eff_pot_pw.size(); }

#ifdef __LCAO
    // calculate the local occupation number matrix
    void cal_occup_m_k(const int iter, 
                       const UnitCell& ucell,
                       const std::vector<std::vector<std::complex<double>>>& dm_k, 
                       const K_Vectors& kv, 
                       const double& mixing_beta, 
                       hamilt::Hamilt<std::complex<double>>* p_ham);
    void cal_occup_m_gamma(const int iter, 
                           const UnitCell& ucell,
                           const std::vector<std::vector<double>>& dm_gamma, 
                           const double& mixing_beta, 
                           hamilt::Hamilt<double>* p_ham);
#endif

    // dftu can be calculated only after locale has been initialed
    bool initialed_locale = false;

  private:
    void copy_locale(const UnitCell& ucell);
    void zero_locale(const UnitCell& ucell);
    void mix_locale(const UnitCell& ucell,const double& mixing_beta);

    std::vector<std::complex<double>> eff_pot_pw;
    std::vector<int> eff_pot_pw_index;

public:
    // local occupancy matrix of the correlated subspace
    // locale: the out put local occupation number matrix of correlated electrons in the current electronic step
    // locale_save: the input local occupation number matrix of correlated electrons in the current electronic step
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>> locale; // locale[iat][l][n][spin](m1,m2)
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>> locale_save; // locale_save[iat][l][n][spin](m1,m2)
#ifdef __LCAO
private:
    //=============================================================
    // In dftu_tools.cpp
    // For calculating onsite potential, which is used
    // for both Hamiltonian and force/stress
    //=============================================================

    void cal_VU_pot_mat_complex(const int spin, const bool newlocale, std::complex<double>* VU);
    void cal_VU_pot_mat_real(const int spin, const bool newlocale, double* VU);

    double get_onebody_eff_pot(const int T,
                               const int iat,
                               const int L,
                               const int N,
                               const int spin,
                               const int m0,
                               const int m1,
                               const bool newlocale);

    //=============================================================
    // In dftu_folding.cpp
    // Subroutines for folding S and dS matrix
    //=============================================================

    void fold_dSR_gamma(const UnitCell& ucell,
                        const Parallel_Orbitals& pv,
                        const Grid_Driver* gd,
                        double* dsloc_x,
                        double* dsloc_y,
                        double* dsloc_z,
                        double* dh_r,
                        const int dim1,
                        const int dim2,
                        double* dSR_gamma);

    // dim = 0 : S, for Hamiltonian
    // dim = 1-3 : dS, for force
    // dim = 4-6 : dS * dR, for stress

    void folding_matrix_k(const UnitCell& ucell,
                          const Grid_Driver& gd,
                          ForceStressArrays& fsr,
                          const Parallel_Orbitals& pv,
                          const int ik,
                          const int dim1,
                          const int dim2,
                          std::complex<double>* mat_k,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d);

    /**
     * @brief new function of folding_S_matrix
     * only for Hamiltonian now, for force and stress will be developed later
     * use HContainer as input and output in mat_k
    */
   void folding_matrix_k_new(const int ik, hamilt::Hamilt<std::complex<double>>* p_ham);

    //=============================================================
    // In dftu_force.cpp
    // For calculating force and stress fomr DFT+U
    //=============================================================
 public:
   void force_stress(const UnitCell& ucell,
                     const Grid_Driver& gd,
                     const elecstate::ElecState* pelec,
                     const Parallel_Orbitals& pv,
                     ForceStressArrays& fsr,
                     ModuleBase::matrix& force_dftu,
                     ModuleBase::matrix& stress_dftu,
                     const K_Vectors& kv);

 private:
   void cal_force_k(const UnitCell& ucell,
                    const Grid_Driver& gd,
                    ForceStressArrays& fsr,
                    const Parallel_Orbitals& pv,
                    const int ik,
                    const std::complex<double>* rho_VU,
                    ModuleBase::matrix& force_dftu,
                    const std::vector<ModuleBase::Vector3<double>>& kvec_d);

   void cal_stress_k(const UnitCell& ucell,
                     const Grid_Driver& gd,
                     ForceStressArrays& fsr,
                     const Parallel_Orbitals& pv,
                     const int ik,
                     const std::complex<double>* rho_VU,
                     ModuleBase::matrix& stress_dftu,
                     const std::vector<ModuleBase::Vector3<double>>& kvec_d);

   void cal_force_gamma(const UnitCell& ucell,
                        const double* rho_VU,
                        const Parallel_Orbitals& pv,
                        double* dsloc_x,
                        double* dsloc_y,
                        double* dsloc_z,
                        ModuleBase::matrix& force_dftu);

   void cal_stress_gamma(const UnitCell& ucell,
                         const Parallel_Orbitals& pv,
                         const Grid_Driver* gd,
                         double* dsloc_x,
                         double* dsloc_y,
                         double* dsloc_z,
                         double* dh_r,
                         const double* rho_VU,
                         ModuleBase::matrix& stress_dftu);
#endif

    //=============================================================
    // In dftu_io.cpp
    // For reading/writing/broadcasting/copying relevant data structures
    //=============================================================
  public:
    void output(const UnitCell& ucell);

  private:
    void write_occup_m(const UnitCell& ucell,
                       std::ofstream& ofs, 
                       bool diag=false);
    void read_occup_m(const UnitCell& ucell,
                      const std::string& fn);
    void local_occup_bcast(const UnitCell& ucell);

    //=============================================================
    // In dftu_yukawa.cpp
    // Relevant for calculating U using Yukawa potential
    //=============================================================

  public:
    bool Yukawa; // 1:use Yukawa potential; 0: do not use Yukawa potential
    void cal_slater_UJ(const UnitCell& ucell, double** rho, const int& nrxx);

  private:
    double lambda; // the parameter in Yukawa potential
    std::vector<std::vector<std::vector<std::vector<double>>>> Fk; // slater integral:Fk[T][L][N][k]
    std::vector<std::vector<std::vector<double>>> U_Yukawa; // U_Yukawa[T][L][N]
    std::vector<std::vector<std::vector<double>>> J_Yukawa; // J_Yukawa[T][L][N]

    void cal_slater_Fk(const UnitCell& ucell,const int L, const int T); // L:angular momnet, T:atom type
    void cal_yukawa_lambda(double** rho, const int& nrxx);

    double spherical_Bessel(const int k, const double r, const double lambda);
    double spherical_Hankel(const int k, const double r, const double lambda);

#ifdef __LCAO
  public:
    /**
     * @brief get the density matrix of target spin
     * nspin = 1 and 4 : ispin should be 0
     * nspin = 2 : ispin should be 0/1
    */
    const hamilt::HContainer<double>* get_dmr(int ispin) const;
    /**
     * @brief set the density matrix for DFT+U calculation
     * if the density matrix is not set or set to nullptr, the DFT+U calculation will not be performed
    */
    void set_dmr(const elecstate::DensityMatrix<double, double>* dm_in_dftu_d);
    void set_dmr(const elecstate::DensityMatrix<std::complex<double>, double>* dm_in_dftu_cd);
  
  private:
    const UnitCell* ucell = nullptr;
    const elecstate::DensityMatrix<double, double>* dm_in_dftu_d = nullptr;
    const elecstate::DensityMatrix<std::complex<double>, double>* dm_in_dftu_cd = nullptr;
#endif
};

#ifdef __LCAO
template <typename T>
void dftu_cal_occup_m(const int iter,
                      const UnitCell& ucell,
                      const std::vector<std::vector<T>>& dm,
                      const K_Vectors& kv,
                      const double& mixing_beta,
                      hamilt::Hamilt<T>* p_ham);
#endif

} // namespace ModuleDFTU

namespace GlobalC
{
	extern ModuleDFTU::DFTU dftu;
}
#endif

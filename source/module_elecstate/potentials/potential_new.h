#ifndef POTENTIALNEW_H
#define POTENTIALNEW_H

#include "module_base/complexmatrix.h"
#include "module_hamilt_general/module_surchem/surchem.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "pot_base.h"

#include <vector>

namespace elecstate
{
/**
 * Potential is the main class of potentials module, it contains:
 * 1. Constructors and deconstructor
 * 2. Func init_pot()
 *     a. need istep for update_for_tddft();
 *     b. need Charge for update_from_charge();
 *     c. it will reset fixed_done to false, v_effective_fixed will be calculated;
 *     d. it should be called after Charge is initialized;
 *     e. it can only be called once in one SCF loop
 * 3. Func pot_register() and components
 *     a. need vector<string> for choose target potentials
 *     b. "local", PotLocal introduces local pseudopotential part of potentials;
 *     c. "hartree", PotHartree introduces Coulombic interaction of electrons part of potentials;
 *     d. "xc", PotXC introduces exchange-correlation including meta-gga part of potentials;
 *     e. "surchem", PotSurChem introduces surface chemistry part of potentials;
 *     f. "efield", PotEfield introduces electronic field including dipole correction part of potentials;
 *     g. "gatefield", PotGate introduces gate field part of potentials;
 * 4. Func update_from_charge()
 *     a. regenerate v_effective
 *     b. if Meta-GGA is choosed, it will regenerate vofk_effective
 * 5. Func update_for_tddft()
 *     a. in principle, it should be added to components, but it related to real time(istep)
 *     b. it should be called after update_from_charge() as a compensation;
 * 6. Func get_vnew()
 *     a. this function is designed for a special demand:
 *         1. update etxc and vtxc when SCF converged and
 *         2. use the final delta_V_eff for calculating force correction
 * 7. Func write_potential()
 * 8. Func write_elecstat_pot()
 * 9. interfaces for v_effective_fixed/v_effective/vofk_effective
 * 10. Func interpolate_vrs()
 *    a. interpolate v_effective on the smooth mesh
 */
class Potential : public PotBase
{
  public:
    // default constructor for UT
    Potential(){};
    // In constructor, size of every potential components should be allocated
    // rho_basis_in is the dense grids, rho_basis_smooth_in is the smooth grids in USPP
    // charge density and potential are defined on dense grids,
    // but effective potential needs to be interpolated on smooth grids in order to compute Veff|psi>
    // Note: rho_basis_in and rho_basis_smooth_in are the same in NCPP
    Potential(const ModulePW::PW_Basis* rho_basis_in,
              const ModulePW::PW_Basis* rho_basis_smooth_in,
              const UnitCell* ucell_in,
              const ModuleBase::matrix* vloc_in,
              Structure_Factor* structure_factors_in,
              surchem* solvent_in,
              double* etxc_in,
              double* vtxc_in);
    ~Potential();

    // initialize potential when SCF begin
    void init_pot(int istep, const Charge*const chg);
    // initialize potential components before SCF
    void pot_register(const std::vector<std::string>& components_list);
    // update potential from current charge
    void update_from_charge(const Charge*const chg, const UnitCell*const ucell);
    // interface for SCF-converged, etxc vtxc for Energy, vnew for force_scc
    void get_vnew(const Charge* chg, ModuleBase::matrix& vnew);

    PotBase* get_pot_type(const std::string& pot_type);

    // interfaces to get values
    ModuleBase::matrix& get_effective_v()
    {
        return this->v_effective;
    }
    const ModuleBase::matrix& get_effective_v() const
    {
        return this->v_effective;
    }

    double* get_effective_v(int is)
    {
        if (this->v_effective.nc > 0)
        {
            return &(this->v_effective(is, 0));
        }
        else
        {
            return nullptr;
        }
    }
    const double* get_effective_v(int is) const
    {
        if (this->v_effective.nc > 0)
        {
            return &(this->v_effective(is, 0));
        }
        else
        {
            return nullptr;
        }
    }
    ModuleBase::matrix& get_effective_vofk()
    {
        return this->vofk_effective;
    }
    const ModuleBase::matrix& get_effective_vofk() const
    {
        return this->vofk_effective;
    }
    double* get_effective_vofk(int is)
    {
        if (this->vofk_effective.nc > 0)
        {
            return &(this->vofk_effective(is, 0));
        }
        else
        {
            return nullptr;
        }
    }
    const double* get_effective_vofk(int is) const
    {
        if (this->vofk_effective.nc > 0)
        {
            return &(this->vofk_effective(is, 0));
        }
        else
        {
            return nullptr;
        }
    }

    ModuleBase::matrix& get_veff_smooth()
    {
        return this->veff_smooth;
    }
    const ModuleBase::matrix& get_veff_smooth() const
    {
        return this->veff_smooth;
    }

    ModuleBase::matrix& get_vofk_smooth()
    {
        return this->vofk_smooth;
    }
    const ModuleBase::matrix& get_vofk_smooth() const
    {
        return this->vofk_smooth;
    }

    template <typename FPTYPE>
    FPTYPE* get_veff_smooth_data();

    template <typename FPTYPE>
    FPTYPE* get_vofk_smooth_data();

    double* get_fixed_v()
    {
        return this->v_effective_fixed.data();
    }
    const double* get_fixed_v() const
    {
        return this->v_effective_fixed.data();
    }


    /// @brief get the value of vloc at G=0;
    /// @return vl(0)
    double get_vl_of_0() const
    {
        return this->vl_of_0;
    }

  private:
    void cal_v_eff(const Charge*const chg, const UnitCell*const ucell, ModuleBase::matrix& v_eff) override;
    void cal_fixed_v(double* vl_pseudo) override;
    // interpolate potential on the smooth mesh if necessary
    void interpolate_vrs();

    void allocate();

    std::vector<double> v_effective_fixed;
    ModuleBase::matrix v_effective;

    ModuleBase::matrix veff_smooth; // used in uspp liuyu 2023-10-12
    ModuleBase::matrix vofk_smooth; // used in uspp liuyu 2023-10-12

    ModuleBase::matrix v_xc; // if PAW is used, vxc must be stored separately

    float *s_veff_smooth = nullptr, *s_vofk_smooth = nullptr;
    double *d_veff_smooth = nullptr, *d_vofk_smooth = nullptr;

    ModuleBase::matrix vofk_effective;

    bool fixed_done = false;

    // gather etxc and vtxc in Potential, will be used in ESolver
    double* etxc_ = nullptr;
    double* vtxc_ = nullptr;

    double vl_of_0 = 0.0;

    std::vector<PotBase*> components;

    const UnitCell* ucell_ = nullptr;
    const ModuleBase::matrix* vloc_ = nullptr;
    Structure_Factor* structure_factors_ = nullptr;
    surchem* solvent_ = nullptr;
};

} // namespace elecstate

#endif
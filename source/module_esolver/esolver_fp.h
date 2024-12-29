#ifndef ESOLVER_FP_H
#define ESOLVER_FP_H

#include "esolver.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_elecstate/elecstate.h"
#include "module_elecstate/module_charge/charge_extra.h"
#include "module_hamilt_general/module_surchem/surchem.h"
#include "module_hamilt_pw/hamilt_pwdft/VL_in_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"

#include <fstream>

//! The First-Principles (FP) Energy Solver Class
/**
 * This class represents components that needed in
 * first-principles energy solver, such as the plane
 * wave basis, the structure factors, and the k points.
 *
 */

namespace ModuleESolver
{
class ESolver_FP : public ESolver
{
  public:
    //! Constructor
    ESolver_FP();

    //! Deconstructor
    virtual ~ESolver_FP();

    //! Initialize of the first-principels energy solver
    virtual void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

  protected:
    //! Something to do before SCF iterations.
    virtual void before_scf(UnitCell& ucell, const int istep);

    //! Something to do after SCF iterations when SCF is converged or comes to the max iter step.
    virtual void after_scf(UnitCell& ucell, const int istep);

    //! ------------------------------------------------------------------------------
    //! These pointers will be deleted in the free_pointers() function every ion step.
    //! ------------------------------------------------------------------------------
    elecstate::ElecState* pelec = nullptr; ///< Electronic states

    //! ------------------------------------------------------------------------------

    //! Electorn charge density
    Charge chr;

    //! Structure factors that used with plane-wave basis set
    Structure_Factor sf;

    //! K points in Brillouin zone
    K_Vectors kv;

    //! Plane-wave basis set for charge density
    ModulePW::PW_Basis* pw_rho;

    //! parallel for rho grid
    Parallel_Grid Pgrid;

    //! pointer to local pseudopotential
    pseudopot_cell_vl locpp;

    /**
     * @brief same as pw_rho for ncpp. Here 'd' stands for 'dense'
     * dense grid for for uspp, used for ultrasoft augmented charge density.
     * charge density and potential are defined on dense grids,
     * but effective potential needs to be interpolated on smooth grids in order to compute Veff|psi>
     */
    ModulePW::PW_Basis* pw_rhod;
    ModulePW::PW_Basis_Big* pw_big; ///< [temp] pw_basis_big class

    //! Charge extrapolation
    Charge_Extra CE;

    // solvent model
    surchem solvent;
};
} // namespace ModuleESolver

#endif

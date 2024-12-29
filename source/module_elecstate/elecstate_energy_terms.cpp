#include "elecstate.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"

namespace elecstate
{

double ElecState::get_hartree_energy()
{
    return H_Hartree_pw::hartree_energy;
}

double ElecState::get_etot_efield()
{
    return Efield::etotefield;
}

double ElecState::get_etot_gatefield()
{
    return Gatefield::etotgatefield;
}

double ElecState::get_solvent_model_Ael()
{
    return surchem::Ael;
}

double ElecState::get_solvent_model_Acav()
{
    return surchem::Acav;
}

double ElecState::get_dftu_energy()
{
    return GlobalC::dftu.get_energy();
}

#ifdef __DEEPKS
double ElecState::get_deepks_E_delta()
{
    return GlobalC::ld.E_delta;
}
double ElecState::get_deepks_E_delta_band()
{
    return GlobalC::ld.e_delta_band;
}
#endif

} // namespace elecstate
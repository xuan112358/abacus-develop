#include "esolver_ks_lcaopw.h"

#include "module_hamilt_pw/hamilt_pwdft/elecond.h"
#include "module_io/input_conv.h"
#include "module_io/nscf_band.h"
#include "module_io/output_log.h"
#include "module_io/write_dos_pw.h"
#include "module_io/write_istate_info.h"
#include "module_io/write_wfc_pw.h"

#include <iostream>

//--------------temporary----------------------------
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"
//-----force-------------------
#include "module_hamilt_pw/hamilt_pwdft/forces.h"
//-----stress------------------
#include "module_hamilt_pw/hamilt_pwdft/stress_pw.h"
//---------------------------------------------------
#include "module_base/memory.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/hsolver_lcaopw.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_io/berryphase.h"
#include "module_io/numerical_basis.h"
#include "module_io/numerical_descriptor.h"
#include "module_io/rho_io.h"
#include "module_io/to_wannier90_pw.h"
#include "module_io/winput.h"
#include "module_io/write_pot.h"
#include "module_io/write_wfc_r.h"

#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>

#include <sys/time.h>
#include "module_hamilt_pw/hamilt_pwdft/hamilt_lcaopw.h"
#include "module_parameter/parameter.h"
#ifdef __LCAO
#include "module_io/write_Vxc_lip.hpp"
#endif

namespace ModuleESolver
{

    template <typename T>
    ESolver_KS_LIP<T>::ESolver_KS_LIP()
    {
        this->classname = "ESolver_KS_LIP";
        this->basisname = "LIP";
    }

    template <typename T>
    void ESolver_KS_LIP<T>::allocate_hsolver()
    {
        this->phsol = new hsolver::HSolverLIP<T>(this->pw_wfc);
    }
    template <typename T>
    void ESolver_KS_LIP<T>::deallocate_hsolver()
    {
        if (this->phsol != nullptr)
        {
            delete reinterpret_cast<hsolver::HSolverLIP<T>*>(this->phsol);
            this->phsol = nullptr;
        }
    }
    template <typename T>
    void ESolver_KS_LIP<T>::allocate_hamilt()
    {
        this->p_hamilt = new hamilt::HamiltLIP<T>(this->pelec->pot, this->pw_wfc, &this->kv
#ifdef __EXX
            , *this->exx_lip
#endif
        );
    }

    template <typename T>
    void ESolver_KS_LIP<T>::before_all_runners(Input& inp, UnitCell& cell)
    {
        ESolver_KS_PW<T>::before_all_runners(inp, cell);
#ifdef __EXX
        if (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "relax"
            || GlobalV::CALCULATION == "cell-relax"
            || GlobalV::CALCULATION == "md")
            if (GlobalC::exx_info.info_global.cal_exx)
            {
                XC_Functional::set_xc_first_loop(cell);
                this->exx_lip = std::unique_ptr<Exx_Lip<T>>(new Exx_Lip<T>(GlobalC::exx_info.info_lip,
                    cell.symm, &this->kv, this->p_wf_init, this->kspw_psi, this->pw_wfc, this->pw_rho, this->sf, &cell, this->pelec));
                // this->exx_lip.init(GlobalC::exx_info.info_lip, cell.symm, &this->kv, this->p_wf_init, this->kspw_psi, this->pw_wfc, this->pw_rho, this->sf, &cell, this->pelec);
            }
#endif
    }

    template <typename T>
    void ESolver_KS_LIP<T>::iter_init(const int istep, const int iter)
    {
        ESolver_KS_PW<T>::iter_init(istep, iter);
#ifdef __EXX
        if (GlobalC::exx_info.info_global.cal_exx && !GlobalC::exx_info.info_global.separate_loop && this->two_level_step)
            this->exx_lip->cal_exx();
#endif
    }

    template <typename T>
    void ESolver_KS_LIP<T>::hamilt2density(const int istep, const int iter, const double ethr)
    {
        ModuleBase::TITLE("ESolver_KS_LIP", "hamilt2density");
        ModuleBase::timer::tick("ESolver_KS_LIP", "hamilt2density");

        if (this->phsol != nullptr)
        {
            // reset energy
            this->pelec->f_en.eband = 0.0;
            this->pelec->f_en.demet = 0.0;
            // choose if psi should be diag in subspace
            // be careful that istep start from 0 and iter start from 1
            // if (iter == 1)
            hsolver::DiagoIterAssist<T>::need_subspace = ((istep == 0 || istep == 1) && iter == 1) ? false : true;
            hsolver::DiagoIterAssist<T>::SCF_ITER = iter;
            hsolver::DiagoIterAssist<T>::PW_DIAG_THR = ethr;
            hsolver::DiagoIterAssist<T>::PW_DIAG_NMAX = GlobalV::PW_DIAG_NMAX;

            // It is not a good choice to overload another solve function here, this will spoil the concept of
            // multiple inheritance and polymorphism. But for now, we just do it in this way.
            // In the future, there will be a series of class ESolver_KS_LCAO_PW, HSolver_LCAO_PW and so on.
            std::weak_ptr<psi::Psi<T>> psig = this->p_wf_init->get_psig();

            if (psig.expired())
            {
                ModuleBase::WARNING_QUIT("ESolver_KS_PW::hamilt2density", "psig lifetime is expired");
            }

            // from HSolverLIP
            this->phsol->solve(this->p_hamilt,        // hamilt::Hamilt<T>* pHamilt,
                this->kspw_psi[0],     // psi::Psi<T>& psi,
                this->pelec,           // elecstate::ElecState<T>* pelec,
                psig.lock().get()[0]); // psi::Psi<T>& transform,

            if (GlobalV::out_bandgap)
            {
                if (!GlobalV::TWO_EFERMI)
                {
                    this->pelec->cal_bandgap();
                }
                else
                {
                    this->pelec->cal_bandgap_updw();
                }
            }
        }
        else
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LIP", "HSolver has not been allocated.");
        }
        // add exx
#ifdef __EXX
        if (GlobalC::exx_info.info_global.cal_exx)
            this->pelec->set_exx(this->exx_lip->get_exx_energy()); // Peize Lin add 2019-03-09
#endif

        // calculate the delta_harris energy
        // according to new charge density.
        // mohan add 2009-01-23
        this->pelec->cal_energies(1);

        Symmetry_rho srho;
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            srho.begin(is, *(this->pelec->charge), this->pw_rhod, GlobalC::Pgrid, GlobalC::ucell.symm);
        }

        // compute magnetization, only for LSDA(spin==2)
        GlobalC::ucell.magnet.compute_magnetization(this->pelec->charge->nrxx,
            this->pelec->charge->nxyz,
            this->pelec->charge->rho,
            this->pelec->nelec_spin.data());

        // deband is calculated from "output" charge density calculated
        // in sum_band
        // need 'rho(out)' and 'vr (v_h(in) and v_xc(in))'
        this->pelec->f_en.deband = this->pelec->cal_delta_eband();

        ModuleBase::timer::tick("ESolver_KS_LIP", "hamilt2density");
    }

#ifdef __EXX
    template <typename T>
    bool ESolver_KS_LIP<T>::do_after_converge(int& iter)
    {
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            // no separate_loop case
            if (!GlobalC::exx_info.info_global.separate_loop)
            {
                GlobalC::exx_info.info_global.hybrid_step = 1;

                // in no_separate_loop case, scf loop only did twice
                // in first scf loop, exx updated once in beginning,
                // in second scf loop, exx updated every iter

                if (this->two_level_step)
                    return true;
                else
                {
                    // update exx and redo scf
                    XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].ncpp.xc_func);
                    iter = 0;
                    std::cout << " Entering 2nd SCF, where EXX is updated" << std::endl;
                    this->two_level_step++;
                    return false;
                }
            }
            // has separate_loop case
            // exx converged or get max exx steps
            else if (this->two_level_step == GlobalC::exx_info.info_global.hybrid_step
                || (iter == 1 && this->two_level_step != 0))
                return true;
            else
            {
                // update exx and redo scf
                if (this->two_level_step == 0)
                {
                    XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].ncpp.xc_func);
                }

                std::cout << " Updating EXX " << std::flush;
                timeval t_start;       gettimeofday(&t_start, NULL);

                this->exx_lip->cal_exx();
                iter = 0;
                this->two_level_step++;

                timeval t_end;       gettimeofday(&t_end, NULL);
                std::cout << "and rerun SCF\t"
                    << std::setprecision(3) << std::setiosflags(std::ios::scientific)
                    << (double)(t_end.tv_sec - t_start.tv_sec) + (double)(t_end.tv_usec - t_start.tv_usec) / 1000000.0
                    << std::defaultfloat << " (s)" << std::endl;
                return false;
            }
        }
        else { return true; }
    }
#endif

    template <typename T>
    void ESolver_KS_LIP<T>::after_all_runners()
    {
        ESolver_KS_PW<T>::after_all_runners();
#ifdef __LCAO
        if (PARAM.inp.out_mat_xc)
        {
            ModuleIO::write_Vxc(GlobalV::NSPIN, GlobalV::NLOCAL,
                GlobalV::DRANK, *this->kspw_psi, GlobalC::ucell, this->sf,
                *this->pw_wfc, *this->pw_rho, *this->pw_rhod,
                GlobalC::ppcell.vloc, *this->pelec->charge, this->kv, this->pelec->wg
#ifdef __EXX
                , *this->exx_lip
#endif
            );
        }
#endif
    }
    template class ESolver_KS_LIP<std::complex<float>>;
    template class ESolver_KS_LIP<std::complex<double>>;
    // LIP is not supported on GPU yet.
} // namespace ModuleESolver
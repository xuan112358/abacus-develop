#include "esolver_of.h"

#include "module_parameter/parameter.h"
#include "module_io/cube_io.h"
#include "module_io/output_log.h"
#include "module_io/write_elecstat_pot.h"
//-----------temporary-------------------------
#include "module_base/global_function.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"
#include "module_elecstate/cal_ux.h"
//-----force-------------------
#include "module_hamilt_pw/hamilt_pwdft/forces.h"
//-----stress------------------
#include "module_hamilt_pw/hamilt_ofdft/of_stress_pw.h"

namespace ModuleESolver
{

ESolver_OF::ESolver_OF()
{
    this->classname = "ESolver_OF";
    this->task_ = new char[60];
}

ESolver_OF::~ESolver_OF()
{
    delete psi_;
    delete[] this->pphi_;

    for (int i = 0; i < PARAM.inp.nspin; ++i)
    {
        delete[] this->pdirect_[i];
        delete[] this->pdLdphi_[i];
        delete[] this->pdEdphi_[i];
        delete[] this->precip_dir_[i];
    }
    delete[] this->pdirect_;
    delete[] this->pdLdphi_;
    delete[] this->pdEdphi_;
    delete[] this->precip_dir_;

    delete[] this->nelec_;
    delete[] this->theta_;
    delete[] this->task_;
    delete this->ptemp_rho_;

    delete this->tf_;
    delete this->vw_;
    delete this->wt_;
    delete this->lkt_;

    delete this->opt_cg_;
    delete this->opt_tn_;
    delete this->opt_dcsrch_;
    delete this->opt_cg_mag_;
}

void ESolver_OF::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    ESolver_FP::before_all_runners(ucell, inp);

    // save necessary parameters
    this->of_kinetic_ = inp.of_kinetic;
    this->of_method_ = inp.of_method;
    this->of_conv_ = inp.of_conv;
    this->of_tole_ = inp.of_tole;
    this->of_tolp_ = inp.of_tolp;
    this->max_iter_ = inp.scf_nmax;
    this->dV_ = ucell.omega / this->pw_rho->nxyz;
    this->bound_cal_potential_
        = std::bind(&ESolver_OF::cal_potential, this, std::placeholders::_1, std::placeholders::_2, std::ref(ucell));

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
    int func_type = XC_Functional::get_func_type();
    if (func_type > 2)
    {
        ModuleBase::WARNING_QUIT("esolver_of", "meta-GGA and Hybrid functionals are not supported by OFDFT.");
    }

    // symmetry analysis should be performed every time the cell is changed
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        ucell.symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    // Setup the k points according to symmetry.
    kv.set(ucell,ucell.symm, PARAM.inp.kpoint_file, PARAM.inp.nspin, ucell.G, ucell.latvec, GlobalV::ofs_running);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

    // print information
    // mohan add 2021-01-30
    ModuleIO::setup_parameters(ucell, kv);

    // initialize the real-space uniform grid for FFT and parallel
    // distribution of plane waves
    Pgrid.init(pw_rho->nx,
                        pw_rho->ny,
                        pw_rho->nz,
                        pw_rho->nplane,
                        pw_rho->nrxx,
                        pw_big->nbz,
                        pw_big->bz); // mohan add 2010-07-22, update 2011-05-04
    // Calculate Structure factor
    sf.setup_structure_factor(&ucell, Pgrid, pw_rho);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT BASIS");

    // initialize local pseudopotential
    this->locpp.init_vloc(ucell,pw_rho);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");


    // initialize elecstate, including potential
    this->init_elecstate(ucell);

    // calculate the total local pseudopotential in real space
    this->pelec->init_scf(0, ucell, Pgrid, sf.strucFac, locpp.numeric, ucell.symm); // atomic_rho, v_of_rho, set_vrs

    // liuyu move here 2023-10-09
    // D in uspp need vloc, thus behind init_scf()
    // calculate the effective coefficient matrix for non-local pseudopotential projectors
    ModuleBase::matrix veff = this->pelec->pot->get_effective_v();

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT POTENTIAL");

    // Initialize KEDF
    // Calculate electron numbers, which will be used to initialize WT KEDF
    this->nelec_ = new double[PARAM.inp.nspin];
    if (PARAM.inp.nspin == 1)
    {
        this->nelec_[0] = PARAM.inp.nelec;
    }
    else if (PARAM.inp.nspin == 2)
    {
        // in fact, nelec_spin will not be used anymore
        this->pelec->init_nelec_spin();
        this->nelec_[0] = this->pelec->nelec_spin[0];
        this->nelec_[1] = this->pelec->nelec_spin[1];
    }
    this->init_kedf(inp);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT KEDF");

    // Initialize optimization methods
    this->init_opt();
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT OPTIMIZATION");

    this->allocate_array();
}

void ESolver_OF::runner(UnitCell& ucell, const int istep)
{
    ModuleBase::timer::tick("ESolver_OF", "runner");
    // get Ewald energy, initial rho and phi if necessary
    this->before_opt(istep, ucell);
    this->iter_ = 0;

    while (true)
    {
        // once we get a new rho and phi, update potential
        this->update_potential(ucell);

        // calculate the energy of new rho and phi
        this->energy_llast_ = this->energy_last_;
        this->energy_last_ = this->energy_current_;
        this->energy_current_ = this->cal_energy();

        // check if the job is done
        if (this->check_exit())
        {
            break;
        }

        // find the optimization direction and step lenghth theta according to the potential
        this->optimize(ucell);

        // update the rho and phi based on the direction and theta
        this->update_rho();

        this->iter_++;
    }

    this->after_opt(istep, ucell);

    ModuleBase::timer::tick("ESolver_OF", "runner");
}

/**
 * @brief Prepare to optimize the charge density,
 * update elecstate, kedf, and opts if needed
 * calculate ewald energy, initialize the rho, phi, theta
 *
 * @param istep
 * @param ucell
 */
void ESolver_OF::before_opt(const int istep, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_OF", "before_opt");
    ModuleBase::timer::tick("ESolver_OF", "before_opt");

    //! 1) call before_scf() of ESolver_FP
    ESolver_FP::before_scf(ucell, istep);

    if (ucell.cell_parameter_updated)
    {
        this->dV_ = ucell.omega / this->pw_rho->nxyz;

        // initialize elecstate, including potential
        this->init_elecstate(ucell);

        // Initialize KEDF
        this->init_kedf(PARAM.inp);

        // Initialize optimization methods
        this->init_opt();

        // Refresh the arrays
        delete this->psi_;
        this->psi_ = new psi::Psi<double>(1, PARAM.inp.nspin, this->pw_rho->nrxx);
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            this->pphi_[is] = this->psi_->get_pointer(is);
        }

        delete this->ptemp_rho_;
        this->ptemp_rho_ = new Charge();
        this->ptemp_rho_->set_rhopw(this->pw_rho);
        this->ptemp_rho_->allocate(PARAM.inp.nspin);

        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            delete[] this->pdLdphi_[is];
            delete[] this->pdEdphi_[is];
            delete[] this->pdirect_[is];
            delete[] this->precip_dir_[is];
            this->pdLdphi_[is] = new double[this->pw_rho->nrxx];
            this->pdEdphi_[is] = new double[this->pw_rho->nrxx];
            this->pdirect_[is] = new double[this->pw_rho->nrxx];
            this->precip_dir_[is] = new std::complex<double>[pw_rho->npw];
        }
    }
    if (ucell.ionic_position_updated)
    {
        CE.update_all_dis(ucell);
        CE.extrapolate_charge(&Pgrid, ucell, pelec->charge, &sf, GlobalV::ofs_running, GlobalV::ofs_warning);
    }

    this->pelec->init_scf(istep, ucell, Pgrid, sf.strucFac, locpp.numeric, ucell.symm);

    // calculate ewald energy
    this->pelec->f_en.ewald_energy = H_Ewald_pw::compute_ewald(ucell, this->pw_rho, sf.strucFac);

    Symmetry_rho srho;
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        srho.begin(is, *(pelec->charge), this->pw_rho, ucell.symm);
    }

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        if (PARAM.inp.init_chg != "file")
        {
            for (int ibs = 0; ibs < this->pw_rho->nrxx; ++ibs)
            {
                // Here we initialize rho to be uniform,
                // because the rho got by pot.init_pot -> Charge::atomic_rho may contain minus elements.
                pelec->charge->rho[is][ibs] = this->nelec_[is] / this->pelec->omega;
                this->pphi_[is][ibs] = sqrt(pelec->charge->rho[is][ibs]);
            }
        }
        else
        {
            for (int ibs = 0; ibs < this->pw_rho->nrxx; ++ibs)
            {
                this->pphi_[is][ibs] = sqrt(pelec->charge->rho[is][ibs]);
            }
        }
    }

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        this->pelec->eferm.get_ef(is) = 0.;
        this->theta_[is] = 0.;
        ModuleBase::GlobalFunc::ZEROS(this->pdLdphi_[is], this->pw_rho->nrxx);
        ModuleBase::GlobalFunc::ZEROS(this->pdEdphi_[is], this->pw_rho->nrxx);
        ModuleBase::GlobalFunc::ZEROS(this->pdirect_[is], this->pw_rho->nrxx);
    }
    if (PARAM.inp.nspin == 1)
    {
        this->theta_[0] = 0.2;
    }

    ModuleBase::timer::tick("ESolver_OF", "before_opt");
}

/**
 * @brief Get dL/dphi = dL/drho * drho/dphi = (dE/drho - mu) * 2 * phi,
 * as well as normdLdphi = sqrt{<dL/dphi|dL/dphi>}
 *
 * @param ucell
 */
void ESolver_OF::update_potential(UnitCell& ucell)
{
    // (1) get dL/dphi
    elecstate::cal_ux(ucell);

    this->pelec->pot->update_from_charge(pelec->charge, &ucell); // Hartree + XC + external
    this->kinetic_potential(pelec->charge->rho,
                            this->pphi_,
                            this->pelec->pot->get_effective_v()); // (kinetic + Hartree + XC + external) * 2 * phi
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        const double* vr_eff = this->pelec->pot->get_effective_v(is);
        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            this->pdEdphi_[is][ir] = vr_eff[ir];
        }
        this->pelec->eferm.get_ef(is) = this->cal_mu(this->pphi_[is], this->pdEdphi_[is], this->nelec_[is]);

        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            this->pdLdphi_[is][ir]
                = this->pdEdphi_[is][ir] - 2. * this->pelec->eferm.get_efval(is) * this->pphi_[is][ir];
        }
    }

    // (2) get the norm of dLdphi
    // ===== temporary solution of potential convergence when of_full_pw = 0 =====
    this->normdLdphi_llast_ = this->normdLdphi_last_;
    this->normdLdphi_last_ = this->normdLdphi_;
    // ===========================================================================
    this->normdLdphi_ = 0.;

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        this->normdLdphi_ += this->inner_product(this->pdLdphi_[is], this->pdLdphi_[is], this->pw_rho->nrxx, 1.0);
    }
    Parallel_Reduce::reduce_all(this->normdLdphi_);
    this->normdLdphi_ = sqrt(this->normdLdphi_ / this->pw_rho->nxyz / PARAM.inp.nspin);
}

/**
 * @brief Get the optimization direction (this->pdirection_) and the step length (this->theta)
 *
 * @param ucell
 */
void ESolver_OF::optimize(UnitCell& ucell)
{
    // (1) get |d0> with optimization algorithm
    this->get_direction(ucell);
    // initialize temp_phi and temp_rho used in line search
    double** ptemp_phi = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        ptemp_phi[is] = new double[this->pw_rho->nrxx];
        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            ptemp_phi[is][ir] = this->pphi_[is][ir];
            this->ptemp_rho_->rho[is][ir] = ptemp_phi[is][ir] * ptemp_phi[is][ir];
        }
    }

    // (2) rotate and renormalize the direction
    this->adjust_direction();

    // (3) make sure that dEdtheta<0 at theta = 0
    double* dEdtheta = new double[PARAM.inp.nspin]; // dE/dtheta of tempPhi
    ModuleBase::GlobalFunc::ZEROS(dEdtheta, PARAM.inp.nspin);

    this->check_direction(dEdtheta, ptemp_phi, ucell);
    // this->test_direction(dEdtheta, ptemp_phi, ucell);

    // (4) call line search to find the best theta (step length)
    this->get_step_length(dEdtheta, ptemp_phi, ucell);

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] ptemp_phi[is];
    }
    delete[] ptemp_phi;
    delete[] dEdtheta;
}

/**
 * @brief Update the charge density and "wavefunction" (phi) after one step of optimization
 * phi = cos(theta) * phi + sin(theta) * direction,
 * rho = phi^2
 */
void ESolver_OF::update_rho()
{
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            this->pphi_[is][ir]
                = this->pphi_[is][ir] * cos(this->theta_[is]) + this->pdirect_[is][ir] * sin(this->theta_[is]);
            pelec->charge->rho[is][ir] = this->pphi_[is][ir] * this->pphi_[is][ir];
        }
    }
    // // ------------ turn on symmetry may cause instability in optimization ------------
    // if (ModuleSymmetry::Symmetry::symm_flag == 1)
    // {
    //     Symmetry_rho srho;
    //     for (int is = 0; is < PARAM.inp.nspin; is++)
    //     {
    //         srho.begin(is, *(pelec->charge), this->pw_rho, Pgrid, ucell.symm);
    //         for (int ibs = 0; ibs < this->pw_rho->nrxx; ++ibs)
    //         {
    //             this->pphi_[is][ibs] = sqrt(pelec->charge->rho[is][ibs]);
    //         }
    //     }
    // }
    // // --------------------------------------------------------------------------------
}

/**
 * @brief Check convergence, return ture if converge or iter >= max_iter_,
 * and print the necessary information
 *
 * @return exit or not
 */
bool ESolver_OF::check_exit()
{
    this->conv_esolver = false;
    bool potConv = false;
    bool potHold = false; // if normdLdphi nearly remains unchanged
    bool energyConv = false;

    if (this->normdLdphi_ < this->of_tolp_)
    {
        potConv = true;
    }
    if (this->iter_ >= 3 && std::abs(this->normdLdphi_ - this->normdLdphi_last_) < 1e-10
        && std::abs(this->normdLdphi_ - this->normdLdphi_llast_) < 1e-10)
    {
        potHold = true;
    }

    if (this->iter_ >= 3 && std::abs(this->energy_current_ - this->energy_last_) < this->of_tole_
        && std::abs(this->energy_current_ - this->energy_llast_) < this->of_tole_)
    {
        energyConv = true;
    }

    this->conv_esolver = (this->of_conv_ == "energy" && energyConv) || (this->of_conv_ == "potential" && potConv)
                         || (this->of_conv_ == "both" && potConv && energyConv);

    this->print_info();

    if (this->conv_esolver || this->iter_ >= this->max_iter_)
    {
        return true;
    }
    // ============ temporary solution of potential convergence ===========
    else if (this->of_conv_ == "potential" && potHold)
    {
        GlobalV::ofs_warning << "ESolver_OF WARNING: "
                             << "The convergence of potential has not been reached, but the norm of potential nearly "
                                "remains unchanged, set of_full_pw = 1 may work."
                             << std::endl;
        return true;
    }
    // ====================================================================
    else
    {
        return false;
    }
}

/**
 * @brief After optimization, output the charge density, effective potential, ..., if needed.
 *
 * @param istep
 * @param ucell
 */
void ESolver_OF::after_opt(const int istep, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_OF", "after_opt");
    ModuleBase::timer::tick("ESolver_OF", "after_opt");

    // 1) calculate the kinetic energy density
    if (PARAM.inp.out_elf[0] > 0)
    {
        this->kinetic_energy_density(this->pelec->charge->rho, this->pphi_, this->pelec->charge->kin_r);
    }

    // 2) call after_scf() of ESolver_FP
    ESolver_FP::after_scf(ucell, istep);

    ModuleBase::timer::tick("ESolver_OF", "after_opt");
}

/**
 * @brief Output the FINAL_ETOT
 */
void ESolver_OF::after_all_runners(UnitCell& ucell)
{

    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << this->pelec->f_en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;
}

/**
 * @brief Calculate the total energy.
 * NOTE THIS FUNCTION SHOULD BE CALLEDD AFTER POTENTIAL HAS BEEN UPDATED
 *
 * @return total energy
 */
double ESolver_OF::cal_energy()
{
    this->pelec->cal_energies(2);
    double kinetic_energy = this->kinetic_energy(); // kinetic energy
    double pseudopot_energy = 0.;                   // electron-ion interaction energy
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        pseudopot_energy += this->inner_product(this->pelec->pot->get_fixed_v(),
                                                pelec->charge->rho[is],
                                                this->pw_rho->nrxx,
                                                this->dV_);
    }
    Parallel_Reduce::reduce_all(pseudopot_energy);
    this->pelec->f_en.ekinetic = kinetic_energy;
    this->pelec->f_en.eion_elec = pseudopot_energy;
    this->pelec->f_en.etot += kinetic_energy + pseudopot_energy;
    return this->pelec->f_en.etot;
}

/**
 * @brief Calculate the force
 *
 * @param [out] force
 */
void ESolver_OF::cal_force(UnitCell& ucell, ModuleBase::matrix& force)
{
    Forces<double> ff(ucell.nat);
    ff.cal_force(ucell, force, *pelec, this->pw_rho, &ucell.symm, &sf, this->solvent, &this->locpp);
}

/**
 * @brief Calculate the stress
 *
 * @param [out] stress
 */
void ESolver_OF::cal_stress(UnitCell& ucell, ModuleBase::matrix& stress)
{
    ModuleBase::matrix kinetic_stress_;
    kinetic_stress_.create(3, 3);
    this->kinetic_stress(kinetic_stress_);

    OF_Stress_PW ss(this->pelec, this->pw_rho);
    ss.cal_stress(stress, kinetic_stress_, ucell, &ucell.symm, this->locpp, &sf, &kv);
}
} // namespace ModuleESolver

#include "dftu.h"
#include "module_base/timer.h"
#include "module_parameter/parameter.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#ifdef __LCAO
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#endif

extern "C"
{
  //I'm not sure what's happenig here, but the interface in scalapack_connecter.h
  //does not seem to work, so I'll use this one here
  void pzgemm_(
		const char *transa, const char *transb,
		const int *M, const int *N, const int *K,
		const std::complex<double> *alpha,
		const std::complex<double> *A, const int *IA, const int *JA, const int *DESCA,
		const std::complex<double> *B, const int *IB, const int *JB, const int *DESCB,
		const std::complex<double> *beta,
		std::complex<double> *C, const int *IC, const int *JC, const int *DESCC);
  
  void pdgemm_(
		const char *transa, const char *transb,
		const int *M, const int *N, const int *K,
		const double *alpha,
		const double *A, const int *IA, const int *JA, const int *DESCA,
		const double *B, const int *IB, const int *JB, const int *DESCB,
		const double *beta,
		double *C, const int *IC, const int *JC, const int *DESCC);
}

namespace ModuleDFTU
{
void DFTU::copy_locale(const UnitCell& ucell)
{
    ModuleBase::TITLE("DFTU", "copy_locale");
    ModuleBase::timer::tick("DFTU", "copy_locale");

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (orbital_corr[T] == -1) {
            continue;
}

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);

            for (int l = 0; l < ucell.atoms[T].nwl + 1; l++)
            {
                const int N = ucell.atoms[T].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    if (PARAM.inp.nspin == 4)
                    {
                        locale_save[iat][l][n][0] = locale[iat][l][n][0];
                    }
                    else if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 2)
                    {
                        locale_save[iat][l][n][0] = locale[iat][l][n][0];
                        locale_save[iat][l][n][1] = locale[iat][l][n][1];
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("DFTU", "copy_locale");
}

void DFTU::zero_locale(const UnitCell& ucell)
{
    ModuleBase::TITLE("DFTU", "zero_locale");
    ModuleBase::timer::tick("DFTU", "zero_locale");

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (orbital_corr[T] == -1) { continue;
}

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);

            for (int l = 0; l < ucell.atoms[T].nwl + 1; l++)
            {
                const int N = ucell.atoms[T].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    if (PARAM.inp.nspin == 4)
                    {
                        locale[iat][l][n][0].zero_out();
                    }
                    else if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 2)
                    {
                        locale[iat][l][n][0].zero_out();
                        locale[iat][l][n][1].zero_out();
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("DFTU", "zero_locale");
}

void DFTU::mix_locale(const UnitCell& ucell,
                      const double& mixing_beta)
{
    ModuleBase::TITLE("DFTU", "mix_locale");
    ModuleBase::timer::tick("DFTU", "mix_locale");

    double beta = mixing_beta;

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (orbital_corr[T] == -1) {
            continue;
}

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);

            for (int l = 0; l < ucell.atoms[T].nwl + 1; l++)
            {
                const int N = ucell.atoms[T].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    if (PARAM.inp.nspin == 4)
                    {
                        locale[iat][l][n][0] = locale[iat][l][n][0]*beta + locale_save[iat][l][n][0]*(1.0-beta);
                    }
                    else if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 2)
                    {
                        locale[iat][l][n][0] = locale[iat][l][n][0] * beta + locale_save[iat][l][n][0] * (1.0-beta);
                        locale[iat][l][n][1] = locale[iat][l][n][1] * beta + locale_save[iat][l][n][1] * (1.0-beta);
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("DFTU", "mix_locale");
}

#ifdef __LCAO

void DFTU::cal_occup_m_k(const int iter, 
                         const UnitCell& ucell,
                         const std::vector<std::vector<std::complex<double>>>& dm_k,
                         const K_Vectors& kv,
                         const double& mixing_beta,
                         hamilt::Hamilt<std::complex<double>>* p_ham)
{
    ModuleBase::TITLE("DFTU", "cal_occup_m_k");
    ModuleBase::timer::tick("DFTU", "cal_occup_m_k");

    this->copy_locale(ucell);
    this->zero_locale(ucell);

    //=================Part 1======================
    // call SCALAPACK routine to calculate the product of the S and density matrix
    const char transN = 'N', transT = 'T';
    const int one_int = 1;
    const std::complex<double> beta(0.0,0.0), alpha(1.0,0.0);

    std::vector<std::complex<double>> srho(this->paraV->nloc);

    for (int ik = 0; ik < kv.get_nks(); ik++)
    {
        // srho(mu,nu) = \sum_{iw} S(mu,iw)*dm_k(iw,nu)
        this->folding_matrix_k_new(ik, p_ham);
        std::complex<double>* s_k_pointer = nullptr;
        if(PARAM.inp.nspin != 4)
        {
            s_k_pointer = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham)->getSk();
        }
        else
        {
            s_k_pointer = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham)->getSk();
        }

#ifdef __MPI
        pzgemm_(&transN,
                &transT,
                &PARAM.globalv.nlocal,
                &PARAM.globalv.nlocal,
                &PARAM.globalv.nlocal,
                &alpha,
                s_k_pointer,
                &one_int,
                &one_int,
                this->paraV->desc,
                dm_k[ik].data(),
                //dm_k[ik].c,
                &one_int,
                &one_int,
                this->paraV->desc,
                &beta,
                &srho[0],
                &one_int,
                &one_int,
                this->paraV->desc);
#endif

        const int spin = kv.isk[ik];
        for (int it = 0; it < ucell.ntype; it++)
        {
            const int NL = ucell.atoms[it].nwl + 1;
            const int LC = orbital_corr[it];

            if (LC == -1) {
                continue;
}

            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                const int iat = ucell.itia2iat(it, ia);

                for (int l = 0; l < NL; l++)
                {
                    if (l != orbital_corr[it]) {
                        continue;
}

                    const int N = ucell.atoms[it].l_nchi[l];

                    for (int n = 0; n < N; n++)
                    {
                        // if(!Yukawa && n!=0) continue;
                        if (n != 0) {
                            continue;
}

                        // Calculate the local occupation number matrix
                        for (int m0 = 0; m0 < 2 * l + 1; m0++)
                        {
                            for (int ipol0 = 0; ipol0 < PARAM.globalv.npol; ipol0++)
                            {
                                const int iwt0 = this->iatlnmipol2iwt[iat][l][n][m0][ipol0];
                                const int mu = this->paraV->global2local_row(iwt0);
                                const int mu_prime = this->paraV->global2local_col(iwt0);

                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    for (int ipol1 = 0; ipol1 < PARAM.globalv.npol; ipol1++)
                                    {
                                        const int iwt1 = this->iatlnmipol2iwt[iat][l][n][m1][ipol1];
                                        const int nu = this->paraV->global2local_col(iwt1);
                                        const int nu_prime = this->paraV->global2local_row(iwt1);

                                        const int irc = nu * this->paraV->nrow + mu;
                                        const int irc_prime = mu_prime * this->paraV->nrow + nu_prime;

                                        const int m0_all = m0 + ipol0 * (2 * l + 1);
                                        const int m1_all = m1 + ipol1 * (2 * l + 1);

                                        if ((nu >= 0) && (mu >= 0))
                                            locale[iat][l][n][spin](m0_all, m1_all) += (srho[irc]).real() / 4.0;

                                        if ((nu_prime >= 0) && (mu_prime >= 0))
                                            locale[iat][l][n][spin](m0_all, m1_all)
                                                += (std::conj(srho[irc_prime])).real() / 4.0;
                                    } // ipol1
                                } // m1
                            } // ipol0
                        } // m0
                    } // end n
                } // end l
            } // end ia
        } // end it
    } // ik

    for (int it = 0; it < ucell.ntype; it++)
    {
        const int NL = ucell.atoms[it].nwl + 1;
        const int LC = orbital_corr[it];

        if (LC == -1) {
            continue;
}

        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const int iat = ucell.itia2iat(it, ia);

            for (int l = 0; l < NL; l++)
            {
                if (l != orbital_corr[it]) {
                    continue;
}

                const int N = ucell.atoms[it].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    // if(!Yukawa && n!=0) continue;
                    if (n != 0) {
                        continue;
}
                        // set the local occupation mumber matrix of spin up and down zeros

#ifdef __MPI
                    if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 4)
                    {
                        ModuleBase::matrix temp(locale[iat][l][n][0]);
                        MPI_Allreduce(&temp(0, 0),
                                      &locale[iat][l][n][0](0, 0),
                                      (2 * l + 1) * PARAM.globalv.npol * (2 * l + 1) * PARAM.globalv.npol,
                                      MPI_DOUBLE,
                                      MPI_SUM,
                                      MPI_COMM_WORLD);
                    }
                    else if (PARAM.inp.nspin == 2)
                    {
                        ModuleBase::matrix temp0(locale[iat][l][n][0]);
                        MPI_Allreduce(&temp0(0, 0),
                                      &locale[iat][l][n][0](0, 0),
                                      (2 * l + 1) * (2 * l + 1),
                                      MPI_DOUBLE,
                                      MPI_SUM,
                                      MPI_COMM_WORLD);

                        ModuleBase::matrix temp1(locale[iat][l][n][1]);
                        MPI_Allreduce(&temp1(0, 0),
                                      &locale[iat][l][n][1](0, 0),
                                      (2 * l + 1) * (2 * l + 1),
                                      MPI_DOUBLE,
                                      MPI_SUM,
                                      MPI_COMM_WORLD);
                    }
#endif

                    // for the case spin independent calculation
                    switch (PARAM.inp.nspin)
                    {
                    case 1:
                        locale[iat][l][n][0] += transpose(locale[iat][l][n][0]);
                        locale[iat][l][n][0] *= 0.5;
                        locale[iat][l][n][1] += locale[iat][l][n][0];
                        break;

                    case 2:
                        for (int is = 0; is < PARAM.inp.nspin; is++)
                            locale[iat][l][n][is] += transpose(locale[iat][l][n][is]);
                        break;

                    case 4: // SOC
                        locale[iat][l][n][0] += transpose(locale[iat][l][n][0]);
                        break;

                    default:
                        std::cout << "Not supported NSPIN parameter" << std::endl;
                        exit(0);
                    }
                } // end n
            } // end l
        } // end ia
    } // end it

    if(mixing_dftu && initialed_locale)
    {
        this->mix_locale(ucell,mixing_beta);
    }

    this->initialed_locale = true;
    ModuleBase::timer::tick("DFTU", "cal_occup_m_k");
    return;
}

void DFTU::cal_occup_m_gamma(const int iter,
                             const UnitCell &ucell,
                             const std::vector<std::vector<double>> &dm_gamma, 
                             const double& mixing_beta, 
                             hamilt::Hamilt<double>* p_ham)
{
    ModuleBase::TITLE("DFTU", "cal_occup_m_gamma");
    ModuleBase::timer::tick("DFTU", "cal_occup_m_gamma");
    this->copy_locale(ucell);
    this->zero_locale(ucell);

    //=================Part 1======================
    // call PBLAS routine to calculate the product of the S and density matrix
    const char transN = 'N', transT = 'T';
    const int one_int = 1;
    const double alpha = 1.0, beta = 0.0;

    std::vector<double> srho(this->paraV->nloc);
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        // srho(mu,nu) = \sum_{iw} S(mu,iw)*dm_gamma(iw,nu)
        double* s_gamma_pointer = dynamic_cast<hamilt::HamiltLCAO<double, double>*>(p_ham)->getSk();

#ifdef __MPI
        pdgemm_(&transN,
                &transT,
                &PARAM.globalv.nlocal,
                &PARAM.globalv.nlocal,
                &PARAM.globalv.nlocal,
                &alpha,
                s_gamma_pointer,
                &one_int,
                &one_int,
                this->paraV->desc,
                dm_gamma[is].data(),
                //dm_gamma[is].c,
                &one_int,
                &one_int,
                this->paraV->desc,
                &beta,
                &srho[0],
                &one_int,
                &one_int,
                this->paraV->desc);
#endif

        for (int it = 0; it < ucell.ntype; it++)
        {
            const int NL = ucell.atoms[it].nwl + 1;
            if (orbital_corr[it] == -1) {
                continue;
}
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                const int iat = ucell.itia2iat(it, ia);

                for (int l = 0; l < NL; l++)
                {
                    if (l != orbital_corr[it]) {
                        continue;
}

                    const int N = ucell.atoms[it].l_nchi[l];

                    for (int n = 0; n < N; n++)
                    {
                        if (n != 0) {
                            continue;
}

                        // Calculate the local occupation number matrix
                        for (int m0 = 0; m0 < 2 * l + 1; m0++)
                        {
                            for (int ipol0 = 0; ipol0 < PARAM.globalv.npol; ipol0++)
                            {
                                const int iwt0 = this->iatlnmipol2iwt[iat][l][n][m0][ipol0];
                                const int mu = this->paraV->global2local_row(iwt0);
                                const int mu_prime = this->paraV->global2local_col(iwt0);

                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    for (int ipol1 = 0; ipol1 < PARAM.globalv.npol; ipol1++)
                                    {
                                        const int iwt1 = this->iatlnmipol2iwt[iat][l][n][m1][ipol1];
                                        const int nu = this->paraV->global2local_col(iwt1);
                                        const int nu_prime = this->paraV->global2local_row(iwt1);

                                        const int irc = nu * this->paraV->nrow + mu;
                                        const int irc_prime = mu_prime * this->paraV->nrow + nu_prime;

                                        if ((nu >= 0) && (mu >= 0))
                                        {
                                            int m0_all = m0 + (2 * l + 1) * ipol0;
                                            int m1_all = m0 + (2 * l + 1) * ipol1;

                                            locale[iat][l][n][is](m0, m1) += srho[irc] / 4.0;
                                        }

                                        if ((nu_prime >= 0) && (mu_prime >= 0))
                                        {
                                            int m0_all = m0 + (2 * l + 1) * ipol0;
                                            int m1_all = m0 + (2 * l + 1) * ipol1;

                                            locale[iat][l][n][is](m0, m1) += srho[irc_prime] / 4.0;
                                        }
                                    }
                                }
                            }
                        }

                        ModuleBase::matrix temp(locale[iat][l][n][is]);

#ifdef __MPI
                        MPI_Allreduce(&temp(0, 0),
                                      &locale[iat][l][n][is](0, 0),
                                      (2 * l + 1) * PARAM.globalv.npol * (2 * l + 1) * PARAM.globalv.npol,
                                      MPI_DOUBLE,
                                      MPI_SUM,
                                      MPI_COMM_WORLD);
#endif

                        // for the case spin independent calculation
                        switch (PARAM.inp.nspin)
                        {
                        case 1:
                            locale[iat][l][n][0] += transpose(locale[iat][l][n][0]);
                            locale[iat][l][n][0] *= 0.5;
                            locale[iat][l][n][1] += locale[iat][l][n][0];
                            break;

                        case 2:
                            locale[iat][l][n][is] += transpose(locale[iat][l][n][is]);
                            break;

                        default:
                            std::cout << "Not supported NSPIN parameter" << std::endl;
                            exit(0);
                        }

                    } // end for(n)
                } // L
            } // ia
        } // it
    } // is

    if(mixing_dftu && initialed_locale)
    {
        this->mix_locale(ucell,mixing_beta);
    }

    this->initialed_locale = true;
    ModuleBase::timer::tick("DFTU", "cal_occup_m_gamma");
    return;
}
#endif
} // namespace ModuleDFTU
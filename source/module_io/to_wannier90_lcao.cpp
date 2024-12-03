#include "to_wannier90_lcao.h"

#include "module_parameter/parameter.h"
#include "fR_overlap.h"
#include "module_base/math_integral.h"
#include "module_base/math_polyint.h"
#include "module_base/math_sphbes.h"
#include "module_base/math_ylmreal.h"
#include "module_base/parallel_reduce.h"
#include "module_base/scalapack_connector.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#include <fstream>
#include <functional>

#ifdef __LCAO
toWannier90_LCAO::toWannier90_LCAO(const bool& out_wannier_mmn,
                                   const bool& out_wannier_amn,
                                   const bool& out_wannier_unk,
                                   const bool& out_wannier_eig,
                                   const bool& out_wannier_wvfn_formatted,
                                   const std::string& nnkpfile,
                                   const std::string& wannier_spin,
                                   const LCAO_Orbitals& orb)
    : toWannier90(out_wannier_mmn,
                  out_wannier_amn,
                  out_wannier_unk,
                  out_wannier_eig,
                  out_wannier_wvfn_formatted,
                  nnkpfile,
                  wannier_spin),
    orb_(orb)
{
}

toWannier90_LCAO::~toWannier90_LCAO()
{
}

void toWannier90_LCAO::calculate(const ModuleBase::matrix& ekb,
                                 const K_Vectors& kv,
                                 const psi::Psi<std::complex<double>>& psi,
                                 const Parallel_Orbitals* pv)
{
    this->ParaV = pv;

    read_nnkp(kv);

    if (PARAM.inp.nspin == 2)
    {
        if (wannier_spin == "up")
        {
            start_k_index = 0;
        }
        else if (wannier_spin == "down")
        {
            start_k_index = num_kpts / 2;
        }
        else
        {
            ModuleBase::WARNING_QUIT("toWannier90::calculate", "Error wannier_spin set,is not \"up\" or \"down\" ");
        }
    }

    if (out_wannier_mmn || out_wannier_amn)
    {
        iw2it.resize(PARAM.globalv.nlocal);
        iw2ia.resize(PARAM.globalv.nlocal);
        iw2iL.resize(PARAM.globalv.nlocal);
        iw2iN.resize(PARAM.globalv.nlocal);
        iw2im.resize(PARAM.globalv.nlocal);
        iw2iorb.resize(PARAM.globalv.nlocal);

        std::map<size_t, std::map<size_t, std::map<size_t, size_t>>> temp_orb_index;
        int count = 0;
        for (int it = 0; it < GlobalC::ucell.ntype; it++)
        {
            for (int iL = 0; iL < GlobalC::ucell.atoms[it].nwl + 1; iL++)
            {
                for (int iN = 0; iN < GlobalC::ucell.atoms[it].l_nchi[iL]; iN++)
                {
                    temp_orb_index[it][iL][iN] = count;
                    count++;
                }
            }
        }

        int iw = 0;
        for (int it = 0; it < GlobalC::ucell.ntype; it++)
        {
            for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
            {
                for (int iL = 0; iL < GlobalC::ucell.atoms[it].nwl + 1; iL++)
                {
                    for (int iN = 0; iN < GlobalC::ucell.atoms[it].l_nchi[iL]; iN++)
                    {
                        for (int im = 0; im < (2 * iL + 1); im++)
                        {
                            iw2it[iw] = it;
                            iw2ia[iw] = ia;
                            iw2iL[iw] = iL;
                            iw2iN[iw] = iN;
                            iw2im[iw] = im;
                            iw2iorb[iw] = temp_orb_index[it][iL][iN];

                            iw++;
                        }
                    }
                }
            }
        }

        initialize_orb_table();
        produce_basis_orb();
        set_R_coor();
        count_delta_k(kv);
    }

    if (out_wannier_eig)
    {
        out_eig(ekb);
    }

    if (out_wannier_mmn)
    {
        int dk_size = delta_k_all.size();
        this->FR.resize(dk_size);
        std::vector<std::function<std::complex<double>(ModuleBase::Vector3<double>)>> fr_ptr(dk_size);
        for (int i = 0; i < dk_size; i++)
        {
            ModuleBase::Vector3<double> delta_k(delta_k_all[i].x, delta_k_all[i].y, delta_k_all[i].z);

            fr_ptr[i] = [delta_k](ModuleBase::Vector3<double> r) -> std::complex<double> {
                double phase = delta_k * r;
                std::complex<double> exp_idkr = std::exp(-1.0 * ModuleBase::IMAG_UNIT * phase);
                return exp_idkr;
            };

            FR[i].set_parameters(fr_ptr[i], &GlobalC::ucell, &orb_, &GlobalC::GridD, ParaV, 140, 110);
            FR[i].calculate_FR();
        }

        cal_Mmn(kv, psi);
    }

    if (out_wannier_amn)
    {
        cal_Amn(kv, psi);
    }

    if (out_wannier_unk)
    {
        out_unk(psi);
    }
}

void toWannier90_LCAO::cal_Mmn(const K_Vectors& kv, const psi::Psi<std::complex<double>>& psi)
{
    // write .mmn file
    std::ofstream mmn_file;

    if (GlobalV::MY_RANK == 0)
    {
        std::string fileaddress = PARAM.globalv.global_out_dir + wannier_file_name + ".mmn";
        mmn_file.open(fileaddress.c_str(), std::ios::out);

        time_t time_now = time(nullptr);
        mmn_file << " Created on " << ctime(&time_now);
        mmn_file << std::setw(12) << num_bands << std::setw(12) << cal_num_kpts << std::setw(12) << nntot << std::endl;
    }

    for (int ik = 0; ik < cal_num_kpts; ik++)
    {
        for (int ib = 0; ib < nntot; ib++)
        {
            int ikb = nnlist[ik][ib];
            ModuleBase::Vector3<double> phase_G = nncell[ik][ib];
            ModuleBase::ComplexMatrix Mmn;

            int cal_ik = ik + start_k_index;
            int cal_ikb = ikb + start_k_index;
            unkdotkb(kv, psi, cal_ik, cal_ikb, phase_G, Mmn);

            if (GlobalV::MY_RANK == 0)
            {
                mmn_file << std::setw(5) << ik + 1 << std::setw(5) << ikb + 1 << std::setw(5) << int(phase_G.x)
                         << std::setw(5) << int(phase_G.y) << std::setw(5) << int(phase_G.z) << std::endl;

                for (int n = 0; n < num_bands; n++)
                {
                    for (int m = 0; m < num_bands; m++)
                    {
                        mmn_file << std::setw(18) << std::setprecision(12) << std::showpoint << std::fixed
                                 << Mmn(m, n).real() << std::setw(18) << std::setprecision(12) << std::showpoint
                                 << std::fixed
                                 << Mmn(m, n).imag()
                                 // jingan test
                                 // << "    " << std::setw(12) << std::setprecision(9) << std::abs(Mmn(m, n))
                                 << std::endl;
                    }
                }
            }
        }
    }

    if (GlobalV::MY_RANK == 0) {
        mmn_file.close();
}
}

void toWannier90_LCAO::cal_Amn(const K_Vectors& kv, const psi::Psi<std::complex<double>>& psi)
{
    produce_trial_in_lcao();
    construct_overlap_table_project();
    cal_orbA_overlap_R();

    // write .amn file
    std::ofstream Amn_file;

    if (GlobalV::MY_RANK == 0)
    {
        time_t time_now = time(nullptr);
        std::string fileaddress = PARAM.globalv.global_out_dir + wannier_file_name + ".amn";
        Amn_file.open(fileaddress.c_str(), std::ios::out);
        Amn_file << " Created on " << ctime(&time_now);
        Amn_file << std::setw(12) << num_bands << std::setw(12) << cal_num_kpts << std::setw(12) << num_wannier
                 << std::endl;
    }

    for (int ik = start_k_index; ik < (cal_num_kpts + start_k_index); ik++)
    {
        ModuleBase::ComplexMatrix Amn;
        unkdotA(kv, psi, ik, Amn);

        if (GlobalV::MY_RANK == 0)
        {
            for (int iw = 0; iw < num_wannier; iw++)
            {
                for (int ib = 0; ib < num_bands; ib++)
                {
                    Amn_file << std::setw(5) << ib + 1 << std::setw(5) << iw + 1 << std::setw(5)
                             << ik + 1 - start_k_index << std::setw(18) << std::showpoint << std::fixed
                             << std::setprecision(12) << Amn(ib, iw).real() << std::setw(18) << std::showpoint
                             << std::fixed << std::setprecision(12)
                             << Amn(ib, iw).imag()
                             // jingan test
                             //<< "   " << std::setw(18) << std::setprecision(13) << std::abs(Amn(ib, iw))
                             << std::endl;
                }
            }
        }
    }

    if (GlobalV::MY_RANK == 0) {
        Amn_file.close();
}
}

void toWannier90_LCAO::out_unk(const psi::Psi<std::complex<double>>& psi)
{
}

void toWannier90_LCAO::initialize_orb_table()
{
    int Lmax_used = 0;
    int Lmax = 0;
    int exx_lmax = 0;
#ifdef __EXX
    exx_lmax = GlobalC::exx_info.info_ri.abfs_Lmax;
#endif

#ifdef __LCAO
    const int ntype = orb_.get_ntype();
    int lmax_orb = -1, lmax_beta = -1;
    for (int it = 0; it < ntype; it++)
    {
        lmax_orb = std::max(lmax_orb, orb_.Phi[it].getLmax());
        lmax_beta = std::max(lmax_beta, GlobalC::ucell.infoNL.Beta[it].getLmax());
    }
    const double dr = orb_.get_dR();
    const double dk = orb_.get_dk();
    const int kmesh = orb_.get_kmesh() * 4 + 1;
    int Rmesh = static_cast<int>(orb_.get_Rmax() / dr) + 4;
    Rmesh += 1 - Rmesh % 2;

    Center2_Orb::init_Table_Spherical_Bessel(2,
                                             3,
                                             Lmax_used,
                                             Lmax,
                                             exx_lmax,
                                             lmax_orb,
                                             lmax_beta,
                                             dr,
                                             dk,
                                             kmesh,
                                             Rmesh,
                                             psb_);
    ModuleBase::Ylm::set_coefficients();
    MGT.init_Gaunt_CH(Lmax);
    MGT.init_Gaunt(Lmax);
#endif
}

void toWannier90_LCAO::set_R_coor()
{
    int R_minX = int(-GlobalC::GridD.getTrueCellX());
    int R_minY = int(-GlobalC::GridD.getTrueCellY());
    int R_minZ = int(-GlobalC::GridD.getTrueCellZ());

    int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

    int R_num = R_x * R_y * R_z;
    R_coor_car.resize(R_num);

    int count = 0;
    for (int ix = 0; ix < R_x; ix++)
    {
        for (int iy = 0; iy < R_y; iy++)
        {
            for (int iz = 0; iz < R_z; iz++)
            {
                ModuleBase::Vector3<double> tmpR(ix + R_minX, iy + R_minY, iz + R_minZ);
                R_coor_car[count] = tmpR * GlobalC::ucell.latvec;
                count++;
            }
        }
    }
}

void toWannier90_LCAO::count_delta_k(const K_Vectors& kv)
{
    std::set<Coordinate_3D> delta_k_all_tmp;
    for (int ik = 0; ik < cal_num_kpts; ik++)
    {
        for (int ib = 0; ib < nntot; ib++)
        {
            int ikb = nnlist[ik][ib];
            ModuleBase::Vector3<double> G = nncell[ik][ib];

            int cal_ik = ik + start_k_index;
            int cal_ikb = ikb + start_k_index;

            ModuleBase::Vector3<double> ik_car = kv.kvec_c[ik];
            ModuleBase::Vector3<double> ikb_car = kv.kvec_c[ikb] + G * GlobalC::ucell.G;
            Abfs::Vector3_Order<double> dk = (ikb_car - ik_car) * GlobalC::ucell.tpiba;
            Coordinate_3D temp_dk(dk.x, dk.y, dk.z);
            delta_k_all_tmp.insert(temp_dk);
        }
    }

    delta_k_all.resize(delta_k_all_tmp.size());

    int index = 0;
    for (auto& delta_k: delta_k_all_tmp)
    {
        delta_k_all_index[delta_k] = index;
        delta_k_all[index] = delta_k;
        index++;
    }
}

void toWannier90_LCAO::unkdotkb(const K_Vectors& kv,
                                const psi::Psi<std::complex<double>>& psi_in,
                                const int& ik,
                                const int& ikb,
                                const ModuleBase::Vector3<double> G,
                                ModuleBase::ComplexMatrix& Mmn)
{
    Mmn.create(num_bands, num_bands);

    int row = this->ParaV->get_row_size();
    int col = this->ParaV->get_col_size();
    int nloc = row * col;
    std::complex<double>* midmatrix = new std::complex<double>[nloc];
    ModuleBase::GlobalFunc::ZEROS(midmatrix, nloc);

    int R_num = R_coor_car.size();
    ModuleBase::Vector3<double> ik_car = kv.kvec_c[ik];
    ModuleBase::Vector3<double> ikb_car = kv.kvec_c[ikb] + G * GlobalC::ucell.G;
    Abfs::Vector3_Order<double> dk = (ikb_car - ik_car) * GlobalC::ucell.tpiba;
    Coordinate_3D temp_dk(dk.x, dk.y, dk.z);
    int delta_k_index = delta_k_all_index[temp_dk];

    hamilt::HContainer<std::complex<double>>* tmp_FR_container = FR[delta_k_index].get_FR_pointer();
    auto row_indexes = ParaV->get_indexes_row();
    auto col_indexes = ParaV->get_indexes_col();

    for (int iap = 0; iap < tmp_FR_container->size_atom_pairs(); ++iap)
    {
        int atom_i = tmp_FR_container->get_atom_pair(iap).get_atom_i();
        int atom_j = tmp_FR_container->get_atom_pair(iap).get_atom_j();
        int start_i = ParaV->atom_begin_row[atom_i];
        int start_j = ParaV->atom_begin_col[atom_j];
        int row_size = ParaV->get_row_size(atom_i);
        int col_size = ParaV->get_col_size(atom_j);
        for (int iR = 0; iR < tmp_FR_container->get_atom_pair(iap).get_R_size(); ++iR)
        {
            auto& matrix = tmp_FR_container->get_atom_pair(iap).get_HR_values(iR);
            const ModuleBase::Vector3<int> r_index = tmp_FR_container->get_atom_pair(iap).get_R_index(iR);
            ModuleBase::Vector3<double> dR
                = ModuleBase::Vector3<double>(r_index.x, r_index.y, r_index.z) * GlobalC::ucell.latvec;
            double phase = ikb_car * dR * ModuleBase::TWO_PI;
            std::complex<double> kRn_phase = std::exp(ModuleBase::IMAG_UNIT * phase);
            for (int i = 0; i < row_size; ++i)
            {
                int mu = row_indexes[start_i + i];
                int ir = ParaV->global2local_row(mu);
                for (int j = 0; j < col_size; ++j)
                {
                    int nu = col_indexes[start_j + j];
                    int ic = ParaV->global2local_col(nu);
                    int index = ic * row + ir;
                    midmatrix[index] += kRn_phase * matrix.get_value(i, j);
                }
            }
        }
    }

    char transa = 'C';
    char transb = 'N';
    int Bands = PARAM.inp.nbands;
    int nlocal = PARAM.globalv.nlocal;
    std::complex<double> alpha = {1.0, 0.0}, beta = {0.0, 0.0};
    int one = 1;

    std::complex<double>* C_matrix = new std::complex<double>[nloc];
    std::complex<double>* out_matrix = new std::complex<double>[nloc];
    ModuleBase::GlobalFunc::ZEROS(C_matrix, nloc);
    ModuleBase::GlobalFunc::ZEROS(out_matrix, nloc);

#ifdef __MPI
    pzgemm_(&transa,
            &transb,
            &Bands,
            &nlocal,
            &nlocal,
            &alpha,
            &psi_in(ik, 0, 0),
            &one,
            &one,
            this->ParaV->desc,
            midmatrix,
            &one,
            &one,
            this->ParaV->desc,
            &beta,
            C_matrix,
            &one,
            &one,
            this->ParaV->desc);

    pzgemm_(&transb,
            &transb,
            &Bands,
            &Bands,
            &nlocal,
            &alpha,
            C_matrix,
            &one,
            &one,
            this->ParaV->desc,
            &psi_in(ikb, 0, 0),
            &one,
            &one,
            this->ParaV->desc,
            &beta,
            out_matrix,
            &one,
            &one,
            this->ParaV->desc);
#endif

    int count_m = -1;
    for (int m = 0; m < PARAM.inp.nbands; m++)
    {
        if (exclude_bands.count(m)) {
            continue;
}
        count_m++;

        int ir = this->ParaV->global2local_row(m);
        if (ir >= 0)
        {
            int count_n = -1;
            for (int n = 0; n < PARAM.inp.nbands; n++)
            {
                if (exclude_bands.count(n)) {
                    continue;
}
                count_n++;

                int ic = this->ParaV->global2local_col(n);
                if (ic >= 0)
                {
                    int index = ic * row + ir;
                    Mmn(count_m, count_n) = out_matrix[index];
                }
            }
        }
    }

#ifdef __MPI
    Parallel_Reduce::reduce_all(Mmn.c, num_bands * num_bands);
#endif

    delete[] midmatrix;
    delete[] C_matrix;
    delete[] out_matrix;
}

void toWannier90_LCAO::produce_basis_orb()
{
    int mat_Nr = orb_.Phi[0].PhiLN(0, 0).getNr();
    int count_Nr = 0;

    orbs.resize(orb_.get_ntype());
    for (int T = 0; T < orb_.get_ntype(); ++T)
    {
        count_Nr = orb_.Phi[T].PhiLN(0, 0).getNr();
        if (count_Nr > mat_Nr)
        {
            mat_Nr = count_Nr;
            orb_r_ntype = T;
        }

        orbs[T].resize(orb_.Phi[T].getLmax() + 1);
        for (int L = 0; L <= orb_.Phi[T].getLmax(); ++L)
        {
            orbs[T][L].resize(orb_.Phi[T].getNchi(L));
            for (int N = 0; N < orb_.Phi[T].getNchi(L); ++N)
            {
                const auto& orb_origin = orb_.Phi[T].PhiLN(L, N);
                orbs[T][L][N].set_orbital_info(orb_origin.getLabel(),
                                               orb_origin.getType(),
                                               orb_origin.getL(),
                                               orb_origin.getChi(),
                                               orb_origin.getNr(),
                                               orb_origin.getRab(),
                                               orb_origin.getRadial(),
                                               Numerical_Orbital_Lm::Psi_Type::Psi,
                                               orb_origin.getPsi(),
                                               static_cast<int>(orb_origin.getNk() * kmesh_times) | 1,
                                               orb_origin.getDk(),
                                               orb_origin.getDruniform(),
                                               false,
                                               true,
                                               PARAM.inp.cal_force);
            }
        }
    }
}

void toWannier90_LCAO::produce_trial_in_lcao()
{
    A_orbs.resize(num_wannier);

    double* r = new double[mesh_r];
    double* rab = new double[mesh_r];

    for (int ir = 0; ir < mesh_r; ir++)
    {
        rab[ir] = dr;
        r[ir] = ir * dr;
    }

    const auto& orb_origin = orb_.Phi[orb_r_ntype].PhiLN(0, 0);

    double* psi = new double[mesh_r];
    double* psir = new double[mesh_r];
    double* inner = new double[mesh_r];
    for (int i = 0; i < num_wannier; i++)
    {
        double alfa32 = pow(alfa[i], 3.0 / 2.0);
        double alfa_new = alfa[i];
        int wannier_index = i;

        if (rvalue[i] == 1)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = 2.0 * alfa32 * exp(-alfa_new * r[ir]);
            }
        }

        if (rvalue[i] == 2)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = 1.0 / sqrt(8.0) * alfa32 * (2.0 - alfa_new * r[ir]) * exp(-alfa_new * r[ir] * 0.5);
            }
        }

        if (rvalue[i] == 3)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = sqrt(4.0 / 27.0) * alfa32
                          * (1.0 - 2.0 / 3.0 * alfa_new * r[ir] + 2.0 / 27.0 * pow(alfa_new, 2.0) * r[ir] * r[ir])
                          * exp(-alfa_new * r[ir] * 1.0 / 3.0);
            }
        }

        for (int ir = 0; ir < mesh_r; ir++)
        {
            psir[ir] = psi[ir] * r[ir];
        }

        // renormalize radial wave functions
        for (int ir = 0; ir < mesh_r; ir++)
        {
            inner[ir] = psir[ir] * psir[ir];
        }
        double unit = 0.0;
        ModuleBase::Integral::Simpson_Integral(mesh_r, inner, rab, unit);

        for (int ir = 0; ir < mesh_r; ir++)
        {
            psi[ir] /= sqrt(unit);
        }

        if (L[i] >= 0)
        {
            A_orbs[i].resize(1);

            A_orbs[i][0].set_orbital_info(orb_origin.getLabel(),
                                          orb_origin.getType(),
                                          L[i],
                                          1,
                                          mesh_r,
                                          rab,
                                          r,
                                          Numerical_Orbital_Lm::Psi_Type::Psi,
                                          psi,
                                          static_cast<int>(orb_origin.getNk() * kmesh_times) | 1,
                                          orb_origin.getDk(),
                                          orb_origin.getDruniform(),
                                          false,
                                          true,
                                          PARAM.inp.cal_force);
        }
        else
        {
            int tmp_size = 0;

            if (L[i] == -1 || L[i] == -2 || L[i] == -3) {
                tmp_size = 2;
}

            if (L[i] == -4 || L[i] == -5) {
                tmp_size = 3;
}

            A_orbs[i].resize(tmp_size);

            for (int tmp_L = 0; tmp_L < tmp_size; tmp_L++)
            {
                A_orbs[i][tmp_L].set_orbital_info(orb_origin.getLabel(),
                                                  orb_origin.getType(),
                                                  tmp_L,
                                                  1,
                                                  mesh_r,
                                                  rab,
                                                  r,
                                                  Numerical_Orbital_Lm::Psi_Type::Psi,
                                                  psi,
                                                  static_cast<int>(orb_origin.getNk() * kmesh_times) | 1,
                                                  orb_origin.getDk(),
                                                  orb_origin.getDruniform(),
                                                  false,
                                                  true,
                                                  PARAM.inp.cal_force);
            }
        }
    }

    delete[] r;
    delete[] rab;
    delete[] psi;
    delete[] psir;
    delete[] inner;
}

void toWannier90_LCAO::construct_overlap_table_project()
{
    int row = this->ParaV->get_row_size();
    int global_ir = 0;

    for (int ir = 0; ir < row; ir += PARAM.globalv.npol)
    {
        global_ir = ParaV->local2global_row(ir);
        int orb_index_row = global_ir / PARAM.globalv.npol;

        for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
        {
            if (L[wannier_index] >= 0)
            {
                center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].insert(std::make_pair(
                    0,
                    Center2_Orb::Orb11(orbs[iw2it[orb_index_row]][iw2iL[orb_index_row]][iw2iN[orb_index_row]],
                                       A_orbs[wannier_index][0],
                                       psb_,
                                       MGT)));
            }
            else
            {
                int tmp_size = 0;

                if (L[wannier_index] == -1 || L[wannier_index] == -2 || L[wannier_index] == -3) {
                    tmp_size = 2;
}

                if (L[wannier_index] == -4 || L[wannier_index] == -5) {
                    tmp_size = 3;
}

                for (int tmp_L = 0; tmp_L < tmp_size; tmp_L++)
                {
                    center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].insert(std::make_pair(
                        tmp_L,
                        Center2_Orb::Orb11(orbs[iw2it[orb_index_row]][iw2iL[orb_index_row]][iw2iN[orb_index_row]],
                                           A_orbs[wannier_index][tmp_L],
                                           psb_,
                                           MGT)));
                }
            }
        }
    }

    for (auto& co1: center2_orb11_A)
    {
        for (auto& co2: co1.second)
        {
            for (auto& co3: co2.second)
            {
                co3.second.init_radial_table();
            }
        }
    }
}

void toWannier90_LCAO::cal_orbA_overlap_R()
{
    int row = this->ParaV->get_row_size();
    int R_num = R_coor_car.size();
    int global_ir = 0;

    psi_psiA_R.resize(row);
    for (int ir = 0; ir < row; ir++)
    {
        psi_psiA_R[ir].resize(num_wannier);

        for (int ic = 0; ic < num_wannier; ++ic)
        {
            psi_psiA_R[ir][ic].resize(R_num);
        }
    }

    double bs2, bs3, bs6, bs12;
    bs2 = 1.0 / sqrt(2.0);
    bs3 = 1.0 / sqrt(3.0);
    bs6 = 1.0 / sqrt(6.0);
    bs12 = 1.0 / sqrt(12.0);

    for (int ir = 0; ir < row; ir++)
    {
        global_ir = ParaV->local2global_row(ir);
        int orb_index_row = global_ir / PARAM.globalv.npol;

        int it1 = iw2it[orb_index_row];
        int ia1 = iw2ia[orb_index_row];
        int im1 = iw2im[orb_index_row];

        for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
        {
            if (L[wannier_index] >= 0)
            {
                for (int iR = 0; iR < R_num; iR++)
                {
                    ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                    ModuleBase::Vector3<double> orb_center
                        = (GlobalC::ucell.atoms[it1].tau[ia1] + R_car) * GlobalC::ucell.lat0;
                    ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * GlobalC::ucell.lat0;

                    double overlap_o
                        = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(orb_center,
                                                                                                   project_orb_center,
                                                                                                   im1,
                                                                                                   m[wannier_index]);
                    psi_psiA_R[ir][wannier_index][iR] = overlap_o;
                }
            }
            else
            {
                if (L[wannier_index] == -1)
                {
                    double tmp_bs2 = 0;
                    if (m[wannier_index] == 0) {
                        tmp_bs2 = bs2;
}
                    if (m[wannier_index] == -1) {
                        tmp_bs2 = -bs2;
}

                    for (int iR = 0; iR < R_num; iR++)
                    {
                        ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                        ModuleBase::Vector3<double> orb_center
                            = (GlobalC::ucell.atoms[it1].tau[ia1] + R_car) * GlobalC::ucell.lat0;
                        ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * GlobalC::ucell.lat0;

                        double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                            orb_center,
                            project_orb_center,
                            im1,
                            0);
                        double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap(
                            orb_center,
                            project_orb_center,
                            im1,
                            1);
                        psi_psiA_R[ir][wannier_index][iR] = bs2 * overlap_s + tmp_bs2 * overlap_px;
                    }
                }

                if (L[wannier_index] == -2)
                {
                    if (m[wannier_index] == 0 || m[wannier_index] == 1)
                    {
                        double tmp_bs2 = bs2;
                        if (m[wannier_index] == -1) {
                            tmp_bs2 = -bs2;
}

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (GlobalC::ucell.atoms[it1].tau[ia1] + R_car) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                                orb_center,
                                project_orb_center,
                                im1,
                                0);
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 1);
                            double overlap_py = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 2);
                            psi_psiA_R[ir][wannier_index][iR]
                                = bs3 * overlap_s - bs6 * overlap_px + tmp_bs2 * overlap_py;
                        }
                    }
                    else if (m[wannier_index] == 2)
                    {
                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (GlobalC::ucell.atoms[it1].tau[ia1] + R_car) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                                orb_center,
                                project_orb_center,
                                im1,
                                0);
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 1);
                            psi_psiA_R[ir][wannier_index][iR] = bs3 * overlap_s + 2.0 * bs6 * overlap_px;
                        }
                    }
                }

                if (L[wannier_index] == -3)
                {
                    double m_px = 1.0;
                    double m_py = 1.0;
                    double m_pz = 1.0;

                    if (m[wannier_index] == 1)
                    {
                        m_py = -1.0;
                        m_pz = -1.0;
                    }
                    else if (m[wannier_index] == 2)
                    {
                        m_px = -1.0;
                        m_pz = -1.0;
                    }
                    else if (m[wannier_index] == 3)
                    {
                        m_px = -1.0;
                        m_py = -1.0;
                    }

                    for (int iR = 0; iR < R_num; iR++)
                    {
                        ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                        ModuleBase::Vector3<double> orb_center
                            = (GlobalC::ucell.atoms[it1].tau[ia1] + R_car) * GlobalC::ucell.lat0;
                        ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * GlobalC::ucell.lat0;

                        double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                            orb_center,
                            project_orb_center,
                            im1,
                            0);
                        double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap(
                            orb_center,
                            project_orb_center,
                            im1,
                            1);
                        double overlap_py = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap(
                            orb_center,
                            project_orb_center,
                            im1,
                            2);
                        double overlap_pz = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap(
                            orb_center,
                            project_orb_center,
                            im1,
                            0);
                        psi_psiA_R[ir][wannier_index][iR]
                            = 0.5 * (overlap_s + m_px * overlap_px + m_py * overlap_py + m_pz * overlap_pz);
                    }
                }

                if (L[wannier_index] == -4)
                {
                    if (m[wannier_index] == 0 || m[wannier_index] == 1)
                    {
                        double tmp_bs2 = bs2;
                        if (m[wannier_index] == -1) {
                            tmp_bs2 = -bs2;
}

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (GlobalC::ucell.atoms[it1].tau[ia1] + R_car) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                                orb_center,
                                project_orb_center,
                                im1,
                                0);
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 1);
                            double overlap_py = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 2);
                            psi_psiA_R[ir][wannier_index][iR]
                                = bs3 * overlap_s - bs6 * overlap_px + tmp_bs2 * overlap_py;
                        }
                    }
                    else if (m[wannier_index] == 2)
                    {
                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (GlobalC::ucell.atoms[it1].tau[ia1] + R_car) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                                orb_center,
                                project_orb_center,
                                im1,
                                0);
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 1);
                            psi_psiA_R[ir][wannier_index][iR] = bs3 * overlap_s + 2.0 * bs6 * overlap_px;
                        }
                    }
                    else if (m[wannier_index] == 3 || m[wannier_index] == 4)
                    {
                        double m_pz = 1.0;
                        if (m[wannier_index] == 4) {
                            m_pz = -1.0;
}

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (GlobalC::ucell.atoms[it1].tau[ia1] + R_car) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_pz = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 0);
                            double overlap_dz2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                     .at(2)
                                                     .cal_overlap(orb_center, project_orb_center, im1, 0);
                            psi_psiA_R[ir][wannier_index][iR] = bs2 * (m_pz * overlap_pz + overlap_dz2);
                        }
                    }
                }

                if (L[wannier_index] == -5)
                {
                    if (m[wannier_index] == 0 || m[wannier_index] == 1)
                    {
                        double tmp_bs2 = -bs2;
                        double tmp_bs12 = -bs12;
                        double tmp_d = 0.5;

                        if (m[wannier_index] == 1)
                        {
                            tmp_bs2 = bs2;
                        }

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (GlobalC::ucell.atoms[it1].tau[ia1] + R_car) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                                orb_center,
                                project_orb_center,
                                im1,
                                0);
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 1);
                            double overlap_dz2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                     .at(2)
                                                     .cal_overlap(orb_center, project_orb_center, im1, 0);
                            double overlap_dx2_y2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                        .at(2)
                                                        .cal_overlap(orb_center, project_orb_center, im1, 3);
                            psi_psiA_R[ir][wannier_index][iR] = bs6 * overlap_s + tmp_bs2 * overlap_px
                                                                + tmp_bs12 * overlap_dz2 + tmp_d * overlap_dx2_y2;
                        }
                    }
                    else if (m[wannier_index] == 2 || m[wannier_index] == 3)
                    {
                        double tmp_bs2 = -bs2;
                        double tmp_bs12 = -bs12;
                        double tmp_d = -0.5;

                        if (m[wannier_index] == 3)
                        {
                            tmp_bs2 = bs2;
                        }

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (GlobalC::ucell.atoms[it1].tau[ia1] + R_car) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                                orb_center,
                                project_orb_center,
                                im1,
                                0);
                            double overlap_py = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 2);
                            double overlap_dz2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                     .at(2)
                                                     .cal_overlap(orb_center, project_orb_center, im1, 0);
                            double overlap_dx2_y2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                        .at(2)
                                                        .cal_overlap(orb_center, project_orb_center, im1, 3);
                            psi_psiA_R[ir][wannier_index][iR] = bs6 * overlap_s + tmp_bs2 * overlap_py
                                                                + tmp_bs12 * overlap_dz2 + tmp_d * overlap_dx2_y2;
                        }
                    }
                    else if (m[wannier_index] == 4 || m[wannier_index] == 5)
                    {
                        double tmp_pz = -1.0;

                        if (m[wannier_index] == 5) {
                            tmp_pz = 1.0;
}

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center
                                = (GlobalC::ucell.atoms[it1].tau[ia1] + R_car) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center
                                = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap(
                                orb_center,
                                project_orb_center,
                                im1,
                                0);
                            double overlap_pz = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                    .at(1)
                                                    .cal_overlap(orb_center, project_orb_center, im1, 0);
                            double overlap_dz2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index]
                                                     .at(2)
                                                     .cal_overlap(orb_center, project_orb_center, im1, 0);
                            psi_psiA_R[ir][wannier_index][iR]
                                = bs6 * overlap_s + tmp_pz * bs2 * overlap_pz + bs3 * overlap_dz2;
                        }
                    }
                }
            }
        }
    }
}

void toWannier90_LCAO::unkdotA(const K_Vectors& kv,
                               const psi::Psi<std::complex<double>>& psi_in,
                               const int& ik,
                               ModuleBase::ComplexMatrix& Amn)
{
    Amn.create(num_bands, num_wannier);

    int row = this->ParaV->get_row_size();
    int index_band = -1;
    int R_num = R_coor_car.size();
    if (PARAM.inp.nspin != 4)
    {
        for (int ib = 0; ib < PARAM.inp.nbands; ib++)
        {
            if (exclude_bands.count(ib)) {
                continue;
}
            index_band++;

            int ic = this->ParaV->global2local_col(ib);
            if (ic >= 0)
            {
                for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
                {
                    for (int iR = 0; iR < R_num; iR++)
                    {
                        ModuleBase::Vector3<double> R = R_coor_car[iR];
                        double kRn = -1.0 * kv.kvec_c[ik] * R * ModuleBase::TWO_PI;
                        std::complex<double> kRn_phase(cos(kRn), sin(kRn));
                        std::complex<double> tmp(0.0, 0.0);

                        for (int ir = 0; ir < row; ir++)
                        {
                            tmp += std::conj(psi_in(ik, ic, ir)) * psi_psiA_R[ir][wannier_index][iR];
                        }

                        Amn(index_band, wannier_index) += kRn_phase * tmp;
                    }
                }
            }
        }
    }
    else
    {
        for (int ib = 0; ib < PARAM.inp.nbands; ib++)
        {
            if (exclude_bands.count(ib)) {
                continue;
}
            index_band++;

            int ic = this->ParaV->global2local_col(ib);
            if (ic >= 0)
            {
                for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
                {
                    for (int iR = 0; iR < R_num; iR++)
                    {
                        ModuleBase::Vector3<double> R = R_coor_car[iR];
                        double kRn = kv.kvec_c[ik] * R * ModuleBase::TWO_PI;
                        std::complex<double> kRn_phase(cos(kRn), sin(kRn));
                        std::complex<double> tmp(0.0, 0.0);
                        int half_row = row / 2;
                        for (int ir = 0; ir < half_row; ir++)
                        {
                            tmp += up_con[wannier_index] * std::conj(psi_in(ik, ic, 2 * ir))
                                       * psi_psiA_R[2 * ir][wannier_index][iR]
                                   + dn_con[wannier_index] * std::conj(psi_in(ik, ic, 2 * ir + 1))
                                         * psi_psiA_R[2 * ir + 1][wannier_index][iR];
                        }

                        Amn(index_band, wannier_index) += kRn_phase * tmp;
                    }
                }
            }
        }
    }

#ifdef __MPI
    Parallel_Reduce::reduce_all(Amn.c, Amn.size);
#endif
}
#endif

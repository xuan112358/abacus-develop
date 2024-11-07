#include "LCAO_deepks_test.h"
#define private public
#include "module_parameter/parameter.h"
#undef private
namespace Test_Deepks
{
Grid_Driver GridD(PARAM.input.test_deconstructor, PARAM.input.test_grid);
}

test_deepks::test_deepks()
{
}

test_deepks::~test_deepks()
{
}

void test_deepks::check_dstable()
{
    // OGT.talpha.print_Table_DSR(ORB);
    // this->compare_with_ref("S_I_mu_alpha.dat","S_I_mu_alpha_ref.dat");
}

void test_deepks::check_psialpha()
{
    std::vector<int> na;
    na.resize(ucell.ntype);
    for (int it = 0; it < ucell.ntype; it++)
    {
        na[it] = ucell.atoms[it].na;
    }
    ld.init(ORB, ucell.nat, ucell.ntype, ParaO, na);

    ld.build_psialpha(PARAM.input.cal_force, ucell, ORB, Test_Deepks::GridD, overlap_orb_alpha_);

    ld.check_psialpha(PARAM.input.cal_force, ucell, ORB, Test_Deepks::GridD);

    this->compare_with_ref("psialpha.dat", "psialpha_ref.dat");
    this->compare_with_ref("dpsialpha_x.dat", "dpsialpha_x_ref.dat");
    this->compare_with_ref("dpsialpha_y.dat", "dpsialpha_y_ref.dat");
    this->compare_with_ref("dpsialpha_z.dat", "dpsialpha_z_ref.dat");
}

void test_deepks::read_dm()
{
    std::ifstream ifs("dm");
    dm.resize(1);
    dm[0].create(PARAM.sys.nlocal, PARAM.sys.nlocal);

    for (int mu = 0; mu < PARAM.sys.nlocal; mu++)
    {
        for (int nu = 0; nu < PARAM.sys.nlocal; nu++)
        {
            double c;
            ifs >> c;
            dm[0](mu, nu) = c;
        }
    }
}

void test_deepks::read_dm_k(const int nks)
{
    dm_k.resize(nks);
    std::stringstream ss;
    for (int ik = 0; ik < nks; ik++)
    {
        ss.str("");
        ss << "dm_" << ik;
        std::ifstream ifs(ss.str().c_str());
        dm_k[ik].create(PARAM.sys.nlocal, PARAM.sys.nlocal);

        for (int mu = 0; mu < PARAM.sys.nlocal; mu++)
        {
            for (int nu = 0; nu < PARAM.sys.nlocal; nu++)
            {
                std::complex<double> c;
                ifs >> c;
                dm_k[ik](mu, nu) = c;
            }
        }
    }
}

void test_deepks::set_dm_new()
{
    // dm_gamma
    dm_new.resize(dm.size());
    for (int i = 0; i < dm.size(); i++)
    {
        dm_new[i].resize(dm[i].nr * dm[i].nc);
        dm_new[i].assign(dm[i].c, dm[i].c + dm[i].nr * dm[i].nc);
    }
}

void test_deepks::set_dm_k_new()
{
    // dm_k
    dm_k_new.resize(dm_k.size());
    for (int i = 0; i < dm_k.size(); i++)
    {
        dm_k_new[i].resize(dm_k[i].nr * dm_k[i].nc);
        dm_k_new[i].assign(dm_k[i].c, dm_k[i].c + dm_k[i].nr * dm_k[i].nc);
    }
}

void test_deepks::set_p_elec_DM()
{
    // gamma
    int nspin=PARAM.inp.nspin;
    this->p_elec_DM=new elecstate::DensityMatrix<double, double>(&ParaO,nspin);
    p_elec_DM->init_DMR(&Test_Deepks::GridD, &ucell);
    for (int ik = 0; ik < nspin; ik++)
    {
        p_elec_DM->set_DMK_pointer(ik, dm_new[ik].data());
    }
    p_elec_DM->cal_DMR();
}

void test_deepks::set_p_elec_DM_k()
{
    // multi k
    this->p_elec_DM_k=new elecstate::DensityMatrix<std::complex<double>, double>(&ParaO, PARAM.inp.nspin, kv.kvec_d, kv.nkstot / PARAM.inp.nspin);
    p_elec_DM_k->init_DMR(&Test_Deepks::GridD, &ucell);
    for (int ik = 0; ik < kv.nkstot; ++ik)
    {
        p_elec_DM_k->set_DMK_pointer(ik, dm_k_new[ik].data());
    }
    p_elec_DM_k->cal_DMR();
}

void test_deepks::check_pdm()
{
    if (PARAM.sys.gamma_only_local)
    {
        this->read_dm();
        this->set_dm_new();
        this->set_p_elec_DM();
        this->ld.cal_projected_DM(p_elec_DM, ucell, ORB, Test_Deepks::GridD);
    }
    else
    {
        this->read_dm_k(kv.nkstot);
        this->set_dm_k_new();
        this->set_p_elec_DM_k();
        this->ld.cal_projected_DM_k(p_elec_DM_k, ucell, ORB, Test_Deepks::GridD);
    }
    this->ld.check_projected_dm();
    this->compare_with_ref("deepks_projdm.dat", "pdm_ref.dat");
}

void test_deepks::check_gdmx()
{
    this->ld.init_gdmx(ucell.nat);
    if (PARAM.sys.gamma_only_local)
    {
        this->ld.cal_gdmx(dm_new[0], ucell, ORB, Test_Deepks::GridD, 0);
    }
    else
    {
        this->ld.cal_gdmx_k(dm_k_new, ucell, ORB, Test_Deepks::GridD, kv.nkstot, kv.kvec_d, 0);
    }
    this->ld.check_gdmx(ucell.nat);

    for (int ia = 0; ia < ucell.nat; ia++)
    {
        std::stringstream ss;
        std::stringstream ss1;
        ss.str("");
        ss << "gdmx_" << ia << ".dat";
        ss1.str("");
        ss1 << "gdmx_" << ia << "_ref.dat";

        this->compare_with_ref(ss.str(), ss1.str());

        ss.str("");
        ss << "gdmy_" << ia << ".dat";
        ss1.str("");
        ss1 << "gdmy_" << ia << "_ref.dat";
        this->compare_with_ref(ss.str(), ss1.str());

        ss.str("");
        ss << "gdmz_" << ia << ".dat";
        ss1.str("");
        ss1 << "gdmz_" << ia << "_ref.dat";
        this->compare_with_ref(ss.str(), ss1.str());
    }
}

void test_deepks::check_descriptor()
{
    this->ld.cal_descriptor(ucell.nat);
    this->ld.check_descriptor(ucell,"./");
    this->compare_with_ref("deepks_desc.dat", "descriptor_ref.dat");
}

void test_deepks::check_gvx()
{
    this->ld.cal_gvx(ucell.nat);
    this->ld.check_gvx(ucell.nat);

    for (int ia = 0; ia < ucell.nat; ia++)
    {
        std::stringstream ss;
        std::stringstream ss1;
        ss.str("");
        ss << "gvx_" << ia << ".dat";
        ss1.str("");
        ss1 << "gvx_" << ia << "_ref.dat";
        this->compare_with_ref(ss.str(), ss1.str());

        ss.str("");
        ss << "gvy_" << ia << ".dat";
        ss1.str("");
        ss1 << "gvy_" << ia << "_ref.dat";
        this->compare_with_ref(ss.str(), ss1.str());

        ss.str("");
        ss << "gvz_" << ia << ".dat";
        ss1.str("");
        ss1 << "gvz_" << ia << "_ref.dat";
        this->compare_with_ref(ss.str(), ss1.str());
    }
}

void test_deepks::check_edelta()
{
    this->ld.load_model("model.ptg");
    if (PARAM.sys.gamma_only_local)
    {
        this->ld.allocate_V_delta(ucell.nat);
    }
    else
    {
        this->ld.allocate_V_delta(ucell.nat, kv.nkstot);
    }
    this->ld.cal_gedm(ucell.nat);

    std::ofstream ofs("E_delta.dat");
    ofs << std::setprecision(10) << this->ld.E_delta << std::endl;
    ofs.close();
    this->compare_with_ref("E_delta.dat", "E_delta_ref.dat");

    this->ld.check_gedm();
    this->compare_with_ref("gedm.dat", "gedm_ref.dat");
}

void test_deepks::check_e_deltabands()
{
    if (PARAM.sys.gamma_only_local)
    {
        this->ld.cal_e_delta_band(dm_new);
    }
    else
    {
        this->folding_nnr(kv);
        this->ld.cal_e_delta_band_k(dm_k_new, kv.nkstot);
    }

    std::ofstream ofs("E_delta_bands.dat");
    ofs << std::setprecision(10) << this->ld.e_delta_band << std::endl;
    ofs.close();
    this->compare_with_ref("E_delta_bands.dat", "E_delta_bands_ref.dat");
}

void test_deepks::check_f_delta()
{
    ModuleBase::matrix svnl_dalpha;
    svnl_dalpha.create(3, 3);
    if (PARAM.sys.gamma_only_local)
    {
        DeePKS_domain::cal_f_delta_gamma(dm_new, ucell, ORB, Test_Deepks::GridD, ParaO, this->ld.lmaxd, this->ld.nlm_save, this->ld.gedm, this->ld.inl_index, this->ld.F_delta, 1, svnl_dalpha);
    }
    else
    {
        DeePKS_domain::cal_f_delta_k(dm_k_new, ucell, ORB, Test_Deepks::GridD, ParaO, this->ld.lmaxd, kv.nkstot, kv.kvec_d, this->ld.nlm_save_k, this->ld.gedm, this->ld.inl_index, this->ld.F_delta, 1, svnl_dalpha);
    }
    DeePKS_domain::check_f_delta(ucell.nat, this->ld.F_delta, svnl_dalpha);

    this->compare_with_ref("F_delta.dat", "F_delta_ref.dat");
}

void test_deepks::compare_with_ref(const std::string f1, const std::string f2)
{
    this->total_check += 1;
    std::ifstream file1(f1.c_str());
    std::ifstream file2(f2.c_str());
    double test_thr = 1e-8;

    std::string word1;
    std::string word2;
    while (file1 >> word1)
    {
        file2 >> word2;
        if ((word1[0] - '0' >= 0 && word1[0] - '0' < 10) || word1[0] == '-')
        {
            double num1 = std::stof(word1);
            double num2 = std::stof(word2);
            if (std::abs(num1 - num2) > test_thr)
            {
                this->failed_check += 1;
                std::cout << "\e[1;31m [  FAILED  ] \e[0m" << f1.c_str() << " inconsistent!" << std::endl;
                return;
            }
        }
        else
        {
            if (word1 != word2)
            {
                this->failed_check += 1;
                return;
            }
        }
    }
    return;
}

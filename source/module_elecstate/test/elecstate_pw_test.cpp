#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "module_parameter/parameter.h"
#undef private
#define protected public
#include "module_elecstate/elecstate_pw.h"
#undef protected
// mock functions for testing
namespace elecstate
{
double get_ucell_omega()
{
    return 500.0;
}
double get_ucell_tpiba()
{
    return 2.0;
}
int tmp_xc_func_type = 1;
int get_xc_func_type()
{
    return tmp_xc_func_type;
}
void Potential::init_pot(int, Charge const*)
{
}
void Potential::cal_v_eff(const Charge* chg, const UnitCell* ucell, ModuleBase::matrix& v_eff)
{
}
void Potential::cal_fixed_v(double* vl_pseudo)
{
}
Potential::~Potential()
{
}
} // namespace elecstate
Charge::Charge()
{
}
Charge::~Charge()
{
}
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
ORB_gaunt_table::ORB_gaunt_table()
{
}
ORB_gaunt_table::~ORB_gaunt_table()
{
}
#endif
pseudopot_cell_vl::pseudopot_cell_vl()
{
}
pseudopot_cell_vl::~pseudopot_cell_vl()
{
}
pseudopot_cell_vnl::pseudopot_cell_vnl()
{
}
pseudopot_cell_vnl::~pseudopot_cell_vnl()
{
}
template <>
void pseudopot_cell_vnl::radial_fft_q<float, base_device::DEVICE_CPU>(base_device::DEVICE_CPU* ctx,
                                                                      const int ng,
                                                                      const int ih,
                                                                      const int jh,
                                                                      const int itype,
                                                                      const float* qnorm,
                                                                      const float* ylm,
                                                                      std::complex<float>* qg) const
{
}
template <>
void pseudopot_cell_vnl::radial_fft_q<double, base_device::DEVICE_CPU>(base_device::DEVICE_CPU* ctx,
                                                                       const int ng,
                                                                       const int ih,
                                                                       const int jh,
                                                                       const int itype,
                                                                       const double* qnorm,
                                                                       const double* ylm,
                                                                       std::complex<double>* qg) const
{
}
template <>
std::complex<float>* pseudopot_cell_vnl::get_vkb_data<float>() const
{
    return nullptr;
}
template <>
std::complex<double>* pseudopot_cell_vnl::get_vkb_data<double>() const
{
    return nullptr;
}
template <>
void pseudopot_cell_vnl::getvnl<float, base_device::DEVICE_CPU>(base_device::DEVICE_CPU*,
                                                                int const&,
                                                                std::complex<float>*) const
{
}
template <>
void pseudopot_cell_vnl::getvnl<double, base_device::DEVICE_CPU>(base_device::DEVICE_CPU*,
                                                                 int const&,
                                                                 std::complex<double>*) const
{
}
Soc::~Soc()
{
}
Fcoef::~Fcoef()
{
}
#include "module_cell/klist.h"
K_Vectors::K_Vectors()
{
}
K_Vectors::~K_Vectors()
{
}

void Charge::set_rho_core(ModuleBase::ComplexMatrix const&, const bool*)
{
}
void Charge::set_rho_core_paw()
{
}
void Charge::init_rho(elecstate::efermi&,
                      ModuleBase::ComplexMatrix const&,
                      ModuleSymmetry::Symmetry& symm,
                      const void*,
                      const void*)
{
}
void Charge::set_rhopw(ModulePW::PW_Basis*)
{
}
void Charge::renormalize_rho()
{
}
void Charge::check_rho()
{
}

void Set_GlobalV_Default()
{
    PARAM.input.device = "cpu";
    PARAM.input.precision = "double";
    PARAM.sys.domag = false;
    PARAM.sys.domag_z = false;
    // Base class dependent
    PARAM.input.nspin = 1;
    PARAM.input.nelec = 10.0;
    PARAM.input.nupdown  = 0.0;
    PARAM.sys.two_fermi = false;
    PARAM.input.nbands = 6;
    PARAM.sys.nlocal = 6;
    PARAM.input.esolver_type = "ksdft";
    PARAM.input.lspinorb = false;
    PARAM.input.basis_type = "pw";
    GlobalV::KPAR = 1;
    GlobalV::NPROC_IN_POOL = 1;
    PARAM.sys.use_uspp = false;
}

/************************************************
 *  unit test of elecstate_pw.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Constructor: elecstate::ElecStatePW constructor and destructor
 *      - including double and single precision versions
 *   - InitRhoData: elecstate::ElecStatePW::init_rho_data()
 *      - get rho and kin_r for ElecStatePW
 *   - ParallelK: elecstate::ElecStatePW::parallelK()
 *      - trivial call due to removing of __MPI
 *   - todo: psiToRho: elecstate::ElecStatePW::psiToRho()
 */

class ElecStatePWTest : public ::testing::Test
{
  protected:
    elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_CPU>* elecstate_pw_d = nullptr;
    elecstate::ElecStatePW<std::complex<float>, base_device::DEVICE_CPU>* elecstate_pw_s = nullptr;
    ModulePW::PW_Basis_K* wfcpw = nullptr;
    Charge* chg = nullptr;
    K_Vectors* klist = nullptr;
    UnitCell* ucell = nullptr;
    pseudopot_cell_vnl* ppcell = nullptr;
    ModulePW::PW_Basis* rhodpw = nullptr;
    ModulePW::PW_Basis* rhopw = nullptr;
    ModulePW::PW_Basis_Big* bigpw = nullptr;
    void SetUp() override
    {
        Set_GlobalV_Default();
        wfcpw = new ModulePW::PW_Basis_K;
        chg = new Charge;
        klist = new K_Vectors;
        klist->set_nks(5);
        ucell = new UnitCell;
        ppcell = new pseudopot_cell_vnl;
        rhodpw = new ModulePW::PW_Basis;
        rhopw = new ModulePW::PW_Basis;
        bigpw = new ModulePW::PW_Basis_Big;
    }

    void TearDown() override
    {
        delete wfcpw;
        delete chg;
        delete klist;
        delete ucell;
        delete ppcell;
        delete rhodpw;
        delete rhopw;
        if (elecstate_pw_d != nullptr)
        {
            delete elecstate_pw_d;
        }
        if (elecstate_pw_s != nullptr)
        {
            delete elecstate_pw_s;
        }
    }
};

TEST_F(ElecStatePWTest, ConstructorDouble)
{
    elecstate_pw_d = new elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_CPU>(wfcpw,
                                                                                               chg,
                                                                                               klist,
                                                                                               ucell,
                                                                                               ppcell,
                                                                                               rhodpw,
                                                                                               rhopw,
                                                                                               bigpw);
    EXPECT_EQ(elecstate_pw_d->classname, "ElecStatePW");
    EXPECT_EQ(elecstate_pw_d->charge, chg);
    EXPECT_EQ(elecstate_pw_d->klist, klist);
    EXPECT_EQ(elecstate_pw_d->bigpw, bigpw);
}

TEST_F(ElecStatePWTest, ConstructorSingle)
{
    elecstate_pw_s = new elecstate::ElecStatePW<std::complex<float>, base_device::DEVICE_CPU>(wfcpw,
                                                                                              chg,
                                                                                              klist,
                                                                                              ucell,
                                                                                              ppcell,
                                                                                              rhodpw,
                                                                                              rhopw,
                                                                                              bigpw);
    EXPECT_EQ(elecstate_pw_s->classname, "ElecStatePW");
    EXPECT_EQ(elecstate_pw_s->charge, chg);
    EXPECT_EQ(elecstate_pw_s->klist, klist);
    EXPECT_EQ(elecstate_pw_s->bigpw, bigpw);
}

TEST_F(ElecStatePWTest, InitRhoDataDouble)
{
    elecstate::tmp_xc_func_type = 3;
    chg->nrxx = 1000;
    elecstate_pw_d = new elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_CPU>(wfcpw,
                                                                                               chg,
                                                                                               klist,
                                                                                               ucell,
                                                                                               ppcell,
                                                                                               rhodpw,
                                                                                               rhopw,
                                                                                               bigpw);
    elecstate_pw_d->init_rho_data();
    EXPECT_EQ(elecstate_pw_d->init_rho, true);
    EXPECT_EQ(elecstate_pw_d->rho, chg->rho);
    EXPECT_EQ(elecstate_pw_d->kin_r, chg->kin_r);
}

TEST_F(ElecStatePWTest, InitRhoDataSingle)
{
    PARAM.input.precision = "single";
    elecstate::tmp_xc_func_type = 3;
    chg->nspin = PARAM.input.nspin;
    chg->nrxx = 1000;
    elecstate_pw_s = new elecstate::ElecStatePW<std::complex<float>, base_device::DEVICE_CPU>(wfcpw,
                                                                                              chg,
                                                                                              klist,
                                                                                              ucell,
                                                                                              ppcell,
                                                                                              rhodpw,
                                                                                              rhopw,
                                                                                              bigpw);
    elecstate_pw_s->init_rho_data();
    EXPECT_EQ(elecstate_pw_s->init_rho, true);
    EXPECT_NE(elecstate_pw_s->rho, nullptr);
    EXPECT_NE(elecstate_pw_s->kin_r, nullptr);
}

TEST_F(ElecStatePWTest, ParallelKDouble)
{
    //this is a trivial call due to removing of __MPI
    elecstate_pw_d = new elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_CPU>(wfcpw,
                                                                                               chg,
                                                                                               klist,
                                                                                               ucell,
                                                                                               ppcell,
                                                                                               rhodpw,
                                                                                               rhopw,
                                                                                               bigpw);
    EXPECT_NO_THROW(elecstate_pw_d->parallelK());
}

TEST_F(ElecStatePWTest, ParallelKSingle)
{
    //this is a trivial call due to removing of __MPI
    elecstate_pw_s = new elecstate::ElecStatePW<std::complex<float>, base_device::DEVICE_CPU>(wfcpw,
                                                                                              chg,
                                                                                              klist,
                                                                                              ucell,
                                                                                              ppcell,
                                                                                              rhodpw,
                                                                                              rhopw,
                                                                                              bigpw);
    EXPECT_NO_THROW(elecstate_pw_s->parallelK());
}

#undef protected
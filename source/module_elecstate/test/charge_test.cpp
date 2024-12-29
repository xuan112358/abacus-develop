#include "gtest/gtest.h"
#include "gmock/gmock.h"

#define private public
#define protected public
#include "module_parameter/parameter.h"
#include "module_cell/unitcell.h"
#include "module_elecstate/module_charge/charge.h"
#include "prepare_unitcell.h"
#undef protected
#undef private

// mock functions for UnitCell
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
#endif
Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}

// mock functions for Charge
namespace elecstate
{
int tmp_xc_func_type = 1;
int get_xc_func_type()
{
    return tmp_xc_func_type;
}
double tmp_ucell_omega = 500.0;
double tmp_gridecut = 80.0;
void Set_GlobalV_Default()
{
    PARAM.input.nspin = 1;
    PARAM.input.test_charge = 0;
    PARAM.input.nelec = 8;
}
} // namespace elecstate

/************************************************
 *  unit test of module_charge/charge.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Constructor: Charge::Charge() and Charge::~Charge()
 *     - this is a trivial test
 *   - Allocate: Charge::set_rhopw(), Charge::allocate(), Charge::destroy()
 *     - allocate rho, rhog, rho_save, rhog_save, kin_r, kin_r_save
 *     - using rhopw and PARAM.input.nspin
 *   - SumRho: Charge::sum_rho()
 *     - calculate \sum_{is}^nspin \sum_{ir}^nrxx rho[is][ir]
 *   - RenormalizeRho: Charge::renormalize_rho()
 *     - renormalize rho so as to ensure the sum of rho equals to total number of electrons
 *   - CheckNe: Charge::cal_rho2ne()
 *     - check the total number of electrons summed from rho[is]
 *   - SaveRhoBeforeSumBand: Charge::save_rho_before_sum_band()
 *     - meaning as the function name
 *   - InitFinalScf:: Charge::init_final_scf()
 *     - similar to Charge::allocate(), but for final scf
 */

class ChargeTest : public ::testing::Test
{
  protected:
    UcellTestPrepare utp = UcellTestLib["Si"];
    std::unique_ptr<UnitCell> ucell;
    Charge* charge;
    ModulePW::PW_Basis* rhopw;
    std::string output;
    void SetUp() override
    {
        elecstate::Set_GlobalV_Default();
        ucell = utp.SetUcellInfo();
        charge = new Charge;
        rhopw = new ModulePW::PW_Basis;
        rhopw->initgrids(ucell->lat0, ucell->latvec, elecstate::tmp_gridecut);
        rhopw->distribute_r();
        rhopw->initparameters(false, elecstate::tmp_gridecut);
        rhopw->distribute_g();
    }
    void TearDown() override
    {
        delete charge;
        delete rhopw;
    }
};

TEST_F(ChargeTest, Constructor)
{
    EXPECT_FALSE(charge->allocate_rho);
    EXPECT_FALSE(charge->allocate_rho_final_scf);
}

TEST_F(ChargeTest, Allocate)
{
    // ucell info
    EXPECT_DOUBLE_EQ(ucell->omega, 265.302);
    // rhopw info
    EXPECT_DOUBLE_EQ(rhopw->lat0, 10.2);
    EXPECT_EQ(rhopw->nx, 24);
    EXPECT_EQ(rhopw->ny, 24);
    EXPECT_EQ(rhopw->nz, 24);
    EXPECT_EQ(rhopw->nxyz, 13824);
    EXPECT_EQ(rhopw->nrxx, 13824);
    EXPECT_EQ(rhopw->npw, 3143);
    EXPECT_EQ(rhopw->npwtot, 3143);
    // call Charge::allocate()
    PARAM.input.test_charge = 2;
    elecstate::tmp_xc_func_type = 3;
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->allocate_rho);
    charge->allocate(PARAM.input.nspin);
    EXPECT_TRUE(charge->allocate_rho);
    // test if Charge::allocate() be called twice
    EXPECT_NO_THROW(charge->allocate(PARAM.input.nspin));
    EXPECT_TRUE(charge->allocate_rho);
}

TEST_F(ChargeTest, SumRho)
{
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->allocate_rho);
    charge->allocate(PARAM.input.nspin);
    EXPECT_TRUE(charge->allocate_rho);
    int nspin = (PARAM.input.nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin; is++)
    {
        for (int ir = 0; ir < rhopw->nrxx; ir++)
        {
            charge->rho[is][ir] = 0.1;
        }
    }
    charge->set_omega(&ucell->omega);;
    EXPECT_NEAR(charge->sum_rho(), 0.1 * nspin * rhopw->nrxx * ucell->omega / rhopw->nxyz, 1E-10);
}

TEST_F(ChargeTest, RenormalizeRho)
{
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->allocate_rho);
    charge->allocate(PARAM.input.nspin);
    EXPECT_TRUE(charge->allocate_rho);
    int nspin = (PARAM.input.nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin; is++)
    {
        for (int ir = 0; ir < rhopw->nrxx; ir++)
        {
            charge->rho[is][ir] = 0.1;
        }
    }
    EXPECT_EQ(PARAM.input.nelec, 8);
    charge->set_omega(&ucell->omega);;
    charge->renormalize_rho();
    EXPECT_NEAR(charge->sum_rho(), 8.0, 1e-10);
}

TEST_F(ChargeTest, CheckNe)
{
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->allocate_rho);
    charge->allocate(PARAM.input.nspin);
    EXPECT_TRUE(charge->allocate_rho);
    int nspin = (PARAM.input.nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin; is++)
    {
        for (int ir = 0; ir < rhopw->nrxx; ir++)
        {
            charge->rho[is][ir] = 0.1;
        }
    }
    EXPECT_EQ(PARAM.input.nelec, 8);
    charge->set_omega(&ucell->omega);;
    charge->renormalize_rho();
    EXPECT_NEAR(charge->sum_rho(), 8.0, 1e-10);
    EXPECT_NEAR(charge->cal_rho2ne(charge->rho[0]), 8.0, 1e-10);
}

TEST_F(ChargeTest, SaveRhoBeforeSumBand)
{
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->allocate_rho);
    charge->allocate(PARAM.input.nspin);
    EXPECT_TRUE(charge->allocate_rho);
    int nspin = (PARAM.input.nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin; is++)
    {
        for (int ir = 0; ir < rhopw->nrxx; ir++)
        {
            charge->rho[is][ir] = 0.1;
        }
    }
    EXPECT_EQ(PARAM.input.nelec, 8);
    elecstate::tmp_xc_func_type = 3;
    charge->set_omega(&ucell->omega);;
    charge->renormalize_rho();
    charge->save_rho_before_sum_band();
    EXPECT_NEAR(charge->cal_rho2ne(charge->rho_save[0]), 8.0, 1e-10);
}

TEST_F(ChargeTest, InitFinalScf)
{
    charge->set_rhopw(rhopw);
    elecstate::tmp_xc_func_type = 1;
    PARAM.input.test_charge = 2;
    charge->init_final_scf();
    EXPECT_TRUE(charge->allocate_rho_final_scf);
}


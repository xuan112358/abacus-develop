#include "../sltk_atom_input.h"

#define private public
#include "module_parameter/parameter.h"
#undef private
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "prepare_unitcell.h"

#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
LCAO_Orbitals::LCAO_Orbitals()
{
}
LCAO_Orbitals::~LCAO_Orbitals()
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

/************************************************
 * unit test of sltk_atom_input.cpp
 * *********************************************/

/**
 * - Tested Functions:
 *  - Constructor:
 *    - the constructor of Atom_input do almost everything
 *      - (1) it finds the boundary of atomic coordinates (x_min, x_max, etc.),
 *      - (2) expands the unitcell according to the boundary condition
 *      and the searching radius, which is (2 * rcutmax_Phi + 0.01) in
 *      gamma_only calculation and (2 * (rcutmax_Phi +rcutmax_Beta) + 0.01)
 *      in multi-k calculation (see atom_arrange::set_sr_NL). And determine
 *      the number of positive and negative layers in each direction in
 *      Check_Expand_Condition()
 *      - (3) records the amount of atoms after expansion (d_amount_expand) and
 *      sets the expanded lattice grids, including their coordinates and their
 *      atomic coordinates inside in Expand_Grid()
 *      - (4) and calculate the number of unitcells in x, y, z directions
 *      in calculate_cells()
 *   - Getters:
 *     - get the values obtained in constructor, for example
 *       - Clength0() is to get the length in x direction of the expanded big cell in
 *       unit of lat0
 *       - minX() will return x_min if no expanding, or -glayerX_minus in case of
 *       expansion
 *       - getCellX() will return the number of unitcells in x direction
 *       - getCellXLength() will return search radius if no expanding, or 1 in case of
 *       expansion
 *       - getRadius() will return the search radius
 *       - getLatNow() will return lat0 in Bohr
 *       - getAmount() will return the total number of lattice grids after expansion
 */

void SetGlobalV()
{
    PARAM.input.test_grid = 0;
}

class SltkAtomInputTest : public ::testing::Test
{
  protected:
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    std::ofstream ofs;
    std::ifstream ifs;
    bool pbc = true;
    double radius = ((8 + 5.01) * 2.0 + 0.01) / 10.2;
    int test_atom_in = 0;
    std::string output;
    void SetUp()
    {
        SetGlobalV();
        ucell = utp.SetUcellInfo();
    }
    void TearDown()
    {
        delete ucell;
    }
};

using SltkAtomInputDeathTest = SltkAtomInputTest;

TEST_F(SltkAtomInputTest, Constructor)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 2;
    PARAM.input.test_grid = 1;
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    ofs.close();
    ifs.open("test.out");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("ntype = 1"));
    EXPECT_THAT(str, testing::HasSubstr("Amount(atom number) = 2"));
    EXPECT_THAT(str, testing::HasSubstr("Periodic_boundary = 1"));
    EXPECT_THAT(str, testing::HasSubstr("Searching radius(lat0) = 2.55"));
    // EXPECT_THAT(str, testing::HasSubstr("CellLength(unit: lat0) = [ 0.707107, 0.707107, 0.707107 ]"));
    EXPECT_THAT(str, testing::HasSubstr("min_tau = [ -0.75, 0, 0 ]"));
    EXPECT_THAT(str, testing::HasSubstr("max_tau = [ 0, 0.75, 0.75 ]"));
    EXPECT_THAT(str, testing::HasSubstr("glayer+ = [ 6, 6, 6 ]"));
    EXPECT_THAT(str, testing::HasSubstr("glayer- = [ 5, 5, 5 ]"));
    EXPECT_THAT(str, testing::HasSubstr("CellDim = [ 11, 11, 11 ]"));
    ifs.close();
    remove("test.out");
}

TEST_F(SltkAtomInputTest, Getters)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 2;
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    EXPECT_TRUE(Atom_inp.getExpandFlag());
    EXPECT_EQ(Atom_inp.getBoundary(), 1);
    EXPECT_NEAR(Atom_inp.getRadius(), 2.55196, 1e-5);
    EXPECT_DOUBLE_EQ(Atom_inp.minX(), -0.75);
    EXPECT_DOUBLE_EQ(Atom_inp.minY(), 0);
    EXPECT_DOUBLE_EQ(Atom_inp.minZ(), 0);
    EXPECT_EQ(Atom_inp.getCell_nX(), 11);
    EXPECT_EQ(Atom_inp.getCell_nY(), 11);
    EXPECT_EQ(Atom_inp.getCell_nZ(), 11);
    EXPECT_EQ(Atom_inp.getGrid_layerX(), 6);
    EXPECT_EQ(Atom_inp.getGrid_layerY(), 6);
    EXPECT_EQ(Atom_inp.getGrid_layerZ(), 6);
    EXPECT_EQ(Atom_inp.getGrid_layerX_minus(), 5);
    EXPECT_EQ(Atom_inp.getGrid_layerY_minus(), 5);
    EXPECT_EQ(Atom_inp.getGrid_layerZ_minus(), 5);
    ofs.close();
    remove("test.out");
}

TEST_F(SltkAtomInputDeathTest, ConstructorWarning1)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 1;
    radius = -1;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in),
                ::testing::ExitedWithCode(1),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("search radius < 0,forbidden"));
    ofs.close();
    remove("test.out");
}
/* 
TEST_F(SltkAtomInputDeathTest, ConstructorWarning2)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 1;
    ucell->atoms[0].taud[1].x = -0.25;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in),
                ::testing::ExitedWithCode(1),
                "");
    // output = testing::internal::GetCapturedStdout();
    // EXPECT_THAT(output, testing::HasSubstr("dminX<0.0"));
    ofs.close();
    remove("test.out");
} 

TEST_F(SltkAtomInputDeathTest, ConstructorWarning3)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 1;
    ucell->atoms[0].taud[1].y = -0.25;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in),
                ::testing::ExitedWithCode(1),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("dminY<0.0"));
    ofs.close();
    remove("test.out");
}

TEST_F(SltkAtomInputDeathTest, ConstructorWarning4)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 1;
    ucell->atoms[0].taud[1].z = -0.25;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in),
                ::testing::ExitedWithCode(1),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("dminZ<0.0"));
    ofs.close();
    remove("test.out");
}*/

TEST_F(SltkAtomInputTest, ConstructorNoExpand)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 1;
    PARAM.input.test_grid = 1;
    // this is a bug if radius is too small
    // because the expand_flag will be false!
    radius = 0;
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    EXPECT_FALSE(Atom_inp.getExpandFlag());
    ofs.close();
    ifs.open("test.out");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("ntype = 1"));
    EXPECT_THAT(str, testing::HasSubstr("Amount(atom number) = 2"));
    EXPECT_THAT(str, testing::HasSubstr("Periodic_boundary = 1"));
    EXPECT_THAT(str, testing::HasSubstr("Searching radius(lat0) = 0"));
    // EXPECT_THAT(str, testing::HasSubstr("CellLength(unit: lat0) = [ 0.707107, 0.707107, 0.707107 ]"));
    EXPECT_THAT(str, testing::HasSubstr("min_tau = [ -0.75, 0, 0 ]"));
    EXPECT_THAT(str, testing::HasSubstr("max_tau = [ 0, 0.75, 0.75 ]"));
    EXPECT_THAT(str, testing::HasSubstr("glayer+ = [ 2, 2, 2 ]"));
    EXPECT_THAT(str, testing::HasSubstr("glayer- = [ 0, 0, 0 ]"));
    ifs.close();
    remove("test.out");
}

TEST_F(SltkAtomInputTest, ConstructorSmallSearchRadius)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 1;
    PARAM.input.test_grid = 1;
    radius = 0.5;
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    EXPECT_TRUE(Atom_inp.getExpandFlag());
    ofs.close();
    ifs.open("test.out");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("ntype = 1"));
    EXPECT_THAT(str, testing::HasSubstr("Amount(atom number) = 2"));
    EXPECT_THAT(str, testing::HasSubstr("Periodic_boundary = 1"));
    EXPECT_THAT(str, testing::HasSubstr("Searching radius(lat0) = 0.5"));
    // EXPECT_THAT(str, testing::HasSubstr("CellLength(unit: lat0) = [ 0.707107, 0.707107, 0.707107 ]"));
    EXPECT_THAT(str, testing::HasSubstr("min_tau = [ -0.75, 0, 0 ]"));
    EXPECT_THAT(str, testing::HasSubstr("max_tau = [ 0, 0.75, 0.75 ]"));
    EXPECT_THAT(str, testing::HasSubstr("glayer+ = [ 2, 2, 2 ]"));
    EXPECT_THAT(str, testing::HasSubstr("glayer- = [ 2, 2, 2 ]"));
    EXPECT_THAT(str, testing::HasSubstr("CellDim = [ 4, 4, 4 ]"));
    ifs.close();
    remove("test.out");
}
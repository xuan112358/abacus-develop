#include "../sltk_atom_arrange.h"

#define private public
#include "module_parameter/parameter.h"
#undef private
#include <iostream>
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
 *  unit test of atom_arrange
 ***********************************************/

/**
 * - Tested Functions:
 *   - atom_arrange::delete_vector(void)
 *     - delete vector
 *   - atom_arrange::set_sr_NL
 * 	   - set the sr: search radius including nonlocal beta
 *   - filter_adjs function
 *     - filter AdjacentAtomInfo to the minimized adjacent atoms
 */

void SetGlobalV()
{
    PARAM.input.test_grid = false;
}

class SltkAtomArrangeTest : public testing::Test
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

TEST_F(SltkAtomArrangeTest, setsrNL)
{
    atom_arrange test;
    const std::string teststring = "m";
    double rcutmax_Phi = 1;
    double rcutmax_Beta = 2;
    bool gamma_only_local = true;
    double test_sr = 0;
    std::ofstream ofs;
    ofs.open("./to_test_arrange.txt");
    test_sr = test.set_sr_NL(ofs, teststring, rcutmax_Phi, rcutmax_Beta, gamma_only_local);
    EXPECT_DOUBLE_EQ(test_sr, 2.001);

    gamma_only_local = false;
    test_sr = test.set_sr_NL(ofs, teststring, rcutmax_Phi, rcutmax_Beta, gamma_only_local);
    EXPECT_DOUBLE_EQ(test_sr, 6.001);

    const std::string teststring2 = "no";
    test_sr = test.set_sr_NL(ofs, teststring2, rcutmax_Phi, rcutmax_Beta, gamma_only_local);
    ofs.close();
    std::ifstream ifs;
    std::string test2, s;
    ifs.open("./to_test_arrange.txt");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("longest orb rcut (Bohr) = 1"));
    EXPECT_THAT(str, testing::HasSubstr("longest nonlocal projector rcut (Bohr) = 2"));
    ifs.close();
    remove("./to_test_arrange");
}

TEST_F(SltkAtomArrangeTest, Search)
{
    ucell->check_dtau();
    Grid_Driver grid_d(PARAM.input.test_deconstructor, PARAM.input.test_grid);
    ofs.open("test.out");
    bool test_only = true;
    atom_arrange::search(pbc, ofs, grid_d, *ucell, radius, test_atom_in, test_only);
    EXPECT_EQ(grid_d.getType(0),0);
    EXPECT_EQ(grid_d.getNatom(0), 1); // adjacent atom is 1
    ofs.close();
    ifs.open("test.out");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("search neighboring atoms done."));
    remove("test.out");
}

TEST_F(SltkAtomArrangeTest, Filteradjs)
{
    ucell->check_dtau();
    Grid_Driver grid_d(PARAM.input.test_deconstructor, PARAM.input.test_grid);
    ofs.open("test.out");
    bool test_only = true;
    atom_arrange::search(pbc, ofs, grid_d, *ucell, radius, test_atom_in, test_only);
    EXPECT_EQ(grid_d.getType(0),0);
    EXPECT_EQ(grid_d.getNatom(0), 1); // adjacent atom is 1
    ofs.close();
    ifs.open("test.out");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("search neighboring atoms done."));
    remove("test.out");

    AdjacentAtomInfo adjs;
    grid_d.Find_atom(*ucell, ucell->atoms[0].tau[0], 0, 0, &adjs);
    EXPECT_EQ(adjs.adj_num, 0);
    // add one adjacent atom
    adjs.adj_num++;
    adjs.adjacent_tau.push_back(ModuleBase::Vector3<double>(0,0,0));
    adjs.box.push_back(ModuleBase::Vector3<int>(0,0,0));
    adjs.natom.push_back(1);
    adjs.ntype.push_back(0);
    EXPECT_EQ(adjs.adj_num, 1);
    // filter adjs to no adjacent status
    std::vector<bool> is_adjs(adjs.adj_num + 1, false);
    is_adjs[0] = true;
    filter_adjs(is_adjs, adjs);
    EXPECT_EQ(adjs.adj_num, 0);
}
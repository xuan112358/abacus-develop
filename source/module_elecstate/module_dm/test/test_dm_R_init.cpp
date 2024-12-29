#include <chrono>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_cell/klist.h"
#undef private
K_Vectors::K_Vectors()
{
}

K_Vectors::~K_Vectors()
{
}
/************************************************
 *  unit test of DensityMatrix constructor
 ***********************************************/

/**
 * This unit test construct a DensityMatrix object
 */

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 10;
int test_nw = 26;

class DMTest : public testing::Test
{
  protected:
    Parallel_Orbitals* paraV;
    int dsize;
    int my_rank = 0;
    UnitCell ucell;
    void SetUp() override
    {
#ifdef __MPI
        // MPI parallel settings
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

        // set up a unitcell, with one element and test_size atoms, each atom has test_nw orbitals
        ucell.ntype = 1;
        ucell.nat = test_size;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it = new int[ucell.nat];
        ucell.iat2ia = new int[ucell.nat];
        ucell.atoms[0].tau.resize(ucell.nat);
        ucell.itia2iat.create(ucell.ntype, ucell.nat);
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
            ucell.iat2ia[iat] = iat;
            ucell.atoms[0].tau[iat] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
            ucell.itia2iat(0, iat) = iat;
        }
        ucell.atoms[0].na = test_size;
        ucell.atoms[0].nw = test_nw;
        ucell.atoms[0].iw2l.resize(test_nw);
        ucell.atoms[0].iw2m.resize(test_nw);
        ucell.atoms[0].iw2n.resize(test_nw);
        for (int iw = 0; iw < test_nw; ++iw)
        {
            ucell.atoms[0].iw2l[iw] = 0;
            ucell.atoms[0].iw2m[iw] = 0;
            ucell.atoms[0].iw2n[iw] = 0;
        }
        ucell.set_iat2iwt(1);
        init_parav();
        // set paraV
        init_parav();
    }

    void TearDown() override
    {
        delete paraV;
        delete[] ucell.atoms;
    }

#ifdef __MPI
    void init_parav()
    {
        int nb = 2;
        int global_row = test_size * test_nw;
        int global_col = test_size * test_nw;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->init(global_row, global_col, nb, MPI_COMM_WORLD);
        paraV->set_atomic_trace(ucell.get_iat2iwt(), test_size, global_row);
    }
#else
    void init_parav()
    {
    }
#endif
};

// test for construct DMR from GridD and UnitCell
TEST_F(DMTest, DMInit1)
{
    // initalize a kvectors
    K_Vectors* kv = nullptr;
    int nspin = 1;
    int nks = 2; // since nspin = 1
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    kv->kvec_d[1].x = 0.5;
    // construct DM
    std::cout << "dim0: " << paraV->dim0 << "    dim1:" << paraV->dim1 << std::endl;
    std::cout << "nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    elecstate::DensityMatrix<double, double> DM(paraV, nspin, kv->kvec_d, nks);
    // initialize this->_DMR
    Grid_Driver gd(0,0);
    DM.init_DMR(&gd, &ucell);
    // compare
    EXPECT_EQ(DM.get_DMR_pointer(1)->size_atom_pairs(), test_size * test_size);
    EXPECT_EQ(DM.get_DMR_pointer(1)->get_atom_pair(2, 2).get_atom_i(), 2);
    EXPECT_EQ(DM.get_DMR_pointer(1)->get_atom_pair(2, 2).get_atom_j(), 2);
    EXPECT_EQ(DM.get_DMR_pointer(1)->get_atom_pair(2, 2).get_row_size(), paraV->get_row_size(2));
    EXPECT_EQ(DM.get_DMR_pointer(1)->get_atom_pair(2, 2).get_col_size(), paraV->get_col_size(2));
    delete kv;
}

// test for construct DMR from RA and UnitCell
TEST_F(DMTest, DMInit2)
{
    // initalize a kvectors
    K_Vectors* kv = nullptr;
    int nspin = 1;
    int nks = 2; // since nspin = 1
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    kv->kvec_d[1].x = 0.5;
    // construct DM
    std::cout << "dim0: " << paraV->dim0 << "    dim1:" << paraV->dim1 << std::endl;
    std::cout << "nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    elecstate::DensityMatrix<double, double> DM(paraV, nspin, kv->kvec_d, nks);
    // initialize Record_adj using Grid_Driver
    Grid_Driver gd(0,0);
    Record_adj ra;
    ra.na_each = new int[ucell.nat];
    ra.info = new int**[ucell.nat];
    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        auto tau1 = ucell.get_tau(iat1);
        int T1, I1;
        ucell.iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        gd.Find_atom(ucell, tau1, T1, I1, &adjs);
        ra.na_each[iat1] = adjs.adj_num + 1;
        ra.info[iat1] = new int*[ra.na_each[iat1]];
        for (int ad = 0; ad < ra.na_each[iat1]; ++ad)
        {
            ra.info[iat1][ad] = new int[5];
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            ra.info[iat1][ad][3] = T2;
            ra.info[iat1][ad][4] = I2;
            ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            ra.info[iat1][ad][0] = R_index.x;
            ra.info[iat1][ad][1] = R_index.y;
            ra.info[iat1][ad][2] = R_index.z;
            ra.info[iat1][ad][3] = T2;
            ra.info[iat1][ad][4] = I2;
        }
    }
    DM.init_DMR(ra, &ucell);
    // compare
    EXPECT_EQ(DM.get_DMR_pointer(1)->size_atom_pairs(), test_size * test_size);
    EXPECT_EQ(DM.get_DMR_pointer(1)->get_atom_pair(2, 2).get_atom_i(), 2);
    EXPECT_EQ(DM.get_DMR_pointer(1)->get_atom_pair(2, 2).get_atom_j(), 2);
    EXPECT_EQ(DM.get_DMR_pointer(1)->get_atom_pair(2, 2).get_row_size(), paraV->get_row_size(2));
    EXPECT_EQ(DM.get_DMR_pointer(1)->get_atom_pair(2, 2).get_col_size(), paraV->get_col_size(2));
    // release memory
    delete kv;
    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        for (int ad = 0; ad < ra.na_each[iat1]; ++ad)
        {
            delete[] ra.info[iat1][ad];
        }
        delete[] ra.info[iat1];
    }
    delete[] ra.info;
}

// test for construct DMR from another HContainer<double>
TEST_F(DMTest, DMInit3)
{
    // initalize a kvectors
    K_Vectors* kv = nullptr;
    int nspin = 2;
    int nks = 4; // since nspin = 2
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    kv->kvec_d[1].x = 0.5;
    kv->kvec_d[3].x = 0.5;
    // construct a DM
    elecstate::DensityMatrix<std::complex<double>, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    Grid_Driver gd(0, 0);
    DM.init_DMR(&gd, &ucell);
    std::cout << "dim0: " << paraV->dim0 << "    dim1:" << paraV->dim1 << std::endl;
    // construct another DM
    elecstate::DensityMatrix<std::complex<double>, double> DM1(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    DM1.init_DMR(*DM.get_DMR_pointer(1));
    // compare
    EXPECT_EQ(DM1.get_DMR_pointer(2)->size_atom_pairs(), test_size * test_size);
    EXPECT_EQ(DM1.get_DMR_pointer(2)->get_atom_pair(2, 2).get_atom_i(), 2);
    EXPECT_EQ(DM1.get_DMR_pointer(1)->get_atom_pair(2, 2).get_atom_j(), 2);
    EXPECT_EQ(DM1.get_DMR_pointer(1)->get_atom_pair(2, 2).get_row_size(), paraV->get_row_size(2));
    EXPECT_EQ(DM1.get_DMR_pointer(2)->get_atom_pair(2, 2).get_col_size(), paraV->get_col_size(2));
    //
    delete kv;
}

// test for construct DMR from another HContainer<complex<double>>
TEST_F(DMTest, DMInit4)
{
    // initalize a kvectors
    K_Vectors* kv = nullptr;
    int nspin = 2;
    int nks = 4; // since nspin = 2
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    kv->kvec_d[1].x = 0.5;
    kv->kvec_d[3].x = 0.5;
    // construct a new HContainer
    Grid_Driver gd(0, 0);
    hamilt::HContainer<std::complex<double>>* tmp_DMR;
    tmp_DMR = new hamilt::HContainer<std::complex<double>>(paraV);
    // set up a HContainer
    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        auto tau1 = ucell.get_tau(iat1);
        int T1, I1;
        ucell.iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        gd.Find_atom(ucell, tau1, T1, I1, &adjs);
        // std::cout << "adjs.adj_num: " <<adjs.adj_num << std::endl;
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            int iat2 = ucell.itia2iat(T2, I2);
            if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }
            ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            // std::cout << "R_index: " << R_index.x << " " << R_index.y << " " << R_index.z << std::endl;
            hamilt::AtomPair<std::complex<double>> tmp_ap(iat1, iat2, R_index.x, R_index.y, R_index.z, paraV);
            tmp_DMR->insert_pair(tmp_ap);
        }
    }
    // construct a DM from this HContainer
    elecstate::DensityMatrix<std::complex<double>, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    DM.init_DMR(*tmp_DMR);
    std::cout << "dim0: " << paraV->dim0 << "    dim1:" << paraV->dim1 << std::endl;
    // compare
    EXPECT_EQ(DM.get_DMR_pointer(2)->size_atom_pairs(), test_size * test_size);
    EXPECT_EQ(DM.get_DMR_pointer(2)->get_atom_pair(2, 2).get_atom_i(), 2);
    EXPECT_EQ(DM.get_DMR_pointer(1)->get_atom_pair(2, 2).get_atom_j(), 2);
    EXPECT_EQ(DM.get_DMR_pointer(1)->get_atom_pair(2, 2).get_row_size(), paraV->get_row_size(2));
    EXPECT_EQ(DM.get_DMR_pointer(2)->get_atom_pair(2, 2).get_col_size(), paraV->get_col_size(2));
    //
    delete kv;
}

// test for save_DMR
TEST_F(DMTest, saveDMR)
{
    // initalize a kvectors
    K_Vectors* kv = nullptr;
    int nspin = 2;
    int nks = 4; // since nspin = 2
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    kv->kvec_d[1].x = 0.5;
    kv->kvec_d[3].x = 0.5;
    // construct a DM
    elecstate::DensityMatrix<std::complex<double>, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    Grid_Driver gd(0, 0);
    DM.init_DMR(&gd, &ucell);
    // construct another DM
    elecstate::DensityMatrix<std::complex<double>, double> DM_test(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    DM_test.init_DMR(*DM.get_DMR_pointer(1));
    DM_test.save_DMR();
    EXPECT_EQ(DM_test.get_DMR_pointer(1)->get_nnr(), DM.get_DMR_pointer(1)->get_nnr());
    EXPECT_EQ(DM_test.get_DMR_pointer(1)->get_nnr(), DM_test._DMR_save[0].size());
    // add a new AtomPair, act as a relaxation
    hamilt::AtomPair<double> tmp_ap(9, 9, 1, 0, 0, paraV);
    DM_test.get_DMR_pointer(1)->insert_pair(tmp_ap);
    DM_test.get_DMR_pointer(1)->allocate();
    // update DMR_save
    DM_test.save_DMR();
    EXPECT_EQ(DM_test.get_DMR_pointer(1)->get_nnr(), DM_test._DMR_save[0].size());
    // delete 
    delete kv;   
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}

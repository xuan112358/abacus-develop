#include <chrono>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_cell/klist.h"

/************************************************
 *  unit test of DensityMatrix constructor
 ***********************************************/

/**
 * This unit test cal_dmk_psi
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
        ucell.atoms[0].tau = new ModuleBase::Vector3<double>[ucell.nat];
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
        ucell.atoms[0].iw2l = new int[test_nw];
        ucell.atoms[0].iw2m = new int[test_nw];
        ucell.atoms[0].iw2n = new int[test_nw];
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

    void TearDown()
    {
        delete paraV;
        delete[] ucell.atoms[0].tau;
        delete[] ucell.atoms[0].iw2l;
        delete[] ucell.atoms[0].iw2m;
        delete[] ucell.atoms[0].iw2n;
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

TEST_F(DMTest, cal_dmk_psi_nspin1)
{
    // initalize a kvectors
    K_Vectors* kv = nullptr;
    int nks = 2;
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    kv->kvec_d[1].x = 0.5;
    // construct DM
    std::cout << "dim0: " << paraV->dim0 << "    dim1:" << paraV->dim1 << std::endl;
    std::cout << "nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    int nspin = 1;
    elecstate::DensityMatrix<double, double> DM(kv, paraV, nspin);
    // compare
    EXPECT_EQ(DM.get_DMK_nks(), kv->get_nks());
    EXPECT_EQ(DM.get_DMK_nrow(), paraV->nrow);
    EXPECT_EQ(DM.get_DMK_ncol(), paraV->ncol);

    // set elements of DMK
    for (int is = 1; is <= nspin; is++)
    {
        for (int ik = 0; ik < kv->get_nks() / nspin; ik++)
        {
            for (int i = 0; i < paraV->nrow; i++)
            {
                for (int j = 0; j < paraV->ncol; j++)
                {
                    DM.set_DMK(is, ik, i, j, is + ik * i + j);
                }
            }
        }
    }
    // compare
    for (int is = 1; is <= nspin; is++)
    {
        for (int ik = 0; ik < kv->get_nks() / nspin; ik++)
        {
            for (int i = 0; i < paraV->nrow; i++)
            {
                for (int j = 0; j < paraV->ncol; j++)
                {
                    EXPECT_EQ(DM.get_DMK(is, ik, i, j), is + ik * i + j);
                }
            }
        }
    }
    // test for get_DMK_pointer
    for (int is = 1; is <= nspin; is++)
    {
        int ik_begin = (is - 1) * kv->get_nks() / nspin;
        for (int ik = 0; ik < kv->get_nks() / nspin; ik++)
        {
            double* ptr = DM.get_DMK_pointer(ik + ik_begin);
            for (int i = 0; i < paraV->nrow; i++)
            {
                for (int j = 0; j < paraV->ncol; j++)
                {
                    // std::cout << ptr[i*paraV->ncol+j] << " ";
                    EXPECT_EQ(ptr[i * paraV->ncol + j], is + ik * i + j);
                }
            }
        }
    }
    // delete kv
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

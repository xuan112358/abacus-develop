#include <chrono>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_cell/klist.h"

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
int test_nw = 10;

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
        delete[] ucell.iat2it;
        delete[] ucell.iat2ia;
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

TEST_F(DMTest, cal_DMR_full)
{
    // get my rank of this process
    int my_rank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
    // output dim and nrow, ncol
    if (my_rank == 0)
    {
        std::cout << "my rank: " << my_rank << " dim0: " << paraV->dim0 << "    dim1:" << paraV->dim1 << std::endl;
        std::cout << "my rank: " << my_rank << " nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    }
    else
    {
        std::cout << "my rank: " << my_rank << " nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    }
    // initalize a kvectors, Gamma-only
    K_Vectors* kv = nullptr;
    int nspin = 4;
    int nks = 2; // since nspin = 2
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    // construct DM
    elecstate::DensityMatrix<std::complex<double>, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks());
    // set this->_DMK
    for (int is = 1; is <= nspin; is++)
    {
        for (int ik = 0; ik < kv->get_nks(); ik++)
        {
            for (int i = 0; i < paraV->nrow; i++)
            {
                for (int j = 0; j < paraV->ncol; j++)
                {
                    DM.set_DMK(is, ik, i, j, std::complex<double>(0.77, 0.77));
                }
            }
        }
    }
    // initialize dmR_full
    hamilt::HContainer<std::complex<double>> dmR_full(ucell, paraV);
    // calculate this->_DMR
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    DM.cal_DMR_full(&dmR_full);
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "my rank: " << my_rank << " elapsed time blas: " << elapsed_time.count() << std::endl;
    // compare the result
    for (int i = 0; i < dmR_full.size_atom_pairs(); i++)
    {
        const std::complex<double>* ptr1 = dmR_full.get_atom_pair(i).get_HR_values(0, 0, 0).get_pointer();
        //
        for (int j = 0; j < dmR_full.get_atom_pair(i).get_size(); j++)
        {
            //std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j].real(), 1.54, 1e-10);
            EXPECT_NEAR(ptr1[j].imag(), 1.54, 1e-10);
        }
    }
    delete kv;
}

TEST_F(DMTest, cal_DMR_blas_double)
{
    // get my rank of this process
    int my_rank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
    // output dim and nrow, ncol
    if (my_rank == 0)
    {
        std::cout << "my rank: " << my_rank << " dim0: " << paraV->dim0 << "    dim1:" << paraV->dim1 << std::endl;
        std::cout << "my rank: " << my_rank << " nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    }
    else
    {
        std::cout << "my rank: " << my_rank << " nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    }
    // initalize a kvectors, Gamma-only
    K_Vectors* kv = nullptr;
    int nspin = 2;
    int nks = 2; // since nspin = 2
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    // construct DM
    elecstate::DensityMatrix<double, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    // set this->_DMK
    for (int is = 1; is <= nspin; is++)
    {
        for (int ik = 0; ik < kv->get_nks() / nspin; ik++)
        {
            for (int i = 0; i < paraV->nrow; i++)
            {
                for (int j = 0; j < paraV->ncol; j++)
                {
                    DM.set_DMK(is, ik, i, j, 0.77);
                }
            }
        }
    }
    // initialize this->_DMR
    Grid_Driver gd(0, 0);
    DM.init_DMR(&gd, &ucell);
    // set Gamma-only
    for (int is = 1; is <= nspin; is++)
    {
        DM.get_DMR_pointer(is)->fix_gamma();
    }
    // calculate this->_DMR
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    DM.cal_DMR();
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "my rank: " << my_rank << " elapsed time blas: " << elapsed_time.count() << std::endl;
    // compare the result
    for (int i = 0; i < DM.get_DMR_pointer(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_DMR_pointer(1)->get_atom_pair(i).get_HR_values(0, 0, 0).get_pointer();
        //
        for (int j = 0; j < DM.get_DMR_pointer(1)->get_atom_pair(i).get_size(); j++)
        {
            // std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], 0.77, 1e-10);
        }
    }
    delete kv;
}

TEST_F(DMTest, cal_DMR_blas_complex)
{
    // get my rank of this process
    int my_rank = 0;
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
    // output dim and nrow, ncol
    if (my_rank == 0)
    {
        std::cout << "my rank: " << my_rank << " dim0: " << paraV->dim0 << "    dim1:" << paraV->dim1 << std::endl;
        std::cout << "my rank: " << my_rank << " nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    }
    else
    {
        std::cout << "my rank: " << my_rank << " nrow: " << paraV->nrow << "    ncol:" << paraV->ncol << std::endl;
    }
    // initalize a kvectors
    K_Vectors* kv = nullptr;
    int nspin = 2;
    int nks = 4; // since nspin = 2
    kv = new K_Vectors;
    kv->set_nks(nks);
    kv->kvec_d.resize(nks);
    kv->kvec_d[1].x = 0.5;
    kv->kvec_d[3].x = 0.5;
    // construct DM
    elecstate::DensityMatrix<std::complex<double>, double> DM(paraV, nspin, kv->kvec_d, kv->get_nks() / nspin);
    // set this->_DMK
    for (int is = 1; is <= nspin; is++)
    {
        for (int ik = 0; ik < kv->get_nks() / nspin; ik++)
        {
            for (int i = 0; i < paraV->nrow; i++)
            {
                for (int j = 0; j < paraV->ncol; j++)
                {
                    DM.set_DMK(is, ik, i, j, is * 0.77 * (ik + 1));
                }
            }
        }
    }
    // initialize this->_DMR
    Grid_Driver gd(0, 0);
    DM.init_DMR(&gd, &ucell);
    // calculate this->_DMR
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    DM.cal_DMR();
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "my rank: " << my_rank << " elapsed time blas: " << elapsed_time.count() << std::endl;
    // compare the result for spin-up
    for (int i = 0; i < DM.get_DMR_pointer(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_DMR_pointer(1)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_DMR_pointer(1)->get_atom_pair(i).get_size(); j++)
        {
            // std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], -0.77, 1e-10);
        }
    }
    // compare the result for spin-down
    for (int i = 0; i < DM.get_DMR_pointer(2)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_DMR_pointer(2)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_DMR_pointer(2)->get_atom_pair(i).get_size(); j++)
        {
            // std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], -0.77 * 2, 1e-10);
        }
    }
    // calculate DMR_total
    DM.switch_dmr(1);
    // compare the result for spin-up after sum
    for (int i = 0; i < DM.get_DMR_pointer(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_DMR_pointer(1)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_DMR_pointer(1)->get_atom_pair(i).get_size(); j++)
        {
            //std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], -0.77 * 3, 1e-10);
        }
    }
    // restore to normal DMR 
    DM.switch_dmr(0);
    for (int i = 0; i < DM.get_DMR_pointer(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_DMR_pointer(1)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_DMR_pointer(1)->get_atom_pair(i).get_size(); j++)
        {
            //std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], -0.77, 1e-10);
        }
    }
    // calculate DMR_differenct
    DM.switch_dmr(2);
    for (int i = 0; i < DM.get_DMR_pointer(1)->size_atom_pairs(); i++)
    {
        double* ptr1 = DM.get_DMR_pointer(1)->get_atom_pair(i).get_HR_values(1, 1, 1).get_pointer();
        //
        for (int j = 0; j < DM.get_DMR_pointer(1)->get_atom_pair(i).get_size(); j++)
        {
            //std::cout << "my rank: " << my_rank << " i: " << i << " j: " << j << " value: " << ptr1[j] << std::endl;
            EXPECT_NEAR(ptr1[j], 0.77, 1e-10);
        }
    }
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

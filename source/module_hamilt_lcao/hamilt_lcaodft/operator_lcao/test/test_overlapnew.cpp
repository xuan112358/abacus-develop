#include "../overlap_new.h"

#include "gtest/gtest.h"

//---------------------------------------
// Unit test of OverlapNew class
// OverlapNew is a derivative class of Operator, it is used to calculate the overlap matrix
// It use HContainer to store the real space SR matrix
// In this test, we test the correctness and time consuming of 3 functions in OverlapNew class
// - initialize_SR() called in constructor
// - contributeHR()
// - contributeHk()
// - SR(double) and SK(complex<double>) are tested in constructHRd2cd
// - SR(double) and SK(double) are tested in constructHRd2d
//---------------------------------------

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 10;
int test_nw = 10;
class OverlapNewTest : public ::testing::Test
{
  protected:
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
        // set up a HContainer with ucell
        SR = new hamilt::HContainer<double>(paraV);
    }

    void TearDown() override
    {
        delete SR;
        delete paraV;
        delete[] ucell.atoms;
    }

#ifdef __MPI
    void init_parav()
    {
        int nb = 10;
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

    UnitCell ucell;
    hamilt::HContainer<double>* SR;
    Parallel_Orbitals* paraV;
    TwoCenterIntegrator intor_;

    int dsize;
    int my_rank = 0;
};

// using TEST_F to test OverlapNew
TEST_F(OverlapNewTest, constructHRd2d)
{
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);
    hamilt::OverlapNew<hamilt::OperatorLCAO<double, double>>
        op(&hsk, kvec_d_in, nullptr, SR, &ucell, {1.0}, &gd, &intor_);
    op.contributeHR();
    // check the value of SR
    for (int iap = 0; iap < SR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = SR->get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        auto indexes1 = paraV->get_indexes_row(iat1);
        auto indexes2 = paraV->get_indexes_col(iat2);
        int nwt = indexes1.size() * indexes2.size();
        for (int i = 0; i < nwt; ++i)
        {
            EXPECT_EQ(tmp.get_pointer(0)[i], 1.0);
        }
    }
    // calculate SK
    op.contributeHk(0);
    // check the value of SK
    double* sk = hsk.get_sk();
    for (int i = 0; i < hsk.get_size(); ++i)
    {
        EXPECT_EQ(sk[i], 1.0);
    }
}

TEST_F(OverlapNewTest, constructHRd2cd)
{
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(2, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    kvec_d_in[1] = ModuleBase::Vector3<double>(0.1, 0.2, 0.3);
    hamilt::HS_Matrix_K<std::complex<double>> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);
    hamilt::OverlapNew<hamilt::OperatorLCAO<std::complex<double>, double>>
        op(&hsk, kvec_d_in, nullptr, SR, &ucell, {1.0}, &gd, &intor_);
    op.contributeHR();
    // check the value of SR
    for (int iap = 0; iap < SR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = SR->get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        auto indexes1 = paraV->get_indexes_row(iat1);
        auto indexes2 = paraV->get_indexes_col(iat2);
        int nwt = indexes1.size() * indexes2.size();
        for (int i = 0; i < nwt; ++i)
        {
            EXPECT_EQ(tmp.get_pointer(0)[i], 1.0);
        }
    }
    // calculate SK for gamma point
    op.contributeHk(0);
    // check the value of SK of gamma point
    auto* sk = hsk.get_sk();
    for (int i = 0; i < paraV->get_row_size() * paraV->get_col_size(); ++i)
    {
        EXPECT_EQ(sk[i].real(), 1.0);
        EXPECT_EQ(sk[i].imag(), 0.0);
    }
    // calculate SK for k point
    hsk.set_zero_sk();
    op.contributeHk(1);
    // check the value of SK
    for (int i = 0; i < paraV->get_row_size() * paraV->get_col_size(); ++i)
    {
        EXPECT_NEAR(sk[i].real(), -0.80901699437494723, 1e-10);
        EXPECT_NEAR(sk[i].imag(), -0.58778525229247336, 1e-10);
    }
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

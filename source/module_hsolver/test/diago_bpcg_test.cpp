#include "module_base/inverse_matrix.h"
#include "module_base/lapack_connector.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_psi/psi.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "../diago_iter_assist.h"
#include "../diago_bpcg.h"
#include "diago_mock.h"
#include "mpi.h"
#include "module_basis/module_pw/test/test_tool.h"

#include <gtest/gtest.h>
#include <complex>
#include <random>

/************************************************
 *  unit test of functions in Diago_BPCG
 ***********************************************/

/**
 * Class Diago_BPCG is an approach for eigenvalue problems
 * This unittest test the function Diago_BPCG::diag() for FPTYPE=double, Device=cpu
 * with different examples.
 *  - the Hermite matrices (npw=500,1000) produced using random numbers and with sparsity of 0%, 60%, 80%
 *  - the Hamiltonian matrix read from "data-H", produced by using out_hs in INPUT of a LCAO calculation
 *  - a 2x2 Hermite matrix for learning and checking
 *
 * Note:
 * The test is passed when the eignvalues are closed to these calculated by LAPACK.
 * It is used together with a header file diago_mock.h.
 * The default Hermite matrix generated here is real symmetric, one can add an imaginary part
 * by changing two commented out lines in diago_mock.h.
 *
 */

// call lapack in order to compare to bpcg
void lapackEigen(int &npw, std::vector<std::complex<double>> &hm, double *e, bool outtime = false)
{
    clock_t start, end;
    start = clock();
    int lwork = 2 * npw;
    std::complex<double> *work2 = new std::complex<double>[lwork];
    double *rwork = new double[3 * npw - 2];
    int info = 0;
    char tmp_c1 = 'V', tmp_c2 = 'U';
    zheev_(&tmp_c1, &tmp_c2, &npw, hm.data(), &npw, e, work2, &lwork, rwork, &info);
    end = clock();
    if (outtime)
        std::cout << "Lapack Run time: " << (double)(end - start) / CLOCKS_PER_SEC << " S" << std::endl;
    delete[] rwork;
    delete[] work2;
}

class DiagoBPCGPrepare
{
  public:
    DiagoBPCGPrepare(int nband, int npw, int sparsity, bool reorder, double eps, int maxiter, double threshold)
        : nband(nband), npw(npw), sparsity(sparsity), reorder(reorder), eps(eps), maxiter(maxiter),
          threshold(threshold)
    {
#ifdef __MPI	
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
#endif	
    }

    int nband, npw, sparsity, maxiter, notconv;
    // eps is the convergence threshold within cg_diago
    double eps, avg_iter;
    bool reorder;
    double threshold;
    int nprocs=1, mypnum=0;
    // threshold is the comparison standard between bpcg and lapack

    void CompareEigen(double *precondition)
    {
        // calculate eigenvalues by LAPACK;
        double *e_lapack = new double[npw];
        auto ev = DIAGOTEST::hmatrix;
        if(mypnum == 0)  lapackEigen(npw, ev, e_lapack, false);
        // initial guess of psi by perturbing lapack psi
        ModuleBase::ComplexMatrix psiguess(nband, npw);
        std::default_random_engine p(1);
        std::uniform_int_distribution<unsigned> u(1, 10);
        for (int i = 0; i < nband; i++)
        {
            for (int j = 0; j < npw; j++)
            {
		        double rand = static_cast<double>(u(p))/10.;
                // psiguess(i,j) = ev(j,i)*(1+rand);
                psiguess(i, j) = ev[j * DIAGOTEST::h_nc + i] * rand;
            }
        }
        // run bpcg
	//======================================================================
        double *en = new double[npw];
        int ik = 1;
	    hamilt::Hamilt<std::complex<double>>* ha;
	    ha =new hamilt::HamiltPW<std::complex<double>>(nullptr, nullptr, nullptr, nullptr,nullptr);
	    int* ngk = new int [1];
	    //psi::Psi<std::complex<double>> psi(ngk,ik,nband,npw);
	    psi::Psi<std::complex<double>> psi;
	    psi.resize(ik,nband,npw);
	    //psi.fix_k(0);
        for (int i = 0; i < nband; i++)
        {
            for (int j = 0; j < npw; j++)
            {
	            psi(i,j)=psiguess(i,j);
	        }
	    }	

        psi::Psi<std::complex<double>> psi_local;
        double* precondition_local;
        DIAGOTEST::npw_local = new int[nprocs];
#ifdef __MPI				
	    DIAGOTEST::cal_division(DIAGOTEST::npw);
        DIAGOTEST::divide_hpsi(psi, psi_local, DIAGOTEST::hmatrix, DIAGOTEST::hmatrix_local); //will distribute psi and Hmatrix to each process
	    precondition_local = new double[DIAGOTEST::npw_local[mypnum]];
	    DIAGOTEST::divide_psi<double>(precondition,precondition_local);	
#else
	    DIAGOTEST::hmatrix_local = DIAGOTEST::hmatrix;
	    DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
	    psi_local = psi;
	    precondition_local = new double[DIAGOTEST::npw];
	    for(int i=0;i<DIAGOTEST::npw;i++) precondition_local[i] = precondition[i];
#endif
        hsolver::DiagoBPCG<std::complex<double>> bpcg(precondition_local);
        psi_local.fix_k(0);
        double start, end;
        start = MPI_Wtime();
        using T = std::complex<double>;
        const int dim = DIAGOTEST::npw;
        const std::vector<T> &h_mat = DIAGOTEST::hmatrix_local;
        auto hpsi_func = [h_mat, dim](T *psi_in, T *hpsi_out,
                                const int ld_psi, const int nvec) {
            auto one = std::make_unique<T>(1.0);
            auto zero = std::make_unique<T>(0.0);
            const T *one_ = one.get();
            const T *zero_ = zero.get();

            base_device::DEVICE_CPU *ctx = {};
            // hpsi_out(dim * nvec) = h_mat(dim * dim) * psi_in(dim * nvec)
            hsolver::gemm_op<T, base_device::DEVICE_CPU>()(
                ctx, 'N', 'N',
                dim, nvec, dim,
                one_,
                h_mat.data(), dim,
                psi_in, ld_psi,
                zero_,
                hpsi_out, ld_psi);
        };
        bpcg.init_iter(nband, npw);
        std::vector<double> ethr_band(nband, 1e-5);
        bpcg.diag(hpsi_func, psi_local.get_pointer(), en, ethr_band);
        bpcg.diag(hpsi_func, psi_local.get_pointer(), en, ethr_band);
        bpcg.diag(hpsi_func, psi_local.get_pointer(), en, ethr_band);
        bpcg.diag(hpsi_func, psi_local.get_pointer(), en, ethr_band);
        end = MPI_Wtime();
        //if(mypnum == 0) printf("diago time:%7.3f\n",end-start);
        delete [] DIAGOTEST::npw_local;
	    delete [] precondition_local;
	    //======================================================================
        for (int i = 0; i < nband; i++)
        {
            EXPECT_NEAR(en[i], e_lapack[i], threshold);
        }

        delete[] en;
        delete[] e_lapack;
        delete ha;
    }
};

class DiagoBPCGTest : public ::testing::TestWithParam<DiagoBPCGPrepare>
{
};

TEST_P(DiagoBPCGTest, RandomHamilt)
{
    DiagoBPCGPrepare dcp = GetParam();
    //std::cout << "npw=" << dcp.npw << ", nband=" << dcp.nband << ", sparsity="
    //		  << dcp.sparsity << ", eps=" << dcp.eps << std::endl;
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = dcp.eps;
    //std::cout<<"maxiter "<<hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX<<std::endl;
    //std::cout<<"eps "<<hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR<<std::endl;
    HPsi<std::complex<double>> hpsi(dcp.nband, dcp.npw, dcp.sparsity);
    DIAGOTEST::hmatrix = hpsi.hamilt();

    DIAGOTEST::npw = dcp.npw;
    // ModuleBase::ComplexMatrix psi = hpsi.psi();
    dcp.CompareEigen(hpsi.precond());
}

INSTANTIATE_TEST_SUITE_P(VerifyCG,
                         DiagoBPCGTest,
                         ::testing::Values(
                             // nband, npw, sparsity, reorder, eps, maxiter, threshold
                             DiagoBPCGPrepare(10, 500, 0, true, 1e-5, 300, 5e-2),
                             DiagoBPCGPrepare(20, 500, 6, true, 1e-5, 300, 5e-2),
                             DiagoBPCGPrepare(20, 1000, 8, true, 1e-5, 300, 5e-2),
                             DiagoBPCGPrepare(40, 1000, 8, true, 1e-6, 300, 5e-2))); 
                            //DiagoBPCGPrepare(40, 2000, 8, true, 1e-5, 500, 1e-2))); 
			    // the last one is passed but time-consumming.

// check that the mock class HPsi work well
// in generating a Hermite matrix
TEST(DiagoBPCGTest, Hamilt)
{
    int dim = 2;
    int nbnd = 2;
    HPsi<std::complex<double>> hpsi(nbnd, dim);
    std::vector<std::complex<double>> hm = hpsi.hamilt();
    EXPECT_EQ(DIAGOTEST::h_nr, 2);
    EXPECT_EQ(DIAGOTEST::h_nc, 2);
    EXPECT_EQ(hm[0].imag(), 0.0);
    EXPECT_EQ(hm[DIAGOTEST::h_nc + 1].imag(), 0.0);
    EXPECT_EQ(conj(hm[DIAGOTEST::h_nc]).real(), hm[1].real());
    EXPECT_EQ(conj(hm[DIAGOTEST::h_nc]).imag(), hm[1].imag());
}

// check that lapack work well
// for an eigenvalue problem
/*TEST(DiagoBPCGTest, ZHEEV)
{
    int dim = 100;
    int nbnd = 2;
    HPsi hpsi(nbnd, dim);
    std::vector<std::complex<double>> hm = hpsi.hamilt();
    std::vector<std::complex<double>> hm_backup = hm;
    ModuleBase::ComplexMatrix eig(dim, dim);
    double e[dim];
    // using zheev to do a direct test
    lapackEigen(dim, hm, e);
    eig = transpose(hm, true) * hm_backup * hm;
    // for (int i=0;i<dim;i++) std::cout<< " e[i] "<<e[i]<<std::endl;
    for (int i = 0; i < dim; i++)
    {
        EXPECT_NEAR(e[i], eig(i, i).real(), 1e-10);
    }
}*/


TEST(DiagoBPCGTest, readH)
{
    // read Hamilt matrix from file data-H
    std::vector<std::complex<double>> hm;
    std::ifstream ifs;
    ifs.open("H-KPoints-Si64.dat");
    DIAGOTEST::readh(ifs, hm);
    ifs.close();
    int dim = DIAGOTEST::npw;
    int nband = 10; // not nband < dim, here dim = 26 in data-H
    // nband, npw, sub, sparsity, reorder, eps, maxiter, threshold
    DiagoBPCGPrepare dcp(nband, dim, 0, true, 1e-5, 500, 1e-1);
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR = dcp.eps;
    hsolver::DiagoIterAssist<std::complex<double>>::SCF_ITER = 1;
    HPsi<std::complex<double>> hpsi;
    hpsi.create(nband, dim);
    DIAGOTEST::hmatrix = hpsi.hamilt();
    DIAGOTEST::npw = dim;
    dcp.CompareEigen(hpsi.precond());
}

int main(int argc, char **argv)
{
	int nproc = 1, myrank = 0;

#ifdef __MPI
	int nproc_in_pool, kpar=1, mypool, rank_in_pool;
    setupmpi(argc,argv,nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
    GlobalV::NPROC_IN_POOL = nproc;
#else
	MPI_Init(&argc, &argv);	
#endif

    testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners &listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (myrank != 0) delete listeners.Release(listeners.default_result_printer());

    int result = RUN_ALL_TESTS();
    if (myrank == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
        return result;
	}

    MPI_Finalize();
	return 0;
}

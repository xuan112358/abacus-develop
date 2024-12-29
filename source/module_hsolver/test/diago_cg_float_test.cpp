#include "gtest/gtest.h"
#define private public
#include "module_parameter/parameter.h"
#undef private
#include "module_base/inverse_matrix.h"
#include "module_base/lapack_connector.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_psi/psi.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "../diago_cg.h"
#include "../diago_iter_assist.h"
#include "diago_mock.h"
#include "mpi.h"
#include "module_basis/module_pw/test/test_tool.h"
#include <complex>

#include <random>

#include <ATen/core/tensor_map.h>

/************************************************
 *  unit test of functions in Diago_CG
 ***********************************************/

/**
 * Class Diago_CG is an approach for eigenvalue problems
 * This unittest test the function Diago_CG::diag() for FPTYPE=float, Device=cpu
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

// call lapack in order to compare to cg
void lapackEigen(int &npw, std::vector<std::complex<float>> &hm, float *e, bool outtime = false)
{
    clock_t start, end;
    start = clock();
    int lwork = 2 * npw;
    std::complex<float> *work2 = new std::complex<float>[lwork];
    float *rwork = new float[3 * npw - 2];
    int info = 0;
    char tmp_c1 = 'V', tmp_c2 = 'U';
    cheev_(&tmp_c1, &tmp_c2, &npw, hm.data(), &npw, e, work2, &lwork, rwork, &info);
    end = clock();
    if (outtime) {
        std::cout << "Lapack Run time: " << (float)(end - start) / CLOCKS_PER_SEC << " S" << std::endl;
}
    delete[] rwork;
    delete[] work2;
}

class DiagoCGPrepare
{
  public:
    DiagoCGPrepare(int nband, int npw, int sparsity, bool reorder, float eps, int maxiter, float threshold)
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
    float eps, avg_iter;
    bool reorder;
    float threshold;
    int nprocs=1, mypnum=0;
    // threshold is the comparison standard between cg and lapack

    void CompareEigen(float *precondition)
    {
        // calculate eigenvalues by LAPACK;
        float *e_lapack = new float[npw];
        auto ev = DIAGOTEST::hmatrix_f;

        if(mypnum == 0) {  lapackEigen(npw, ev, e_lapack, false);
}

        // initial guess of psi by perturbing lapack psi
        std::vector<std::complex<float>> psiguess(nband * npw);
        std::default_random_engine p(1);
        std::uniform_int_distribution<unsigned> u(1, 10);

        for (int i = 0; i < nband; i++)
        {
            for (int j = 0; j < npw; j++)
            {
		        float rand = static_cast<float>(u(p))/10.;
                // psiguess(i,j) = ev(j,i)*(1+rand);
                psiguess[i * npw + j] = ev[j * DIAGOTEST::h_nc + i] * rand;
            }
        }

        // run cg
	//======================================================================
        float *en = new float[npw];
        int ik = 1;
	    hamilt::Hamilt<std::complex<float>>* ha;
	    ha =new hamilt::HamiltPW<std::complex<float>>(nullptr, nullptr, nullptr, nullptr,nullptr);
	    int* ngk = new int [1];
	    //psi::Psi<std::complex<float>> psi(ngk,ik,nband,npw);
	    psi::Psi<std::complex<float>> psi;
	    psi.resize(ik,nband,npw);
	    //psi.fix_k(0);
        for (int i = 0; i < nband; i++)
        {
            for (int j = 0; j < npw; j++)
            {
	            psi(i,j)=psiguess[i * npw + j];
	        }
	    }	

        psi::Psi<std::complex<float>> psi_local;
        float* precondition_local;
        DIAGOTEST::npw_local = new int[nprocs];
#ifdef __MPI				
	    DIAGOTEST::cal_division(DIAGOTEST::npw);
        DIAGOTEST::divide_hpsi(psi, psi_local, DIAGOTEST::hmatrix_f, DIAGOTEST::hmatrix_local_f); //will distribute psi and Hmatrix to each process
	    precondition_local = new float[DIAGOTEST::npw_local[mypnum]];
        DIAGOTEST::divide_psi<float>(precondition, precondition_local);
#else
	    DIAGOTEST::hmatrix_local_f = DIAGOTEST::hmatrix_f;
	    DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
	    psi_local = psi;
	    precondition_local = new float[DIAGOTEST::npw];
	    for(int i=0;i<DIAGOTEST::npw;i++) precondition_local[i] = precondition[i];
#endif
        // hsolver::DiagoCG<std::complex<float>> cg(precondition_local);
        psi_local.fix_k(0);
        // cg.diag(ha,psi_local,en); 
        /**************************************************************/
        //  New interface of cg method
        /**************************************************************/
        // warp the subspace_func into a lambda function
        auto subspace_func = [ha](const ct::Tensor& psi_in, ct::Tensor& psi_out) { /*do nothing*/ };
        hsolver::DiagoCG<std::complex<float>> cg(
            PARAM.input.basis_type,
            PARAM.input.calculation,
            hsolver::DiagoIterAssist<std::complex<float>>::need_subspace,
            subspace_func,
            hsolver::DiagoIterAssist<std::complex<float>>::PW_DIAG_THR,
            hsolver::DiagoIterAssist<std::complex<float>>::PW_DIAG_NMAX,
            GlobalV::NPROC_IN_POOL);
        // hsolver::DiagoCG<std::complex<float>> cg(precondition_local);
        psi_local.fix_k(0);
        float start, end;
        start = MPI_Wtime();

        auto hpsi_func = [ha](const ct::Tensor& psi_in, ct::Tensor& hpsi_out) {
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");
            auto psi_wrapper = psi::Psi<std::complex<float>>(
                psi_in.data<std::complex<float>>(), 1, 
                ndim == 1 ? 1 : psi_in.shape().dim_size(0), 
                ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1), true);
            psi::Range all_bands_range(true, psi_wrapper.get_current_k(), 0, psi_wrapper.get_nbands() - 1);
            using hpsi_info = typename hamilt::Operator<std::complex<float>>::hpsi_info;
            hpsi_info info(&psi_wrapper, all_bands_range, hpsi_out.data<std::complex<float>>());
            ha->ops->hPsi(info);
        };
        auto spsi_func = [ha](const ct::Tensor& psi_in, ct::Tensor& spsi_out) {
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");
            ha->sPsi(psi_in.data<std::complex<float>>(), spsi_out.data<std::complex<float>>(), 
                ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1), 
                ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1), 
                ndim == 1 ? 1 : psi_in.shape().dim_size(0));
        };
        auto psi_tensor = ct::TensorMap(
            psi_local.get_pointer(), 
            ct::DataType::DT_COMPLEX, 
            ct::DeviceType::CpuDevice,
            ct::TensorShape({psi_local.get_nbands(), psi_local.get_nbasis()})).slice({0, 0}, {psi_local.get_nbands(), psi_local.get_current_nbas()});
        auto eigen_tensor = ct::TensorMap(
            en,
            ct::DataType::DT_FLOAT,
            ct::DeviceType::CpuDevice,
            ct::TensorShape({psi_local.get_nbands()}));
        auto prec_tensor = ct::TensorMap(
            precondition_local,
            ct::DataType::DT_FLOAT, 
            ct::DeviceType::CpuDevice,
            ct::TensorShape({static_cast<int>(psi_local.get_current_nbas())})).slice({0}, {psi_local.get_current_nbas()});

        std::vector<double> ethr_band(nband, 1e-5);
        cg.diag(hpsi_func, spsi_func, psi_tensor, eigen_tensor, ethr_band, prec_tensor);
        // TODO: Double check tensormap's potential problem
        ct::TensorMap(psi_local.get_pointer(), psi_tensor, {psi_local.get_nbands(), psi_local.get_nbasis()}).sync(psi_tensor);
        /**************************************************************/

        end = MPI_Wtime();
        //if(mypnum == 0) printf("diago time:%7.3f\n",end-start);
        delete [] DIAGOTEST::npw_local;
	    delete [] precondition_local;
	    //======================================================================
        for (int i = 0; i < nband; i++)
        {
            EXPECT_NEAR(en[i], e_lapack[i], threshold)<<i;
        }

        delete[] en;
        delete[] e_lapack;
        delete ha;
    }
};

class DiagoCGFloatTest : public ::testing::TestWithParam<DiagoCGPrepare>
{
};

TEST_P(DiagoCGFloatTest, RandomHamilt)
{
    DiagoCGPrepare dcp = GetParam();
    //std::cout << "npw=" << dcp.npw << ", nband=" << dcp.nband << ", sparsity="
    //		  << dcp.sparsity << ", eps=" << dcp.eps << std::endl;
    hsolver::DiagoIterAssist<std::complex<float>>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<std::complex<float>>::PW_DIAG_THR = dcp.eps;
    //std::cout<<"maxiter "<<hsolver::DiagoIterAssist<std::complex<float>>::PW_DIAG_NMAX<<std::endl;
    //std::cout<<"eps "<<hsolver::DiagoIterAssist<std::complex<float>>::PW_DIAG_THR<<std::endl;
    HPsi<std::complex<float>> hpsi(dcp.nband, dcp.npw, dcp.sparsity);
    DIAGOTEST::hmatrix_f = hpsi.hamilt();

    DIAGOTEST::npw = dcp.npw;
    // ModuleBase::ComplexMatrix psi = hpsi.psi();
    dcp.CompareEigen(hpsi.precond());
}

INSTANTIATE_TEST_SUITE_P(VerifyCG,
                         DiagoCGFloatTest,
                         ::testing::Values(
                             // nband, npw, sparsity, reorder, eps, maxiter, threshold
                             DiagoCGPrepare(10, 200, 0, true, 1e-6, 300, 1e-0),
                             DiagoCGPrepare(10, 200, 6, true, 1e-6, 300, 1e-0),
                             DiagoCGPrepare(10, 400, 8, true, 1e-6, 300, 1e-0),
                             DiagoCGPrepare(10, 600, 8, true, 1e-6, 300, 1e-0))); 
                            //DiagoCGPrepare(40, 2000, 8, true, 1e-5, 500, 1e-2))); 
			    // the last one is passed but time-consumming.

// check that the mock class HPsi work well
// in generating a Hermite matrix
TEST(DiagoCGFloatTest, Hamilt)
{
    int dim = 2;
    int nbnd = 2;
    HPsi<std::complex<float>> hpsi(nbnd, dim);
    std::vector<std::complex<float>> hm = hpsi.hamilt();
    EXPECT_EQ(DIAGOTEST::h_nr, 2);
    EXPECT_EQ(DIAGOTEST::h_nc, 2);
    EXPECT_EQ(hm[0].imag(), 0.0);
    EXPECT_EQ(hm[DIAGOTEST::h_nc + 1].imag(), 0.0);
    EXPECT_EQ(conj(hm[DIAGOTEST::h_nc]).real(), hm[1].real());
    EXPECT_EQ(conj(hm[DIAGOTEST::h_nc]).imag(), hm[1].imag());
}

// check that lapack work well
// for an eigenvalue problem
/*TEST(DiagoCGFloatTest, ZHEEV)
{
    int dim = 100;
    int nbnd = 2;
    HPsi hpsi(nbnd, dim);
    std::vector<std::complex<float>> hm = hpsi.hamilt();
    std::vector<std::complex<float>> hm_backup = hm;
    ModuleBase::ComplexMatrix eig(dim, dim);
    float e[dim];
    // using zheev to do a direct test
    lapackEigen(dim, hm, e);
    eig = transpose(hm, true) * hm_backup * hm;
    // for (int i=0;i<dim;i++) std::cout<< " e[i] "<<e[i]<<std::endl;
    for (int i = 0; i < dim; i++)
    {
        EXPECT_NEAR(e[i], eig(i, i).real(), 1e-10);
    }
}*/

// cg for a 2x2 matrix
#ifdef __MPI
#else
TEST(DiagoCGFloatTest, TwoByTwo)
{
    int dim = 2;
    int nband = 2;
    ModuleBase::ComplexMatrix hm(2, 2);
    hm(0, 0) = std::complex<float>{4.0, 0.0};
    hm(0, 1) = std::complex<float>{1.0, 0.0};
    hm(1, 0) = std::complex<float>{1.0, 0.0};
    hm(1, 1) = std::complex<float>{3.0, 0.0};
    // nband, npw, sub, sparsity, reorder, eps, maxiter, threshold
    DiagoCGPrepare dcp(nband, dim, 0, true, 1e-4, 50, 1e-0);
    hsolver::DiagoIterAssist<std::complex<float>>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<std::complex<float>>::PW_DIAG_THR = dcp.eps;
    HPsi hpsi;
    hpsi.create(nband, dim);
    DIAGOTEST::hmatrix_f = hm;
    DIAGOTEST::npw = dim;
    dcp.CompareEigen(hpsi.precond());
}
#endif

TEST(DiagoCGFloatTest, readH)
{
    // read Hamilt matrix from file data-H
    std::vector<std::complex<float>> hm;
    std::ifstream ifs;
    ifs.open("H-KPoints-Si64.dat");
    DIAGOTEST::readh(ifs, hm);
    ifs.close();
    int dim = DIAGOTEST::npw;
    int nband = 10; // not nband < dim, here dim = 26 in data-H
    // nband, npw, sparsity, reorder, eps, maxiter, threshold
    DiagoCGPrepare dcp(nband, dim, 0, true, 1e-5, 500, 1e-0);
    hsolver::DiagoIterAssist<std::complex<float>>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<std::complex<float>>::PW_DIAG_THR = dcp.eps;
    HPsi<std::complex<float>> hpsi;
    hpsi.create(nband, dim);
    DIAGOTEST::hmatrix_f = hpsi.hamilt();
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
    if (myrank != 0) { delete listeners.Release(listeners.default_result_printer());
}

    int result = RUN_ALL_TESTS();
    if (myrank == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
        return result;
	}

    MPI_Finalize();
	return 0;
}

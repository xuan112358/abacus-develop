#include "blas_connector.h"

#ifdef __DSP
#include "module_base/kernels/dsp/dsp_connector.h"
#include "module_base/global_variable.h"
#endif

void BlasConnector::axpy( const int n, const float alpha, const float *X, const int incX, float *Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		saxpy_(&n, &alpha, X, &incX, Y, &incY);
}
}

void BlasConnector::axpy( const int n, const double alpha, const double *X, const int incX, double *Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		daxpy_(&n, &alpha, X, &incX, Y, &incY);
}
}

void BlasConnector::axpy( const int n, const std::complex<float> alpha, const std::complex<float> *X, const int incX, std::complex<float> *Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		caxpy_(&n, &alpha, X, &incX, Y, &incY);
}
}

void BlasConnector::axpy( const int n, const std::complex<double> alpha, const std::complex<double> *X, const int incX, std::complex<double> *Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		zaxpy_(&n, &alpha, X, &incX, Y, &incY);
}
}


// x=a*x
void BlasConnector::scal( const int n,  const float alpha, float *X, const int incX, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		sscal_(&n, &alpha, X, &incX);
}
}

void BlasConnector::scal( const int n, const double alpha, double *X, const int incX, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		dscal_(&n, &alpha, X, &incX);
}
}

void BlasConnector::scal( const int n, const std::complex<float> alpha, std::complex<float> *X, const int incX, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		cscal_(&n, &alpha, X, &incX);
}
}

void BlasConnector::scal( const int n, const std::complex<double> alpha, std::complex<double> *X, const int incX, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		zscal_(&n, &alpha, X, &incX);
}
}


// d=x*y
float BlasConnector::dot( const int n, const float *X, const int incX, const float *Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		return sdot_(&n, X, &incX, Y, &incY);
	}
	return sdot_(&n, X, &incX, Y, &incY);
}

double BlasConnector::dot( const int n, const double *X, const int incX, const double *Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		return ddot_(&n, X, &incX, Y, &incY);
	}
	return ddot_(&n, X, &incX, Y, &incY);
}

// C = a * A.? * B.? + b * C
// Row-Major part
void BlasConnector::gemm(const char transa, const char transb, const int m, const int n, const int k,
	const float alpha, const float *a, const int lda, const float *b, const int ldb,
	const float beta, float *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		sgemm_(&transb, &transa, &n, &m, &k,
		&alpha, b, &ldb, a, &lda,
		&beta, c, &ldc);
	}
	#ifdef __DSP
	else if (device_type == base_device::AbacusDevice_t::DspDevice){
		sgemm_mth_(&transb, &transa, &n, &m, &k,
		&alpha, b, &ldb, a, &lda,
		&beta, c, &ldc, GlobalV::MY_RANK);
	}
	#endif
}

void BlasConnector::gemm(const char transa, const char transb, const int m, const int n, const int k,
	const double alpha, const double *a, const int lda, const double *b, const int ldb,
	const double beta, double *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		dgemm_(&transb, &transa, &n, &m, &k,
		&alpha, b, &ldb, a, &lda,
		&beta, c, &ldc);
	}
	#ifdef __DSP
	else if (device_type == base_device::AbacusDevice_t::DspDevice){
		dgemm_mth_(&transb, &transa, &n, &m, &k,
		&alpha, b, &ldb, a, &lda,
		&beta, c, &ldc, GlobalV::MY_RANK);
	}
	#endif
}

void BlasConnector::gemm(const char transa, const char transb, const int m, const int n, const int k,
    const std::complex<float> alpha, const std::complex<float> *a, const int lda, const std::complex<float> *b, const int ldb,
    const std::complex<float> beta, std::complex<float> *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
    	cgemm_(&transb, &transa, &n, &m, &k,
        &alpha, b, &ldb, a, &lda,
        &beta, c, &ldc);
	}
	#ifdef __DSP
	else if (device_type == base_device::AbacusDevice_t::DspDevice) {
    	cgemm_mth_(&transb, &transa, &n, &m, &k,
        &alpha, b, &ldb, a, &lda,
        &beta, c, &ldc, GlobalV::MY_RANK);
	}
	#endif
}

void BlasConnector::gemm(const char transa, const char transb, const int m, const int n, const int k,
	const std::complex<double> alpha, const std::complex<double> *a, const int lda, const std::complex<double> *b, const int ldb,
	const std::complex<double> beta, std::complex<double> *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		zgemm_(&transb, &transa, &n, &m, &k,
		&alpha, b, &ldb, a, &lda,
		&beta, c, &ldc);
	}
	#ifdef __DSP
	else if (device_type == base_device::AbacusDevice_t::DspDevice) {
    	zgemm_mth_(&transb, &transa, &n, &m, &k,
        &alpha, b, &ldb, a, &lda,
        &beta, c, &ldc, GlobalV::MY_RANK);
	}
	#endif
}

// Col-Major part
void BlasConnector::gemm_cm(const char transa, const char transb, const int m, const int n, const int k,
	const float alpha, const float *a, const int lda, const float *b, const int ldb,
	const float beta, float *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		sgemm_(&transa, &transb, &m, &n, &k,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
	#ifdef __DSP
	else if (device_type == base_device::AbacusDevice_t::DspDevice){
		sgemm_mth_(&transb, &transa, &m, &n, &k,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc, GlobalV::MY_RANK);
	}
	#endif
}

void BlasConnector::gemm_cm(const char transa, const char transb, const int m, const int n, const int k,
	const double alpha, const double *a, const int lda, const double *b, const int ldb,
	const double beta, double *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		dgemm_(&transa, &transb, &m, &n, &k,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
	#ifdef __DSP
	else if (device_type == base_device::AbacusDevice_t::DspDevice){
		dgemm_mth_(&transa, &transb, &m, &n, &k,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc, GlobalV::MY_RANK);
	}
	#endif
}

void BlasConnector::gemm_cm(const char transa, const char transb, const int m, const int n, const int k,
    const std::complex<float> alpha, const std::complex<float> *a, const int lda, const std::complex<float> *b, const int ldb,
    const std::complex<float> beta, std::complex<float> *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
    	cgemm_(&transa, &transb, &m, &n, &k,
        &alpha, a, &lda, b, &ldb,
        &beta, c, &ldc);
	}
	#ifdef __DSP
	else if (device_type == base_device::AbacusDevice_t::DspDevice) {
    	cgemm_mth_(&transa, &transb, &m, &n, &k,
        &alpha, a, &lda, b, &ldb,
        &beta, c, &ldc, GlobalV::MY_RANK);
	}
	#endif
}

void BlasConnector::gemm_cm(const char transa, const char transb, const int m, const int n, const int k,
	const std::complex<double> alpha, const std::complex<double> *a, const int lda, const std::complex<double> *b, const int ldb,
	const std::complex<double> beta, std::complex<double> *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		zgemm_(&transa, &transb, &m, &n, &k,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
	#ifdef __DSP
	else if (device_type == base_device::AbacusDevice_t::DspDevice) {
    	zgemm_mth_(&transa, &transb, &m, &n, &k,
        &alpha, a, &lda, b, &ldb,
        &beta, c, &ldc, GlobalV::MY_RANK);
	}
	#endif
}

// Symm and Hemm part. Only col-major is supported.

void BlasConnector::symm_cm(const char side, const char uplo, const int m, const int n,
	const float alpha, const float *a, const int lda, const float *b, const int ldb,
	const float beta, float *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		ssymm_(&side, &uplo, &m, &n,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
}

void BlasConnector::symm_cm(const char side, const char uplo, const int m, const int n,
	const double alpha, const double *a, const int lda, const double *b, const int ldb,
	const double beta, double *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		dsymm_(&side, &uplo, &m, &n,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
}

void BlasConnector::symm_cm(const char side, const char uplo, const int m, const int n,
    const std::complex<float> alpha, const std::complex<float> *a, const int lda, const std::complex<float> *b, const int ldb,
    const std::complex<float> beta, std::complex<float> *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
    	csymm_(&side, &uplo, &m, &n,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
}

void BlasConnector::symm_cm(const char side, const char uplo, const int m, const int n,
	const std::complex<double> alpha, const std::complex<double> *a, const int lda, const std::complex<double> *b, const int ldb,
	const std::complex<double> beta, std::complex<double> *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		zsymm_(&side, &uplo, &m, &n,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
}

void BlasConnector::hemm_cm(char side, char uplo, int m, int n,
    std::complex<float> alpha, std::complex<float> *a, int lda, std::complex<float> *b, int ldb,
    std::complex<float> beta, std::complex<float> *c, int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
    	chemm_(&side, &uplo, &m, &n,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
}

void BlasConnector::hemm_cm(char side, char uplo, int m, int n,
	std::complex<double> alpha, std::complex<double> *a, int lda, std::complex<double> *b, int ldb,
	std::complex<double> beta, std::complex<double> *c, int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		zhemm_(&side, &uplo, &m, &n,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
}

void BlasConnector::gemv(const char trans, const int m, const int n,
    const float alpha, const float* A, const int lda, const float* X, const int incx,
    const float beta, float* Y, const int incy, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
    	sgemv_(&trans, &m, &n, &alpha, A, &lda, X, &incx, &beta, Y, &incy);
}
}

void BlasConnector::gemv(const char trans, const int m, const int n,
    const double alpha, const double* A, const int lda, const double* X, const int incx,
    const double beta, double* Y, const int incy, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
    	dgemv_(&trans, &m, &n, &alpha, A, &lda, X, &incx, &beta, Y, &incy);
}
}

void BlasConnector::gemv(const char trans, const int m, const int n,
    const std::complex<float> alpha, const std::complex<float> *A, const int lda, const std::complex<float> *X, const int incx,
    const std::complex<float> beta, std::complex<float> *Y, const int incy, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
    	cgemv_(&trans, &m, &n, &alpha, A, &lda, X, &incx, &beta, Y, &incy);
}
}

void BlasConnector::gemv(const char trans, const int m, const int n,
    const std::complex<double> alpha, const std::complex<double> *A, const int lda, const std::complex<double> *X, const int incx,
    const std::complex<double> beta, std::complex<double> *Y, const int incy, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
    	zgemv_(&trans, &m, &n, &alpha, A, &lda, X, &incx, &beta, Y, &incy);
}
}

// out = ||x||_2
float BlasConnector::nrm2( const int n, const float *X, const int incX, base_device::AbacusDevice_t device_type )
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		return snrm2_( &n, X, &incX );
	}
	return snrm2_( &n, X, &incX );
}


double BlasConnector::nrm2( const int n, const double *X, const int incX, base_device::AbacusDevice_t device_type )
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		return dnrm2_( &n, X, &incX );
	}
	return dnrm2_( &n, X, &incX );
}


double BlasConnector::nrm2( const int n, const std::complex<double> *X, const int incX, base_device::AbacusDevice_t device_type )
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		return dznrm2_( &n, X, &incX );
	}
	return dznrm2_( &n, X, &incX );
}

// copies a into b
void BlasConnector::copy(const long n, const double *a, const int incx, double *b, const int incy, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		dcopy_(&n, a, &incx, b, &incy);
	}
}

void BlasConnector::copy(const long n, const std::complex<double> *a, const int incx, std::complex<double> *b, const int incy, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		zcopy_(&n, a, &incx, b, &incy);
	}
}
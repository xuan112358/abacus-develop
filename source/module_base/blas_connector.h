#ifndef BLAS_CONNECTOR_H
#define BLAS_CONNECTOR_H

#include <complex>
#include "module_base/module_device/types.h"

// These still need to be linked in the header file
// Because quite a lot of code will directly use the original cblas kernels.

extern "C"
{
	// level 1: std::vector-std::vector operations, O(n) data and O(n) work.

	// Peize Lin add ?scal 2016-08-04, to compute x=a*x
	void sscal_(const int *N, const float *alpha, float *X, const int *incX);
	void dscal_(const int *N, const double *alpha, double *X, const int *incX);
	void cscal_(const int *N, const std::complex<float> *alpha, std::complex<float> *X, const int *incX);
	void zscal_(const int *N, const std::complex<double> *alpha, std::complex<double> *X, const int *incX);

	// Peize Lin add ?axpy 2016-08-04, to compute y=a*x+y
	void saxpy_(const int *N, const float *alpha, const float *X, const int *incX, float *Y, const int *incY);
	void daxpy_(const int *N, const double *alpha, const double *X, const int *incX, double *Y, const int *incY);
	void caxpy_(const int *N, const std::complex<float> *alpha, const std::complex<float> *X, const int *incX, std::complex<float> *Y, const int *incY);
	void zaxpy_(const int *N, const std::complex<double> *alpha, const std::complex<double> *X, const int *incX, std::complex<double> *Y, const int *incY);

	void dcopy_(long const *n, const double *a, int const *incx, double *b, int const *incy);
	void zcopy_(long const *n, const std::complex<double> *a, int const *incx, std::complex<double> *b, int const *incy);

	//reason for passing results as argument instead of returning it:
	//see https://www.numbercrunch.de/blog/2014/07/lost-in-translation/
	// void zdotc_(std::complex<double> *result, const int *n, const std::complex<double> *zx,
	// 	const int *incx, const std::complex<double> *zy, const int *incy);
	// Peize Lin add ?dot 2017-10-27, to compute d=x*y
	float sdot_(const int *N, const float *X, const int *incX, const float *Y, const int *incY);
	double ddot_(const int *N, const double *X, const int *incX, const double *Y, const int *incY);

	// Peize Lin add ?nrm2 2018-06-12, to compute out = ||x||_2 = \sqrt{ \sum_i x_i**2 }
	float snrm2_( const int *n, const float *X, const int *incX );
	double dnrm2_( const int *n, const double *X, const int *incX );
	double dznrm2_( const int *n, const std::complex<double> *X, const int *incX );

    // symmetric rank-k update
    void dsyrk_(
        const char* uplo,
        const char* trans,
        const int* n,
        const int* k,
        const double* alpha,
        const double* a,
        const int* lda,
        const double* beta,
        double* c,
        const int* ldc
    );

	// level 2: matrix-std::vector operations, O(n^2) data and O(n^2) work.
	void sgemv_(const char*const transa, const int*const m, const int*const n,
		const float*const alpha, const float*const a, const int*const lda, const float*const x, const int*const incx,
		const float*const beta, float*const y, const int*const incy);
	void dgemv_(const char*const transa, const int*const m, const int*const n,
		const double*const alpha, const double*const a, const int*const lda, const double*const x, const int*const incx,
		const double*const beta, double*const y, const int*const incy);

	void cgemv_(const char *trans, const int *m, const int *n, const std::complex<float> *alpha,
			const std::complex<float> *a, const int *lda, const std::complex<float> *x, const int *incx,
			const std::complex<float> *beta, std::complex<float> *y, const int *incy);
		
	void zgemv_(const char *trans, const int *m, const int *n, const std::complex<double> *alpha,
			const std::complex<double> *a, const int *lda, const std::complex<double> *x, const int *incx,
			const std::complex<double> *beta, std::complex<double> *y, const int *incy);

	void dsymv_(const char *uplo, const int *n,
		const double *alpha, const double *a, const int *lda,
		const double *x, const int *incx,
		const double *beta, double *y, const int *incy);

    // A := alpha x * y.T + A
    void dger_(const int* m,
               const int* n,
               const double* alpha,
               const double* x,
               const int* incx,
               const double* y,
               const int* incy,
               double* a,
               const int* lda);
    void zgerc_(const int* m,
                const int* n,
                const std::complex<double>* alpha,
                const std::complex<double>* x,
                const int* incx,
                const std::complex<double>* y,
                const int* incy,
                std::complex<double>* a,
                const int* lda);

    // level 3: matrix-matrix operations, O(n^2) data and O(n^3) work.

	// Peize Lin add ?gemm 2017-10-27, to compute C = a * A.? * B.? + b * C
	// A is general
	void sgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
		const float *alpha, const float *a, const int *lda, const float *b, const int *ldb,
		const float *beta, float *c, const int *ldc);
	void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
		const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
		const double *beta, double *c, const int *ldc);
	void cgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
		const std::complex<float> *alpha, const std::complex<float> *a, const int *lda, const std::complex<float> *b, const int *ldb,
		const std::complex<float> *beta, std::complex<float> *c, const int *ldc);
	void zgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
		const std::complex<double> *alpha, const std::complex<double> *a, const int *lda, const std::complex<double> *b, const int *ldb,
		const std::complex<double> *beta, std::complex<double> *c, const int *ldc);

	//a is symmetric
	void dsymm_(const char *side, const char *uplo, const int *m, const int *n,
		const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
		const double *beta, double *c, const int *ldc);
	//a is hermitian
	void zhemm_(char *side, char *uplo, int *m, int *n,std::complex<double> *alpha,
		std::complex<double> *a,  int *lda,  std::complex<double> *b, int *ldb, std::complex<double> *beta, std::complex<double> *c, int *ldc);

	//solving triangular matrix with multiple right hand sides
	void dtrsm_(char *side, char* uplo, char *transa, char *diag, int *m, int *n,
		double* alpha, double* a, int *lda, double*b, int *ldb);
	void ztrsm_(char *side, char* uplo, char *transa, char *diag, int *m, int *n,
	std::complex<double>* alpha, std::complex<double>* a, int *lda, std::complex<double>*b, int *ldb);

}

// Class BlasConnector provide the connector to fortran lapack routine.
// The entire function in this class are static and inline function.
// Usage example:	BlasConnector::functionname(parameter list).
class BlasConnector
{
public:

	// Peize Lin add 2016-08-04
	// y=a*x+y
	static
	void axpy( const int n, const float alpha, const float *X, const int incX, float *Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void axpy( const int n, const double alpha, const double *X, const int incX, double *Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void axpy( const int n, const std::complex<float> alpha, const std::complex<float> *X, const int incX, std::complex<float> *Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void axpy( const int n, const std::complex<double> alpha, const std::complex<double> *X, const int incX, std::complex<double> *Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);


	// Peize Lin add 2016-08-04
	// x=a*x
	static
	void scal( const int n,  const float alpha, float *X, const int incX, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void scal( const int n, const double alpha, double *X, const int incX, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void scal( const int n, const std::complex<float> alpha, std::complex<float> *X, const int incX, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void scal( const int n, const std::complex<double> alpha, std::complex<double> *X, const int incX, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);


	// Peize Lin add 2017-10-27
	// d=x*y
	static
	float dot( const int n, const float *X, const int incX, const float *Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	double dot( const int n, const double *X, const int incX, const double *Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);


	// Peize Lin add 2017-10-27, fix bug trans 2019-01-17
	// C = a * A.? * B.? + b * C
	static
	void gemm(const char transa, const char transb, const int m, const int n, const int k,
		const float alpha, const float *a, const int lda, const float *b, const int ldb,
		const float beta, float *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void gemm(const char transa, const char transb, const int m, const int n, const int k,
		const double alpha, const double *a, const int lda, const double *b, const int ldb,
		const double beta, double *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

    static
    void gemm(const char transa, const char transb, const int m, const int n, const int k,
              const std::complex<float> alpha, const std::complex<float> *a, const int lda, const std::complex<float> *b, const int ldb,
              const std::complex<float> beta, std::complex<float> *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void gemm(const char transa, const char transb, const int m, const int n, const int k,
		const std::complex<double> alpha, const std::complex<double> *a, const int lda, const std::complex<double> *b, const int ldb,
		const std::complex<double> beta, std::complex<double> *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void gemv(const char trans, const int m, const int n,
        const float alpha, const float* A, const int lda, const float* X, const int incx,
        const float beta, float* Y, const int incy, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

    static
    void gemv(const char trans, const int m, const int n,
        const double alpha, const double* A, const int lda, const double* X, const int incx,
        const double beta, double* Y, const int incy, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

    static
    void gemv(const char trans, const int m, const int n,
          const std::complex<float> alpha, const std::complex<float> *A, const int lda, const std::complex<float> *X, const int incx,
          const std::complex<float> beta, std::complex<float> *Y, const int incy, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

    static
    void gemv(const char trans, const int m, const int n,
              const std::complex<double> alpha, const std::complex<double> *A, const int lda, const std::complex<double> *X, const int incx,
              const std::complex<double> beta, std::complex<double> *Y, const int incy, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);
 

	// Peize Lin add 2018-06-12
	// out = ||x||_2
	static
	float nrm2( const int n, const float *X, const int, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice );

	static
	double nrm2( const int n, const double *X, const int incX, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice );

	static
	double nrm2( const int n, const std::complex<double> *X, const int incX, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice );


	// copies a into b
	static
	void copy(const long n, const double *a, const int incx, double *b, const int incy, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void copy(const long n, const std::complex<double> *a, const int incx, std::complex<double> *b, const int incy, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);
};

// If GATHER_INFO is defined, the original function is replaced with a "i" suffix,
// preventing changes on the original code.
// The real function call is at gather_math_lib_info.cpp
#ifdef GATHER_INFO

#define zgemm_ zgemm_i
void zgemm_i(const char *transa,
             const char *transb,
             const int *m,
             const int *n,
             const int *k,
             const std::complex<double> *alpha,
             const std::complex<double> *a,
             const int *lda,
             const std::complex<double> *b,
             const int *ldb,
             const std::complex<double> *beta,
             std::complex<double> *c,
             const int *ldc);

#define zaxpy_  zaxpy_i
void zaxpy_i(const int *N,
            const std::complex<double> *alpha,
            const std::complex<double> *X,
            const int *incX,
            std::complex<double> *Y,
            const int *incY);

/*
#define zgemv_ zgemv_i

void zgemv_i(const char *trans,
             const int *m,
             const int *n,
             const std::complex<double> *alpha,
             const std::complex<double> *a,
             const int *lda,
             const std::complex<double> *x,
             const int *incx,
             const std::complex<double> *beta,
             std::complex<double> *y,
             const int *incy);
*/

#endif // GATHER_INFO
#endif // BLAS_CONNECTOR_H

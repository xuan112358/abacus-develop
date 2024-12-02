#ifndef OPT_TN_H
#define OPT_TN_H

#include "opt_CG.h"

#include <limits>

namespace ModuleBase
{
/**
 * @brief A class designed to deal with optimization problems with Truncated-Newton (TN) method.
 * At each step, the optimization direction d is determined by roughly
 * solving linear equation Hd=-g, where H is Hessian matrix of f(x), and g is the gradient df(x)/dx.
 * In this section, we get Hd by interpolation, and solve Hd=-g with CG method.
 * We set three truncated conditions:
 * 1) the residual of CG method has decreased more than 90%;
 * 2) the number of CG iteration is more than 50;
 * 3) the residual of CG method doesn't decrease and the CG iteration is larger than 9.
 * @author sunliang
 */
class Opt_TN
{
  public:
    Opt_TN()
    {
        this->mach_prec_ = std::numeric_limits<double>::epsilon(); // get machine precise
    }
    ~Opt_TN() {};

    /**
     * @brief Allocate the space for the arrays in cg_.
     *
     * @param nx length of the solution array x
     */
    void allocate(int nx)
    {
        this->nx_ = nx;
        this->cg_.allocate(this->nx_);
    }

    void set_para(double dV)
    {
        this->dV_ = dV;
        this->cg_.set_para(this->dV_);
    }

    /**
     * @brief Refresh the class.
     * If nx changes, reallocate space in cg_.
     *
     * @param nx_new length of new x, default 0 means the length doesn't change
     */
    void refresh(int nx_new = 0)
    {
        this->iter_ = 0;
        if (nx_new != 0)
        {
            this->nx_ = nx_new;
        }
        this->cg_.refresh(nx_new);
    }

    template <class T>
    void next_direct(
        double* px,        // current x
        double* pgradient, // df(x)/dx
        int& flag,       // record which truncated condition was triggered, 0 for cond.1, 1 for cond.2, and 2 for cond.3
        double* rdirect, // next optimization direction
        T* t,            // point of class T, which contains the gradient function
        void (T::*p_calGradient)(
            double* ptemp_x,
            double* rtemp_gradient) // a function point, which calculates the gradient at provided x
    );

    int get_iter()
    {
        return this->iter_;
    }

  private:
    Opt_CG cg_;
    double dV_ = 1.;

    int nx_ = 0;            // length of the solution array x
    int iter_ = 0;          // number of the iteration
    double mach_prec_ = 0.; // machine precision

    double inner_product(double* pa, double* pb, int length)
    {
        double innerproduct = BlasConnector::dot(length, pa, 1, pb, 1);
        innerproduct *= this->dV_;
        return innerproduct;
    }

    /**
     * @brief Get epsilon used in interpolation.
     * epsilon = 2*sqrt(mach_prec_) * (1+|x|) / |d|.
     * || means modulu.
     * @param px x
     * @param pcg_direction the direction of cg_
     * @return epsilon
     */
    double get_epsilon(double* px, double* pcg_direction)
    {
        double epsilon = 0.;
        double xx = this->inner_product(px, px, this->nx_);
        Parallel_Reduce::reduce_all(xx);
        double dd = this->inner_product(pcg_direction, pcg_direction, this->nx_);
        Parallel_Reduce::reduce_all(dd);
        epsilon = 2 * sqrt(this->mach_prec_) * (1 + sqrt(xx)) / sqrt(dd);
        // epsilon = 2 * sqrt(this->mach_prec_) * (1 + sqrt(this->inner_product(px, px, this->nx_)))
        //         / sqrt(this->inner_product(pcg_direction, pcg_direction, this->nx_));
        return epsilon;
    }
};

/**
 * @brief Get the next direction d with TN method.
 *
 * @tparam T for example, ESolver_OF
 * @param [in] px current x.
 * @param [in] pgradient df(x)/dx.
 * @param [out] flag record which truncated condition was triggered, 0 for cond.1, 1 for cond.2, and 2 for cond.3.
 * @param [out] rdirect next optimization direction.
 * @param [in] t pointer of class T, which contains the gradient function.
 * @param [in] p_calGradient a function pointer, which calculates gradient at provided x.
 */
template <class T>
void Opt_TN::next_direct(double* px,
                         double* pgradient,
                         int& flag,
                         double* rdirect,
                         T* t,
                         void (T::*p_calGradient)(double* px, double* rgradient))
{
    // initialize arrays and parameters
    ModuleBase::GlobalFunc::ZEROS(rdirect, this->nx_); // very important

    double* minus_gradient = new double[this->nx_]; // b=-g, which will be used in CG
    double* temp_x = new double[this->nx_];         // temp_x = x + step * cg_direct, used in interpolation
    double* temp_gradient = new double[this->nx_];  // df(temp_x)/dx
    double* cg_direct = new double[this->nx_];      // rdirect += cg_alpha * cg_direct at each step
    double* temp_Hcgd = new double[this->nx_];      // Hessian * cg_direct
    for (int i = 0; i < this->nx_; ++i)
    {
        minus_gradient[i] = -pgradient[i];
    }
    ModuleBase::GlobalFunc::ZEROS(cg_direct, this->nx_);
    ModuleBase::GlobalFunc::ZEROS(temp_x, this->nx_);
    ModuleBase::GlobalFunc::ZEROS(temp_gradient, this->nx_);
    ModuleBase::GlobalFunc::ZEROS(temp_Hcgd, this->nx_);

    cg_.refresh(0, minus_gradient);
    int cg_iter = 0;
    int cg_ifPD = 0;

    double epsilon = 0.;       // step length in interpolation
    double cg_alpha = 0.;      // step length got by CG
    double init_residual = 0.; // initial residual of CG
    double last_residual = 0.; // last residual of CG
    double curr_residual = 0.; // current residual of CG

    while (true)
    {
        cg_.next_direct(temp_Hcgd, 0, cg_direct);

        // get temp_Hcgd with interpolation
        // Hcgd = (df(temp_x)/dx - df(x)/x) / epsilon, where temp_x = x + step * cg_direct
        epsilon = this->get_epsilon(px, cg_direct);
        // epsilon = 1e-9;
        for (int i = 0; i < this->nx_; ++i)
        {
            temp_x[i] = px[i] + epsilon * cg_direct[i];
        }
        (t->*p_calGradient)(temp_x, temp_gradient);
        for (int i = 0; i < this->nx_; ++i)
        {
            temp_Hcgd[i] = (temp_gradient[i] - pgradient[i]) / epsilon;
        }

        // get CG step length and update rdirect
        cg_alpha = cg_.step_length(temp_Hcgd, cg_direct, cg_ifPD);
        if (cg_ifPD == -1) // Hessian is not positive definite, and cgiter = 1.
        {
            for (int i = 0; i < this->nx_; ++i)
            {
                rdirect[i] += cg_alpha * cg_direct[i];
            }
            flag = -1;
            break;
        }
        else if (cg_ifPD == -2) // Hessian is not positive definite, and cgiter > 1.
        {
            flag = -2;
            break;
        }

        for (int i = 0; i < this->nx_; ++i)
        {
            rdirect[i] += cg_alpha * cg_direct[i];
        }

        // store residuals used in truncated conditions
        last_residual = curr_residual;
        curr_residual = cg_.get_residual();
        cg_iter = cg_.get_iter();
        if (cg_iter == 1)
        {
            init_residual = curr_residual;
        }

        // check truncated conditions
        // if (curr_residual < 1e-12)
        if (curr_residual < 0.1 * init_residual)
        {
            flag = 0;
            // std::cout << "cg_ iter_ = " << cg_iter << "\n";
            break;
        }
        else if (cg_iter > 50)
        {
            flag = 1;
            break;
        }
        else if ((fabs(curr_residual - last_residual) / curr_residual) < 0.01 && cg_iter > 9)
        {
            flag = 2;
            break;
        }
    }
    this->iter_++;
    delete[] minus_gradient;
    delete[] temp_gradient;
    delete[] temp_x;
    delete[] temp_Hcgd;
    delete[] cg_direct;
}
} // namespace ModuleBase
#endif
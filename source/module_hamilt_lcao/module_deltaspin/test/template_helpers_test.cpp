#include <algorithm>
#include <string>

#include "../spin_constrain.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

/************************************************
 *  unit test of functions in template_helpers.cpp
 ***********************************************/

/**
 * - Tested functions:
 *  - Functions in template_helpers.cpp are defined by template specialization.
 *    but they are not used in the code.
 *    So, we just test if they can be called without error.
 */
#include "module_cell/klist.h"
K_Vectors::K_Vectors()
{
}
K_Vectors::~K_Vectors()
{
}

class SpinConstrainTest : public testing::Test
{
  protected:
    spinconstrain::SpinConstrain<double>& sc
        = spinconstrain::SpinConstrain<double>::getScInstance();
};

TEST_F(SpinConstrainTest, TemplatHelpers)
{
    // this is a trivial test as the double version is not used
    std::vector<std::complex<double>> Sloc2;
    EXPECT_NO_THROW(sc.cal_mw_from_lambda(0));
    EXPECT_NO_THROW(sc.cal_mi_lcao(0,false));
    EXPECT_NO_THROW(sc.run_lambda_loop(0));
    EXPECT_FALSE(sc.check_rms_stop(0, 0, 0.0, 0.0, 0.0));
    EXPECT_NO_THROW(sc.print_termination());
    EXPECT_NO_THROW(sc.print_header());
    std::vector<ModuleBase::Vector3<double>> new_spin, old_spin, new_delta_lambda, old_delta_lambda;
    EXPECT_FALSE(sc.check_gradient_decay(new_spin, old_spin, new_delta_lambda, old_delta_lambda, true));
    double alpha = 0.0;
    EXPECT_NO_THROW(sc.check_restriction(new_spin, alpha));
    EXPECT_EQ(sc.cal_alpha_opt(new_spin, old_spin, alpha), 0.0);
}
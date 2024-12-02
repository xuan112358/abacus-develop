#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "module_parameter/parameter.h"
#undef private
#define private public
#define protected public
#include "module_esolver/esolver_lj.h"
#include "module_md/md_func.h"
#include "setcell.h"

#define doublethreshold 1e-12
/************************************************
 *  unit test of functions in md_func.h
 ***********************************************/

/**
 * - Tested Function
 *   - MD_func::gaussrand
 *     - genarate Gaussian random number
 *
 *   - MD_func::rand_vel
 *     - initialize atomic velocity randomly
 *
 *   - MD_func::get_mass_mbl
 *     - initialize atomic mass and degree of freedom
 *
 *   - MD_func::read_vel
 *     - read atomic velocity from STRU
 *
 *   - MD_func::init_vel
 *     - initialize the atomic velocities
 *
 *   - MD_func::rescale_vel
 *     - rescale the velocity to the target temperature
 *
 *   - MD_func::compute_stress
 *     - calculate the contribution of classical kinetic energy of atoms to stress
 *
 *   - MD_func::dump_info
 *     - output MD dump information
 *
 *   - MD_func::print_stress
 *     - output stress
 *
 *   - MD_func::current_md_info
 *     - test the current_md_info function with the correct file path
 *
 *   - MD_func::current_step_warning
 *     - test the current_md_info function with an incorrect file path
 */

class MD_func_test : public testing::Test
{
  protected:
    UnitCell ucell;
    double* allmass;                    // atom mass
    ModuleBase::Vector3<double>* pos;   // atom position
    ModuleBase::Vector3<double>* vel;   // atom velocity
    ModuleBase::Vector3<int>* ionmbl;   // atom is frozen or not
    ModuleBase::Vector3<double>* force; // atom force
    ModuleBase::matrix virial;          // virial for this lattice
    ModuleBase::matrix stress;          // stress for this lattice
    double potential;                   // potential energy
    int natom;                          // atom number
    double temperature;                 // temperature
    int frozen_freedom;                 // frozen_freedom
    Parameter param_in;

    void SetUp()
    {
        Setcell::setupcell(ucell);
        Setcell::parameters(param_in.input);
        natom = ucell.nat;
        allmass = new double[natom];
        pos = new ModuleBase::Vector3<double>[natom];
        ionmbl = new ModuleBase::Vector3<int>[natom];
        vel = new ModuleBase::Vector3<double>[natom];
        force = new ModuleBase::Vector3<double>[natom];
        stress.create(3, 3);
        virial.create(3, 3);
    }

    void TearDown()
    {
        delete[] allmass;
        delete[] pos;
        delete[] vel;
        delete[] ionmbl;
        delete[] force;
    }
};

TEST_F(MD_func_test, gaussrand)
{
    EXPECT_DOUBLE_EQ(MD_func::gaussrand(), 1.1122716058967226);
    EXPECT_DOUBLE_EQ(MD_func::gaussrand(), -0.34532367182326629);
    EXPECT_DOUBLE_EQ(MD_func::gaussrand(), 0.60805637857480721);
}

TEST_F(MD_func_test, randomvel)
{
    ucell.init_vel = 0;
    temperature = 300 / ModuleBase::Hartree_to_K;
    MD_func::init_vel(ucell, GlobalV::MY_RANK, false, temperature, allmass, frozen_freedom, ionmbl, vel);

    EXPECT_NEAR(vel[0].x, 9.9105892783200826e-06, doublethreshold);
    EXPECT_NEAR(vel[0].y, -3.343699576563167e-05, doublethreshold);
    EXPECT_NEAR(vel[0].z, 9.385130426808701e-05, doublethreshold);
    EXPECT_NEAR(vel[1].x, -0.00017919300771203808, doublethreshold);
    EXPECT_NEAR(vel[1].y, 5.7074002254799079e-05, doublethreshold);
    EXPECT_NEAR(vel[1].z, -3.1088136026582953e-05, doublethreshold);
    EXPECT_NEAR(vel[2].x, 0.000141316492668737, doublethreshold);
    EXPECT_NEAR(vel[2].y, -0.00015841124290501442, doublethreshold);
    EXPECT_NEAR(vel[2].z, 1.900921882689748e-05, doublethreshold);
    EXPECT_NEAR(vel[3].x, 2.7965925764981002e-05, doublethreshold);
    EXPECT_NEAR(vel[3].y, 0.00013477423641584702, doublethreshold);
    EXPECT_NEAR(vel[3].z, -8.177238706840153e-05, doublethreshold);
}

TEST_F(MD_func_test, getmassmbl)
{
    ucell.init_vel = 0;
    temperature = 300 / ModuleBase::Hartree_to_K;
    MD_func::init_vel(ucell, GlobalV::MY_RANK, false, temperature, allmass, frozen_freedom, ionmbl, vel);

    for (int i = 0; i < natom; ++i)
    {
        EXPECT_DOUBLE_EQ(allmass[i], 39.948 / ModuleBase::AU_to_MASS);
        EXPECT_TRUE(ionmbl[i].x == 1);
        EXPECT_TRUE(ionmbl[i].y == 1);
        EXPECT_TRUE(ionmbl[i].z == 1);
    }

    EXPECT_TRUE(frozen_freedom == 3);
}

TEST_F(MD_func_test, readvel)
{
    MD_func::read_vel(ucell, vel);

    EXPECT_DOUBLE_EQ(vel[0].x, -0.0001320807363640);
    EXPECT_DOUBLE_EQ(vel[0].y, 7.13429429835e-05);
    EXPECT_DOUBLE_EQ(vel[0].z, -1.40179977966e-05);
    EXPECT_DOUBLE_EQ(vel[1].x, 0.000153039878532);
    EXPECT_DOUBLE_EQ(vel[1].y, -0.000146533266608);
    EXPECT_DOUBLE_EQ(vel[1].z, 9.64491480698e-05);
    EXPECT_DOUBLE_EQ(vel[2].x, -0.000133789480226);
    EXPECT_DOUBLE_EQ(vel[2].y, -3.0451038112e-06);
    EXPECT_DOUBLE_EQ(vel[2].z, -5.40998380137e-05);
    EXPECT_DOUBLE_EQ(vel[3].x, 0.000112830338059);
    EXPECT_DOUBLE_EQ(vel[3].y, 7.82354274358e-05);
    EXPECT_DOUBLE_EQ(vel[3].z, -2.83313122596e-05);
}

TEST_F(MD_func_test, RescaleVel)
{
    int frozen_freedom = 3;
    for (int i = 0; i < natom; ++i)
    {
        allmass[i] = 39.948 / ModuleBase::AU_to_MASS;
        vel[i].x = 0.1;
        vel[i].y = 0.2;
        vel[i].z = 0.3;
    }
    temperature = 300.0 / ModuleBase::Hartree_to_K;
    MD_func::rescale_vel(natom, temperature, allmass, frozen_freedom, vel);

    EXPECT_DOUBLE_EQ(vel[0].x, 4.579010775069125e-05);
    EXPECT_DOUBLE_EQ(vel[0].y, 9.15802155013825e-05);
    EXPECT_DOUBLE_EQ(vel[0].z, 0.00013737032325207373);
    EXPECT_DOUBLE_EQ(vel[1].x, 4.579010775069125e-05);
    EXPECT_DOUBLE_EQ(vel[1].y, 9.15802155013825e-05);
    EXPECT_DOUBLE_EQ(vel[1].z, 0.00013737032325207373);
    EXPECT_DOUBLE_EQ(vel[2].x, 4.579010775069125e-05);
    EXPECT_DOUBLE_EQ(vel[2].y, 9.15802155013825e-05);
    EXPECT_DOUBLE_EQ(vel[2].z, 0.00013737032325207373);
    EXPECT_DOUBLE_EQ(vel[3].x, 4.579010775069125e-05);
    EXPECT_DOUBLE_EQ(vel[3].y, 9.15802155013825e-05);
    EXPECT_DOUBLE_EQ(vel[3].z, 0.00013737032325207373);
}

TEST_F(MD_func_test, InitVelCase1)
{
    ucell.init_vel = 1;
    temperature = -1.0;
    MD_func::init_vel(ucell, GlobalV::MY_RANK, false, temperature, allmass, frozen_freedom, ionmbl, vel);

    EXPECT_NEAR(temperature, 300.0 / ModuleBase::Hartree_to_K, doublethreshold);
}

TEST_F(MD_func_test, InitVelCase2)
{
    ucell.init_vel = 1;
    temperature = 300.0 / ModuleBase::Hartree_to_K;
    MD_func::init_vel(ucell, GlobalV::MY_RANK, false, temperature, allmass, frozen_freedom, ionmbl, vel);

    EXPECT_DOUBLE_EQ(temperature, 300.0 / ModuleBase::Hartree_to_K);
}

TEST_F(MD_func_test, InitVelCase3)
{
    ucell.init_vel = 1;
    temperature = 310.0 / ModuleBase::Hartree_to_K;

    EXPECT_DOUBLE_EQ(temperature, 310.0 / ModuleBase::Hartree_to_K);
}

TEST_F(MD_func_test, InitVelCase4)
{
    ucell.init_vel = 0;
    temperature = 300.0 / ModuleBase::Hartree_to_K;
    MD_func::init_vel(ucell, GlobalV::MY_RANK, false, temperature, allmass, frozen_freedom, ionmbl, vel);

    EXPECT_DOUBLE_EQ(temperature, 300.0 / ModuleBase::Hartree_to_K);
}

TEST_F(MD_func_test, InitVelCase5)
{
    ucell.init_vel = 1;
    temperature = 330.0 / ModuleBase::Hartree_to_K;
    MD_func::init_vel(ucell, GlobalV::MY_RANK, true, temperature, allmass, frozen_freedom, ionmbl, vel);

    EXPECT_DOUBLE_EQ(temperature, 330.0 / ModuleBase::Hartree_to_K);
}

TEST_F(MD_func_test, compute_stress)
{
    temperature = 300.0 / ModuleBase::Hartree_to_K;
    MD_func::init_vel(ucell, GlobalV::MY_RANK, false, temperature, allmass, frozen_freedom, ionmbl, vel);
    MD_func::compute_stress(ucell, vel, allmass, true, virial, stress);
    EXPECT_DOUBLE_EQ(stress(0, 0), 5.2064533063674207e-06);
    EXPECT_DOUBLE_EQ(stress(0, 1), -1.6467487572445666e-06);
    EXPECT_DOUBLE_EQ(stress(0, 2), 1.5039983732220917e-06);
    EXPECT_DOUBLE_EQ(stress(1, 0), -1.6467487572445661e-06);
    EXPECT_DOUBLE_EQ(stress(1, 1), 2.380646437613151e-06);
    EXPECT_DOUBLE_EQ(stress(1, 2), -1.2514149065904968e-06);
    EXPECT_DOUBLE_EQ(stress(2, 0), 1.5039983732220917e-06);
    EXPECT_DOUBLE_EQ(stress(2, 1), -1.2514149065904968e-06);
    EXPECT_DOUBLE_EQ(stress(2, 2), 9.6330189688583664e-07);
}

TEST_F(MD_func_test, dump_info)
{
    MD_func::dump_info(0, PARAM.sys.global_out_dir, ucell, param_in, virial, force, vel);
    std::ifstream ifs("MD_dump");
    std::string output_str;
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("MDSTEP:  0"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("LATTICE_CONSTANT: 0.529177000000 Angstrom"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("LATTICE_VECTORS"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  10.000000000000  0.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  10.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  10.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("VIRIAL (kbar)"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(
        output_str,
        testing::HasSubstr("INDEX    LABEL    POSITION (Angstrom)    FORCE (eV/Angstrom)    VELOCITY (Angstrom/fs)"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  0  Ar  0.000000000000  0.000000000000  0.000000000000  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  1  Ar  2.751720222021  2.751720222021  0.000000000000  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  2  Ar  2.698802525444  0.000000000000  2.645884828867  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  3  Ar  0.000000000000  2.804637918599  2.645884828867  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    ifs.close();

    // append
    MD_func::dump_info(1, PARAM.sys.global_out_dir, ucell, param_in, virial, force, vel);
    std::ifstream ifs2("MD_dump");
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("MDSTEP:  0"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("LATTICE_CONSTANT: 0.529177000000 Angstrom"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("LATTICE_VECTORS"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  10.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  10.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  10.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("VIRIAL (kbar)"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(
        output_str,
        testing::HasSubstr("INDEX    LABEL    POSITION (Angstrom)    FORCE (eV/Angstrom)    VELOCITY (Angstrom/fs)"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  0  Ar  0.000000000000  0.000000000000  0.000000000000  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  1  Ar  2.751720222021  2.751720222021  0.000000000000  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  2  Ar  2.698802525444  0.000000000000  2.645884828867  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  3  Ar  0.000000000000  2.804637918599  2.645884828867  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    getline(ifs2, output_str);
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("MDSTEP:  1"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("LATTICE_CONSTANT: 0.529177000000 Angstrom"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("LATTICE_VECTORS"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  10.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  10.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  10.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("VIRIAL (kbar)"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(
        output_str,
        testing::HasSubstr("INDEX    LABEL    POSITION (Angstrom)    FORCE (eV/Angstrom)    VELOCITY (Angstrom/fs)"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  0  Ar  0.000000000000  0.000000000000  0.000000000000  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  1  Ar  2.751720222021  2.751720222021  0.000000000000  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  2  Ar  2.698802525444  0.000000000000  2.645884828867  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  3  Ar  0.000000000000  2.804637918599  2.645884828867  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    ifs2.close();

    remove("MD_dump");
}

TEST_F(MD_func_test, print_stress)
{
    GlobalV::ofs_running.open("running.log");
    MD_func::print_stress(GlobalV::ofs_running, virial, stress);

    std::ifstream ifs("running.log");
    std::string output_str;
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("Virtual Pressure is 0 kbar "));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("Virial Term is 0 kbar "));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("Kinetic Term is 0 kbar "));
    getline(ifs, output_str);
    getline(ifs, output_str);
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr(" ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><"));
    getline(ifs, output_str);
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr(" MD STRESS (kbar)"));
    getline(ifs, output_str);
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr(" ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><"));
    getline(ifs, output_str);
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("              0              0              0"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("              0              0              0"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("              0              0              0"));

    ifs.close();
    remove("running.log");
}

TEST_F(MD_func_test, current_md_info)
{
    // Set up the file directory and create the Restart_md.dat file
    std::string file_dir = "./";
    std::ofstream file(file_dir + "Restart_md.dat");
    file << 123;
    file.close();
    int istep = -1;
    double temperature = 0.0;
    MD_func::current_md_info(0, file_dir, istep, temperature);

    // Call the function with the correct file path and check the result
    EXPECT_EQ(istep, 123);
    EXPECT_DOUBLE_EQ(temperature, 0.0);
    remove("Restart_md.dat");
}

TEST_F(MD_func_test, current_step_warning)
{
    // Call the function and check that it outputs a warning and quits
    std::string file_dir = "./";
    int istep = 0;
    double temperature = 0.0;
    EXPECT_EXIT(MD_func::current_md_info(0, file_dir, istep, temperature), ::testing::ExitedWithCode(1), "");
}
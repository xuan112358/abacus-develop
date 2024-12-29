#include "gtest/gtest.h"
#include "gmock/gmock.h"
#define private public
#include "module_parameter/parameter.h"
#undef private
#include "module_io/write_orb_info.h"
#include "module_cell/unitcell.h"
#include "prepare_unitcell.h"
#include "module_elecstate/read_pseudo.h"

#ifdef __LCAO
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
LCAO_Orbitals::LCAO_Orbitals(){}
LCAO_Orbitals::~LCAO_Orbitals(){}
#endif
Magnetism::Magnetism()
{
	this->tot_magnetization = 0.0;
	this->abs_magnetization = 0.0;
	this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
	delete[] this->start_magnetization;
}

/************************************************
 *  unit test of write_orb_info
 ***********************************************/

/**
 * - Tested Functions:
 *   - write_orb_info()
 */


TEST(OrbInfo,WriteOrbInfo)
{
    UnitCell* ucell = new UnitCell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();
    ucell->set_iat2itia();
    std::string pp_dir = "./support/";
    std::ofstream ofs;
    ofs.open("running.log");
    PARAM.sys.global_out_dir = "./";
	PARAM.input.pseudo_rcut = 15.0;
    PARAM.input.lspinorb = false;
	PARAM.input.nspin = 1;
    PARAM.input.basis_type = "pw";
    PARAM.input.dft_functional = "default";
    PARAM.sys.nlocal = 18;
    elecstate::read_cell_pseudopots(pp_dir,ofs,*ucell);
    ucell->cal_nwfc(ofs);
    ModuleIO::write_orb_info(ucell);
    ofs.close();
    std::ifstream ifs("Orbital");
    std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("#io    spec    l    m    z  sym"));
    EXPECT_THAT(str, testing::HasSubstr("0      Si    0    0    1              s"));
    EXPECT_THAT(str, testing::HasSubstr("0      Si    2    4    1            dxy"));
    EXPECT_THAT(str, testing::HasSubstr("#sym =Symmetry name of real orbital"));
    ifs.close();
    delete ucell;
    std::remove("Orbital");
    std::remove("running.log");
}

#ifdef __MPI
#include "mpi.h"
int main(int argc, char **argv)
{


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
#endif



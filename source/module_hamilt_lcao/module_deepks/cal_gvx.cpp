/// 1. cal_gvx : gvx is used for training with force label, which is gradient of descriptors,
///      calculated by d(des)/dX = d(pdm)/dX * d(des)/d(pdm) = gdmx * gvdm
///      using einsum
/// 2. check_gvx : prints gvx into gvx.dat for checking

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "module_base/blas_connector.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_parameter/parameter.h"

// calculates gradient of descriptors from gradient of projected density
// matrices
void LCAO_Deepks::cal_gvx(const int nat, const std::vector<torch::Tensor>& gevdm)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_gvx");

    if (!gdmr_vector.empty())
    {
        gdmr_vector.erase(gdmr_vector.begin(), gdmr_vector.end());
    }

    // gdmr_vector : nat(derivative) * 3 * inl(projector) * nm * nm
    if (GlobalV::MY_RANK == 0)
    {
        // make gdmx as tensor
        int nlmax = this->inlmax / nat;
        for (int nl = 0; nl < nlmax; ++nl)
        {
            std::vector<torch::Tensor> bmmv;
            for (int ibt = 0; ibt < nat; ++ibt)
            {
                std::vector<torch::Tensor> xmmv;
                for (int i = 0; i < 3; ++i)
                {
                    std::vector<torch::Tensor> ammv;
                    for (int iat = 0; iat < nat; ++iat)
                    {
                        int inl = iat * nlmax + nl;
                        int nm = 2 * this->inl_l[inl] + 1;
                        std::vector<double> mmv;
                        for (int m1 = 0; m1 < nm; ++m1)
                        {
                            for (int m2 = 0; m2 < nm; ++m2)
                            {
                                if (i == 0)
                                {
                                    mmv.push_back(this->gdmx[ibt][inl][m1 * nm + m2]);
                                }
                                if (i == 1)
                                {
                                    mmv.push_back(this->gdmy[ibt][inl][m1 * nm + m2]);
                                }
                                if (i == 2)
                                {
                                    mmv.push_back(this->gdmz[ibt][inl][m1 * nm + m2]);
                                }
                            }
                        } // nm^2
                        torch::Tensor mm = torch::tensor(mmv, torch::TensorOptions().dtype(torch::kFloat64))
                                               .reshape({nm, nm}); // nm*nm
                        ammv.push_back(mm);
                    }
                    torch::Tensor amm = torch::stack(ammv, 0); // nat*nm*nm
                    xmmv.push_back(amm);
                }
                torch::Tensor bmm = torch::stack(xmmv, 0); // 3*nat*nm*nm
                bmmv.push_back(bmm);
            }
            this->gdmr_vector.push_back(torch::stack(bmmv, 0)); // nbt*3*nat*nm*nm
        }

        assert(this->gdmr_vector.size() == nlmax);

        // einsum for each inl:
        // gdmr_vector : b:nat(derivative) * x:3 * a:inl(projector) * m:nm *
        // n:nm gevdm : a:inl * v:nm (descriptor) * m:nm (pdm, dim1) *
        // n:nm (pdm, dim2) gvx_vector : b:nat(derivative) * x:3 *
        // a:inl(projector) * m:nm(descriptor)
        std::vector<torch::Tensor> gvx_vector;
        for (int nl = 0; nl < nlmax; ++nl)
        {
            gvx_vector.push_back(at::einsum("bxamn, avmn->bxav", {this->gdmr_vector[nl], gevdm[nl]}));
        }

        // cat nv-> \sum_nl(nv) = \sum_nl(nm_nl)=des_per_atom
        // concatenate index a(inl) and m(nm)
        this->gvx_tensor = torch::cat(gvx_vector, -1);

        assert(this->gvx_tensor.size(0) == nat);
        assert(this->gvx_tensor.size(1) == 3);
        assert(this->gvx_tensor.size(2) == nat);
        assert(this->gvx_tensor.size(3) == this->des_per_atom);
    }

    return;
}

void LCAO_Deepks::check_gvx(const int nat)
{
    std::stringstream ss;
    std::ofstream ofs_x;
    std::ofstream ofs_y;
    std::ofstream ofs_z;

    ofs_x << std::setprecision(12);
    ofs_y << std::setprecision(12);
    ofs_z << std::setprecision(12);

    for (int ia = 0; ia < nat; ia++)
    {
        ss.str("");
        ss << "gvx_" << ia << ".dat";
        ofs_x.open(ss.str().c_str());
        ss.str("");
        ss << "gvy_" << ia << ".dat";
        ofs_y.open(ss.str().c_str());
        ss.str("");
        ss << "gvz_" << ia << ".dat";
        ofs_z.open(ss.str().c_str());

        ofs_x << std::setprecision(10);
        ofs_y << std::setprecision(10);
        ofs_z << std::setprecision(10);

        for (int ib = 0; ib < nat; ib++)
        {
            for (int inl = 0; inl < inlmax / nat; inl++)
            {
                int nm = 2 * inl_l[inl] + 1;
                {
                    const int ind = ib * inlmax / nat + inl;
                    ofs_x << gvx_tensor.index({ia, 0, ib, inl}).item().toDouble() << " ";
                    ofs_y << gvx_tensor.index({ia, 1, ib, inl}).item().toDouble() << " ";
                    ofs_z << gvx_tensor.index({ia, 2, ib, inl}).item().toDouble() << " ";
                }
            }
            ofs_x << std::endl;
            ofs_y << std::endl;
            ofs_z << std::endl;
        }
        ofs_x.close();
        ofs_y.close();
        ofs_z.close();
    }
}

#endif

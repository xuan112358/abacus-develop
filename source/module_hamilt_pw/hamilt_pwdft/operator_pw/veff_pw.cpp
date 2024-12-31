#include "veff_pw.h"

#include "module_base/timer.h"
#include "module_base/tool_quit.h"

namespace hamilt {

template<typename T, typename Device>
Veff<OperatorPW<T, Device>>::Veff(const int* isk_in,
                                       const Real* veff_in,
                                       const int veff_row,
                                       const int veff_col,
                                       const ModulePW::PW_Basis_K* wfcpw_in)
{
    if (isk_in == nullptr || wfcpw_in == nullptr) {
        ModuleBase::WARNING_QUIT("VeffPW", "Constuctor of Operator::VeffPW is failed, please check your code!");
    }
    this->classname = "Veff";
    this->cal_type = calculation_type::pw_veff;
    this->isk = isk_in;
    this->veff = veff_in;
    //note: "veff = nullptr" means that this core does not treat potential but still treats wf. 
    this->veff_row = veff_row;
    this->veff_col = veff_col;
    this->wfcpw = wfcpw_in;
    resmem_complex_op()(this->ctx, this->porter, this->wfcpw->nmaxgr, "Veff<PW>::porter");
    resmem_complex_op()(this->ctx, this->porter1, this->wfcpw->nmaxgr, "Veff<PW>::porter1");

}

template<typename T, typename Device>
Veff<OperatorPW<T, Device>>::~Veff()
{
    delmem_complex_op()(this->ctx, this->porter);
    delmem_complex_op()(this->ctx, this->porter1);
}

template<typename T, typename Device>
void Veff<OperatorPW<T, Device>>::act(
    const int nbands,
    const int nbasis,
    const int npol,
    const T* tmpsi_in,
    T* tmhpsi,
    const int ngk_ik,
    const bool is_first_node)const
{
    ModuleBase::timer::tick("Operator", "VeffPW");
    if(is_first_node)
    {
        setmem_complex_op()(this->ctx, tmhpsi, 0, nbasis*nbands/npol);
    }

    int max_npw = nbasis / npol;
    const int current_spin = this->isk[this->ik];
    
    // T *porter = new T[wfcpw->nmaxgr];
    for (int ib = 0; ib < nbands; ib += npol)
    {
        if (npol == 1)
        {
            // wfcpw->recip2real(tmpsi_in, porter, this->ik);
            wfcpw->recip_to_real(this->ctx, tmpsi_in, this->porter, this->ik);
            // NOTICE: when MPI threads are larger than number of Z grids
            // veff would contain nothing, and nothing should be done in real space
            // but the 3DFFT can not be skipped, it will cause hanging
            if(this->veff_col != 0)
            {
                veff_op()(this->ctx, this->veff_col, this->porter, this->veff + current_spin * this->veff_col);
                // const Real* current_veff = &(this->veff[0](current_spin, 0));
                // for (int ir = 0; ir < this->veff->nc; ++ir)
                // {
                //     porter[ir] *= current_veff[ir];
                // }
            }
            // wfcpw->real2recip(porter, tmhpsi, this->ik, true);
            wfcpw->real_to_recip(this->ctx, this->porter, tmhpsi, this->ik, true);
        }
        else
        {
            // T *porter1 = new T[wfcpw->nmaxgr];
            // fft to real space and doing things.
            wfcpw->recip_to_real(this->ctx, tmpsi_in, this->porter, this->ik);
            wfcpw->recip_to_real(this->ctx, tmpsi_in + max_npw, this->porter1, this->ik);
            if(this->veff_col != 0)
            {
                /// denghui added at 20221109
                const Real* current_veff[4];
                for(int is = 0; is < 4; is++) {
                    current_veff[is] = this->veff + is * this->veff_col ; // for CPU device
                }
                veff_op()(this->ctx, this->veff_col, this->porter, this->porter1, current_veff);
                // T sup, sdown;
                // for (int ir = 0; ir < this->veff_col; ir++) {
                //     sup = this->porter[ir] * (current_veff[0][ir] + current_veff[3][ir])
                //         + this->porter1[ir]
                //                 * (current_veff[1][ir]
                //                 - T(0.0, 1.0) * current_veff[2][ir]);
                //     sdown = this->porter1[ir] * (current_veff[0][ir] - current_veff[3][ir])
                //             + this->porter[ir]
                //                 * (current_veff[1][ir]
                //                     + T(0.0, 1.0) * current_veff[2][ir]);
                //     this->porter[ir] = sup;
                //     this->porter1[ir] = sdown;
                // }
            }
            // (3) fft back to G space.
            wfcpw->real_to_recip(this->ctx, this->porter, tmhpsi, this->ik, true);
            wfcpw->real_to_recip(this->ctx, this->porter1, tmhpsi + max_npw, this->ik, true);
        }
        tmhpsi += max_npw * npol;
        tmpsi_in += max_npw * npol;
    }
    ModuleBase::timer::tick("Operator", "VeffPW");
}

template<typename T, typename Device>
template<typename T_in, typename Device_in>
hamilt::Veff<OperatorPW<T, Device>>::Veff(const Veff<OperatorPW<T_in, Device_in>> *veff) {
    this->classname = "Veff";
    this->cal_type = calculation_type::pw_veff;
    this->ik = veff->get_ik();
    this->isk = veff->get_isk();
    this->veff_col = veff->get_veff_col();
    this->veff_row = veff->get_veff_row();
    this->wfcpw = veff->get_wfcpw();
    resmem_complex_op()(this->ctx, this->porter, this->wfcpw->nmaxgr);
    resmem_complex_op()(this->ctx, this->porter1, this->wfcpw->nmaxgr);
    this->veff = veff->get_veff();
    if (this->isk == nullptr || this->veff == nullptr || this->wfcpw == nullptr) {
        ModuleBase::WARNING_QUIT("VeffPW", "Constuctor of Operator::VeffPW is failed, please check your code!");
    }
}

template class Veff<OperatorPW<std::complex<float>, base_device::DEVICE_CPU>>;
template class Veff<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>>;
// template Veff<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>>::Veff(const
// Veff<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>> *veff);
#if ((defined __CUDA) || (defined __ROCM))
template class Veff<OperatorPW<std::complex<float>, base_device::DEVICE_GPU>>;
template class Veff<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>>;
// template Veff<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>>::Veff(const
// Veff<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>> *veff);
#endif
} // namespace hamilt
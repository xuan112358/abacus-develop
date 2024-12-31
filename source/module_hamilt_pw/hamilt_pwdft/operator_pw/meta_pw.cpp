#include "meta_pw.h"

#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_base/tool_quit.h"

namespace hamilt {

template<typename T, typename Device>
Meta<OperatorPW<T, Device>>::Meta(Real tpiba_in,
                                       const int* isk_in,
                                       const Real* vk_in,
                                       const int vk_row,
                                       const int vk_col,
                                       const ModulePW::PW_Basis_K* wfcpw_in)
{
    if(isk_in == nullptr || tpiba_in < 1e-10 || wfcpw_in == nullptr)
    {
        ModuleBase::WARNING_QUIT("MetaPW", "Constuctor of Operator::MetaPW is failed, please check your code!");
    }
    this->classname = "Meta";
    this->cal_type = calculation_type::pw_meta;
    this->isk = isk_in;
    this->tpiba = tpiba_in;
    this->vk = vk_in;
    this->vk_row = vk_row;
    this->vk_col = vk_col;
    this->wfcpw = wfcpw_in;
    resmem_complex_op()(this->ctx, this->porter, this->wfcpw->nmaxgr, "Meta<PW>::porter");

}

template<typename T, typename Device>
Meta<OperatorPW<T, Device>>::~Meta()
{
    delmem_complex_op()(this->ctx, this->porter);
}

template<typename T, typename Device>
void Meta<OperatorPW<T, Device>>::act(
    const int nbands,
    const int nbasis,
    const int npol,
    const T* tmpsi_in,
    T* tmhpsi,
    const int ngk_ik,
    const bool is_first_node)const
{
    if (XC_Functional::get_func_type() != 3)
    {
        return;
    }

    ModuleBase::timer::tick("Operator", "MetaPW");
    if(is_first_node)
    {
        setmem_complex_op()(this->ctx, tmhpsi, 0, nbasis*nbands/npol);
    }

    const int current_spin = this->isk[this->ik];
    int max_npw = nbasis / npol;
    //npol == 2 case has not been considered

    for (int ib = 0; ib < nbands; ++ib)
    {
        for (int j = 0; j < 3; j++)
        {
            meta_op()(this->ctx, this->ik, j, ngk_ik, this->wfcpw->npwk_max, this->tpiba, wfcpw->get_gcar_data<Real>(), wfcpw->get_kvec_c_data<Real>(), tmpsi_in, this->porter);
            wfcpw->recip_to_real(this->ctx, this->porter, this->porter, this->ik);

            if(this->vk_col != 0) {
                vector_mul_vector_op()(this->ctx, this->vk_col, this->porter, this->porter, this->vk + current_spin * this->vk_col);
            }

            wfcpw->real_to_recip(this->ctx, this->porter, this->porter, this->ik);
            meta_op()(this->ctx, this->ik, j, ngk_ik, this->wfcpw->npwk_max, this->tpiba, wfcpw->get_gcar_data<Real>(), wfcpw->get_kvec_c_data<Real>(), this->porter, tmhpsi, true);

        } // x,y,z directions
        tmhpsi += max_npw;
        tmpsi_in += max_npw;
    }
    ModuleBase::timer::tick("Operator", "MetaPW");
}

template<typename T, typename Device>
template<typename T_in, typename Device_in>
Meta<OperatorPW<T, Device>>::Meta(const Meta<OperatorPW<T_in, Device_in>> *meta) {
    this->classname = "Meta";
    this->cal_type = calculation_type::pw_meta;
    this->ik = meta->get_ik();
    this->isk = meta->get_isk();
    this->tpiba = meta->get_tpiba();
    this->vk = meta->get_vk();
    this->vk_row = meta->get_vk_row();
    this->vk_col = meta->get_vk_col();
    this->wfcpw = meta->get_wfcpw();
    if(this->isk == nullptr || this->tpiba < 1e-10 || this->vk == nullptr || this->wfcpw == nullptr)
    {
        ModuleBase::WARNING_QUIT("MetaPW", "Constuctor of Operator::MetaPW is failed, please check your code!");
    }
}

template class Meta<OperatorPW<std::complex<float>, base_device::DEVICE_CPU>>;
template class Meta<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>>;
// template Meta<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>>::Meta(const
// Meta<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>> *meta);
#if ((defined __CUDA) || (defined __ROCM))
template class Meta<OperatorPW<std::complex<float>, base_device::DEVICE_GPU>>;
template class Meta<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>>;
// template Meta<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>>::Meta(const
// Meta<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>> *meta); template
// Meta<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>>::Meta(const Meta<OperatorPW<std::complex<double>,
// base_device::DEVICE_CPU>> *meta); template Meta<OperatorPW<std::complex<double>,
// base_device::DEVICE_GPU>>::Meta(const Meta<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>> *meta);
#endif
} // namespace hamilt
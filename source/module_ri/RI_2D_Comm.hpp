//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef RI_2D_COMM_HPP
#define RI_2D_COMM_HPP

#include "RI_2D_Comm.h"
#include "RI_Util.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/tool_title.h"
#include "module_base/timer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_parameter/parameter.h"
#include <RI/global/Global_Func-2.h>

#include <cmath>
#include <string>
#include <stdexcept>

inline RI::Tensor<double> tensor_conj(const RI::Tensor<double>& t) { return t; }
inline RI::Tensor<std::complex<double>> tensor_conj(const RI::Tensor<std::complex<double>>& t)
{
    RI::Tensor<std::complex<double>> r(t.shape);
    for (int i = 0; i < t.data->size(); ++i) {
        (*r.data)[i] = std::conj((*t.data)[i]);
    }
    return r;
}
template<typename Tdata, typename Tmatrix>
auto RI_2D_Comm::split_m2D_ktoR(const UnitCell& ucell,
                                const K_Vectors & kv, 
                                const std::vector<const Tmatrix*>&mks_2D, 
                                const Parallel_2D & pv, 
                                const int nspin, 
                                const bool spgsym)
-> std::vector<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>>
{
	ModuleBase::TITLE("RI_2D_Comm","split_m2D_ktoR");
	ModuleBase::timer::tick("RI_2D_Comm", "split_m2D_ktoR");

	const TC period = RI_Util::get_Born_vonKarmen_period(kv);
	const std::map<int,int> nspin_k = {{1,1}, {2,2}, {4,1}};
    const double SPIN_multiple = std::map<int, double>{ {1,0.5}, {2,1}, {4,1} }.at(nspin);							// why?

    std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>> mRs_a2D(nspin);
    for (int is_k = 0; is_k < nspin_k.at(nspin); ++is_k)
	{
		const std::vector<int> ik_list = RI_2D_Comm::get_ik_list(kv, is_k);
		for(const TC &cell : RI_Util::get_Born_von_Karmen_cells(period))
		{
            RI::Tensor<Tdata> mR_2D;
            int ik_full = 0;
            for (const int ik : ik_list)
            {
                auto set_mR_2D = [&mR_2D](auto&& mk_frac) {
                    if (mR_2D.empty()) {
                        mR_2D = RI::Global_Func::convert<Tdata>(mk_frac);
                    } else {
                        mR_2D
                            = mR_2D + RI::Global_Func::convert<Tdata>(mk_frac);
                    }
                };
                using Tdata_m = typename Tmatrix::value_type;
                if (!spgsym)
                {
                    RI::Tensor<Tdata_m> mk_2D = RI_Util::Vector_to_Tensor<Tdata_m>(*mks_2D[ik], pv.get_col_size(), pv.get_row_size());
                    const Tdata_m frac = SPIN_multiple
                        * RI::Global_Func::convert<Tdata_m>(std::exp(
                            -ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT * (kv.kvec_c[ik] * (RI_Util::array3_to_Vector3(cell) * ucell.latvec))));
                    if (static_cast<int>(std::round(SPIN_multiple * kv.wk[ik] * kv.get_nkstot_full())) == 2)
                        { set_mR_2D(mk_2D * (frac * 0.5) + tensor_conj(mk_2D * (frac * 0.5))); }
                    else { set_mR_2D(mk_2D * frac); }
                }
                else
                { // traverse kstar, ik means ik_ibz
                    for (auto& isym_kvd : kv.kstars[ik % ik_list.size()])
                    {
                        RI::Tensor<Tdata_m> mk_2D = RI_Util::Vector_to_Tensor<Tdata_m>(*mks_2D[ik_full + is_k * kv.get_nkstot_full()], pv.get_col_size(), pv.get_row_size());
                        const Tdata_m frac = SPIN_multiple
                            * RI::Global_Func::convert<Tdata_m>(std::exp(
                                -ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT * ((isym_kvd.second * ucell.G) * (RI_Util::array3_to_Vector3(cell) * ucell.latvec))));
                        set_mR_2D(mk_2D * frac);
                        ++ik_full;
                    }
                }
            }
			for(int iwt0_2D=0; iwt0_2D!=mR_2D.shape[0]; ++iwt0_2D)
			{
				const int iwt0 =ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver)
                    ? pv.local2global_col(iwt0_2D)
                    : pv.local2global_row(iwt0_2D);
				int iat0, iw0_b, is0_b;
				std::tie(iat0,iw0_b,is0_b) = RI_2D_Comm::get_iat_iw_is_block(ucell,iwt0);
				const int it0 = ucell.iat2it[iat0];
				for(int iwt1_2D=0; iwt1_2D!=mR_2D.shape[1]; ++iwt1_2D)
				{
					const int iwt1 =ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver)
                        ? pv.local2global_row(iwt1_2D)
                        : pv.local2global_col(iwt1_2D);
					int iat1, iw1_b, is1_b;
					std::tie(iat1,iw1_b,is1_b) = RI_2D_Comm::get_iat_iw_is_block(ucell,iwt1);
					const int it1 = ucell.iat2it[iat1];

					const int is_b = RI_2D_Comm::get_is_block(is_k, is0_b, is1_b);
					RI::Tensor<Tdata> &mR_a2D = mRs_a2D[is_b][iat0][{iat1,cell}];
                    if (mR_a2D.empty()) {
                        mR_a2D = RI::Tensor<Tdata>(
                            {static_cast<size_t>(ucell.atoms[it0].nw),
                             static_cast<size_t>(
                                 ucell.atoms[it1].nw)});
                    }
                    mR_a2D(iw0_b,iw1_b) = mR_2D(iwt0_2D, iwt1_2D);
				}
			}
        }
    }
	ModuleBase::timer::tick("RI_2D_Comm", "split_m2D_ktoR");
	return mRs_a2D;
}


template<typename Tdata, typename TK>
void RI_2D_Comm::add_Hexx(
    const UnitCell &ucell,
	const K_Vectors &kv,
	const int ik,
	const double alpha,
	const std::vector<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>> &Hs,
    const Parallel_Orbitals& pv,
    TK* hk)
{
	ModuleBase::TITLE("RI_2D_Comm","add_Hexx");
	ModuleBase::timer::tick("RI_2D_Comm", "add_Hexx");

	const std::map<int, std::vector<int>> is_list = {{1,{0}}, {2,{kv.isk[ik]}}, {4,{0,1,2,3}}};
	for(const int is_b : is_list.at(PARAM.inp.nspin))
	{
		int is0_b, is1_b;
		std::tie(is0_b,is1_b) = RI_2D_Comm::split_is_block(is_b);
		for(const auto &Hs_tmpA : Hs[is_b])
		{
			const TA &iat0 = Hs_tmpA.first;
			for(const auto &Hs_tmpB : Hs_tmpA.second)
			{
				const TA &iat1 = Hs_tmpB.first.first;
				const TC &cell1 = Hs_tmpB.first.second;
				const std::complex<double> frac = alpha
					* std::exp( ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT * (kv.kvec_c[ik] * (RI_Util::array3_to_Vector3(cell1)*ucell.latvec)) );
				const RI::Tensor<Tdata> &H = Hs_tmpB.second;
				for(size_t iw0_b=0; iw0_b<H.shape[0]; ++iw0_b)
				{
					const int iwt0 = RI_2D_Comm::get_iwt(ucell,iat0, iw0_b, is0_b);
                    if (pv.global2local_row(iwt0) < 0) {
                        continue;
                    }
                    for(size_t iw1_b=0; iw1_b<H.shape[1]; ++iw1_b)
					{
						const int iwt1 = RI_2D_Comm::get_iwt(ucell,iat1, iw1_b, is1_b);
                        if (pv.global2local_col(iwt1) < 0) {
                            continue;
                        }
                        LCAO_domain::set_mat2d(iwt0, iwt1, RI::Global_Func::convert<TK>(H(iw0_b, iw1_b)) * RI::Global_Func::convert<TK>(frac), pv, hk);
					}
				}
			}
		}
	}
	ModuleBase::timer::tick("RI_2D_Comm", "add_Hexx");
}

std::tuple<int,int,int>
RI_2D_Comm::get_iat_iw_is_block(const UnitCell& ucell,const int& iwt)
{
	const int iat = ucell.iwt2iat[iwt];
	const int iw = ucell.iwt2iw[iwt];
	switch(PARAM.inp.nspin)
	{
		case 1: case 2:
			return std::make_tuple(iat, iw, 0);
		case 4:
			return std::make_tuple(iat, iw/2, iw%2);
		default:
			throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
}

int RI_2D_Comm::get_is_block(const int is_k, const int is_row_b, const int is_col_b)
{
	switch(PARAM.inp.nspin)
	{
		case 1:		return 0;
		case 2:		return is_k;
		case 4:		return is_row_b*2+is_col_b;
		default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
}

std::tuple<int,int>
RI_2D_Comm::split_is_block(const int is_b)
{
	switch(PARAM.inp.nspin)
	{
		case 1:	case 2:
			return std::make_tuple(0, 0);
		case 4:
			return std::make_tuple(is_b/2, is_b%2);
		default:
			throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
}



int RI_2D_Comm::get_iwt(const UnitCell& ucell,
                        const int iat, 
                        const int iw_b, 
                        const int is_b)
{
	const int it = ucell.iat2it[iat];
	const int ia = ucell.iat2ia[iat];
	int iw=-1;
	switch(PARAM.inp.nspin)
	{
		case 1: case 2:
			iw = iw_b;			break;
		case 4:
			iw = iw_b*2+is_b;	break;
		default:
			throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
	const int iwt = ucell.itiaiw2iwt(it,ia,iw);
	return iwt;
}

template<typename Tdata, typename TR>
void RI_2D_Comm::add_HexxR(
    const int current_spin,
    const double alpha,
    const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>& Hs,
    const Parallel_Orbitals& pv,
    const int npol,
    hamilt::HContainer<TR>& hR,
    const RI::Cell_Nearest<int, int, 3, double, 3>* const cell_nearest)
{
    ModuleBase::TITLE("RI_2D_Comm", "add_HexxR");
    ModuleBase::timer::tick("RI_2D_Comm", "add_HexxR");
    const std::map<int, std::vector<int>> is_list = { {1,{0}}, {2,{current_spin}}, {4,{0,1,2,3}} };
    for (const int is_hs : is_list.at(PARAM.inp.nspin))
    {
        int is0_b = 0, is1_b = 0;
        std::tie(is0_b, is1_b) = RI_2D_Comm::split_is_block(is_hs);
        for (const auto& Hs_tmpA : Hs[is_hs])
        {
            const TA& iat0 = Hs_tmpA.first;
            for (const auto& Hs_tmpB : Hs_tmpA.second)
            {
                const TA& iat1 = Hs_tmpB.first.first;
                const TC& cell = Hs_tmpB.first.second;
                const Abfs::Vector3_Order<int> R = RI_Util::array3_to_Vector3(
                    (cell_nearest ?
                        cell_nearest->get_cell_nearest_discrete(iat0, iat1, cell)
                        : cell));
                hamilt::BaseMatrix<TR>* HlocR = hR.find_matrix(iat0, iat1, R.x, R.y, R.z);
                if (HlocR == nullptr)
                { // add R to HContainer
                    hamilt::AtomPair<TR> tmp(iat0, iat1, R.x, R.y, R.z, &pv);
                    hR.insert_pair(tmp);
                    HlocR = hR.find_matrix(iat0, iat1, R.x, R.y, R.z);
                }
                auto row_indexes = pv.get_indexes_row(iat0);
                auto col_indexes = pv.get_indexes_col(iat1);
                const RI::Tensor<Tdata>& HexxR = (Tdata)alpha * Hs_tmpB.second;
                for (int lw0_b = 0;lw0_b < row_indexes.size();lw0_b += npol)    // block
                {
                    const int& gw0 = row_indexes[lw0_b] / npol;
                    const int& lw0 = (npol == 2) ? (lw0_b + is0_b) : lw0_b;
                    for (int lw1_b = 0;lw1_b < col_indexes.size();lw1_b += npol)
                    {
                        const int& gw1 = col_indexes[lw1_b] / npol;
                        const int& lw1 = (npol == 2) ? (lw1_b + is1_b) : lw1_b;
                        HlocR->add_element(lw0, lw1, RI::Global_Func::convert<TR>(HexxR(gw0, gw1)));
                    }
                }
            }
        }
    }

    ModuleBase::timer::tick("RI_2D_Comm", "add_HexxR");
}

#endif
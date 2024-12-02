#include "wavefunc_in_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_parameter/parameter.h"
#include <cstring>		// Peize Lin fix bug about strcmp 2016-08-02
#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"
#include "module_base/math_polyint.h"
#include "module_base/math_ylmreal.h"
#include "module_hamilt_pw/hamilt_pwdft/soc.h"

void Wavefunc_in_pw::make_table_q(
	const UnitCell &ucell,
	std::vector<std::string> &fn, 
	ModuleBase::realArray &table_local)
{
	ModuleBase::TITLE("Wavefunc_in_pw","make_table_q");

	if( fn.size() != static_cast<size_t>(ucell.ntype) )
	{
		ModuleBase::WARNING_QUIT("Wavefunc_in_pw::make_table_q","maybe NUMERICAL_ORBITAL is not read in, please check.");
	}

	for(int it=0; it<ucell.ntype; it++)
	{
		std::ifstream in(fn[it].c_str());
		if(!in)
		{
			GlobalV::ofs_warning << " File name : " << fn[it] << std::endl;
			ModuleBase::WARNING_QUIT("Wavefunc_in_pw::make_table_q","Can not find file.");
		}
		else
		{
			std::stringstream ss;
			ss << "Orbital of species " << ucell.atoms[it].label;
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,ss.str(),fn[it]);
		}
		in.close();
	}

	table_local.zero_out();
	for(int it=0; it<ucell.ntype; it++)
	{
		int ic=0;
		for(int L=0; L<ucell.atoms[it].nwl+1; L++)
		{
			for(int N=0; N<ucell.atoms[it].l_nchi[L]; N++)
			{
				GlobalV::ofs_running << " L=" << L << " N=" << N;
				std::ifstream in(fn[it].c_str());
				if (!in)
				{
					GlobalV::ofs_warning << " File name : " << fn[it] << std::endl;
					ModuleBase::WARNING_QUIT("Wavefunc_in_pw::make_table_q","Can not find file.");
				}
				int meshr=0;
				double dr=0.0; // only used in uniform grid
				char word[80];     // pengfei Li add 15-1-31
				while (in.good())
				{
					in >> word;
					if (std::strcmp(word , "END") == 0)		// Peize Lin fix bug about strcmp 2016-08-02
					{
						break;
					}
				}

				ModuleBase::CHECK_NAME(in, "Mesh");
				in >> meshr;
				int meshr_read = meshr;
				if(meshr%2==0)
				{
					++meshr;
				}
				GlobalV::ofs_running << " meshr=" << meshr;

				ModuleBase::CHECK_NAME(in, "dr");
				in >> dr;
				GlobalV::ofs_running << " dr=" << dr;

				double* radial = new double[meshr];
				double *psi = new double[meshr];
				double* psir = new double[meshr];
				double* rab = new double[meshr];

				ModuleBase::GlobalFunc::ZEROS(radial, meshr);
				ModuleBase::GlobalFunc::ZEROS(psi, meshr);
				ModuleBase::GlobalFunc::ZEROS(psir, meshr);
				ModuleBase::GlobalFunc::ZEROS(rab, meshr);
				for(int ir=0; ir<meshr; ir++)
				{
					rab[ir] = dr;
					// plus one because we can't read in r = 0 term now.
					radial[ir] = ir*dr;  //mohan modify 2010-04-19
				}
				GlobalV::ofs_running << " Rmax(Angstrom)=" << radial[meshr-1] << std::endl;

				std::string name1;
				std::string name2;
				std::string name3;
				int tmp_it=0;
				int tmp_l=0;
				int tmp_n=0;
				bool find = false;

				while( !find )
				{
					if(in.eof())
					{
						GlobalV::ofs_warning << "\n Can't find l="
						<< L << " n=" << N << " orbital." << std::endl;
						ModuleBase::WARNING_QUIT("Control_Overlap","Read_PAO");
					}
					in >> name1 >> name2 >> name3;
					assert( name1 == "Type" );
					in >> tmp_it >> tmp_l >> tmp_n;
					if( L == tmp_l && N == tmp_n )
					{
						// meshr_read is different from meshr if meshr is even number.
						for(int ir=0; ir<meshr_read; ir++)
						{
							in >> psi[ir];
							//psi[ir] = 1.0; //hahaha
							psir[ir] = psi[ir] * radial[ir];
						}
						find = true;
					}
					else
					{
						 double no_use;
						 for(int ir=0; ir<meshr_read; ir++)
						 {
							 in >> no_use;
						 }
					}
				}
				double* table = new double[PARAM.globalv.nqx];
				Wavefunc_in_pw::integral(ucell,meshr, psir, radial, rab, L, table);
				for(int iq=0; iq<PARAM.globalv.nqx; iq++)
				{
					//double energy_q = pow(iq * PARAM.globalv.dq,2);
					table_local(it,ic,iq) = table[iq];//* Wavefunc_in_pw::smearing(energy_q,150,0.666666);
				}
				delete[] table;
				delete[] radial;
				delete[] rab;
				delete[] psi;
				delete[] psir;
				++ic;
			}// N
		}// L
	}// T


	if(GlobalV::MY_RANK==0)
	{
		for(int it=0; it<ucell.ntype; it++)
		{
			std::stringstream ss;
			ss << PARAM.globalv.global_out_dir << ucell.atoms[it].label << "/LOCAL_G.dat";
			std::ofstream ofs(ss.str().c_str());
			for(int iq=0; iq<PARAM.globalv.nqx; iq++)
			{
				int ic=0;
				double energy_q = pow((double)iq*PARAM.globalv.dq,2);
				ofs << energy_q; // unit (Ry)
				for(int L=0; L<ucell.atoms[it].nwl+1; L++)
				{
					for(int N=0; N<ucell.atoms[it].l_nchi[L]; N++)
					{
						ofs << " " << table_local(it,ic,iq);
						++ic;
					}
				}
				ofs << std::endl;
			}
			ofs.close();
		}
	}

	return;
}

double Wavefunc_in_pw::smearing(const double &energy_x,
                               const double &ecut,
                               const double &beta)
{
    double w = 0.0;
    const double beta_e = beta * ecut ;

    if (beta >= 1.0 || beta<0 )
    {
        ModuleBase::WARNING_QUIT("wavefunc_in_pw::smearing", "beta must between 0 ~ 1 ");
    }

    if (energy_x < beta_e)
    {
        w = 1.0;
    }
    else if (energy_x >= beta_e && energy_x <= ecut)
    {
        const double arg = ModuleBase::PI * (ecut - energy_x) * 0.5 / (1-beta) / ecut ;
        // const double sin_arg = sin(arg);  // gong 2009. 7. 12 , correct
        // w = sin_arg*sin_argi ;
        w = 0.5 * (1 - cos(2.0 * arg));
    }
    else if (energy_x > ecut)
    {
        w = 0.0 ;
    }

    return w ;
}


void Wavefunc_in_pw::integral(const UnitCell& ucell,
							  const int meshr, 
							  const double *psir, 
							  const double *r,
							  const double *rab, 
							  const int &l, 
							  double* table)
{
	const double pref = ModuleBase::FOUR_PI / sqrt(ucell.omega);

	double *inner_part = new double[meshr];
	for(int ir=0; ir<meshr; ir++)
	{
		inner_part[ir] = psir[ir] * psir[ir];
	}

	double unit = 0.0;
	ModuleBase::Integral::Simpson_Integral(meshr, inner_part, rab, unit);
	delete[] inner_part;
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"normalize unit",unit);

	double *aux = new double[meshr];
	double *vchi = new double[meshr];
	for (int iq=0; iq<PARAM.globalv.nqx; iq++)
	{
		const double q = PARAM.globalv.dq * iq;
		ModuleBase::Sphbes::Spherical_Bessel(meshr, r, q, l, aux);
		for (int ir = 0;ir < meshr;ir++)
		{
			vchi[ir] = psir[ir] * aux[ir] * r[ir];
		}

		double vqint = 0.0;
		ModuleBase::Integral::Simpson_Integral(meshr, vchi, rab, vqint);

		table[iq] =  vqint * pref;
	}
	delete[] aux;
	delete[] vchi;
	return;
}

void Wavefunc_in_pw::produce_local_basis_in_pw(const UnitCell& ucell,
											   const int& ik,
                                               const ModulePW::PW_Basis_K* wfc_basis,
                                               const Structure_Factor& sf,
                                               ModuleBase::ComplexMatrix& psi,
                                               const ModuleBase::realArray& table_local)
{
	ModuleBase::TITLE("Wavefunc_in_pw","produce_local_basis_in_pw");
	assert(ik>=0);
	const int npw = wfc_basis->npwk[ik];
	const int total_lm = ( ucell.lmax + 1) * ( ucell.lmax + 1);
	ModuleBase::matrix ylm(total_lm, npw);
	std::complex<double> *aux = new std::complex<double>[npw];
	double *chiaux = nullptr;

	ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
	for(int ig=0;ig<npw;ig++)
	{
		gk[ig] = wfc_basis->getgpluskcar(ik, ig);
	}

	ModuleBase::YlmReal::Ylm_Real(total_lm, npw, gk, ylm);

	//int index = 0;
	double *flq = new double[npw];
	int iwall=0;
	for (int it = 0;it < ucell.ntype;it++)
	{
		for (int ia = 0;ia < ucell.atoms[it].na;ia++)
		{
            std::complex<double>* sk = sf.get_sk(ik, it, ia, wfc_basis);
            int ic = 0;
            for(int L = 0; L < ucell.atoms[it].nwl+1; L++)
			{
				std::complex<double> lphase = pow(ModuleBase::NEG_IMAG_UNIT, L); //mohan 2010-04-19
				for(int N=0; N < ucell.atoms[it].l_nchi[L]; N++)
				{
//					GlobalV::ofs_running << " it=" << it << " ia=" << ia << " L=" << L << " N=" << N << std::endl;

					for(int ig=0; ig<npw; ig++)
					{
						flq[ig] = ModuleBase::PolyInt::Polynomial_Interpolation(table_local,
						it, ic, PARAM.globalv.nqx, PARAM.globalv.dq, gk[ig].norm() * ucell.tpiba );
					}

					if(PARAM.inp.nspin==4)
					{
/*						for(int is_N = 0; is_N < 2; is_N++)*/  //for rotate base
						for(int is_N = 0; is_N < 1; is_N++)
						{
							if(L==0 && is_N==1) { continue;}
							if(ucell.atoms[it].ncpp.has_so)
							{
								const double j = std::abs(double(L+is_N) - 0.5);
								if (!(PARAM.globalv.domag||PARAM.globalv.domag_z))
								{//atomic_wfc_so
									for(int m=0; m<2*L+1; m++)
									{
										std::cout<<"iwall: "<<iwall<<std::endl;
										const int lm = L*L+m;
										for(int ig=0; ig<npw; ig++)
										{
											//if(is_N==0)
											psi(iwall, ig) =
											lphase * sk[ig] * ylm(lm, ig) * flq[ig];
											//else
                                            psi(iwall + 1, ig + wfc_basis->npwk_max)
                                                = lphase * sk[ig] * ylm(lm, ig) * flq[ig];
                                        }
                                        iwall += 2;
                                    }
								}//if
								else
								{//atomic_wfc_so_mag
									double alpha, gamma;
									std::complex<double> fup,fdown;
                              		//int nc;
                              		//This routine creates two functions only in the case j=l+1/2 or exit in the other case
									if(fabs(j-L+0.5)<1e-4) { continue;
}
									delete[] chiaux;
									chiaux = new double [npw];
                              		//Find the functions j= l- 1/2
									if(L==0) {
									for(int ig=0;ig<npw;ig++){
										chiaux[ig] = flq[ig];
									}
									} else
									{
										/*for(int ib = 0;ib < ucell.atoms[it].nchi;ib++)
										{
											if((ucell.atoms[it].lchi[ib] == L)&&(fabs(ucell.atoms[it].jjj[ib]-L+0.5)<1e-4))
											{
											nc=ib;
											break;
											}
										}*/
										for(int ig=0;ig<npw;ig++)
										{//Average the two functions
											chiaux[ig] =  L *
												ModuleBase::PolyInt::Polynomial_Interpolation(table_local,
												it, ic, PARAM.globalv.nqx, PARAM.globalv.dq, gk[ig].norm() * ucell.tpiba );

											chiaux[ig] += flq[ig] * (L+1.0) ;
											chiaux[ig] *= 1/(2.0*L+1.0);
										}
									}
									//and construct the starting wavefunctions as in the noncollinear case.
									//alpha = ucell.magnet.angle1_[it];
									//gamma = -1 * ucell.magnet.angle2_[it] + 0.5 * ModuleBase::PI;
									alpha = ucell.atoms[it].angle1[ia];
									gamma = -1 * ucell.atoms[it].angle2[ia] + 0.5 * ModuleBase::PI;
									for(int m = 0;m<2*L+1;m++)
									{
										const int lm = L*L +m;

                                        if (iwall + 2 * L + 1 > ucell.natomwfc) 
                                        {
                                            ModuleBase::WARNING_QUIT("this->wf.atomic_wfc()", "error: too many wfcs");
                                        }
                                        for (int ig = 0; ig < npw; ig++)
                                        {
                                            aux[ig] = sk[ig] * ylm(lm,ig) * chiaux[ig];
                                        }
                                        //rotate wfc as needed
										//first rotation with angle alpha around (OX)
										for(int ig = 0;ig<npw;ig++)
										{
											fup = cos(0.5 * alpha) * aux[ig];
											fdown = ModuleBase::IMAG_UNIT * sin(0.5* alpha) * aux[ig];
											//build the orthogonal wfc
											//first rotation with angle (alpha + ModuleBase::PI) around (OX)
											psi(iwall,ig) = (cos(0.5 * gamma) + ModuleBase::IMAG_UNIT * sin(0.5*gamma)) * fup;
                                            psi(iwall, ig + wfc_basis->npwk_max)
                                                = (cos(0.5 * gamma) - ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fdown;
                                            // second rotation with angle gamma around(OZ)
                                            fup = cos(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                            fdown = ModuleBase::IMAG_UNIT * sin(0.5 * (alpha + ModuleBase::PI))*aux[ig];
											psi(iwall+2*L+1,ig) = (cos(0.5*gamma) + ModuleBase::IMAG_UNIT*sin(0.5*gamma))*fup;
                                            psi(iwall + 2 * L + 1, ig + wfc_basis->npwk_max)
                                                = (cos(0.5 * gamma) - ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fdown;
                                        }
                                        iwall++;
                                    }
									iwall += 2*L +1;
								} // end else INPUT.starting_spin_angle || !PARAM.globalv.domag
							} // end if ucell.atoms[it].has_so
							else
							{//atomic_wfc_nc
								double alpha, gamman;
								std::complex<double> fup, fdown;
								//alpha = ucell.magnet.angle1_[it];
								//gamman = -ucell.magnet.angle2_[it] + 0.5*ModuleBase::PI;
								alpha = ucell.atoms[it].angle1[ia];
								gamman = -ucell.atoms[it].angle2[ia] + 0.5*ModuleBase::PI;
								for(int m = 0;m<2*L+1;m++)
								{
									const int lm = L*L +m;
									if (iwall + 2 * L + 1 > ucell.natomwfc)
									{
										ModuleBase::WARNING_QUIT("this->wf.atomic_wfc()", "error: too many wfcs");
									}
                                    for (int ig = 0; ig < npw; ig++)
                                    {
                                        aux[ig] = sk[ig] * ylm(lm,ig) * flq[ig];
                                    }
                                    //rotate function
									//first, rotation with angle alpha around(OX)
									for(int ig = 0;ig<npw;ig++)
									{
										fup = cos(0.5*alpha) * aux[ig];
										fdown = ModuleBase::IMAG_UNIT * sin(0.5* alpha) * aux[ig];
										//build the orthogonal wfc
										//first rotation with angle(alpha+ModuleBase::PI) around(OX)
										psi(iwall,ig) = (cos(0.5 * gamman) + ModuleBase::IMAG_UNIT * sin(0.5*gamman)) * fup;
                                        psi(iwall, ig + wfc_basis->npwk_max)
                                            = (cos(0.5 * gamman) - ModuleBase::IMAG_UNIT * sin(0.5 * gamman)) * fdown;
                                        // second rotation with angle gamma around(OZ)
                                        fup = cos(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                        fdown = ModuleBase::IMAG_UNIT * sin(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
										psi(iwall+2*L+1,ig) = (cos(0.5*gamman) + ModuleBase::IMAG_UNIT*sin(0.5*gamman))*fup;
                                        psi(iwall + 2 * L + 1, ig + wfc_basis->npwk_max)
                                            = (cos(0.5 * gamman) - ModuleBase::IMAG_UNIT * sin(0.5 * gamman)) * fdown;
                                    } // end ig
                                    iwall++;
                                } // end m
								iwall += 2*L+1;
							} // end else ucell.atoms[it].has_so
						} // end for is_N
                    } // end if PARAM.inp.noncolin
					else
					{//LSDA and nomagnet case
						for(int m=0; m<2*L+1; m++)
						{
							const int lm = L*L+m;
							for(int ig=0; ig<npw; ig++)
							{
								psi(iwall, ig) =
								lphase * sk[ig] * ylm(lm, ig) * flq[ig];
							}
							++iwall;
						}
					}
					++ic;
				} // end for N
			} // end for L
			delete[] sk;
		} // end for ia
	} // end for it
	assert(iwall == PARAM.globalv.nlocal);
	delete[] flq;
	delete[] aux;
	delete[] chiaux;
	delete[] gk;
}


# Full List of INPUT Keywords

- [Full List of INPUT Keywords](#full-list-of-input-keywords)
  - [System variables](#system-variables)
    - [suffix](#suffix)
    - [calculation](#calculation)
    - [esolver\_type](#esolver_type)
    - [symmetry](#symmetry)
    - [symmetry\_prec](#symmetry_prec)
    - [symmetry\_autoclose](#symmetry_autoclose)
    - [kpar](#kpar)
    - [bndpar](#bndpar)
    - [latname](#latname)
    - [psi\_initializer](#psi_initializer)
    - [init\_wfc](#init_wfc)
    - [init\_chg](#init_chg)
    - [init\_vel](#init_vel)
    - [mem\_saver](#mem_saver)
    - [diago\_proc](#diago_proc)
    - [nbspline](#nbspline)
    - [kspacing](#kspacing)
    - [min\_dist\_coef](#min_dist_coef)
    - [device](#device)
    - [precision](#precision)
  - [Variables related to input files](#variables-related-to-input-files)
    - [stru\_file](#stru_file)
    - [kpoint\_file](#kpoint_file)
    - [pseudo\_dir](#pseudo_dir)
    - [orbital\_dir](#orbital_dir)
    - [read\_file\_dir](#read_file_dir)
    - [restart\_load](#restart_load)
    - [wannier\_card](#wannier_card)
  - [Plane wave related variables](#plane-wave-related-variables)
    - [ecutwfc](#ecutwfc)
    - [ecutrho](#ecutrho)
    - [nx, ny, nz](#nx-ny-nz)
    - [ndx, ndy, ndz](#ndx-ndy-ndz)
    - [pw\_seed](#pw_seed)
    - [pw\_diag\_thr](#pw_diag_thr)
    - [diago\_smooth\_ethr](#diago_smooth_ethr)
    - [pw\_diag\_nmax](#pw_diag_nmax)
    - [pw\_diag\_ndim](#pw_diag_ndim)
    - [erf\_ecut](#erf_ecut)
    - [fft\_mode](#fft_mode)
    - [erf\_height](#erf_height)
    - [erf\_sigma](#erf_sigma)
  - [Numerical atomic orbitals related variables](#numerical-atomic-orbitals-related-variables)
    - [nb2d](#nb2d)
    - [lmaxmax](#lmaxmax)
    - [lcao\_ecut](#lcao_ecut)
    - [lcao\_dk](#lcao_dk)
    - [lcao\_dr](#lcao_dr)
    - [lcao\_rmax](#lcao_rmax)
    - [search\_radius](#search_radius)
    - [search\_pbc](#search_pbc)
    - [bx, by, bz](#bx-by-bz)
    - [elpa\_num\_thread](#elpa_num_thread)
    - [num\_stream](#num_stream)
  - [Electronic structure](#electronic-structure)
    - [basis\_type](#basis_type)
    - [ks\_solver](#ks_solver)
    - [nbands](#nbands)
    - [nelec](#nelec)
    - [nelec\_delta](#nelec_delta)
    - [nupdown](#nupdown)
    - [dft\_functional](#dft_functional)
    - [xc\_temperature](#xc_temperature)
    - [pseudo\_rcut](#pseudo_rcut)
    - [pseudo\_mesh](#pseudo_mesh)
    - [nspin](#nspin)
    - [smearing\_method](#smearing_method)
    - [smearing\_sigma](#smearing_sigma)
    - [smearing\_sigma\_temp](#smearing_sigma_temp)
    - [mixing\_type](#mixing_type)
    - [mixing\_beta](#mixing_beta)
    - [mixing\_beta\_mag](#mixing_beta_mag)
    - [mixing\_ndim](#mixing_ndim)
    - [mixing\_restart](#mixing_restart)
    - [mixing\_dmr](#mixing_dmr)
    - [mixing\_gg0](#mixing_gg0)
    - [mixing\_gg0\_mag](#mixing_gg0_mag)
    - [mixing\_gg0\_min](#mixing_gg0_min)
    - [mixing\_angle](#mixing_angle)
    - [mixing\_tau](#mixing_tau)
    - [mixing\_dftu](#mixing_dftu)
    - [gamma\_only](#gamma_only)
    - [printe](#printe)
    - [scf\_nmax](#scf_nmax)
    - [scf\_thr](#scf_thr)
    - [scf\_ene\_thr](#scf_ene_thr)
    - [scf\_thr\_type](#scf_thr_type)
    - [scf\_os\_stop](#scf_os_stop)
    - [scf\_os\_thr](#scf_os_thr)
    - [scf\_os\_ndim](#scf_os_ndim)
    - [chg\_extrap](#chg_extrap)
    - [lspinorb](#lspinorb)
    - [noncolin](#noncolin)
    - [soc\_lambda](#soc_lambda)
  - [Electronic structure (SDFT)](#electronic-structure-sdft)
    - [method\_sto](#method_sto)
    - [nbands\_sto](#nbands_sto)
    - [nche\_sto](#nche_sto)
    - [emin\_sto](#emin_sto)
    - [emax\_sto](#emax_sto)
    - [seed\_sto](#seed_sto)
    - [initsto\_ecut](#initsto_ecut)
    - [initsto\_freq](#initsto_freq)
    - [npart\_sto](#npart_sto)
  - [Geometry relaxation](#geometry-relaxation)
    - [relax\_method](#relax_method)
    - [relax\_new](#relax_new)
    - [relax\_scale\_force](#relax_scale_force)
    - [relax\_nmax](#relax_nmax)
    - [relax\_cg\_thr](#relax_cg_thr)
    - [cal\_force](#cal_force)
    - [force\_thr](#force_thr)
    - [force\_thr\_ev](#force_thr_ev)
    - [force\_thr\_ev2](#force_thr_ev2)
    - [relax\_bfgs\_w1](#relax_bfgs_w1)
    - [relax\_bfgs\_w2](#relax_bfgs_w2)
    - [relax\_bfgs\_rmax](#relax_bfgs_rmax)
    - [relax\_bfgs\_rmin](#relax_bfgs_rmin)
    - [relax\_bfgs\_init](#relax_bfgs_init)
    - [cal\_stress](#cal_stress)
    - [stress\_thr](#stress_thr)
    - [press1, press2, press3](#press1-press2-press3)
    - [fixed\_axes](#fixed_axes)
    - [fixed\_ibrav](#fixed_ibrav)
    - [fixed\_atoms](#fixed_atoms)
    - [cell\_factor](#cell_factor)
  - [Variables related to output information](#variables-related-to-output-information)
    - [out\_mul](#out_mul)
    - [out\_freq\_elec](#out_freq_elec)
    - [out\_chg](#out_chg)
    - [out\_pot](#out_pot)
    - [out\_dm](#out_dm)
    - [out\_dm1](#out_dm1)
    - [out\_wfc\_pw](#out_wfc_pw)
    - [out\_wfc\_r](#out_wfc_r)
    - [out\_wfc\_lcao](#out_wfc_lcao)
    - [out\_dos](#out_dos)
    - [out\_band](#out_band)
    - [out\_proj\_band](#out_proj_band)
    - [out\_stru](#out_stru)
    - [out\_bandgap](#out_bandgap)
    - [out\_level](#out_level)
    - [out\_alllog](#out_alllog)
    - [out\_mat\_hs](#out_mat_hs)
    - [out\_mat\_tk](#out_mat_tk)
    - [out\_mat\_r](#out_mat_r)
    - [out\_mat\_hs2](#out_mat_hs2)
    - [out\_mat\_t](#out_mat_t)
    - [out\_mat\_dh](#out_mat_dh)
    - [out\_mat\_xc](#out_mat_xc)
    - [out\_eband\_terms](#out_eband_terms)
    - [out\_hr\_npz/out\_dm\_npz](#out_hr_npzout_dm_npz)
    - [dm\_to\_rho](#dm_to_rho)
    - [out\_app\_flag](#out_app_flag)
    - [out\_ndigits](#out_ndigits)
    - [out\_interval](#out_interval)
    - [out\_element\_info](#out_element_info)
    - [restart\_save](#restart_save)
    - [rpa](#rpa)
    - [nbands\_istate](#nbands_istate)
    - [bands\_to\_print](#bands_to_print)
    - [if\_separate\_k](#if_separate_k)
    - [out\_elf](#out_elf)
  - [Density of states](#density-of-states)
    - [dos\_edelta\_ev](#dos_edelta_ev)
    - [dos\_sigma](#dos_sigma)
    - [dos\_scale](#dos_scale)
    - [dos\_emin\_ev](#dos_emin_ev)
    - [dos\_emax\_ev](#dos_emax_ev)
    - [dos\_nche](#dos_nche)
  - [NAOs](#naos)
    - [bessel\_nao\_ecut](#bessel_nao_ecut)
    - [bessel\_nao\_tolerence](#bessel_nao_tolerence)
    - [bessel\_nao\_rcut](#bessel_nao_rcut)
    - [bessel\_nao\_smooth](#bessel_nao_smooth)
    - [bessel\_nao\_sigma](#bessel_nao_sigma)
  - [DeePKS](#deepks)
    - [deepks\_out\_labels](#deepks_out_labels)
    - [deepks\_scf](#deepks_scf)
    - [deepks\_equiv](#deepks_equiv)
    - [deepks\_model](#deepks_model)
    - [bessel\_descriptor\_lmax](#bessel_descriptor_lmax)
    - [bessel\_descriptor\_ecut](#bessel_descriptor_ecut)
    - [bessel\_descriptor\_tolerence](#bessel_descriptor_tolerence)
    - [bessel\_descriptor\_rcut](#bessel_descriptor_rcut)
    - [bessel\_descriptor\_smooth](#bessel_descriptor_smooth)
    - [bessel\_descriptor\_sigma](#bessel_descriptor_sigma)
    - [deepks\_bandgap](#deepks_bandgap)
    - [deepks\_v\_delta](#deepks_v_delta)
    - [deepks\_out\_unittest](#deepks_out_unittest)
  - [OFDFT: orbital free density functional theory](#ofdft-orbital-free-density-functional-theory)
    - [of\_kinetic](#of_kinetic)
    - [of\_method](#of_method)
    - [of\_conv](#of_conv)
    - [of\_tole](#of_tole)
    - [of\_tolp](#of_tolp)
    - [of\_tf\_weight](#of_tf_weight)
    - [of\_vw\_weight](#of_vw_weight)
    - [of\_wt\_alpha](#of_wt_alpha)
    - [of\_wt\_beta](#of_wt_beta)
    - [of\_wt\_rho0](#of_wt_rho0)
    - [of\_hold\_rho0](#of_hold_rho0)
    - [of\_lkt\_a](#of_lkt_a)
    - [of\_read\_kernel](#of_read_kernel)
    - [of\_kernel\_file](#of_kernel_file)
    - [of\_full\_pw](#of_full_pw)
    - [of\_full\_pw\_dim](#of_full_pw_dim)
  - [Electric field and dipole correction](#electric-field-and-dipole-correction)
    - [efield\_flag](#efield_flag)
    - [dip\_cor\_flag](#dip_cor_flag)
    - [efield\_dir](#efield_dir)
    - [efield\_pos\_max](#efield_pos_max)
    - [efield\_pos\_dec](#efield_pos_dec)
    - [efield\_amp](#efield_amp)
  - [Gate field (compensating charge)](#gate-field-compensating-charge)
    - [gate\_flag](#gate_flag)
    - [zgate](#zgate)
    - [block](#block)
    - [block\_down](#block_down)
    - [block\_up](#block_up)
    - [block\_height](#block_height)
  - [Exact Exchange](#exact-exchange)
    - [exx\_hybrid\_alpha](#exx_hybrid_alpha)
    - [exx\_hse\_omega](#exx_hse_omega)
    - [exx\_separate\_loop](#exx_separate_loop)
    - [exx\_hybrid\_step](#exx_hybrid_step)
    - [exx\_mixing\_beta](#exx_mixing_beta)
    - [exx\_lambda](#exx_lambda)
    - [exx\_pca\_threshold](#exx_pca_threshold)
    - [exx\_c\_threshold](#exx_c_threshold)
    - [exx\_v\_threshold](#exx_v_threshold)
    - [exx\_dm\_threshold](#exx_dm_threshold)
    - [exx\_c\_grad\_threshold](#exx_c_grad_threshold)
    - [exx\_v\_grad\_threshold](#exx_v_grad_threshold)
    - [exx\_schwarz\_threshold](#exx_schwarz_threshold)
    - [exx\_cauchy\_threshold](#exx_cauchy_threshold)
    - [exx\_cauchy\_force\_threshold](#exx_cauchy_force_threshold)
    - [exx\_cauchy\_stress\_threshold](#exx_cauchy_stress_threshold)
    - [exx\_ccp\_threshold](#exx_ccp_threshold)
    - [exx\_ccp\_rmesh\_times](#exx_ccp_rmesh_times)
    - [exx\_distribute\_type](#exx_distribute_type)
    - [exx\_opt\_orb\_lmax](#exx_opt_orb_lmax)
    - [exx\_opt\_orb\_ecut](#exx_opt_orb_ecut)
    - [exx\_opt\_orb\_tolerence](#exx_opt_orb_tolerence)
    - [exx\_real\_number](#exx_real_number)
    - [rpa\_ccp\_rmesh\_times](#rpa_ccp_rmesh_times)
    - [exx\_symmetry\_realspace](#exx_symmetry_realspace)
    - [out\_ri\_cv](#out_ri_cv)
  - [Molecular dynamics](#molecular-dynamics)
    - [md\_type](#md_type)
    - [md\_nstep](#md_nstep)
    - [md\_dt](#md_dt)
    - [md\_thermostat](#md_thermostat)
    - [md\_tfirst, md\_tlast](#md_tfirst-md_tlast)
    - [md\_restart](#md_restart)
    - [md\_restartfreq](#md_restartfreq)
    - [md\_dumpfreq](#md_dumpfreq)
    - [dump\_force](#dump_force)
    - [dump\_vel](#dump_vel)
    - [dump\_virial](#dump_virial)
    - [md\_seed](#md_seed)
    - [md\_tfreq](#md_tfreq)
    - [md\_tchain](#md_tchain)
    - [md\_pmode](#md_pmode)
    - [ref\_cell\_factor](#ref_cell_factor)
    - [md\_pcouple](#md_pcouple)
    - [md\_pfirst, md\_plast](#md_pfirst-md_plast)
    - [md\_pfreq](#md_pfreq)
    - [md\_pchain](#md_pchain)
    - [lj\_rule](#lj_rule)
    - [lj\_eshift](#lj_eshift)
    - [lj\_rcut](#lj_rcut)
    - [lj\_epsilon](#lj_epsilon)
    - [lj\_sigma](#lj_sigma)
    - [pot\_file](#pot_file)
    - [dp\_rescaling](#dp_rescaling)
    - [dp\_fparam](#dp_fparam)
    - [dp\_aparam](#dp_aparam)
    - [msst\_direction](#msst_direction)
    - [msst\_vel](#msst_vel)
    - [msst\_vis](#msst_vis)
    - [msst\_tscale](#msst_tscale)
    - [msst\_qmass](#msst_qmass)
    - [md\_damp](#md_damp)
    - [md\_tolerance](#md_tolerance)
    - [md\_nraise](#md_nraise)
    - [cal\_syns](#cal_syns)
    - [dmax](#dmax)
  - [DFT+*U* correction](#dftu-correction)
    - [dft\_plus\_u](#dft_plus_u)
    - [orbital\_corr](#orbital_corr)
    - [hubbard\_u](#hubbard_u)
    - [yukawa\_potential](#yukawa_potential)
    - [yukawa\_lambda](#yukawa_lambda)
    - [uramping](#uramping)
    - [omc](#omc)
    - [onsite\_radius](#onsite_radius)
  - [vdW correction](#vdw-correction)
    - [vdw\_method](#vdw_method)
    - [vdw\_s6](#vdw_s6)
    - [vdw\_s8](#vdw_s8)
    - [vdw\_a1](#vdw_a1)
    - [vdw\_a2](#vdw_a2)
    - [vdw\_d](#vdw_d)
    - [vdw\_abc](#vdw_abc)
    - [vdw\_C6\_file](#vdw_c6_file)
    - [vdw\_C6\_unit](#vdw_c6_unit)
    - [vdw\_R0\_file](#vdw_r0_file)
    - [vdw\_R0\_unit](#vdw_r0_unit)
    - [vdw\_cutoff\_type](#vdw_cutoff_type)
    - [vdw\_cutoff\_radius](#vdw_cutoff_radius)
    - [vdw\_radius\_unit](#vdw_radius_unit)
    - [vdw\_cutoff\_period](#vdw_cutoff_period)
    - [vdw\_cn\_thr](#vdw_cn_thr)
    - [vdw\_cn\_thr\_unit](#vdw_cn_thr_unit)
  - [Berry phase and wannier90 interface](#berry-phase-and-wannier90-interface)
    - [berry\_phase](#berry_phase)
    - [gdir](#gdir)
    - [towannier90](#towannier90)
    - [nnkpfile](#nnkpfile)
    - [wannier\_method](#wannier_method)
    - [wannier\_spin](#wannier_spin)
    - [out\_wannier\_mmn](#out_wannier_mmn)
    - [out\_wannier\_amn](#out_wannier_amn)
    - [out\_wannier\_eig](#out_wannier_eig)
    - [out\_wannier\_unk](#out_wannier_unk)
    - [out\_wannier\_wvfn\_formatted](#out_wannier_wvfn_formatted)
  - [TDDFT: time dependent density functional theory](#tddft-time-dependent-density-functional-theory)
    - [td\_edm](#td_edm)
    - [td\_print\_eij](#td_print_eij)
    - [td\_propagator](#td_propagator)
    - [td\_vext](#td_vext)
    - [td\_vext\_dire](#td_vext_dire)
    - [td\_stype](#td_stype)
    - [td\_ttype](#td_ttype)
    - [td\_tstart](#td_tstart)
    - [td\_tend](#td_tend)
    - [td\_lcut1](#td_lcut1)
    - [td\_lcut2](#td_lcut2)
    - [td\_gauss\_freq](#td_gauss_freq)
    - [td\_gauss\_phase](#td_gauss_phase)
    - [td\_gauss\_sigma](#td_gauss_sigma)
    - [td\_gauss\_t0](#td_gauss_t0)
    - [td\_gauss\_amp](#td_gauss_amp)
    - [td\_trape\_freq](#td_trape_freq)
    - [td\_trape\_phase](#td_trape_phase)
    - [td\_trape\_t1](#td_trape_t1)
    - [td\_trape\_t2](#td_trape_t2)
    - [td\_trape\_t3](#td_trape_t3)
    - [td\_trape\_amp](#td_trape_amp)
    - [td\_trigo\_freq1](#td_trigo_freq1)
    - [td\_trigo\_freq2](#td_trigo_freq2)
    - [td\_trigo\_phase1](#td_trigo_phase1)
    - [td\_trigo\_phase2](#td_trigo_phase2)
    - [td\_trigo\_amp](#td_trigo_amp)
    - [td\_heavi\_t0](#td_heavi_t0)
    - [td\_heavi\_amp](#td_heavi_amp)
    - [out\_dipole](#out_dipole)
    - [out\_current](#out_current)
    - [out\_current\_k](#out_current_k)
    - [out\_efield](#out_efield)
    - [out\_vecpot](#out_vecpot)
    - [init\_vecpot\_file](#init_vecpot_file)
    - [ocp](#ocp)
    - [ocp\_set](#ocp_set)
  - [Variables useful for debugging](#variables-useful-for-debugging)
    - [t\_in\_h](#t_in_h)
    - [vl\_in\_h](#vl_in_h)
    - [vnl\_in\_h](#vnl_in_h)
    - [vh\_in\_h](#vh_in_h)
    - [vion\_in\_h](#vion_in_h)
    - [test\_force](#test_force)
    - [test\_stress](#test_stress)
    - [test\_skip\_ewald](#test_skip_ewald)
  - [Electronic conductivities](#electronic-conductivities)
    - [cal\_cond](#cal_cond)
    - [cond\_che\_thr](#cond_che_thr)
    - [cond\_dw](#cond_dw)
    - [cond\_wcut](#cond_wcut)
    - [cond\_dt](#cond_dt)
    - [cond\_dtbatch](#cond_dtbatch)
    - [cond\_smear](#cond_smear)
    - [cond\_fwhm](#cond_fwhm)
    - [cond\_nonlocal](#cond_nonlocal)
  - [Implicit solvation model](#implicit-solvation-model)
    - [imp\_sol](#imp_sol)
    - [eb\_k](#eb_k)
    - [tau](#tau)
    - [sigma\_k](#sigma_k)
    - [nc\_k](#nc_k)
  - [Quasiatomic Orbital (QO) analysis](#quasiatomic-orbital-qo-analysis)
    - [qo\_switch](#qo_switch)
    - [qo\_basis](#qo_basis)
    - [qo\_strategy](#qo_strategy)
    - [qo\_screening\_coeff](#qo_screening_coeff)
    - [qo\_thr](#qo_thr)
  - [PEXSI](#pexsi)
    - [pexsi\_npole](#pexsi_npole)
    - [pexsi\_inertia](#pexsi_inertia)
    - [pexsi\_nmax](#pexsi_nmax)
    - [pexsi\_comm](#pexsi_comm)
    - [pexsi\_storage](#pexsi_storage)
    - [pexsi\_ordering](#pexsi_ordering)
    - [pexsi\_row\_ordering](#pexsi_row_ordering)
    - [pexsi\_nproc](#pexsi_nproc)
    - [pexsi\_symm](#pexsi_symm)
    - [pexsi\_trans](#pexsi_trans)
    - [pexsi\_method](#pexsi_method)
    - [pexsi\_nproc\_pole](#pexsi_nproc_pole)
    - [pexsi\_temp](#pexsi_temp)
    - [pexsi\_gap](#pexsi_gap)
    - [pexsi\_delta\_e](#pexsi_delta_e)
    - [pexsi\_mu\_lower](#pexsi_mu_lower)
    - [pexsi\_mu\_upper](#pexsi_mu_upper)
    - [pexsi\_mu](#pexsi_mu)
    - [pexsi\_mu\_thr](#pexsi_mu_thr)
    - [pexsi\_mu\_expand](#pexsi_mu_expand)
    - [pexsi\_mu\_guard](#pexsi_mu_guard)
    - [pexsi\_elec\_thr](#pexsi_elec_thr)
    - [pexsi\_zero\_thr](#pexsi_zero_thr)
  - [Linear Response TDDFT](#linear-response-tddft)
    - [xc\_kernel](#xc_kernel)
    - [lr\_init\_xc\_kernel](#lr_init_xc_kernel)
    - [lr\_solver](#lr_solver)
    - [lr\_thr](#lr_thr)
    - [nocc](#nocc)
    - [nvirt](#nvirt)
    - [lr\_nstates](#lr_nstates)
    - [lr\_unrestricted](#lr_unrestricted)
    - [abs\_wavelen\_range](#abs_wavelen_range)
    - [out\_wfc\_lr](#out_wfc_lr)
    - [abs\_broadening](#abs_broadening)
    - [ri\_hartree\_benchmark](#ri_hartree_benchmark)
    - [aims\_nbasis](#aims_nbasis)
  - [Reduced Density Matrix Functional Theory](#Reduced-Density-Matrix-Functional-Theory)
    - [rdmft](#rdmft)
    - [rdmft\_power\_alpha](#rdmft_power_alpha)

[back to top](#full-list-of-input-keywords)
## System variables

These variables are used to control general system parameters.

### suffix

- **Type**: String
- **Description**: In each run, ABACUS will generate a subdirectory in the working directory. This subdirectory contains all the information of the run. The subdirectory name has the format: OUT.suffix, where the `suffix` is the name you can pick up for your convenience.
- **Default**: ABACUS

### calculation

- **Type**: String
- **Description**: Specify the type of calculation.

  - **scf**: perform self-consistent electronic structure calculations
  - **nscf**: perform non-self-consistent electronic structure calculations. A charge density file is required
  - **relax**: perform structure relaxation calculations, the `relax_nmax` parameter depicts the maximal number of ionic iterations
  - **cell-relax**: perform cell relaxation calculations
  - **md**: perform molecular dynamics simulations
  - **get_pchg**: obtain partial (band-decomposed) charge densities (for LCAO basis only). See `nbands_istate` and `bands_to_print` for more information
  - **get_wf**: obtain wave functions (for LCAO basis only). See `nbands_istate` and `bands_to_print` for more information
  - **get_S** : obtain the overlap matrix formed by localized orbitals (for LCAO basis with multiple k points). the file name is `SR.csr` with file format being the same as that generated by [out_mat_hs2](#out_mat_hs2)
  - **gen_bessel** : generates projectors, i.e., a series of Bessel functions, for the DeePKS method (for LCAO basis only); see also keywords `bessel_descriptor_lmax`, `bessel_descriptor_rcut` and `bessel_descriptor_tolerence`. A file named `jle.orb` will be generated which contains the projectors. An example is provided in examples/H2O-deepks-pw
  - **test_memory** : obtain a rough estimation of memory consuption for the calculation
  - **test_neighbour** : obtain information of neighboring atoms (for LCAO basis only), please specify a positive [search_radius](#search_radius) manually
- **Default**: scf

### esolver_type

- **Type**: String
- **Description**: choose the energy solver.
  - ksdft: Kohn-Sham density functional theory
  - ofdft: orbital-free density functional theory
  - sdft: [stochastic density functional theory](#electronic-structure-sdft)
  - tddft: real-time time-dependent density functional theory (TDDFT)
  - lj: Leonard Jones potential
  - dp: DeeP potential, see details in [md.md](../md.md#dpmd)
  - ks-lr: Kohn-Sham density functional theory + LR-TDDFT
  - lr: LR-TDDFT with given KS orbitals
- **Default**: ksdft

### symmetry

- **Type**: Integer
- **Description**: takes value 1, 0 or -1.
  - -1: No symmetry will be considered. It is recommended to set -1 for non-colinear + soc calculations, where time reversal symmetry is broken sometimes.
  - 0: Only time reversal symmetry would be considered in symmetry operations, which implied k point and -k point would be treated as a single k point with twice the weight.
  - 1: Symmetry analysis will be performed to determine the type of Bravais lattice and associated symmetry operations. (point groups, space groups, primitive cells, and irreducible k-points)
- **Default**:
  - 0:
    - if *[calculation](#calculation)==md/nscf/get_pchg/get_wf/get_S* or *[gamma_only](#gamma_only)==True*;
    - If (*[dft_fuctional](#dft_functional)==hse/hf/pbe0/scan0/opt_orb* or *[rpa](#rpa)==True*). Currently *symmetry==1* is not supported in EXX (exact exchange) calculation.
    - If *[efield_flag](#efield_flag)==1*
  - 1: else

### symmetry_prec

- **Type**: Real
- **Description**: The accuracy for symmetry judgment. Usually the default value is good enough, but if the lattice parameters or atom positions in STRU file is not accurate enough, this value should be enlarged.
  > Note: if *[calculation](#calculation)==cell_relax*, this value can be dynamically changed corresponding to the variation of accuracy of the lattice parameters and atom positions during the relaxation. The new value will be printed in `OUT.${suffix}/running_cell-relax.log` in that case.
- **Default**: 1.0e-6
- **Unit**:  Bohr

### symmetry_autoclose

- **Type**: Boolean
- **Availability**: *[symmetry](#symmetry)==1*
- **Description**: Control how to deal with error in symmetry analysis due to inaccurate lattice parameters or atom positions in STRU file, especially useful when *[calculation](#calculation)==cell-relax*
  - False: quit with an error message
  - True: automatically set symmetry to 0 and continue running without symmetry analysis
- **Default**: True

### kpar

- **Type**: Integer
- **Description**: divide all processors into kpar groups, and k points will be distributed among each group. The value taken should be less than or equal to the number of k points as well as the number of MPI processes.
- **Default**: 1

### bndpar

- **Type**: Integer
- **Description**: divide all processors into bndpar groups, and bands (only stochastic orbitals now) will be distributed among each group. It should be larger than 0.
- **Default**: 1

### latname

- **Type**: String
- **Description**: Specifies the type of Bravias lattice. When set to `none`, the three lattice vectors are supplied explicitly in STRU file. When set to a certain Bravais lattice type, there is no need to provide lattice vector, but a few lattice parameters might be required. For more information regarding this parameter, consult the [page on STRU file](stru.md).

  Available options are (correspondence with ibrav in QE(Quantum Espresso) is given in parenthesis):
  - none: free structure
  - sc: simple cubic (1)
  - fcc: face-centered cubic (2)
  - bcc: body-centered cubic (3)
  - hexagonal: hexagonal (4)
  - trigonal: trigonal (5)
  - st: simple tetragonal (6)
  - bct: body-centered tetragonal (7)
  - so: orthorhombic (8)
  - baco: base-centered orthorhombic (9)
  - fco: face-centered orthorhombic (10)
  - bco: body-centered orthorhombic (11)
  - sm: simple monoclinic (12)
  - bacm: base-centered monoclinic (13)
  - triclinic: triclinic (14)
- **Default**: none

### psi_initializer

- **Type**: Integer
- **Description**: enable the experimental feature psi_initializer, to support use numerical atomic orbitals initialize wavefunction (`basis_type pw` case).

  NOTE: this feature is not well-implemented for `nspin 4` case (closed presently), and cannot use with `calculation nscf`/`esolver_type sdft` cases.
  Available options are:
  - 0: disable psi_initializer
  - 1: enable psi_initializer
- **Default**: 0

### init_wfc

- **Type**: String
- **Description**: Only useful for plane wave basis only now. It is the name of the starting wave functions. In the future. we should also make this variable available for localized orbitals set.

  Available options are:

  - atomic: from atomic pseudo wave functions. If they are not enough, other wave functions are initialized with random numbers.
  - atomic+random: add small random numbers on atomic pseudo-wavefunctions
  - file: from binary files `WAVEFUNC*.dat`, which are output by setting [out_wfc_pw](#out_wfc_pw) to `2`.
  - random: random numbers

  with `psi_initializer 1`, two more options are supported:
  - nao: from numerical atomic orbitals. If they are not enough, other wave functions are initialized with random numbers.
  - nao+random: add small random numbers on numerical atomic orbitals
- **Default**: atomic

### init_chg

- **Type**: String
- **Description**: This variable is used for both plane wave set and localized orbitals set. It indicates the type of starting density.

  - atomic: the density is starting from the summation of the atomic density of single atoms.
  - file: the density will be read in from a binary file `charge-density.dat` first. If it does not exist, the charge density will be read in from cube files. Besides, when you do `nspin=1` calculation, you only need the density file SPIN1_CHG.cube. However, if you do `nspin=2` calculation, you also need the density file SPIN2_CHG.cube. The density file should be output with these names if you set out_chg = 1 in INPUT file.
  - wfc: the density will be calculated by wavefunctions and occupations. Wavefunctions are read in from binary files `WAVEFUNC*.dat` while occupations are read in from file `istate.info`.
  - auto: Abacus first attempts to read the density from a file; if not found, it defaults to using atomic density.
- **Default**: atomic

### init_vel

- **Type**: Boolean
- **Description**:

  - True: read the atom velocity (atomic unit : 1 a.u. = 21.877 Angstrom/fs) from the atom file (`STRU`) and determine the initial temperature [md_tfirst](#md_tfirst-md_tlast).  If [md_tfirst](#md_tfirst-md_tlast) is unset or less than zero, `init_vel` is autoset to be `true`.
  - False: assign value to atom velocity using Gaussian distributed random numbers.
- **Default**: False

### mem_saver

- **Type**: Boolean
- **Description**: Used only for nscf calculations.
  - 0: no memory saving techniques are used.
  - 1: a memory saving technique will be used for many k point calculations.

- **Default**: 0

### diago_proc

- **Type**: Integer
- **Availability**: pw base
- **Description**:
  - 0: it will be set to the number of MPI processes. Normally, it is fine just leave it to the default value.
  - `>0`: it specifies the number of processes used for carrying out diagonalization. Must be less than or equal to total number of MPI processes. Also, when cg diagonalization is used, diago_proc must be the same as the total number of MPI processes.
- **Default**: 0

### nbspline

- **Type**: Integer
- **Description**: If set to a natural number, a Cardinal B-spline interpolation will be used to calculate Structure Factor. `nbspline` represents the order of B-spline basis and a larger one can get more accurate results but cost more.
  It is turned off by default.
- **Default**: -1

### kspacing

- **Type**: Real
- **Description**: Set the smallest allowed spacing between k points, unit in 1/bohr. It should be larger than 0.0, and suggest smaller than 0.25. When you have set this value > 0.0, then the KPT file is unnecessary, and the number of K points nk_i = max(1, int(|b_i|/KSPACING_i)+1), where b_i is the reciprocal lattice vector. The default value 0.0 means that ABACUS will read the applied KPT file.
If only one value is set (such as `kspacing 0.5`), then kspacing values of a/b/c direction are all set to it; and one can also set 3 values to set the kspacing value for a/b/c direction separately (such as: `kspacing 0.5 0.6 0.7`).

  Note: if gamma_only is set to be true, kspacing is invalid.
- **Default**: 0.0

### min_dist_coef

- **Type**: Real
- **Description**: a factor related to the allowed minimum distance between two atoms. At the beginning, ABACUS will check the structure, and if the distance of two atoms is shorter than min_dist_coef*(standard covalent bond length), we think this structure is unreasonable. If you want to calculate some structures in extreme conditions like high pressure, you should set this parameter as a smaller value or even 0.
- **Default**: 0.2

### device

- **Type**: String
- **Description**: Specifies the computing device for ABACUS.

  Available options are:

  - cpu: for CPUs via Intel, AMD, or Other supported CPU devices
  - gpu: for GPUs via CUDA or ROCm.

  Known limitations: `ks_solver` must also be set to the algorithms supported. lcao_in_pw currently does not support `gpu`.

- **Default**: cpu

### precision

- **Type**: String
- **Description**: Specifies the precision of the PW_SCF calculation.

  Available options are:

  - single: single precision
  - double: double precision

  Known limitations:

  - pw basis: required by the `single` precision options
  - cg/bpcg/dav ks_solver: required by the `single` precision options
- **Default**: double

[back to top](#full-list-of-input-keywords)

## Variables related to input files

These variables are used to control parameters related to input files.

### stru_file

- **Type**: String
- **Description**: the name of the structure file
  - Containing various information about atom species, including pseudopotential files, local orbitals files, cell information, atom positions, and whether atoms should be allowed to move.
  - When [calculation](#calculation) is set to `md` and [md_restart](#md_restart) is set to `true`, this keyword will NOT work.
  - Refer to [Doc](https://github.com/deepmodeling/abacus-develop/blob/develop/docs/advanced/input_files/stru.md)
- **Default**: STRU

### kpoint_file

- **Type**: String
- **Description**: the name of the k-points file
  - In atomic orbitals basis with `gamma_only` set to true, the `KPT` file is unnecessary, because a `KPT` file will be generated automatically.
  - When more than one k-points are required, an explicit `KPT` file is mandatory.
  - Refer to [Doc](https://github.com/deepmodeling/abacus-develop/blob/develop/docs/advanced/input_files/kpt.md)
- **Default**: KPT

### pseudo_dir

- **Type**: String
- **Description**: the pseudopotential file directory
  - This parameter is combined with the pseudopotential filenames in the STRU file to form the complete pseudopotential file paths.
  - Example: set pseudo_dir to "../" with "Si.upf" which specified under "ATOMIC_SPECIES" in STRU file, ABACUS will open the pseudopotential file in path "../Si.upf".
- **Default**: ""

### orbital_dir

- **Type**: String
- **Description**: the orbital file directory
  - This parameter is combined with orbital filenames in the STRU file to form the complete orbital file paths.
  - Example: set orbital_dir to "../" with "Si.orb" which specified under "NUMERICAL_ORBITAL" in STRU file, ABACUS will open the orbital file in path "../Si.orb".
- **Default**: ""

### read_file_dir

- **Type**: String
- **Description**: Indicates the location of files, such as electron density (`SPIN1_CHG.cube`), required as a starting point.
  - Example: './' implies the files to be read are located in the working directory.
- **Default**: OUT.$suffix

### restart_load

- **Type**: Boolean
- **Availability**: Numerical atomic orbital basis
- **Description**: If [restart_save](#restart_save) is set to true and an electronic iteration is finished, calculations can be restarted from the charge density file, which are saved in the former calculation. Please ensure [read_file_dir](#read_file_dir) is correct, and  the charge density file exist.

  If EXX(exact exchange) is calculated (i.e. *[dft_fuctional](#dft_functional)==hse/hf/pbe0/scan0/opt_orb* or *[rpa](#rpa)==True*), the Hexx(R) files in the same folder for each processor will also be read.
- **Default**: False

### wannier_card

- **Type**: String
- **Availability**: Using ABACUS with Wannier90.
- **Description**: The name of the input file related to Wannier90.
- **Default**: "none"

[back to top](#full-list-of-input-keywords)

## Plane wave related variables

These variables are used to control the plane wave related parameters.

### ecutwfc

- **Type**: Real
- **Description**: Energy cutoff for plane wave functions, the unit is **Rydberg**. Note that even for localized orbitals basis, you still need to setup an energy cutoff for this system. Because our local pseudopotential parts and the related force are calculated from plane wave basis set, etc. Also, because our orbitals are generated by matching localized orbitals to a chosen set of wave functions from a certain energy cutoff, this set of localize orbitals is most accurate under this same plane wave energy cutoff.
- **Default**: 50 Ry (PW basis), 100 Ry (LCAO basis)

### ecutrho

- **Type**: Real
- **Description**: Energy cutoff for charge density and potential, the unit is **Rydberg**. For norm-conserving pseudopotential you should stick to the default value, you can reduce it by a little but it will introduce noise especially on forces and stress. For ultrasoft pseudopotential a larger value than the default is often desirable (`ecutrho` = 8 to 12 times `ecutwfc`, typically). The use of gradient-corrected functional, especially in cells with vacuum, or for pseudopotential without non-linear core correction, usually requires an higher values of `ecutrho` to be accurately converged.
- **Default**: 4*ecutwfc

### nx, ny, nz

- **Type**: Integer
- **Description**: If set to a positive number, then the three variables specify the numbers of FFT grid points in x, y, z directions, respectively. If set to 0, the number will be calculated from ecutrho.

    Note: You must specify all three dimensions for this setting to be used.
- **Default**: 0

### ndx, ndy, ndz

- **Type**: Integer
- **Description**: If set to a positive number, then the three variables specify the numbers of FFT grid (for the dense part of charge density in ultrasoft pseudopotential) points in x, y, z directions, respectively. If set to 0, the number will be calculated from ecutwfc.

    Note: You must specify all three dimensions for this setting to be used.

    Note: These parameters must be used combined with [nx,ny,nz](#nx-ny-nz). If [nx,ny,nz](#nx-ny-nz) are unset, ndx,ndy,ndz are used as [nx,ny,nz](#nx-ny-nz).
- **Default**: 0

### pw_seed

- **Type**: Integer
- **Description**: Only useful for plane wave basis only now. It is the random seed to initialize wave functions. Only positive integers are available.
- **Default**:0

### pw_diag_thr

- **Type**: Real
- **Description**: Only used when you use `ks_solver = cg/dav/dav_subspace/bpcg`. It indicates the threshold for the first electronic iteration, from the second iteration the pw_diag_thr will be updated automatically. **For nscf calculations with planewave basis set, pw_diag_thr should be <= 1e-3.**
- **Default**: 0.01

### diago_smooth_ethr

- **Type**: bool
- **Description**: If `TRUE`, the smooth threshold strategy, which applies a larger threshold (10e-5) for the empty states, will be implemented in the diagonalization methods. (This strategy should not affect total energy, forces, and other ground-state properties, but computational efficiency will be improved.) If `FALSE`, the smooth threshold strategy will not be applied.
- **Default**: false

### pw_diag_nmax

- **Type**: Integer
- **Description**: Only useful when you use `ks_solver = cg/dav/dav_subspace/bpcg`. It indicates the maximal iteration number for cg/david/dav_subspace/bpcg method.
- **Default**: 40

### pw_diag_ndim

- **Type**: Integer
- **Description**: Only useful when you use `ks_solver = dav` or `ks_solver = dav_subspace`. It indicates dimension of workspace(number of wavefunction packets, at least 2 needed) for the Davidson method. A larger value may yield a smaller number of iterations in the algorithm but uses more memory and more CPU time in subspace diagonalization.
- **Default**: 4

### erf_ecut

- **Type**: Real
- **Description**: Used in variable-cell molecular dynamics (or in stress calculation). See [erf_sigma](#erf_sigma) in detail.
- **Default**: 0.0
- **Unit**: Ry

### fft_mode

- **Type**: Integer
- **Description**: Set the mode of FFTW.
  - 0: FFTW_ESTIMATE
  - 1: FFTW_MEASURE
  - 2: FFTW_PATIENT
  - 3: FFTW_EXHAUSTIVE
- **Default**: 0

### erf_height

- **Type**: Real
- **Description**: Used in variable-cell molecular dynamics (or in stress calculation). See [erf_sigma](#erf_sigma) in detail.
- **Default**: 0.0
- **Unit**: Ry

### erf_sigma

- **Type**: Real
- **Description**: In order to recover the accuracy of a constant energy cutoff calculation, the kinetic functional is modified, which is used in variable-cell molecular dynamics (or in stress calculation).

  [erf_ecut](#erf_ecut) is the value of the constant energy cutoff; [erf_height](#erf_height) and [erf_sigma](#erf_sigma) are the height and the width of the energy step for reciprocal vectors whose square modulus is greater than [erf_ecut](#erf_ecut). In the kinetic energy, G^2 is replaced by G^2 + erf_height * (1 + erf ( (G^2 - erf_ecut)/erf_sigma) )

  See: M. Bernasconi et al., J. Phys. Chem. Solids **56**, 501 (1995), [doi:10.1016/0022-3697(94)00228-2](#https://doi.org/10.1016/0022-3697(94)00228-2)
- **Default**: 0.1
- **Unit**: Ry

[back to top](#full-list-of-input-keywords)

## Numerical atomic orbitals related variables

These variables are used to control the numerical atomic orbitals related parameters.

### nb2d

- **Type**: Integer
- **Description**: In LCAO calculations, we arrange the total number of processors in an 2D array, so that we can partition the wavefunction matrix (number of bands*total size of atomic orbital basis) and distribute them in this 2D array. When the system is large, we group processors into sizes of nb2d, so that multiple processors take care of one row block (a group of atomic orbitals) in the wavefunction matrix. If set to 0, nb2d will be automatically set in the program according to the size of atomic orbital basis:
  - if size <= 500 : nb2d = 1
  - if 500 < size <= 1000 : nb2d = 32
  - if size > 1000 : nb2d = 64;
- **Default**: 0

### lmaxmax

- **Type**: Integer
- **Description**: If not equals to 2, then the maximum l channels on LCAO is set to lmaxmax. If 2, then the number of l channels will be read from the LCAO data sets. Normally no input should be supplied for this variable so that it is kept as its default.
- **Default**: 2.

### lcao_ecut

- **Type**: Real
- **Description**: Energy cutoff (in Ry) for two-center integrals in LCAO. The two-center integration table are obtained via a k space integral whose upper limit is about sqrt(`lcao_ecut`).
- **Default**: `ecutwfc`

### lcao_dk

- **Type**: Real
- **Description**: k spacing (in Bohr${}^{-1}$) for two-center integrals. The two-center integration table are obtained via a k space integral on a uniform grid with spacing `lcao_dk`.
- **Default**: 0.01

### lcao_dr

- **Type**: Real
- **Description**: r spacing (in Bohr) of the integration table of two-center integrals.
- **Default**: 0.01

### lcao_rmax

- **Type**: Real
- **Description**: Maximum distance (in Bohr) for the two-center integration table.
- **Default**: 30

### search_radius

- **Type**: Real
- **Description**: Searching radius in finding the neighbouring atoms. By default the radius will be automatically determined by the cutoffs of orbitals and nonlocal beta projectors.
- **Default**: -1
- **Unit**: Bohr

### search_pbc

- **Type**: Boolean
- **Description**: If True, periodic images will be included in searching for the neighbouring atoms. If False, periodic images will be ignored.
- **Default**: True

### bx, by, bz

- **Type**: Integer
- **Description**: In the matrix operation of grid integral, bx/by/bz grids (in x, y, z directions) are treated as a whole as a matrix element. A different value will affect the calculation speed. The default is 0, which means abacus will automatically calculate these values.
- **Default**: 0

### elpa_num_thread

- **Type**: int
- **Description**: Number of threads used in one elpa calculation. 

  If the number is below 0 or 0 or beyond the max number of threads, all elpa calculation will be using all mpi threads
- **Default**: -1

### num_stream

- **Type** :int
- **Description**: choose the number of streams in GPU when we compute the `LCAO`. According to different devices , we may have different effects.For most devices,the stream is
enough when the number is bigger then 2.
- **Default** : "4" 

[back to top](#full-list-of-input-keywords)

## Electronic structure

These variables are used to control the electronic structure and geometry relaxation
calculations.

### basis_type

- **Type**: String
- **Description**: Choose the basis set.
  - **pw**: Using plane-wave basis set only.
  - **lcao**: Using localized atomic orbital sets.
  - **lcao_in_pw**: Expand the localized atomic set in plane-wave basis, non-self-consistent field calculation not tested.
- **Default**: pw

### ks_solver

- **Type**: String
- **Description**: Choose the diagonalization methods for the Hamiltonian matrix expanded in a certain basis set.

  For plane-wave basis,

  - **cg**: cg method.
  - **bpcg**: bpcg method, which is a block-parallel Conjugate Gradient (CG) method, typically exhibits higher acceleration in a GPU environment.
  - **dav**: the Davidson algorithm.
  - **dav_subspace**: Davidson algorithm without orthogonalization operation, this method is the most recommended for efficiency. `pw_diag_ndim` can be set to 2 for this method.

  For atomic orbitals basis,

  - **lapack**: This method is only avaliable for serial version. For parallel version please use **scalapack_gvx**.
  - **genelpa**: This method should be used if you choose localized orbitals.
  - **scalapack_gvx**: Scalapack can also be used for localized orbitals.
  - **cusolver**: This method needs building with CUDA and at least one gpu is available.
  - **cusolvermp**: This method supports multi-GPU acceleration and needs building with CUDA。 Note that when using cusolvermp, you should set the number of MPI processes to be equal to the number of GPUs.
  - **elpa**: The ELPA solver supports both CPU and GPU. By setting the `device` to GPU, you can launch the ELPA solver with GPU acceleration (provided that you have installed a GPU-supported version of ELPA, which requires you to manually compile and install ELPA, and the ABACUS should be compiled with -DUSE_ELPA=ON and -DUSE_CUDA=ON). The ELPA solver also supports multi-GPU acceleration.

  If you set ks_solver=`genelpa` for basis_type=`pw`, the program will be stopped with an error message:

  ```text
  genelpa can not be used with plane wave basis.
  ```

  Then the user has to correct the input file and restart the calculation.
- **Default**: 
  - **PW basis**: cg.
  - **LCAO basis**:
    - genelpa (if compiling option `USE_ELPA` has been set)
    - lapack (if compiling option `ENABLE_MPI` has not been set)
    - scalapack_gvx (if compiling option `USE_ELPA` has not been set and compiling option `ENABLE_MPI` has been set)
    - cusolver (if compiling option `USE_CUDA` has been set)

### nbands

- **Type**: Integer
- **Description**: The number of Kohn-Sham orbitals to calculate. It is recommended to setup this value, especially when smearing techniques are utilized, more bands should be included.
- **Default**:
  - nspin=1: max(1.2\*occupied_bands, occupied_bands + 10)
  - nspin=2: max(1.2\*nelec_spin, nelec_spin + 10), in which nelec_spin = max(nelec_spin_up, nelec_spin_down)
  - nspin=4: max(1.2\*nelec, nelec + 20)

### nelec

- **Type**: Real
- **Description**:

  - 0.0: the total number of electrons will be calculated by the sum of valence electrons (i.e. assuming neutral system).
  - `>0.0`: this denotes the total number of electrons in the system. Must be less than 2*nbands.
- **Default**: 0.0

### nelec_delta

- **Type**: Real
- **Description**:
 the total number of electrons will be calculated by `nelec`+`nelec_delta`.
- **Default**: 0.0

### nupdown

- **Type**: Real
- **Description**:
  - 0.0: no constrain apply to system.
  - `>0.0`: this denotes the difference number of electrons between spin-up and spin-down in the system. The range of value must in [-nelec ~ nelec]. It is one method of constraint DFT, the fermi energy level will separate to E_Fermi_up and E_Fermi_down.
- **Default**: 0.0

### dft_functional

- **Type**: String
- **Description**: In our package, the XC functional can either be set explicitly using the `dft_functional` keyword in `INPUT` file. If `dft_functional` is not specified, ABACUS will use the xc functional indicated in the pseudopotential file.
  On the other hand, if dft_functional is specified, it will overwrite the functional from pseudopotentials and performs calculation with whichever functional the user prefers. We further offer two ways of supplying exchange-correlation functional. The first is using 'short-hand' names such as 'LDA', 'PBE', 'SCAN'. A complete list of 'short-hand' expressions can be found in [the source code](../../../source/module_hamilt_general/module_xc/xc_functional.cpp). The other way is only available when ***compiling with LIBXC***, and it allows for supplying exchange-correlation functionals as combinations of LIBXC keywords for functional components, joined by a plus sign, for example, dft_functional='LDA_X_1D_EXPONENTIAL+LDA_C_1D_CSC'. The list of LIBXC keywords can be found on its [website](https://libxc.gitlab.io/functionals/). In this way, **we support all the LDA,GGA and mGGA functionals provided by LIBXC**.

  Furthermore, the old INPUT parameter exx_hybrid_type for hybrid functionals has been absorbed into dft_functional. Options are `hf` (pure Hartree-Fock), `pbe0`(PBE0), `hse` (Note: in order to use HSE functional, LIBXC is required). Note also that HSE has been tested while PBE0 has NOT been fully tested yet, and the maximum CPU cores for running exx in parallel is $N(N+1)/2$, with N being the number of atoms. And forces for hybrid functionals are not supported yet.

  If set to `opt_orb`, the program will not perform hybrid functional calculation. Instead, it is going to generate opt-ABFs as discussed in this [article](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.0c00481).
- **Default**: same as UPF file.

### xc_temperature

- **Type**: Real
- **Description**: specifies temperature when using temperature-dependent XC functionals (KSDT and so on).
- **Default** : 0.0
- **Unit**: Ry

### pseudo_rcut

- **Type**: Real
- **Description**: Cut-off of radial integration for pseudopotentials
- **Default**: 15
- **Unit**: Bohr

### pseudo_mesh

- **Type**: Integer
- **Description**:
  - 0: use our own mesh for radial integration of pseudopotentials
  - 1: use the mesh that is consistent with quantum espresso
- **Default**: 0

### nspin

- **Type**: Integer
- **Description**: The number of spin components of wave functions.
  - **1**: Spin degeneracy
  - **2**: Collinear spin polarized.
  - **4**: For the case of [noncollinear polarized](../scf/spin.md#noncollinear-spin-polarized-calculations), nspin will be automatically set to 4 without being specified by the user.
- **Default**: 1

### smearing_method

- **Type**: String
- **Description**: It indicates which occupation and smearing method is used in the calculation.
  - **fixed**: fixed occupations (available for non-coductors only)
  - **gauss** or **gaussian**: Gaussian smearing method.
  - **mp**: methfessel-paxton smearing method; recommended for metals.
  - **mp2**: 2-nd methfessel-paxton smearing method; recommended for metals.
  - **mv** or **cold**: marzari-vanderbilt smearing method.
  - **fd**: Fermi-Dirac smearing method: $f=1/\{1+\exp[(E-\mu)/kT]\}$ and smearing_sigma below is the temperature $T$ (in Ry).
- **Default**: gauss

### smearing_sigma

- **Type**: Real
- **Description**: Energy range for smearing.
- **Default**: 0.015
- **Unit**: Ry

### smearing_sigma_temp

- **Type**: Real
- **Description**: Energy range for smearing, `smearing_sigma` = 1/2 *kB* `smearing_sigma_temp`.
- **Default**: 2 * `smearing_sigma` / kB.
- **Unit**: K

### mixing_type

- **Type**: String
- **Description**: Charge mixing methods.
  - **plain**: Just simple mixing.
  - **pulay**: Standard Pulay method. [P. Pulay Chemical Physics Letters, (1980)](https://www.sciencedirect.com/science/article/abs/pii/0009261480803964)
  - **broyden**: Simplified modified Broyden method. [D.D. Johnson Physical Review B (1988)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.38.12807)

  In general, the convergence of the Broyden method is slightly faster than that of the Pulay method.
- **Default**: broyden

### mixing_beta

- **Type**: Real
- **Description**: In general, the formula of charge mixing can be written as $\rho_{new} = \rho_{old} + \beta * \rho_{update}$, where $\rho_{new}$ represents the new charge density after charge mixing, $\rho_{old}$ represents the charge density in previous step, $\rho_{update}$ is obtained through various mixing methods, and $\beta$ is set by the parameter `mixing_beta`. A lower value of 'mixing_beta' results in less influence of $\rho_{update}$ on $\rho_{new}$, making the self-consistent field (SCF) calculation more stable. However, it may require more steps to achieve convergence.
We recommend the following options:
  - **0.8**: `nspin=1`
  - **0.4**: `nspin=2` and `nspin=4`
  - **0**: keep charge density unchanged, usually used for restarting with `init_chg=file` or testing.
  - **0.1 or less**: if convergence of SCF calculation is difficult to reach, please try `0 < mixing_beta < 0.1`.

  Note: For low-dimensional large systems, the setup of `mixing_beta=0.1`, `mixing_ndim=20`, and `mixing_gg0=1.0` usually works well.

- **Default**: 0.8 for `nspin=1`, 0.4 for `nspin=2` and `nspin=4`.

### mixing_beta_mag

- **Type**: Real
- **Description**: Mixing parameter of magnetic density.
- **Default**: `4*mixing_beta`, but the maximum value is 1.6.

Note that `mixing_beta_mag` is not euqal to `mixing_beta` means that $\rho_{up}$ and $\rho_{down}$ mix independently from each other. This setting will fail for one case where the $\rho_{up}$ and $\rho_{down}$ of the ground state refers to different Kohn-Sham orbitals. For an atomic system, the $\rho_{up}$ and $\rho_{down}$ of the ground state refers to different Kohn-Sham orbitals. We all know Kohn-Sham orbitals are orthogonal to each other. So the mixture of $\rho_{up}$ and $\rho_{down}$ should be exactly independent, otherwise SCF cannot find the ground state forever. To sum up, please make sure `mixing_beta_mag` and `mixing_gg0_mag` exactly euqal to `mixing_beta` and `mixing_gg0` if you calculate an atomic system.

### mixing_ndim

- **Type**: Integer
- **Description**: It indicates the mixing dimensions in Pulay or Broyden. Pulay and Broyden method use the density from previous mixing_ndim steps and do a charge mixing based on this density.

  For systems that are difficult to converge, one could try increasing the value of 'mixing_ndim' to enhance the stability of the self-consistent field (SCF) calculation.
- **Default**: 8

### mixing_restart

- **Type**: double
- **Description**: If the density difference between input and output `drho` is smaller than `mixing_restart`, SCF will restart at next step which means SCF will restart by using output charge density from perivos iteration as input charge density directly, and start a new mixing. Notice that `mixing_restart` will only take effect once in one SCF.

- **Default**: 0

### mixing_dmr

- **Type**: bool
- **Availability**: Only for `mixing_restart>=0.0`
- **Description**: At n-th iteration which is calculated by `drho<mixing_restart`, SCF will start a mixing for real-space density matrix by using the same coefficiences as the mixing of charge density.

- **Default**: false

### mixing_gg0

- **Type**: Real
- **Description**: Whether to perfom Kerker scaling for charge density.
  - **>0**: The high frequency wave vectors will be suppressed by multiplying a scaling factor $\frac{k^2}{k^2+gg0^2}$. Setting `mixing_gg0 = 1.0` is normally a good starting point. Kerker preconditioner will be automatically turned off if `mixing_beta <= 0.1`.
  - **0**: No Kerker scaling is performed.

  For systems that are difficult to converge, particularly metallic systems, enabling Kerker scaling may aid in achieving convergence.
- **Default**: 1.0

### mixing_gg0_mag

- **Type**: Real
- **Description**: Whether to perfom Kerker preconditioner of magnetic density.
  Note: we do not recommand to open Kerker preconditioner of magnetic density unless the system is too hard to converge.
- **Default**: 0.0

### mixing_gg0_min

- **Type**: Real
- **Description**: the minimum kerker coefficient
- **Default**: 0.1

### mixing_angle

- **Type**: Real
- **Availability**: Only relevant for non-colinear calculations `nspin=4`.
- **Description**: Normal broyden mixing can give the converged result for a given magnetic configuration. If one is not interested in the energies of a given magnetic configuration but wants to determine the ground state by relaxing the magnetic moments’ directions, one cannot rely on the standard Broyden mixing algorithm. To enhance the ability to find correct magnetic configuration for non-colinear calculations, ABACUS implements a promising mixing method proposed by J. Phys. Soc. Jpn. 82 (2013) 114706. Here, `mixing_angle` is the angle mixing parameter. In fact, only `mixing_angle=1.0` is implemented currently.
  - **<=0**: Normal broyden mixing for $m_{x}, m_{y}, m_{z}$
  - **>0**: Angle mixing for the modulus $|m|$ with `mixing_angle=1.0`
- **Default**: -10.0

Note: In new angle mixing, you should set `mixing_beta_mag >> mixing_beta`. The setup of `mixing_beta=0.2`, `mixing_beta_mag=1.0` usually works well.

### mixing_tau

- **Type**: Boolean
- **Availability**: Only relevant for meta-GGA calculations.
- **Description**: Whether to mix the kinetic energy density.
  - **True**: The kinetic energy density will also be mixed. It seems for general cases, SCF converges fine even without this mixing. However, if there is difficulty in converging SCF for meta-GGA, it might be helpful to turn this on.
  - **False**: The kinetic energy density will not be mixed.
- **Default**: False

### mixing_dftu

- **Type**: Boolean
- **Availability**: Only relevant for DFT+U calculations.
- **Description**: Whether to mix the occupation matrices.
  - **True**: The occupation matrices will also be mixed by plain mixing. From experience this is not very helpful if the +U calculation does not converge.
  - **False**: The occupation matrices will not be mixed.
- **Default**: False

### gamma_only

- **Type**: Integer
- **Availability**: Only used in localized orbitals set
- **Description**: Whether to use gamma_only algorithm.
  - **0**: more than one k-point is used and the ABACUS is slower compared to the gamma only algorithm.
  - **1**: ABACUS uses gamma only, the algorithm is faster and you don't need to specify the k-points file.

  Note: If gamma_only is set to 1, the KPT file will be overwritten. So make sure to turn off gamma_only for multi-k calculations.

- **Default**: 0

### printe

- **Type**: Integer
- **Description**: Print out energy for each band for every printe step
- **Default**: `scf_nmax`

### scf_nmax

- **Type**: Integer
- **Description**: This variable indicates the maximal iteration number for electronic iterations.
- **Default**: 100

### scf_thr

- **Type**: Real
- **Description**: It's the density threshold for electronic iteration. It represents the charge density error between two sequential densities from electronic iterations. Usually for local orbitals, usually 1e-6 may be accurate enough.
- **Default**: 1.0e-9 (plane-wave basis), or 1.0e-7 (localized atomic orbital basis).

### scf_ene_thr

- **Type**: Real
- **Description**: It's the energy threshold for electronic iteration. It represents the total energy error between two sequential densities from electronic iterations.
- **Default**: -1.0. If the user does not set this parameter, it will not take effect.
- **Unit**: eV

### scf_thr_type

- **Type**: Integer
- **Description**: Choose the calculation method of convergence criterion.
  - **1**: the criterion is defined as $\Delta\rho_G = \frac{1}{2}\iint{\frac{\Delta\rho(r)\Delta\rho(r')}{|r-r'|}d^3r d^3r'}$.
  - **2**: the criterion is defined as $\Delta\rho_R = \frac{1}{N_e}\int{|\Delta\rho(r)|d^3r}$, where $N_e$ is the number of electron.

  Note: This parameter is still under testing and the default setting is usually sufficient.

- **Default**: 1 (plane-wave basis), or 2 (localized atomic orbital basis).

### scf_os_stop

- **Type**: bool
- **Description**: For systems that are difficult to converge, the SCF process may exhibit oscillations in charge density, preventing further progress toward the specified convergence criteria and resulting in continuous oscillation until the maximum number of steps is reached; this greatly wastes computational resources. To address this issue, this function allows ABACUS to terminate the SCF process early upon detecting oscillations, thus reducing subsequent meaningless calculations. The detection of oscillations is based on the slope of the logarithm of historical drho values.. To this end, Least Squares Method is used to calculate the slope of the logarithmically taken drho for the previous `scf_os_ndim` iterations. If the calculated slope is larger than `scf_os_thr`, stop the SCF.

  - **0**: The SCF will continue to run regardless of whether there is oscillation or not. 
  - **1**: If the calculated slope is larger than `scf_os_thr`, stop the SCF.

- **Default**: false

### scf_os_thr

- **Type**: double
- **Description**: The slope threshold to determine if the SCF is stuck in a charge density oscillation. If the calculated slope is larger than `scf_os_thr`, stop the SCF.

- **Default**: -0.01

### scf_os_ndim

- **Type**: int
- **Description**: To determine the number of old iterations' `drho` used in slope calculations.
- **Default**: `mixing_ndim`

### sc_os_ndim

- **Type**: int
- **Description**: To determine the number of old iterations to judge oscillation, it occured,  more accurate lambda with DeltaSpin method would be calculated, only for PW base.
- **Default**: 5

### chg_extrap

- **Type**: String
- **Description**: Methods to do extrapolation of density when ABACUS is doing geometry relaxations or molecular dynamics.
  - **atomic**: atomic extrapolation.
  - **first-order**: first-order extrapolation.
  - **second-order**: second-order extrapolation.
- **Default**: first-order (geometry relaxations), second-order (molecular dynamics), else atomic

### lspinorb

- **Type**: Boolean
- **Description**: Whether to consider spin-orbital coupling effect in the calculation.
  - **True**: Consider spin-orbital coupling effect, and `nspin` is also automatically set to 4.
  - **False**: Do not consider spin-orbital coupling effect.
- **Default**: False

### noncolin

- **Type**: Boolean
- **Description**: Whether to allow non-collinear polarization, in which case the coupling between spin up and spin down will be taken into account.
  - **True**: Allow non-collinear polarization, and `nspin` is also automatically set to 4.
  - **False**: Do not allow non-collinear polarization.
- **Default**: False

### soc_lambda

- **Type**: Real
- **Availability**: Relevant for soc calculations.
- **Description**: Sometimes, for some real materials, both scalar-relativistic and full-relativistic can not describe the exact spin-orbit coupling. Artificial modulation may help in such cases.

  `soc_lambda`, which has value range [0.0, 1.0] , is used for modulate SOC effect.

  In particular, `soc_lambda 0.0` refers to scalar-relativistic case and `soc_lambda 1.0` refers to full-relativistic case.
- **Default**: 1.0

[back to top](#full-list-of-input-keywords)

## Electronic structure (SDFT)

These variables are used to control the parameters of stochastic DFT (SDFT),  mix stochastic-deterministic DFT (MDFT), or complete-basis Chebyshev method (CT). In the following text, stochastic DFT is used to refer to these three methods. We suggest using SDFT to calculate high-temperature systems and we only support [smearing_method](#smearing_method) "fd". Both "scf" and "nscf" [calculation](#calculation) are supported.

### method_sto

- **Type**: Integer
- **Availability**: [esolver_type](#esolver_type) = `sdft`
- **Description**: Different methods to do stochastic DFT
  - 1: Calculate $T_n(\hat{h})\ket{\chi}$ twice, where $T_n(x)$ is the n-th order Chebyshev polynomial and $\hat{h}=\frac{\hat{H}-\bar{E}}{\Delta E}$ owning eigenvalues $\in(-1,1)$. This method cost less memory but is slower.
  - 2: Calculate $T_n(\hat{h})\ket{\chi}$ once but needs much more memory. This method is much faster. Besides, it calculates $N_e$ with $\bra{\chi}\sqrt{\hat f}\sqrt{\hat f}\ket{\chi}$, which needs a smaller [nche_sto](#nche_sto). However, when the memory is not enough, only method 1 can be used.
  - other: use 2
- **Default**: 2

### nbands_sto

- **Type**: Integer or string
- **Availability**: [esolver_type](#esolver_type) = `sdft`
- **Description**: The number of stochastic orbitals
  - \> 0: Perform stochastic DFT.
   Increasing the number of bands improves accuracy and reduces stochastic errors, which scale as $1/\sqrt{N_{\chi}}$;
   To perform mixed stochastic-deterministic DFT, you should set [nbands](#nbands), which represents the number of KS orbitals.
  - 0: Perform Kohn-Sham DFT.
  - all: All complete basis sets are used to replace stochastic orbitals with the Chebyshev method (CT), resulting in the same results as KSDFT without stochastic errors.
- **Default**: 256

### nche_sto

- **Type**: Integer
- **Availability**: [esolver_type](#esolver_type) = `sdft`
- **Description**: Chebyshev expansion orders for stochastic DFT.
- **Default**: 100

### emin_sto

- **Type**: Real
- **Availability**: [esolver_type](#esolver_type) = `sdft`
- **Description**: Trial energy to guess the lower bound of eigen energies of the Hamiltonian Operator $\hat{H}$.
- **Default**: 0.0
- **Unit**: Ry

### emax_sto

- **Type**: Real
- **Availability**: [esolver_type](#esolver_type) = `sdft`
- **Description**: Trial energy to guess the upper bound of eigen energies of the Hamiltonian Operator $\hat{H}$.
- **Default**: 0.0
- **Unit**: Ry

### seed_sto

- **Type**: Integer
- **Availability**: [esolver_type](#esolver_type) = `sdft`
- **Description**: The random seed to generate stochastic orbitals.
  - \>= 0: Stochastic orbitals have the form of $\exp(i2\pi\theta(G))$, where $\theta$ is a uniform distribution in $(0,1)$.
  - 0: the seed is decided by time(NULL).
  - \<= -1: Stochastic orbitals have the form of $\pm1$ with equal probability.
  - -1: the seed is decided by time(NULL).
- **Default**: 0

### initsto_ecut

- **Type**: Real
- **Availability**: [esolver_type](#esolver_type) = `sdft`
- **Description**: Stochastic wave functions are initialized in a large box generated by "4*`initsto_ecut`". `initsto_ecut` should be larger than [ecutwfc](#ecutwfc). In this method, SDFT results are the same when using different cores. Besides, coefficients of the same G are the same when ecutwfc is rising to initsto_ecut. If it is smaller than [ecutwfc](#ecutwfc), it will be turned off.
- **Default**: 0.0
- **Unit**: Ry

### initsto_freq

- **Type**: Integer
- **Availability**: [esolver_type](#esolver_type) = `sdft`
- **Description**: Frequency (once each initsto_freq steps) to generate new stochastic orbitals when running md.
  - positive integer: Update stochastic orbitals
  - 0:                Never change stochastic orbitals.
- **Default**: 0

### npart_sto

- **Type**: Integer
- **Availability**: [method_sto](#method_sto) = `2` and [out_dos](#out_dos) = 1 or [cal_cond](#cal_cond) = `True`
- **Description**: Make memory cost to 1/npart_sto times of the previous one when running the post process of SDFT like DOS or conductivities.
- **Default**: 1

[back to top](#full-list-of-input-keywords)

## Geometry relaxation

These variables are used to control the geometry relaxation.

### relax_method

- **Type**: String
- **Description**: The methods to do geometry optimization.
  - cg: using the conjugate gradient (CG) algorithm. Note that there are two implementations of the conjugate gradient (CG) method, see [relax_new](#relax_new).
  - bfgs: using the Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm.
  - bfgs_trad: using the traditional Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm. 
  - cg_bfgs: using the CG method for the initial steps, and switching to BFGS method when the force convergence is smaller than [relax_cg_thr](#relax_cg_thr).
  - sd: using the steepest descent (SD) algorithm.
  - fire: the Fast Inertial Relaxation Engine method (FIRE), a kind of molecular-dynamics-based relaxation algorithm, is implemented in the molecular dynamics (MD) module. The algorithm can be used by setting [calculation](#calculation) to `md` and [md_type](#md_type) to `fire`. Also ionic velocities should be set in this case. See [fire](../md.md#fire) for more details.
- **Default**: cg

### relax_new

- **Type**: Boolean
- **Description**: At around the end of 2022 we made a new implementation of the Conjugate Gradient (CG) method for `relax` and `cell-relax` calculations. But the old implementation was also kept.
  - True: use the new implementation of CG method for `relax` and `cell-relax` calculations.
  - False: use the old implementation of CG method for `relax` and `cell-relax` calculations.
- **Default**: True

### relax_scale_force

- **Type**: Real
- **Availability**: only used when `relax_new` set to `True`
- **Description**: The paramether controls the size of the first conjugate gradient step. A smaller value means the first step along a new CG direction is smaller. This might be helpful for large systems, where it is safer to take a smaller initial step to prevent the collapse of the whole configuration.
- **Default**: 0.5

### relax_nmax

- **Type**: Integer
- **Description**: The maximal number of ionic iteration steps. If set to 0, the code performs a quick "dry run", stopping just after initialization. This is useful to check for input correctness and to have the summary printed.
- **Default**: 1 for SCF, 50 for relax and cell-relax calcualtions

### relax_cg_thr

- **Type**: Real
- **Description**: When move-method is set to `cg_bfgs`, a mixed algorithm of conjugate gradient (CG) method and Broyden–Fletcher–Goldfarb–Shanno (BFGS) method is used. The ions first move according to CG method, then switched to BFGS method when the maximum of force on atoms is reduced below the CG force threshold, which is set by this parameter.
- **Default**: 0.5
- **Unit**: eV/Angstrom

### cal_force

- **Type**: Boolean
- **Description**:
  - **True** calculate the force at the end of the electronic iteration
  - **False** no force calculation at the end of the electronic iteration
- **Default**: False if `calculation` is set to `scf`, True if `calculation` is set to `cell-relax`, `relax`, or `md`.

### force_thr

- **Type**: Real
- **Description**: Threshold of the force convergence in Ry/Bohr. The threshold is compared with the largest force among all of the atoms. The recommended value for using atomic orbitals is 0.04 eV/Angstrom (0.0016 Ry/Bohr). The parameter is equivalent to [force_thr_ev](#force_thr_ev) except for the unit. You may choose either you like.
- **Default**: 0.001
- **Unit**: Ry/Bohr (25.7112 eV/Angstrom)

### force_thr_ev

- **Type**: Real
- **Description**: Threshold of the force convergence in eV/Angstrom. The threshold is compared with the largest force among all of the atoms. The recommended value for using atomic orbitals is 0.04 eV/Angstrom (0.0016 Ry/Bohr). The parameter is equivalent to [force_thr](#force_thr) except for the unit. You may choose either you like.
- **Default**: 0.0257112
- **Unit**: eV/Angstrom (0.03889 Ry/Bohr)

### force_thr_ev2

- **Type**: Real
- **Description**: The calculated force will be set to 0 when it is smaller than the parameter `force_thr_ev2`.
- **Default**: 0.0
- **Unit**: eV/Angstrom

### relax_bfgs_w1

- **Type**: Real
- **Description**: This variable controls the Wolfe condition for Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm used in geometry relaxation. You can look into the paper Phys.Chem.Chem.Phys.,2000,2,2177 for more information.
- **Default**: 0.01

### relax_bfgs_w2

- **Type**: Real
- **Description**: This variable controls the Wolfe condition for Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm used in geometry relaxation. You can look into the paper Phys.Chem.Chem.Phys.,2000,2,2177 for more information.
- **Default**: 0.5

### relax_bfgs_rmax

- **Type**: Real
- **Description**: This variable is for geometry optimization. It stands for the maximal movement of all the atoms. The sum of the movements from all atoms can be increased during the optimization steps. However, it can not be larger than `relax_bfgs_rmax`
- **Unit**: Bohr
- **Default**: 0.8

### relax_bfgs_rmin

- **Type**: Real
- **Description**: This variable is for geometry optimization. It indicates the minimal movement of all the atoms. When the movement of all the atoms is smaller than relax_bfgs_rmin Bohr, and the force convergence is still not achieved, the calculation will break down.
- **Default**: 1e-5
- **Unit**: Bohr

### relax_bfgs_init

- **Type**: Real
- **Description**: This variable is for geometry optimization. It stands for the sum of initial movements of all of the atoms.
- **Default**: 0.5
- **Unit**: Bohr

### cal_stress

- **Type**: Boolean
- **Description**:
  - **True**: calculate the stress at the end of the electronic iteration
  - **False**: no calculation of the stress at the end of the electronic iteration
- **Default**: True if `calculation` is `cell-relax`, False otherwise.

### stress_thr

- **Type**: Real
- **Description**: The threshold of the stress convergence. The threshold is compared with the largest component of the stress tensor.
- **Default**: 0.5
- **Unit**: kbar

### press1, press2, press3

- **Type**: Real
- **Description**: The external pressures along three axes. Positive input value is taken as compressive stress.
- **Default**: 0
- **Unit**: kbar

### fixed_axes

- **Type**: String
- **Availability**: only used when `calculation` set to `cell-relax`
- **Description**: Axes that are fixed during cell relaxation. Possible choices are:
  - **None**: default; all of the axes can relax
  - **volume**: relaxation with fixed volume
  - **shape**: fix shape but change volume (i.e. only lattice constant changes)
  - **a**: fix a axis during relaxation
  - **b**: fix b axis during relaxation
  - **c**: fix c axis during relaxation
  - **ab**: fix both a and b axes during relaxation
  - **ac**: fix both a and c axes during relaxation
  - **bc**: fix both b and c axes during relaxation

> Note : fixed_axes = "shape" and "volume" are only available for [relax_new](#relax_new) = True

- **Default**: None

### fixed_ibrav

- **Type**: Boolean
- **Availability**: Must be used along with [relax_new](#relax_new) set to True, and a specific [latname](#latname) must be provided
- **Description**:
  - **True**: the lattice type will be preserved during relaxation
  - **False**: No restrictions are exerted during relaxation in terms of lattice type

> Note: it is possible to use `fixed_ibrav` with `fixed_axes`, but please make sure you know what you are doing. For example, if we are doing relaxation of a simple cubic lattice (`latname` = "sc"), and we use `fixed_ibrav` along with `fixed_axes` = "volume", then the cell is never allowed to move and as a result, the relaxation never converges.

- **Default**: False

### fixed_atoms

- **Type**: Boolean
- **Description**:
  - **True**: The direct coordinates of atoms will be preserved during variable-cell relaxation.
  - **False**: No restrictions are exerted on positions of all atoms. However, users can still fix certain components of certain atoms by using the `m` keyword in `STRU` file. For the latter option, check the end of this [instruction](stru.md).
- **Default**: False

### cell_factor

- **Type**: Real
- **Description**: Used in the construction of the pseudopotential tables. It should exceed the maximum linear contraction of the cell during a simulation.
- **Default**: 1.2

[back to top](#full-list-of-input-keywords)

## Variables related to output information

These variables are used to control the output of properties.

### out_mul

- **Type**: Boolean
- **Availability**: Numerical atomic orbital basis
- **Description**: Whether to print the Mulliken population analysis result into `OUT.${suffix}/mulliken.txt`. In molecular dynamics calculations, the output frequency is controlled by [out_interval](#out_interval).
- **Default**: False

### out_freq_elec

- **Type**: Integer
- **Description**: The output frequency of the charge density (controlled by [out_chg](#out_chg)), wavefunction (controlled by [out_wfc_pw](#out_wfc_pw) or [out_wfc_r](#out_wfc_r)), and density matrix of localized orbitals (controlled by [out_dm](#out_dm)).
  - \>0: Output them every `out_freq_elec` iteration numbers in electronic iterations.
  - 0: Output them when the electronic iteration is converged or reaches the maximal iteration number.
- **Default**: 0

### out_chg

- **Type**: Integer \[Integer\](optional)
- **Description**: 
  The first integer controls whether to output the charge density on real space grids:
  - 1: Output the charge density (in Bohr^-3) on real space grids into the density files in the folder `OUT.${suffix}`. The files are named as:
    - nspin = 1: SPIN1_CHG.cube;
    - nspin = 2: SPIN1_CHG.cube, and SPIN2_CHG.cube;
    - nspin = 4: SPIN1_CHG.cube, SPIN2_CHG.cube, SPIN3_CHG.cube, and SPIN4_CHG.cube.
  - 2: On top of 1, also output the initial charge density. The files are named as:
    - nspin = 1: SPIN1_CHG_INI.cube
    - nspin = 2: SPIN1_CHG_INI.cube, and SPIN2_CHG_INI.cube;
    - nspin = 4: SPIN1_CHG_INI.cube, SPIN2_CHG_INI.cube, SPIN3_CHG_INI.cube, and SPIN4_CHG_INI.cube.
  - -1: disable the charge density auto-back-up file `{suffix}-CHARGE-DENSITY.restart`, useful for large systems.
    
  The second integer controls the precision of the charge density output, if not given, will use `3` as default. For purpose restarting from this file and other high-precision involved calculation, recommend to use `10`.

  ---
  The circle order of the charge density on real space grids is: x is the outer loop, then y and finally z (z is moving fastest).

  If EXX(exact exchange) is calculated, (i.e. *[dft_fuctional](#dft_functional)==hse/hf/pbe0/scan0/opt_orb* or *[rpa](#rpa)==True*), the Hexx(R) files will be output in the folder `OUT.${suffix}` too, which can be read in NSCF calculation.

  In molecular dynamics calculations, the output frequency is controlled by [out_interval](#out_interval).
- **Default**: 0 3

### out_pot

- **Type**: Integer
- **Description**:
  - 1: Output the **total local potential** (i.e., local pseudopotential + Hartree potential + XC potential + external electric field (if exists) + dipole correction potential (if exists) + ...) on real space grids (in Ry) into files in the folder `OUT.${suffix}`. The files are named as:
    - nspin = 1: SPIN1_POT.cube;
    - nspin = 2: SPIN1_POT.cube, and SPIN2_POT.cube;
    - nspin = 4: SPIN1_POT.cube, SPIN2_POT.cube, SPIN3_POT.cube, and SPIN4_POT.cube.
  - 2: Output the **electrostatic potential** on real space grids into `OUT.${suffix}/ElecStaticPot.cube`. The Python script named `tools/average_pot/aveElecStatPot.py` can be used to calculate the average electrostatic potential along the z-axis and outputs it into ElecStaticPot_AVE.

    Please note that the total local potential refers to the local component of the self-consistent potential, excluding the non-local pseudopotential. The distinction between the local potential and the electrostatic potential is as follows: local potential = electrostatic potential + XC potential.
  - 3: Apart from 1, also output the **total local potential** of the initial charge density. The files are named as:
    - nspin = 1: SPIN1_POT_INI.cube;
    - nspin = 2: SPIN1_POT_INI.cube, and SPIN2_POT_INI.cube;
    - nspin = 4: SPIN1_POT_INI.cube, SPIN2_POT_INI.cube, SPIN3_POT_INI.cube, and SPIN4_POT_INI.cube.

  In molecular dynamics calculations, the output frequency is controlled by [out_interval](#out_interval).
- **Default**: 0

### out_dm

- **Type**: Boolean
- **Availability**: Numerical atomic orbital basis (gamma-only algorithm)
- **Description**: Whether to output the density matrix of localized orbitals into files in the folder `OUT.${suffix}`. The files are named as:
  - nspin = 1: SPIN1_DM;
  - nspin = 2: SPIN1_DM, and SPIN2_DM.
- **Default**: False

### out_dm1

- **Type**: Boolean
- **Availability**: Numerical atomic orbital basis (multi-k points)
- **Description**: Whether to output the density matrix of localized orbitals into files in the folder `OUT.${suffix}`. The density matrices are written in the format of sparse matrices, as mentioned in [out_mat_hs2](#out_mat_hs2). The files are named as:
  - nspin = 1: data-DMR-sparse_SPIN0.csr;
  - nspin = 2: data-DMR-sparse_SPIN0.csr, and data-DMR-sparse_SPIN1.csr.
- **Default**: False

### out_wfc_pw

- **Type**: Integer
- **Availability**: Plane wave basis or get_wf calculation in numerical atomic orbital basis
- **Description**:
  - 1: Output the coefficients of wave functions into text files named `OUT.${suffix}/WAVEFUNC${K}.txt`, where ${K} is the index of k points.
  - 2: results are stored in binary files named `OUT.${suffix}/WAVEFUNC${K}.dat`.
- **Default**: 0

### out_wfc_r

- **Type**: Boolean
- **Availability**: Plane wave basis or get_wf calculation in numerical atomic orbital basis
- **Description**: Whether to output real-space wave functions into `OUT.suffix/wfc_realspace/wfc_realspace_${K}_${B}`, where `${K}` is the index of k points, `${B}` is the index of bands.
- **Default**: False

### out_wfc_lcao

- **Type**: Integer
- **Availability**: Numerical atomic orbital basis
- **Description**: Whether to output the wavefunction coefficients into files in the folder `OUT.${suffix}`. The files are named as `WFC_{GAMMA|K}{index of K point}`, and if [out_app_flag](#out_app_flag) is false, the file name will also contains `_ION{ION step}`, where `ION step` is the index of ionic step:
  - 0: no output
  - 1: (txt format) 
    - gamma-only: `WFC_NAO_GAMMA1_ION1.txt` or `WFC_NAO_GAMMA1.txt`, ...;
    - non-gamma-only: `WFC_NAO_K1_ION1.txt` or `WFC_NAO_K1.txt`, ...;
  - 2: (binary format)
    - gamma-only: `WFC_NAO_GAMMA1_ION1.dat` or `WFC_NAO_GAMMA1.dat`, ...;
    - non-gamma-only: `WFC_NAO_K1_ION1.dat` or `WFC_NAO_K1.dat`, ....

  The corresponding sequence of the orbitals can be seen in [Basis Set](../pp_orb.md#basis-set).

  Also controled by [out_interval](#out_interval) and [out_app_flag](#out_app_flag).
- **Default**: Flase

### out_dos

- **Type**: Integer
- **Description**: Whether to output the density of states (DOS). For more information, refer to the [dos.md](../elec_properties/dos.md).
  - 0: no output
  - 1: output the density of states (DOS)
  - 2: 
    - lcao-only: output the density of states (DOS) and the projected density of states (PDOS)
- **Default**: 0

### out_band

- **Type**: Boolean \[Integer\](optional)
- **Description**: Whether to output the band structure (in eV), optionally output precision can be set by a second parameter, default is 8. For more information, refer to the [band.md](../elec_properties/band.md)
- **Default**: False

### out_proj_band

- **Type**: Boolean
- **Description**: Whether to output the projected band structure. For more information, refer to the [band.md](../elec_properties/band.md)
- **Default**: False

### out_stru

- **Type**: Boolean
- **Description**: Whether to output structure files per ionic step in geometry relaxation calculations into `OUT.${suffix}/STRU_ION${istep}_D`, where `${istep}` is the ionic step.
- **Default**: False

### out_bandgap

- **Type**: Boolean
- **Description**: Whether to print the bandgap per electronic iteration into `OUT.${suffix}/running_${calculation}.log`. The value of bandgaps can be obtained by searching for the keyword:
  - [nupdown](#nupdown) > 0: `E_bandgap_up` and `E_bandgap_dw`
  - [nupdown](#nupdown) = 0: `E_bandgap`
- **Default**: False

### out_level

- **Type**: String
- **Description**: Control the output level of information in `OUT.${suffix}/running_${calculation}.log`.
  - ie: electronic iteration level, which prints useful information for electronic iterations;
  - i: geometry relaxation level, which prints some information for geometry relaxations additionally;
  - m: molecular dynamics level, which does not print some information for simplicity.

- **Default**: ie

### out_alllog

- **Type**: Boolean
- **Description**: Whether to print information into individual logs from all ranks in an MPI run.
  - True: Information from each rank will be written into individual files named `OUT.${suffix}/running_${calculation}_${rank+1}.log`.
  - False: Information will only be written from rank 0 into a file named `OUT.${suffix}/running_${calculation}.log`.
- **Default**: False

### out_mat_hs

- **Type**: Boolean \[Integer\](optional)
- **Availability**: Numerical atomic orbital basis
- **Description**: Whether to print the upper triangular part of the Hamiltonian matrices (in Ry) and overlap matrices for each k point into files in the directory `OUT.${suffix}`. The second number controls precision. For more information, please refer to [hs_matrix.md](../elec_properties/hs_matrix.md#out_mat_hs). Also controled by [out_interval](#out_interval) and [out_app_flag](#out_app_flag).
- **Default**: False 8

### out_mat_tk

- **Type**: Boolean \[Integer\](optional)
- **Availability**: Numerical atomic orbital basis
- **Description**: Whether to print the upper triangular part of the kinetic matrices (in Ry) for each k point into `OUT.${suffix}/data-i-T`, where i is the index of k points (see `OUT.${suffix}/kpoints`). One may optionally provide a second parameter to specify the precision. 
- **Default**: False \[8\]

### out_mat_r

- **Type**: Boolean
- **Availability**: Numerical atomic orbital basis (not gamma-only algorithm)
- **Description**: Whether to print the matrix representation of the position matrix (in Bohr) into a file named `data-rR-tr` in the directory `OUT.${suffix}`. If [calculation](#calculation) is set to `get_S`, the position matrix can be obtained without scf iterations. For more information, please refer to [position_matrix.md](../elec_properties/position_matrix.md#extracting-position-matrices).
- **Default**: False

### out_mat_hs2

- **Type**: Boolean
- **Availability**: Numerical atomic orbital basis (not gamma-only algorithm)
- **Description**: Whether to print files containing the Hamiltonian matrix $H(R)$ (in Ry) and overlap matrix $S(R)$ into files in the directory `OUT.${suffix}`. For more information, please refer to [hs_matrix.md](../elec_properties/hs_matrix.md#out_mat_hs2).
- **Default**: False

### out_mat_t

- **Type**: Boolean
- **Availability**: Numerical atomic orbital basis (not gamma-only algorithm)
- **Description**: For LCAO calculations, if out_mat_t is set to 1, ABACUS will generate files containing the kinetic energy matrix $T(R)$ (in Ry). The format will be the same as the Hamiltonian matrix $H(R)$ and overlap matrix $S(R)$ as mentioned in [out_mat_hs2](#out_mat_hs2). The name of the files will be `data-TR-sparse_SPIN0.csr` and so on. Also controled by [out_interval](#out_interval) and [out_app_flag](#out_app_flag).
- **Default**: False

### out_mat_dh

- **Type**: Boolean
- **Availability**: Numerical atomic orbital basis (not gamma-only algorithm)
- **Description**: Whether to print files containing the derivatives of the Hamiltonian matrix (in Ry/Bohr). The format will be the same as the Hamiltonian matrix $H(R)$ and overlap matrix $S(R)$ as mentioned in [out_mat_hs2](#out_mat_hs2). The name of the files will be `data-dHRx-sparse_SPIN0.csr` and so on. Also controled by [out_interval](#out_interval) and [out_app_flag](#out_app_flag).
- **Default**: False

### out_mat_xc

- **Type**: Boolean
- **Availability**: Numerical atomic orbital (NAO) and NAO-in-PW basis
- **Description**: Whether to print the upper triangular part of the exchange-correlation matrices in **Kohn-Sham orbital representation** (unit: Ry): $\braket{\psi_i|V_\text{xc}^\text{(semi-)local}+V_\text{exx}+V_\text{DFTU}|\psi_j}$ for each k point into files in the directory `OUT.${suffix}`, which is useful for the subsequent GW calculation. (Note that currently DeePKS term is not included. ) The files are named `k-$k-Vxc`, the meaning of `$k`corresponding to k point and spin  is same as [hs_matrix.md](../elec_properties/hs_matrix.md#out_mat_hs).
The band (KS orbital) energy for each (k-point, spin, band) will be printed in the file `OUT.${suffix}/vxc_out.dat`. If EXX is calculated, the local and EXX part of band energy will also be printed in `OUT.${suffix}/vxc_local_out.dat`and `OUT.${suffix}/vxc_exx_out.dat`, respectively. All the `vxc*_out.dat` files contains 3 integers (nk, nspin, nband) followed by nk\*nspin\*nband lines of energy Hartree and eV.
- **Default**: False

### out_eband_terms

- **Type**: Boolean
- **Availability**: Numerical atomic orbital basis
- **Description**: Whether to print the band energy terms separately in the file `OUT.${suffix}/${term}_out.dat`. The terms include the kinetic, pseudopotential (local + nonlocal), Hartree and exchange-correlation (including exact exchange if calculated).
- **Default**: False

### out_hr_npz/out_dm_npz

- **Type**: Boolean
- **Availability**: Numerical atomic orbital basis
- **Description**: Whether to print Hamiltonian matrices $H(R)$/density matrics $DM(R)$ in npz format. This feature does not work for gamma-only calculations. Currently only intended for internal usage.
- **Default**: False

### dm_to_rho

- **Type**: Boolean
- **Availability**: Numerical atomic orbital basis
- **Description**: Reads density matrix $DM(R)$ in npz format and creates electron density on grids. This feature does not work for gamma-only calculations. Only supports serial calculations. Currently only intended for internal usage.
- **Default**: False

### out_app_flag

- **Type**: Boolean
- **Availability**: Numerical atomic orbital basis (not gamma-only algorithm)
- **Description**: Whether to output $r(R)$, $H(R)$, $S(R)$, $T(R)$, $dH(R)$, $H(k)$, $S(k)$ and $wfc(k)$ matrices in an append manner during molecular dynamics calculations. Check input parameters [out_mat_r](#out_mat_r), [out_mat_hs2](#out_mat_hs2), [out_mat_t](#out_mat_t), [out_mat_dh](#out_mat_dh), [out_mat_hs](#out_mat_hs) and [out_wfc_lcao](#out_wfc_lcao) for more information.
- **Default**: true

### out_ndigits

- **Type**: Integar
- **Availability**: `out_mat_hs 1` case presently.
- **Description**: Controls the length of decimal part of output data, such as charge density, Hamiltonian matrix, Overlap matrix and so on.
- **Default**: 8

### out_interval

- **Type**: Integer
- **Description**: Control the interval for printing charge density, local potential, electrostatic potential, Mulliken population analysis, $r(R)$, $H(R)$, $S(R)$, $T(R)$, $dH(R)$, $H(k)$, $S(k)$ and $\psi(k)$ matrices during molecular dynamics calculations. Check input parameters [out_chg](#out_chg), [out_pot](#out_pot), [out_mul](#out_mul), [out_mat_r](#out_mat_r), [out_mat_hs2](#out_mat_hs2), [out_mat_t](#out_mat_t), [out_mat_dh](#out_mat_dh), [out_mat_hs](#out_mat_hs) and [out_wfc_lcao](#out_wfc_lcao) for more information, respectively.
- **Default**: 1

### out_element_info

- **Type**: Boolean
- **Description**: Whether to print element information into files in the directory `OUT.${suffix}/${element_label}`, including pseudopotential and orbital information of the element (in atomic Ryberg units).
- **Default**: False

### restart_save

- **Type**: Boolean
- **Availability**: Numerical atomic orbital basis
- **Description**: Whether to save charge density files per ionic step, which are used to restart calculations. According to the value of [read_file_dir](#read_file_dir):
  - auto: These files are saved in folder `OUT.${suffix}/restart/`;
  - other: These files are saved in folder `${read_file_dir}/restart/`.

  If EXX(exact exchange) is calculated (i.e. *[dft_fuctional](#dft_functional)==hse/hf/pbe0/scan0/opt_orb* or *[rpa](#rpa)==True*), the Hexx(R) files for each processor will also be saved in the above folder, which can be read in EXX calculation with *[restart_load](#restart_load)==True*.
- **Default**: False

### rpa

- **Type**: Boolean
- **Description**: Generate output files used in rpa calculations.
- **Default**: False

### nbands_istate

- **Type**: Integer
- **Availability**: Only for LCAO, used when `calculation = get_wf` or `calculation = get_pchg`.
- **Description**: The number of bands around the Fermi level you would like to calculate. `get_wf` means to calculate the envelope functions of wave functions $\Psi_{i}=\Sigma_{\mu}C_{i\mu}\Phi_{\mu}$, where $\Psi_{i}$ is the ith wave function with the band index $i$ and $\Phi_{\mu}$ is the localized atomic orbital set. `get_pchg` means to calculate the density of each wave function $|\Psi_{i}|^{2}$. Specifically, suppose we have highest occupied bands at 100th wave functions. And if you set this variable to 5, it will print five wave functions from 96th to 105th. But before all this can be carried out, the wave functions coefficients should be first calculated and written into a file by setting the flag `out_wfc_lcao = 1`.
- **Default**: 5

### bands_to_print

- **Type**: String
- **Availability**: For both PW and LCAO. When `basis_type = lcao`, used when `calculation = get_wf` or `calculation = get_pchg`.
- **Description**: Specifies the bands to calculate the wave functions/charge densities for, using a space-separated string of 0s and 1s, providing a more flexible selection compared to `nbands_istate`. Each digit in the string corresponds to a band, starting from the first band. A `1` indicates that the charge density should be calculated for that band, while a `0` means the band will be ignored. The parameter allows a compact and flexible notation (similar to [`ocp_set`](#ocp_set)), for example the syntax `1 4*0 5*1 0` is used to denote the selection of bands: `1` means calculate for the first band, `4*0` skips the next four bands, `5*1` means calculate for the following five bands, and the final `0` skips the next band. It's essential that the total count of bands does not exceed the total number of bands (`nbands`); otherwise, it results in an error, and the process exits. The input string must contain only numbers and the asterisk (`*`) for repetition, ensuring correct format and intention of band selection.
- **Default**: none

### if_separate_k

- **Type**: Boolean
- **Availability**: Only for LCAO, used only when `calculation = get_pchg` and `gamma_only` is turned off.
- **Description**: Specifies whether to write the partial charge densities for all k-points to individual files or merge them. **Warning**: Enabling symmetry may produce incorrect results due to incorrect k-point weights. Therefore, when calculating partial charge densities, it is strongly recommended to set `symmetry = -1`.
- **Default**: false

### out_elf

- **Type**: Integer \[Integer\](optional)
- **Availability**: Only for Kohn-Sham DFT and Orbital Free DFT.
- **Description**: Whether to output the electron localization function (ELF) in the folder `OUT.${suffix}`. The files are named as 
    - nspin = 1:
      - ELF.cube: ${\rm{ELF}} = \frac{1}{1+\chi^2}$, $\chi = \frac{\frac{1}{2}\sum_{i}{f_i |\nabla\psi_{i}|^2} - \frac{|\nabla\rho|^2}{8\rho}}{\frac{3}{10}(3\pi^2)^{2/3}\rho^{5/3}}$;
    - nspin = 2:
      - ELF_SPIN1.cube, ELF_SPIN2.cube: ${\rm{ELF}}_\sigma = \frac{1}{1+\chi_\sigma^2}$, $\chi_\sigma = \frac{\frac{1}{2}\sum_{i}{f_i |\nabla\psi_{i,\sigma}|^2} - \frac{|\nabla\rho_\sigma|^2}{8\rho_\sigma}}{\frac{3}{10}(6\pi^2)^{2/3}\rho_\sigma^{5/3}}$;
      - ELF.cube: ${\rm{ELF}} = \frac{1}{1+\chi^2}$, $\chi = \frac{\frac{1}{2}\sum_{i,\sigma}{f_i |\nabla\psi_{i,\sigma}|^2} - \sum_{\sigma}{\frac{|\nabla\rho_\sigma|^2}{8\rho_\sigma}}}{\sum_{\sigma}{\frac{3}{10}(6\pi^2)^{2/3}\rho_\sigma^{5/3}}}$;

  The second integer controls the precision of the kinetic energy density output, if not given, will use `3` as default. For purpose restarting from this file and other high-precision involved calculation, recommend to use `10`.

  ---
  In molecular dynamics calculations, the output frequency is controlled by [out_interval](#out_interval).
- **Default**: 0 3

[back to top](#full-list-of-input-keywords)

## Density of states

These variables are used to control the calculation of DOS. [Detailed introduction](https://github.com/deepmodeling/abacus-develop/blob/develop/docs/advanced/elec_properties/dos.md)

### dos_edelta_ev

- **Type**: Real
- **Description**: the step size in writing Density of States (DOS)
- **Default**: 0.01
- **Unit**: eV

### dos_sigma

- **Type**: Real
- **Description**: the width of the Gaussian factor when obtaining smeared Density of States (DOS)
- **Default**: 0.07
- **Unit**: eV

### dos_scale

- **Type**: Real
- **Description**: Defines the energy range of DOS output as (emax-emin)*(1+dos_scale), centered at (emax+emin)/2. This parameter will be used when dos_emin and dos_emax are not set.
- **Default**: 0.01
- **Unit**: eV

### dos_emin_ev

- **Type**: Real
- **Description**: the minimal range for Density of States (DOS)
  - If set, "dos_scale" will be ignored.
- **Default**: Minimal eigenenergy of $\hat{H}$
- **Unit**: eV

### dos_emax_ev

- **Type**: Real
- **Description**: the maximal range for Density of States (DOS)
  - If set, "dos_scale" will be ignored.
- **Default**: Maximal eigenenergy of $\hat{H}$
- **Unit**: eV

### dos_nche

- **Type**: Integer
The order of Chebyshev expansions when using Stochastic Density Functional Theory (SDFT) to calculate DOS.
- **Default**: 100

[back to top](#full-list-of-input-keywords)

## NAOs

These variables are used to control the generation of numerical atomic orbitals (NAOs), whose radial parts are linear combinations of spherical Bessel functions with a node (i.e., evaluate to zero) at the cutoff radius.
In plane-wave-based calculations, necessary information will be printed into `OUT.${suffix}/orb_matrix.${i}.dat`, which serves as an input file for the generation of NAOs. Please check [SIAB package](../../../tools/SIAB/README.md#siab-package-description) for more information.

### bessel_nao_ecut

- **Type**: Real
- **Description**: "Energy cutoff" (in Ry) of spherical Bessel functions. The number of spherical Bessel functions that constitute the radial parts of NAOs is determined by sqrt(`bessel_nao_ecut`)$\times$`bessel_nao_rcut`/$\pi$.
- **Default**: `ecutwfc`

### bessel_nao_tolerence

- **Type**: Real
- **Description**: tolerance when searching for the zeros of spherical Bessel functions.
- **Default**: 1.0e-12

### bessel_nao_rcut

- **Type**: Real
- **Description**: Cutoff radius (in Bohr) and the common node of spherical Bessel functions used to construct the NAOs.
- **Default**: 6.0

### bessel_nao_smooth

- **Type**: Boolean
- **Description**: if True, NAOs will be smoothed near the cutoff radius by $1-\exp\left(-\frac{(r-r_{cut})^2}{2\sigma^2}\right)$. See `bessel_nao_rcut` for $r_{cut}$ and `bessel_nao_sigma` for $\sigma$.
- **Default**: True

### bessel_nao_sigma

- **Type**: Real
- **Description**: Smoothing range (in Bohr). See also `bessel_nao_smooth`.
- **Default**: 0.1

[back to top](#full-list-of-input-keywords)

## DeePKS

These variables are used to control the usage of DeePKS method (a comprehensive data-driven approach to improve the accuracy of DFT).
Warning: this function is not robust enough for the current version. Please try the following variables at your own risk:

### deepks_out_labels

- **Type**: Boolean
- **Availability**: numerical atomic orbital basis
- **Description**: Print labels and descriptors for DeePKS training in OUT.${suffix}. The names of these files start with "deepks".
- **Note**: In `LCAO` calculation, the path of a numerical descriptor (an `orb` file) is needed to be specified under the `NUMERICAL_DESCRIPTOR` tag in the `STRU` file. For example:

  ```text
  NUMERICAL_ORBITAL
  H_gga_8au_60Ry_2s1p.orb
  O_gga_7au_60Ry_2s2p1d.orb

  NUMERICAL_DESCRIPTOR
  jle.orb
  ```

- **Default**: False

### deepks_scf

- **Type**: Boolean
- **Availability**: numerical atomic orbital basis
- **Description**: perform self-consistent field iteration in DeePKS method
- **Note**: A trained, traced model file is needed.
- **Default**: False

### deepks_equiv

- **Type**: Boolean
- **Availability**: numerical atomic orbital basis
- **Description**: whether to use equivariant version of DeePKS
- **Note**: the equivariant version of DeePKS-kit is still under development, so this feature is currently only intended for internal usage.
- **Default**: False

### deepks_model

- **Type**: String
- **Availability**: numerical atomic orbital basis and `deepks_scf` is true
- **Description**: the path of the trained, traced neural network model file generated by [deepks-kit](https://github.com/deepmodeling/deepks-kit)
- **Default**: None

### bessel_descriptor_lmax

- **Type**: Integer
- **Availability**: `gen_bessel` calculation
- **Description**: the maximum angular momentum of the Bessel functions generated as the projectors in DeePKS
- **NOte**: To generate such projectors, set calculation type to `gen_bessel` in ABACUS. See also [calculation](#calculation).
- **Default**: 2

### bessel_descriptor_ecut

- **Type**: Real
- **Availability**: `gen_bessel` calculation
- **Description**: energy cutoff of Bessel functions
- **Default**: same as ecutwfc
- **Unit**: Ry

### bessel_descriptor_tolerence

- **Type**: Real
- **Availability**: `gen_bessel` calculation
- **Description**: tolerance for searching the zeros of Bessel functions
- **Default**: 1.0e-12

### bessel_descriptor_rcut

- **Type**: Real
- **Availability**: `gen_bessel` calculation
- **Description**: cutoff radius of Bessel functions
- **Default**: 6.0
- **Unit**: Bohr

### bessel_descriptor_smooth

- **Type**: Boolean
- **Availability**: `gen_bessel` calculation
- **Description**: smooth the Bessel functions at radius cutoff
- **Default**: False

### bessel_descriptor_sigma

- **Type**: Real
- **Availability**: `gen_bessel` calculation
- **Description**: smooth parameter at the cutoff radius of projectors
- **Default**: 0.1
- **Unit**: Bohr

### deepks_bandgap

- **Type**: Boolean
- **Availability**: numerical atomic orbital basis and `deepks_scf` is true
- **Description**: include bandgap label for DeePKS training
- **Default**: False

### deepks_v_delta

- **Type**: int
- **Availability**: numerical atomic orbital basis
- **Description**: Include V_delta label for DeePKS training. When `deepks_out_labels` is true and `deepks_v_delta` > 0, ABACUS will output h_base.npy, v_delta.npy and h_tot.npy(h_tot=h_base+v_delta). 
  Meanwhile, when `deepks_v_delta` equals 1, ABACUS will also output v_delta_precalc.npy, which is used to calculate V_delta during DeePKS training. However, when the number of atoms grows, the size of v_delta_precalc.npy will be very large. In this case, it's recommended to set `deepks_v_delta` as 2, and ABACUS will output phialpha.npy and grad_evdm.npy but not v_delta_precalc.npy. These two files are small and can be used to calculate v_delta_precalc in the procedure of training DeePKS.
- **Default**: 0

### deepks_out_unittest

- **Type**: Boolean
- **Description**: generate files for constructing DeePKS unit test
- **Note**: Not relevant when running actual calculations. When set to 1, ABACUS needs to be run with only 1 process.
- **Default**: False

[back to top](#full-list-of-input-keywords)

## OFDFT: orbital free density functional theory

### of_kinetic

- **Type**: String
- **Availability**: OFDFT
- **Description**: The type of KEDF (kinetic energy density functional).
  - **wt**: Wang-Teter.
  - **tf**: Thomas-Fermi.
  - **vw**: von Weizsäcker.
  - **tf+**: TF$\rm{\lambda}$vW, the parameter $\rm{\lambda}$ can be set by `of_vw_weight`.
  - **lkt**: Luo-Karasiev-Trickey.
- **Default**: wt

### of_method

- **Type**: String
- **Availability**: OFDFT
- **Description**: The optimization method used in OFDFT.
  - **cg1**: Polak-Ribiere. Standard CG algorithm.
  - **cg2**: Hager-Zhang (generally faster than cg1).
  - **tn**: Truncated Newton algorithm.
- **Default**:tn

### of_conv

- **Type**: String
- **Availability**: OFDFT
- **Description**: Criterion used to check the convergence of OFDFT.
  - **energy**: Ttotal energy changes less than `of_tole`.
  - **potential**: The norm of potential is less than `of_tolp`.
  - **both**: Both energy and potential must satisfy the convergence criterion.
- **Default**: energy

### of_tole

- **Type**: Real
- **Availability**: OFDFT
- **Description**: Tolerance of the energy change for determining the convergence.
- **Default**: 2e-6
- **Unit**: Ry

### of_tolp

- **Type**: Real
- **Availability**: OFDFT
- **Description**: Tolerance of potential for determining the convergence.
- **Default**: 1e-5
- **Unit**: Ry

### of_tf_weight

- **Type**: Real
- **Availability**: OFDFT with `of_kinetic=tf, tf+, wt`
- **Description**: Weight of TF KEDF (kinetic energy density functional).
- **Default**: 1.0

### of_vw_weight

- **Type**: Real
- **Availability**: OFDFT with `of_kinetic=vw, tf+, wt, lkt`
- **Description**: Weight of vW KEDF (kinetic energy density functional).
- **Default**: 1.0

### of_wt_alpha

- **Type**: Real
- **Availability**: OFDFT with `of_kinetic=wt`
- **Description**: Parameter alpha of WT KEDF (kinetic energy density functional).
- **Default**: $5/6$

### of_wt_beta

- **Type**: Real
- **Availability**: OFDFT with `of_kinetic=wt`
- **Description**: Parameter beta of WT KEDF (kinetic energy density functional).
- **Default**: $5/6$

### of_wt_rho0

- **Type**: Real
- **Availability**: OFDFT with `of_kinetic=wt`
- **Description**: The average density of system.
- **Default**: 0.0
- **Unit**: Bohr^-3

### of_hold_rho0

- **Type**: Boolean
- **Availability**: OFDFT with `of_kinetic=wt`
- **Description**: Whether to fix the average density rho0.
  - **True**: rho0 will be fixed even if the volume of system has changed, it will be set to True automatically if `of_wt_rho0` is not zero.
  - **False**: rho0 will change if volume of system has changed.
- **Default**: False

### of_lkt_a

- **Type**: Real
- **Availability**: OFDFT with `of_kinetic=lkt`
- **Description**: Parameter a of LKT KEDF (kinetic energy density functional).
- **Default**: 1.3

### of_read_kernel

- **Type**: Boolean
- **Availability**: OFDFT with `of_kinetic=wt`
- **Description**: Whether to read in the kernel file.
  - **True**: The kernel of WT KEDF (kinetic energy density functional) will be filled from the file specified by `of_kernel_file`.
  - **False**: The kernel of WT KEDF (kinetic energy density functional) will be filled from formula.
- **Default**: False

### of_kernel_file

- **Type**: String
- **Availability**: OFDFT with `of_read_kernel=True`
- **Description**: The name of WT kernel file.
- **Default**: WTkernel.txt

### of_full_pw

- **Type**: Boolean
- **Availability**: OFDFT
- **Description**: Whether to use full planewaves.
  - **True**: Ecut will be ignored while collecting planewaves, so that all planewaves will be used in FFT.
  - **False**: Only use the planewaves inside ecut, the same as KSDFT.
- **Default**: True

### of_full_pw_dim

- **Type**: Integer
- **Availability**: OFDFT with `of_full_pw = True`
- **Description**: Specify the parity of FFT dimensions.
  - **0**: either odd or even.
  - **1**: odd only.
  - **2**: even only.

  Note: Even dimensions may cause slight errors in FFT. It should be ignorable in ofdft calculation, but it may make Cardinal B-spline interpolation unstable, so please set `of_full_pw_dim = 1` if `nbspline != -1`.

- **Default**: 0

[back to top](#full-list-of-input-keywords)

## Electric field and dipole correction

These variables are relevant to electric field and dipole correction

### efield_flag

- **Type**: Boolean
- **Description**: added the electric field.
  - True: A saw-like potential simulating an electric field is added to the bare ionic potential.
  - False: Not added the electric field.
- **Default**: False

### dip_cor_flag

- **Type**: Boolean
- **Availability**: with dip_cor_flag = True and efield_flag = True.
- **Description**: Added a dipole correction to the bare ionic potential.
  - True：A dipole correction is also added to the bare ionic potential.
  - False: A dipole correction is not added to the bare ionic potential.

> Note: If you want no electric field, parameter efield_amp  should be zero. Must be used ONLY in a slab geometry for surface alculations, with the discontinuity FALLING IN THE EMPTY SPACE.

- **Default**: False

### efield_dir

- **Type**: Integer
- **Availability**: with efield_flag = True.
- **Description**: The direction of the electric field or dipole correction is parallel to the reciprocal lattice vector, so the potential is constant in planes defined by FFT grid points, efield_dir can set to 0, 1 or 2.
  - 0: parallel to $b_1=\frac{2\pi(a_2\times a_3)}{a_1\cdot(a_2\times a_3)}$
  - 1: parallel to $b_2=\frac{2\pi(a_3\times a_1)}{a_1\cdot(a_2\times a_3)}$
  - 2: parallel to $b_3=\frac{2\pi(a_1\times a_2)}{a_1\cdot(a_2\times a_3)}$
- **Default**: 2

### efield_pos_max

- **Type**: Real
- **Availability**: with efield_flag = True.
- **Description**: Position of the maximum of the saw-like potential along crystal axis efield_dir, within the  unit cell, 0 <= efield_pos_max < 1.
- **Default**: Autoset to `center of vacuum - width of vacuum / 20`

### efield_pos_dec

- **Type**: Real
- **Availability**: with efield_flag = True.
- **Description**: Zone in the unit cell where the saw-like potential decreases, 0 < efield_pos_dec < 1.
- **Default**: Autoset to `width of vacuum / 10`

### efield_amp

- **Type**: Real
- **Availability**: with efield_flag = True.
- **Description**: Amplitude of the electric field. The saw-like potential increases with slope efield_amp  in the region from  efield_pos_max+efield_pos_dec-1) to (efield_pos_max), then decreases until (efield_pos_max+efield_pos_dec), in units of the crystal vector efield_dir.

> Note: The change of slope of this potential must be located in the empty region, or else unphysical forces will result.

- **Default**: 0.0
- **Unit**: a.u., 1 a.u. = 51.4220632*10^10 V/m.

[back to top](#full-list-of-input-keywords)

## Gate field (compensating charge)

These variables are relevant to gate field (compensating charge) [Detailed introduction](https://github.com/deepmodeling/abacus-develop/blob/develop/docs/advanced/scf/advanced.md#compensating-charge)

### gate_flag

- **Type**: Boolean
- **Description**: Controls the addition of compensating charge by a charged plate for charged cells.
  - true: A charged plate is placed at the **zgate** position to add compensating charge. The direction is determined by **efield_dir**.
  - false: No compensating charge is added.
- **Default**: false

### zgate

- **Type**: Real
- **Description**: position of the charged plate in the unit cell
- **Unit**: Unit cell size
- **Default**: 0.5
- **Constraints**: 0 <= **zgate** < 1

### block

- **Type**: Boolean
- **Description**: Controls the addition of a potential barrier to prevent electron spillover.
  - true: A potential barrier is added from **block_down** to **block_up** with a height of **block_height**. If **dip_cor_flag** is set to true, **efield_pos_dec** is used to smoothly increase and decrease the potential barrier.
  - false: No potential barrier is added.
- **Default**: false

### block_down

- **Type**: Real
- **Description**: lower beginning of the potential barrier
- **Unit**: Unit cell size
- **Default**: 0.45
- **Constraints**: 0 <= **block_down** < **block_up** < 1

### block_up

- **Type**: Real
- **Description**: upper beginning of the potential barrier
- **Unit**: Unit cell size
- **Default**: 0.55
- **Constraints**: 0 <= **block_down** < **block_up** < 1

### block_height

- **Type**: Real
- **Description**: height of the potential barrier
- **Unit**: Rydberg
- **Default**: 0.1

[back to top](#full-list-of-input-keywords)

## Exact Exchange

These variables are relevant when using hybrid functionals.

**Availablity**: *[dft_functional](#dft_functional)==hse/hf/pbe0/scan0/opt_orb* or *[rpa](#rpa)==True*, and *[basis_type](#basis_type)==lcao/lcao_in_pw*

### exx_hybrid_alpha

- **Type**: Real
- **Description**: fraction of Fock exchange in hybrid functionals, so that $E_{X}=\alpha E_{X}+(1-\alpha)E_{X,\text{LDA/GGA}}$
- **Default**:
  - 1: if *[dft_functional](#dft_functional)==hf*
  - 0.25: else

### exx_hse_omega

- **Type**: Real
- **Description**: range-separation parameter in HSE functional, such that $1/r=\text{erfc}(\omega r)/r+\text{erf}(\omega r)/r$
- **Default**: 0.11

### exx_separate_loop

- **Type**: Boolean
- **Description**: There are two types of iterative approaches provided by ABACUS to evaluate Fock exchange.
  - False: Start with a GGA-Loop, and then Hybrid-Loop, in which EXX Hamiltonian $H_{exx}$ is updated with electronic iterations.
  - True: A two-step method is employed, i.e. in the inner iterations, density matrix is updated, while in the outer iterations, $H_{exx}$ is calculated based on density matrix that converges in the inner iteration.
- **Default**: True

### exx_hybrid_step

- **Type**: Integer
- **Availability**: *[exx_separate_loop](#exx_separate_loop)==1*
- **Description**: the maximal iteration number of the outer-loop, where the Fock exchange is calculated
- **Default**: 100

### exx_mixing_beta

- **Type**: Real
- **Availability**: *[exx_separate_loop](#exx_separate_loop)==1*
- **Description**: mixing_beta for densty matrix in each iteration of the outer-loop
- **Default**: 1.0

### exx_lambda

- **Type**: Real
- **Availability**: *[basis_type](#basis_type)==lcao_in_pw*
- **Description**: It is used to compensate for divergence points at G=0 in the evaluation of Fock exchange using *lcao_in_pw* method.
- **Default**: 0.3

### exx_pca_threshold

- **Type**: Real
- **Description**: To accelerate the evaluation of four-center integrals ($ik|jl$), the product of atomic orbitals are expanded in the basis of auxiliary basis functions (ABF): $\Phi_{i}\Phi_{j}\sim C^{k}_{ij}P_{k}$. The size of the ABF (i.e. number of $P_{k}$) is reduced using principal component analysis. When a large PCA threshold is used, the number of ABF will be reduced, hence the calculation becomes faster. However, this comes at the cost of computational accuracy. A relatively safe choice of the value is 1e-4.
- **Default**: 1E-4

### exx_c_threshold

- **Type**: Real
- **Description**: See also the entry [exx_pca_threshold](#exx_pca_threshold). Smaller components (less than exx_c_threshold) of the $C^{k}_{ij}$ matrix are neglected to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-4.
- **Default**: 1E-4

### exx_v_threshold

- **Type**: Real
- **Description**: See also the entry [exx_pca_threshold](#exx_pca_threshold). With the approximation $\Phi_{i}\Phi_{j}\sim C^{k}_{ij}P_{k}$, the four-center integral in Fock exchange is expressed as $(ik|jl)=\Sigma_{a,b}C^{a}_{ij}V_{ab}C^{b}_{kl}$, where $V_{ab}=(P_{a}|P_{b})$ is a double-center integral. Smaller values of the V matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 0, i.e. no truncation.
- **Default**: 1E-1

### exx_dm_threshold

- **Type**: Real
- **Description**: The Fock exchange can be expressed as $\Sigma_{k,l}(ik|jl)D_{kl}$ where D is the density matrix. Smaller values of the density matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-4.
- **Default**: 1E-4

### exx_c_grad_threshold

- **Type**: Real
- **Description**: See also the entry [exx_pca_threshold](#exx_pca_threshold). $\nabla C^{k}_{ij}$ is used in force and stress. Smaller components (less than exx_c_grad_threshold) of the $\nabla C^{k}_{ij}$ matrix are neglected to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-4.
- **Default**: 1E-4

### exx_v_grad_threshold

- **Type**: Real
- **Description**: See also the entry [exx_pca_threshold](#exx_pca_threshold). With the approximation $\Phi_{i}\Phi_{j}\sim C^{k}_{ij}P_{k}$, the four-center integral in Fock exchange is expressed as $(ik|jl)=\Sigma_{a,b}C^{a}_{ij}V_{ab}C^{b}_{kl}$, where $V_{ab}=(P_{a}|P_{b})$ is a double-center integral. $\nabla V_{ab}$ is used in force and stress. Smaller values of the V matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 0, i.e. no truncation.
- **Default**: 1E-1

### exx_schwarz_threshold

- **Type**: Real
- **Description**: In practice the four-center integrals are sparse, and using Cauchy-Schwartz inequality, we can find an upper bound of each integral before carrying out explicit evaluations. Those that are smaller than exx_schwarz_threshold will be truncated. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-5.  (Currently not used)
- **Default**: 0

### exx_cauchy_threshold

- **Type**: Real
- **Description**: In practice the Fock exchange matrix is sparse, and using Cauchy-Schwartz inequality, we can find an upper bound of each matrix element before carrying out explicit evaluations. Those that are smaller than exx_cauchy_threshold will be truncated. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-7.
- **Default**: 1E-7

### exx_cauchy_force_threshold

- **Type**: Real
- **Description**: In practice the Fock exchange matrix in force is sparse, and using Cauchy-Schwartz inequality, we can find an upper bound of each matrix element before carrying out explicit evaluations. Those that are smaller than exx_cauchy_force_threshold will be truncated. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-7.
- **Default**: 1E-7

### exx_cauchy_stress_threshold

- **Type**: Real
- **Description**: In practice the Fock exchange matrix in stress is sparse, and using Cauchy-Schwartz inequality, we can find an upper bound of each matrix element before carrying out explicit evaluations. Those that are smaller than exx_cauchy_stress_threshold will be truncated. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-7.
- **Default**: 1E-7

### exx_ccp_threshold

- **Type**: Real
- **Description**: It is related to the cutoff of on-site Coulomb potentials. (Currently not used)
- **Default**: 1e-8

### exx_ccp_rmesh_times

- **Type**: Real
- **Description**: This parameter determines how many times larger the radial mesh required for calculating Columb potential is to that of atomic orbitals. The value should be at least 1. Reducing this value can effectively increase the speed of self-consistent calculations using hybrid functionals.
- **Default**:
  - 5: if *[dft_functional](#dft_functional)==hf/pbe0/scan0/muller/power/wp22*
  - 1.5: if *[dft_functional](#dft_functional)==hse/cwp22*
  - 1: else

### exx_distribute_type

- **Type**: String
- **Description**: When running in parallel, the evaluation of Fock exchange is done by distributing atom pairs on different processes, then gather the results. exx_distribute_type governs the mechanism of distribution. Available options are `htime`, `order`, `kmean1` and `kmeans2`.
  - `order`: Atom pairs are simply distributed by their orders.
  - `htime`: The balance in time is achieved on each processor, hence if the memory is sufficient, this is the recommended method.
  - `kmeans1` ,   `kmeans2`: Two methods where the k-means clustering method is used to reduce memory requirement. They might be necessary for very large systems. (Currently not used)
- **Default**: `htime`

### exx_opt_orb_lmax

- **Type**: Integer
- **Availability**: *[dft_functional](#dft_functional)==opt_orb*
- **Description**: The maximum l of the spherical Bessel functions, when the radial part of opt-ABFs are generated as linear combinations of spherical Bessel functions. A reasonable choice is 2.
- **Default**: 0

### exx_opt_orb_ecut

- **Type**: Real
- **Availability**: *[dft_functional](#dft_functional)==opt_orb*
- **Description**: The cut-off of plane wave expansion, when the plane wave basis is used to optimize the radial ABFs. A reasonable choice is 60.
- **Default**: 0
- **Unit**: Ry

### exx_opt_orb_tolerence

- **Type**: Real
- **Availability**: *[dft_functional](#dft_functional)==opt_orb*
- **Description**: The threshold when solving for the zeros of spherical Bessel functions. A reasonable choice is 1e-12.
- **Default**: 0

### exx_real_number

- **Type**: Boolean
- **Description**:
  - True: Enforce LibRI to use `double` data type.
  - False: Enforce LibRI to use `complex` data type.
  Setting it to True can effectively improve the speed of self-consistent calculations with hybrid functionals.
- **Default**: depends on the [gamma_only](#gamma_only) option
  - True: if gamma_only
  - False: else

### rpa_ccp_rmesh_times

- **Type**: Real
- **Description**: How many times larger the radial mesh required is to that of atomic orbitals in the postprocess calculation of the **bare** Coulomb matrix for RPA, GW, etc.
- **Default**: 10

### exx_symmetry_realspace

- **Type**: Boolean
- **Availability**: *[symmetry](#symmetry)==1* and exx calculation  (*[dft_fuctional](#dft_functional)==hse/hf/pbe0/scan0/opt_orb* or *[rpa](#rpa)==True*)
- **Description**: 
  - False: only rotate k-space density matrix D(k) from irreducible k-points to accelerate diagonalization
  - True: rotate both D(k) and Hexx(R) to accelerate both diagonalization and EXX calculation
- **Default**: True

### out_ri_cv

- **Type**: Boolean
- **Description**: Whether to output the coefficient tensor C(R) and ABFs-representation Coulomb matrix V(R) for each atom pair and cell in real space.
- **Default**: false

[back to top](#full-list-of-input-keywords)

## Molecular dynamics

These variables are used to control molecular dynamics calculations. For more information, please refer to [md.md](../md.md#molecular-dynamics) in detail.

### md_type

- **Type**: String
- **Description**: Control the algorithm to integrate the equation of motion for molecular dynamics (MD), see [md.md](../md.md#molecular-dynamics) in detail.

  - fire: a MD-based relaxation algorithm, named fast inertial relaxation engine.
  - nve: NVE ensemble with velocity Verlet algorithm.
  - nvt: NVT ensemble, see [md_thermostat](#md_thermostat) in detail.
  - npt: Nose-Hoover style NPT ensemble, see [md_pmode](#md_pmode) in detail.
  - langevin: NVT ensemble with Langevin thermostat, see [md_damp](#md_damp) in detail.
  - msst: MSST method, see [msst_direction](#msst_direction), [msst_vel](#msst_vel), [msst_qmass](#msst_qmass), [msst_vis](#msst_vis), [msst_tscale](#msst_tscale) in detail.

- **Default**: nvt

### md_nstep

- **Type**: Integer
- **Description**: The total number of molecular dynamics steps.
- **Default**: 10

### md_dt

- **Type**: Real
- **Description**: The time step used in molecular dynamics calculations.
- **Default**: 1.0
- **Unit**: fs

### md_thermostat

- **Type**: String
- **Description**: Specify the temperature control method used in NVT ensemble.

  - nhc: Nose-Hoover chain, see [md_tfreq](#md_tfreq) and [md_tchain](#md_tchain) in detail.
  - anderson: Anderson thermostat, see [md_nraise](#md_nraise) in detail.
  - berendsen: Berendsen thermostat, see [md_nraise](#md_nraise) in detail.
  - rescaling: velocity Rescaling method 1, see [md_tolerance](#md_tolerance) in detail.
  - rescale_v: velocity Rescaling method 2, see [md_nraise](#md_nraise) in detail.

- **Default**: nhc

### md_tfirst, md_tlast

- **Type**: Real
- **Description**: The temperature used in molecular dynamics calculations.

  If `md_tfirst` is unset or less than zero, [init_vel](#init_vel) is autoset to be `true`. If [init_vel](#init_vel) is `true`, the initial temperature will be determined by the velocities read from `STRU`. In this case, if velocities are unspecified in `STRU`, the initial temperature is set to zero.

  If `md_tfirst` is set to a positive value and [init_vel](#init_vel) is `true` simultaneously, please make sure they are consistent, otherwise abacus will exit immediately.

  Note that `md_tlast` is only used in NVT/NPT simulations. If `md_tlast` is unset or less than zero, `md_tlast` is set to `md_tfirst`. If `md_tlast` is set to be different from `md_tfirst`, ABACUS will automatically change the temperature from `md_tfirst` to `md_tlast`.
- **Default**: No default
- **Unit**: K

### md_restart

- **Type**: Boolean
- **Description**: Control whether to restart molecular dynamics calculations and time-dependent density functional theory calculations.
  - True: ABACUS will read in `${read_file_dir}/Restart_md.dat` to determine the current step `${md_step}`, then read in the corresponding `STRU_MD_${md_step}` in the folder `OUT.$suffix/STRU/` automatically. For tddft, ABACUS will also read in `WFC_NAO_K${kpoint}` of the last step (You need to set out_wfc_lcao=1 and out_app_flag=0 to obtain this file).
  - False: ABACUS will start molecular dynamics calculations normally from the first step.
- **Default**: False

### md_restartfreq

- **Type**: Integer
- **Description**: The output frequency of `OUT.${suffix}/Restart_md.dat` and structural files in the directory `OUT.${suffix}/STRIU/`, which are used to restart molecular dynamics calculations, see [md_restart](#md_restart) in detail.
- **Default**: 5

### md_dumpfreq

- **Type**: Integer
- **Description**: The output frequency of `OUT.${suffix}/MD_dump` in molecular dynamics calculations, which including the information of lattices and atoms.
- **Default**: 1

### dump_force

- **Type**: Boolean
- **Description**: Whether to output atomic forces into the file `OUT.${suffix}/MD_dump`.
- **Default**: True

### dump_vel

- **Type**: Boolean
- **Description**: Whether to output atomic velocities into the file `OUT.${suffix}/MD_dump`.
- **Default**: True

### dump_virial

- **Type**: Boolean
- **Description**: Whether to output lattice virials into the file `OUT.${suffix}/MD_dump`.
- **Default**: True

### md_seed

- **Type**: Integer
- **Description**: The random seed to initialize random numbers used in molecular dynamics calculations.
  - \< 0: No srand() function is called.
  - \>= 0: The function srand(md_seed) is called.
- **Default**: -1

### md_tfreq

- **Type**: Real
- **Description**: Control the frequency of temperature oscillations during the simulation. If it is too large, the temperature will fluctuate violently; if it is too small, the temperature will take a very long time to equilibrate with the atomic system.

  Note: It is a system-dependent empirical parameter, ranging from 1/(40\*md_dt) to 1/(100\*md_dt). An improper choice might lead to the failure of jobs.

- **Default**: 1/40/md_dt
- **Unit**: $\mathrm{fs^{-1}}$

### md_tchain

- **Type**: Integer
- **Description**: Number of thermostats coupled with the particles in the NVT/NPT ensemble based on the Nose-Hoover style non-Hamiltonian equations of motion.
- **Default**: 1

### md_pmode

- **Type**: String
- **Description**: Specify the cell fluctuation mode in NPT ensemble based on the Nose-Hoover style non-Hamiltonian equations of motion.
  - iso: The three diagonal elements of the lattice are fluctuated isotropically.
  - aniso: The three diagonal elements of the lattice are fluctuated anisotropically.
  - tri: The lattice must be a lower-triangular matrix, and all six freedoms are fluctuated.
- **Default**: iso
- **Relavent**: [md_tfreq](#md_tfreq), [md_tchain](#md_tchain), [md_pcouple](#md_pcouple), [md_pfreq](#md_pfreq), and [md_pchain](#md_pchain).

<!-- ### md_prec_level

- **Type**: Integer
- **Description**: Determine the precision level of variable-cell molecular dynamics calculations.
  - 0: FFT grids do not change, only G vectors and K vectors are changed due to the change of lattice vector. This level is suitable for cases where the variation of the volume and shape is not large, and the efficiency is relatively higher.
  - 2: FFT grids change per step. This level is suitable for cases where the variation of the volume and shape is large, such as the MSST method. However, accuracy comes at the cost of efficiency.

- **Default**: 0 -->

### ref_cell_factor

- **Type**: Real
- **Description**: Construct a reference cell bigger than the initial cell. The reference cell has to be large enough so that the lattice vectors of the fluctuating cell do not exceed the reference lattice vectors during MD. Typically, 1.02 ~ 1.10 is sufficient. However, the cell fluctuations depend on the specific system and thermodynamic conditions. So users must test for a proper choice. This parameters should be used in conjunction with [erf_ecut](#erf_ecut), [erf_height](#erf_height), and [erf_sigma](#erf_sigma).
- **Default**: 1.0

### md_pcouple

- **Type**: String
- **Description**: The coupled lattice vectors will scale proportionally in NPT ensemble based on the Nose-Hoover style non-Hamiltonian equations of motion.
  - none: Three lattice vectors scale independently.
  - xyz: Lattice vectors x, y, and z scale proportionally.
  - xy: Lattice vectors x and y scale proportionally.
  - xz: Lattice vectors x and z scale proportionally.
  - yz: Lattice vectors y and z scale proportionally.
- **Default**: none

### md_pfirst, md_plast

- **Type**: Real
- **Description**: The target pressure used in NPT ensemble simulations, the default value of `md_plast` is `md_pfirst`. If `md_plast` is set to be different from `md_pfirst`, ABACUS will automatically change the target pressure from `md_pfirst` to `md_plast`.
- **Default**: -1.0
- **Unit**: kbar

### md_pfreq

- **Type**: Real
- **Description**: The frequency of pressure oscillations during the NPT ensemble simulation. If it is too large, the pressure will fluctuate violently; if it is too small, the pressure will take a very long time to equilibrate with the atomic system.

  Note: It is a system-dependent empirical parameter. An improper choice might lead to the failure of jobs.
- **Default**: 1/400/md_dt
- **Unit**: $\mathrm{kbar^{-1}}$

### md_pchain

- **Type**: Integer
- **Description**: The number of thermostats coupled with the barostat in the NPT ensemble based on the Nose-Hoover style non-Hamiltonian equations of motion.
- **Default**: 1

### lj_rule

- **Type**: Integer
- **Description**: The Lennard-Jones potential between two atoms equals: 
  $$V_{LJ}(r_{ij})=4\epsilon_{ij}\left(\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12}-\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{6}\right)=\frac{C_{ij}^{(12)}}{{r_{ij}}^{12}}-\frac{C_{ij}^{(6)}}{{r_{ij}}^{6}}.$$ 
  
  The parameters [lj_epsilon](#lj_epsilon) and [lj_sigma](#lj_sigma) should be multiple-component vectors. For example, there are two choices in the calculations of 3 atom species: 

  Supply six-component vectors that describe the interactions between all possible atom pairs. The six-component vectors represent lower triangular symmetric matrixs, and the correspondence between the vector component $\sigma _k$ and the matrix element $\sigma (i,j)$ is
  $$k= i(i+1)/2 +j$$
  
  Supply three-component vectors that describe the interactions between atoms of the same species. In this case, two types of combination rules can be used to construct non-diagonal elements in the parameter matrix. 

  - 1: geometric average:
  $$\begin{array}{rcl}C_{ij}^{(6)}&=&\left(C_{ii}^{(6)}C_{jj}^{(6)}\right)^{1/2}\\C_{ij}^{(12)}&=&\left(C_{ii}^{(12)}C_{jj}^{(12)}\right)^{1/2}\end{array}$$

  - 2: arithmetic average:
  $$\begin{array}{rcl}\sigma_{ij}&=&\frac{1}{2}\left(\sigma_{ii}+\sigma_{jj}\right)\\ \epsilon_{ij}&=&\left(\epsilon_{ii}\epsilon_{jj}\right)^{1/2}\end{array}$$
- **Default**: 2

### lj_eshift

- **Type**: Boolean
- **Description**: It True, the LJ potential is shifted by a constant such that it is zero at the cut-off distance.
- **Default**: False

### lj_rcut

- **Type**: Real
- **Description**: Cut-off radius for Leonard Jones potential, beyond which the interaction will be neglected. It can be a single value, which means that all pairs of atoms types share the same cut-off radius. Otherwise, it should be a multiple-component vector, containing $N(N+1)/2$ values, see details in [lj_rule](#lj_rule).
- **Default**: No default
- **Unit**: Angstrom

### lj_epsilon

- **Type**: Real
- **Description**: The vector representing the $\epsilon$ matrix for Leonard Jones potential. See details in [lj_rule](#lj_rule).
- **Default**: No default
- **Unit**: eV

### lj_sigma

- **Type**: Real
- **Description**: The vector representing the $\sigma$ matrix for Leonard Jones potential. See details in [lj_rule](#lj_rule).
- **Default**: No default
- **Unit**: Angstrom

### pot_file

- **Type**: String
- **Description**: The filename of DP potential files, see [md.md](../md.md#dpmd) in detail.
- **Default**: graph.pb

### dp_rescaling

- **Type**: Real
- **Availability**: [esolver_type](#esolver_type) = `dp`.
- **Description**: Rescaling factor to use a temperature-dependent DP. Energy, stress and force calculated by DP will be multiplied by this factor.
- **Default**: 1.0

### dp_fparam

- **Type**: Real
- **Availability**: [esolver_type](#esolver_type) = `dp`.
- **Description**: The frame parameter for dp potential. The array size is dim_fparam, then all frames are assumed to be provided with the same fparam.
- **Default**: {}

### dp_aparam

- **Type**: Real
- **Availability**: [esolver_type](#esolver_type) = `dp`.
- **Description**: The atomic parameter for dp potential. The array size can be (1) natoms x dim_aparam, then all frames are assumed to be provided with the same aparam; (2) dim_aparam, then all frames and atoms are assumed to be provided with the same aparam.
- **Default**: {}

### msst_direction

- **Type**: Integer
- **Description**: The direction of the shock wave in the MSST method.
  - 0: x direction
  - 1: y direction
  - 2: z direction
- **Default**: 2

### msst_vel

- **Type**: Real
- **Description**: The velocity of the shock wave in the MSST method.
- **Default**: 0.0
- **Unit**: Angstrom/fs

### msst_vis

- **Type**: Real
- **Description**: Artificial viscosity in the MSST method.
- **Default**: 0.0
- **Unit**: g/(mol\*Angstrom\*fs)

### msst_tscale

- **Type**: Real
- **Description**: The reduction percentage of the initial temperature used to compress volume in the MSST method.
- **Default**: 0.01

### msst_qmass

- **Type**: Real
- **Description**: Inertia of the extended system variable. You should set a number larger than 0.
- **Default**: No default
- **Unit**: $\mathrm{g^{2}/(mol^{2}*Angstrom^{4})}$

### md_damp

- **Type**: Real
- **Description**: The damping parameter used to add fictitious force in the Langevin method.
- **Default**: 1.0
- **Unit**: fs

### md_tolerance

- **Type**: Real
- **Description**: Thr temperature tolerance for velocity rescaling. Velocities are rescaled if the current and target temperature differ more than `md_tolerance`.
- **Default**: 100.0
- **Unit**: K

### md_nraise

- **Type**: Integer
- **Description**:
  - Anderson: The "collision frequency" parameter is given as 1/`md_nraise`.
  - Berendsen: The "rise time" parameter is given in units of the time step: tau = `md_nraise`*`md_dt`, so `md_dt`/tau = 1/`md_nraise`.
  - Rescale_v: Every `md_nraise` steps the current temperature is rescaled to the target temperature.
- **Default**: 1

### cal_syns

- **Type**: Boolean
- **Description**: Whether the asynchronous overlap matrix is calculated for Hefei-NAMD.
- **Default**: False

### dmax

- **Type**: Real
- **Description**: The maximum displacement of all atoms in one step. This parameter is useful when [cal_syns](#cal_syns) = True.
- **Default**: 0.01
- **Unit**: bohr

[back to top](#full-list-of-input-keywords)

## DFT+*U* correction

These variables are used to control DFT+U correlated parameters

### dft_plus_u

- **Type**: Integer
- **Description**: Determines whether to calculate the plus U correction, which is especially important for correlated electrons.
  - 1: Calculate plus U correction with radius-adjustable localized projections (with parameter `onsite_radius`).
  - 2: Calculate plus U correction using first zeta of NAOs as projections (this is old method for testing).
  - 0: Do not calculate plus U correction.
- **Default**: 0

### orbital_corr

- **Type**: Integer
- **Description**: Specifies which orbits need plus U correction for each atom type ($l_1,l_2,l_3,\ldots$ for atom type 1, 2, 3, respectively).
  - -1: The plus U correction will not be calculated for this atom.
  - 1: For p-electron orbits, the plus U correction is needed.
  - 2: For d-electron orbits, the plus U correction is needed.
  - 3: For f-electron orbits, the plus U correction is needed.
- **Default**: -1

### hubbard_u

- **Type**: Real
- **Description**: Specifies the Hubbard Coulomb interaction parameter U (eV) in plus U correction, which should be specified for each atom unless the Yukawa potential is used.

> Note: Since only the simplified scheme by Duradev is implemented, the 'U' here is actually U-effective, which is given by Hubbard U minus Hund J.

- **Default**: 0.0

### yukawa_potential

- **Type**: Boolean
- **Description**: Determines whether to use the local screen Coulomb potential method to calculate the values of U and J.
  - True: `hubbard_u` does not need to be specified.
  - False: `hubbard_u` does need to be specified.
- **Default**: False

### yukawa_lambda

- **Type**: Real
- **Availability**: DFT+U with `yukawa_potential` = True.
- **Description**: The screen length of Yukawa potential. If left to default, the screen length will be calculated as an average of the entire system. It's better to stick to the default setting unless there is a very good reason.
- **Default**: Calculated on the fly.

### uramping

- **Type**: Real
- **Unit**: eV
- **Availability**: DFT+U calculations with `mixing_restart > 0`.
- **Description**: Once `uramping` > 0.15 eV. DFT+U calculations will start SCF with U = 0 eV, namely normal LDA/PBE calculations. Once SCF restarts when `drho<mixing_restart`, U value will increase by `uramping` eV. SCF will repeat above calcuations until U values reach target defined in `hubbard_u`. As for `uramping=1.0 eV`, the recommendations of `mixing_restart` is around `5e-4`.
- **Default**: -1.0.

### omc

- **Type**: Integer
- **Description**: The parameter controls the form of occupation matrix control used.
  - 0: No occupation matrix control is performed, and the onsite density matrix will be calculated from wavefunctions in each SCF step.
  - 1: The first SCF step will use an initial density matrix read from a file named `[initial_onsite.dm](http://initial_onsite.dm/)`, but for later steps, the onsite density matrix will be updated.
  - 2: The same onsite density matrix from `initial_onsite.dm` will be used throughout the entire calculation.

> Note : The easiest way to create `initial_onsite.dm` is to run a DFT+U calculation, look for a file named `onsite.dm` in the OUT.prefix directory, and make replacements there. The format of the file is rather straight-forward.

- **Default**: 0

### onsite_radius

- **Type**: Real
- **Availability**: `dft_plus_u` is set to 1
- **Description**:

  - The `Onsite-radius` parameter facilitates modulation of the single-zeta portion of numerical atomic orbitals for projections for DFT+U.
  - The modulation algorithm includes a smooth truncation applied directly to the tail of the original orbital, followed by normalization.  Consider the function:
  $$
  g(r;\sigma)=\begin{cases}
  1-\exp\left(-\frac{(r-r_c)^2}{2\sigma^2}\right), & r < r_c\\
  0, & r \geq r_c
  \end{cases}
  $$
  - where $\sigma$ is a parameter that controls the smoothing interval. A normalized function truncated smoothly at $r_c$ can be represented as:

  $$
  \alpha(r) = \frac{\chi(r)g(r;\sigma)}{\langle\chi(r)g(r;\sigma), \chi(r)g(r;\sigma)\rangle}
  $$

  - To find an appropriate $\sigma$, the optimization process is as follows:

  - Maximizing the overlap integral under a normalization constraint is equivalent to minimizing an error function:

  $$
  \min \langle \chi(r)-\alpha(r), \chi(r)-\alpha(r)\rangle \quad \text{subject to} \quad \langle \alpha(r),\alpha(r)\rangle=1
  $$

  - Similar to the process of generating numerical atomic orbitals, this optimization choice often induces additional oscillations in the outcome. To suppress these oscillations, we may include a derivative term in the objective function ($f'(r)\equiv \mathrm{d}f(r)/\mathrm{d}r$):

  $$
  \min \left[\gamma\langle \chi(r)-\alpha(r), \chi(r)-\alpha(r)\rangle + \langle \chi'(r)-\alpha'(r), \chi'(r)-\alpha'(r)\rangle\right] \quad \text{subject to} \quad \langle \alpha(r),\alpha(r)\rangle=1
  $$

  - where $\gamma$ is a parameter that adjusts the relative weight of the error function to the derivative error function.
- **Unit**: Bohr
- **Default**: 3.0

[back to top](#full-list-of-input-keywords)

## vdW correction

These variables are used to control vdW-corrected related parameters.

### vdw_method

- **Type**: String
- **Description**: Specifies the method used for Van der Waals (VdW) correction. Available options are:
  - `d2`: [Grimme's D2](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.20495) dispersion correction method
  - `d3_0`: [Grimme's DFT-D3(0)](https://aip.scitation.org/doi/10.1063/1.3382344) dispersion correction method (zero-damping)
  - `d3_bj`: [Grimme's DFTD3(BJ)](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.21759) dispersion correction method (BJ-damping)
  - `none`: no vdW correction
- **Default**: none
- **Note**: ABACUS supports automatic setting on DFT-D3 parameters for common functionals after version 3.8.3 (and several develop versions earlier). To benefit from this feature, please specify the parameter `dft_functional` explicitly (for more details on this parameter, please see [dft_functional](#dft_functional)), otherwise the autoset procedure will crash with error message like `cannot find DFT-D3 parameter for XC(***)`. If not satisfied with those in-built parameters, any manually setting on `vdw_s6`, `vdw_s8`, `vdw_a1` and `vdw_a2` will overwrite. 
- **Special**: There are special cases for functional family wB97 (Omega-B97): if want to use the functional wB97X-D3BJ, one needs to specify the `dft_functional` as `HYB_GGA_WB97X_V` and `vdw_method` as `d3_bj`. If want to use the functional wB97X-D3, specify `dft_functional` as `HYB_GGA_WB97X_D3` and `vdw_method` as `d3_0`.

### vdw_s6

- **Type**: Real
- **Availability**: `vdw_method` is set to `d2`, `d3_0`, or `d3_bj`
- **Description**: This scale factor is used to optimize the interaction energy deviations in van der Waals (vdW) corrected calculations. The recommended values of this parameter are dependent on the chosen vdW correction method and the DFT functional being used. For DFT-D2, the recommended values are 0.75 (PBE), 1.2 (BLYP), 1.05 (B-P86), 1.0 (TPSS), and 1.05 (B3LYP). If not set, will use values of PBE functional. For DFT-D3, recommended values with different DFT functionals can be found on the [here](https://github.com/dftd3/simple-dftd3/blob/main/assets/parameters.toml). If not set, will search in ABACUS built-in dataset based on the `dft_functional` keywords. User set value will overwrite the searched value.
- **Default**:
  - 0.75: if `vdw_method` is set to `d2`

### vdw_s8

- **Type**: Real
- **Availability**: `vdw_method` is set to `d3_0` or `d3_bj`
- **Description**: This scale factor is relevant for D3(0) and D3(BJ) van der Waals (vdW) correction methods. The recommended values of this parameter with different DFT functionals can be found on the [webpage](https://github.com/dftd3/simple-dftd3/blob/main/assets/parameters.toml). If not set, will search in ABACUS built-in dataset based on the `dft_functional` keywords. User set value will overwrite the searched value.

### vdw_a1

- **Type**: Real
- **Availability**: `vdw_method` is set to `d3_0` or `d3_bj`
- **Description**: This damping function parameter is relevant for D3(0) and D3(BJ) van der Waals (vdW) correction methods. The recommended values of this parameter with different DFT functionals can be found on the [webpage](https://github.com/dftd3/simple-dftd3/blob/main/assets/parameters.toml). If not set, will search in ABACUS built-in dataset based on the `dft_functional` keywords. User set value will overwrite the searched value.

### vdw_a2

- **Type**: Real
- **Availability**: `vdw_method` is set to `d3_0` or `d3_bj`
- **Description**: This damping function parameter is only relevant for D3(0) and D3(BJ) van der Waals (vdW) correction methods. The recommended values of this parameter with different DFT functionals can be found on the [webpage](https://github.com/dftd3/simple-dftd3/blob/main/assets/parameters.toml). If not set, will search in ABACUS built-in dataset based on the `dft_functional` keywords. User set value will overwrite the searched value.

### vdw_d

- **Type**: Real
- **Availability**: `vdw_method` is set to `d2`
- **Description**: Controls the damping rate of the damping function in the DFT-D2 method.
- **Default**: 20

### vdw_abc

- **Type**: Integer
- **Availability**: `vdw_method` is set to `d3_0` or `d3_bj`
- **Description**: Determines whether three-body terms are calculated for DFT-D3 methods.
  - True: ABACUS will calculate the three-body term.
  - False: The three-body term is not included.
- **Default**: False

### vdw_C6_file

- **Type**: String
- **Availability**: `vdw_method` is set to `d2`
- **Description**: Specifies the name of the file containing $C_6$ parameters for each element when using the D2 method. If not set, ABACUS uses the default $C_6$ parameters (Jnm6/mol) stored in the [program](https://github.com/deepmodeling/abacus-develop/blob/develop/source/module_hamilt_general/module_vdw/vdwd2_parameters.cpp). To manually set the $C_6$ parameters, provide a file containing the parameters. An example is given by:

  ```text
  H  0.1
  Si 9.0
  ```

  Namely, each line contains the element name and the corresponding $C_6$ parameter.
- **Default**: default

### vdw_C6_unit

- **Type**: String
- **Availability**: `vdw_C6_file` is not default
- **Description**: Specifies the unit of the provided $C_6$ parameters in the D2 method. Available options are:
  - `Jnm6/mol` (J·nm^6/mol)
  - `eVA` (eV·Angstrom)
- **Default**: Jnm6/mol

### vdw_R0_file

- **Type**: String
- **Availability**: `vdw_method` is set to `d2`
- **Description**: Specifies the name of the file containing $R_0$ parameters for each element when using the D2 method. If not set, ABACUS uses the default $R_0$ parameters (Angstrom) stored in the [program](https://github.com/deepmodeling/abacus-develop/blob/develop/source/module_hamilt_general/module_vdw/vdwd2_parameters.cpp). To manually set the $R_0$ parameters, provide a file containing the parameters. An example is given by:

  ```text
  Li 1.0
  Cl 2.0
  ```

  Namely, each line contains the element name and the corresponding $R_0$ parameter.
- **Default**: default

### vdw_R0_unit

- **Type**: String
- **Availability**: `vdw_R0_file` is not default
- **Description**: Specifies the unit for the $R_0$ parameters in the D2 method when manually set by the user. Available options are:
  - `A` (Angstrom)
  - `Bohr`
- **Default**: A

### vdw_cutoff_type

- **Type**: String
- **Description**: Determines the method used for specifying the cutoff radius in periodic systems when applying Van der Waals correction. Available options are:
  - `radius`: The supercell is selected within a sphere centered at the origin with a radius defined by `vdw_cutoff_radius`.
  - `period`: The extent of the supercell is explicitly specified using the `vdw_cutoff_period` keyword.
- **Default**: radius

### vdw_cutoff_radius

- **Type**: Real
- **Availability**: `vdw_cutoff_type` is set to `radius`
- **Description**: Defines the radius of the cutoff sphere when `vdw_cutoff_type` is set to `radius`. The default values depend on the chosen `vdw_method`.
- **Default**:
  - 56.6918 if `vdw_method` is set to `d2`
  - 95 if `vdw_method` is set to `d3_0` or `d3_bj`
- **Unit**: defined by `vdw_radius_unit` (default `Bohr`)

### vdw_radius_unit

- **Type**: String
- **Availability**: `vdw_cutoff_type` is set to `radius`
- **Description**: specify the unit of `vdw_cutoff_radius`. Available options are:
  - `A`(Angstrom)
  - `Bohr`
- **Default**: Bohr

### vdw_cutoff_period

- **Type**: Integer Integer Integer
- **Availability**: `vdw_cutoff_type` is set to `period`
- **Description**: The three integers supplied here explicitly specify the extent of the supercell in the directions of the three basis lattice vectors.
- **Default**: 3 3 3

### vdw_cn_thr

- **Type**: Real
- **Availability**: `vdw_method` is set to `d3_0` or `d3_bj`
- **Description**: The cutoff radius when calculating coordination numbers.
- **Default**: 40
- **Unit**: defined by `vdw_cn_thr_unit` (default: `Bohr`)

### vdw_cn_thr_unit

- **Type**: String
- **Description**: Unit of the coordination number cutoff (`vdw_cn_thr`). Available options are:
  - `A`(Angstrom)
  - `Bohr`
- **Default**: Bohr

[back to top](#full-list-of-input-keywords)

## Berry phase and wannier90 interface

These variables are used to control berry phase and wannier90 interface parameters. [Detail introduce](https://github.com/deepmodeling/abacus-develop/blob/develop/docs/advanced/interface/Wannier90.md#wannier90)

### berry_phase

- **Type**: Boolean
- **Description**: controls the calculation of Berry phase
  - true: Calculate Berry phase.
  - false: Do not calculate Berry phase.
- **Default**: false

### gdir

- **Type**: Integer
- **Description**: the direction of the polarization in the lattice vector for Berry phase calculation
  - 1: Calculate the polarization in the direction of the lattice vector a_1 defined in the STRU file.
  - 2: Calculate the polarization in the direction of the lattice vector a_2 defined in the STRU file.
  - 3: Calculate the polarization in the direction of the lattice vector a_3 defined in the STRU file.
- **Default**: 3

### towannier90

- **Type**: Integer
- **Description**: Controls the generation of files for the Wannier90 code.
  - 1: Generate files for the Wannier90 code.
  - 0: Do not generate files for the Wannier90 code.
- **Default**: 0

### nnkpfile

- **Type**: String
- **Description**: the file name generated when running "wannier90 -pp ..." command
- **Default**: seedname.nnkp

### wannier_method

- **Type**: Integer
- **Description**: Only available on LCAO basis, using different methods to generate "\*.mmn" file and "\*.amn" file.
  - 1: Calculated using the `lcao_in_pw` method, the calculation accuracy can be improved by increasing `ecutwfc` to maintain consistency with the pw basis set results.
  - 2: The overlap between atomic orbitals is calculated using grid integration. The radial grid points are generated using the Gauss-Legendre method, while the spherical grid points are generated using the Lebedev-Laikov method.
- **Default**: 1

### wannier_spin

- **Type**: String
- **Description**: the spin direction for the Wannier function calculation when nspin is set to 2
  - `up`: Calculate spin up for the Wannier function.
  - `down`: Calculate spin down for the Wannier function.
- **Default**: `up`

### out_wannier_mmn

- **Type**: Bool
- **Description**: write the "*.mmn" file or not.
  - 0: don't write the "*.mmn" file.
  - 1: write the "*.mmn" file.
- **Default**: 1

### out_wannier_amn

- **Type**: Bool
- **Description**: write the "*.amn" file or not.
  - 0: don't write the "*.amn" file.
  - 1: write the "*.amn" file.
- **Default**: 1

### out_wannier_eig

- **Type**: Bool
- **Description**: write the "*.eig" file or not.
  - 0: don't write the "*.eig" file.
  - 1: write the "*.eig" file.
- **Default**: 1

### out_wannier_unk

- **Type**: Bool
- **Description**: write the "UNK.*" file or not.
  - 0: don't write the "UNK.*" file.
  - 1: write the "UNK.*" file.
- **Default**: 0

### out_wannier_wvfn_formatted

- **Type**: Bool
- **Description**: write the "UNK.*" file in ASCII format or binary format.
  - 0: write the "UNK.*" file in binary format.
  - 1: write the "UNK.*" file in ASCII format (text file format).
- **Default**: 1

[back to top](#full-list-of-input-keywords)

## TDDFT: time dependent density functional theory

### td_edm

- **Type**: Integer
- **Description**: the method to calculate the energy density matrix
  - 0: new method (use the original formula).
  - 1: old method (use the formula for ground state).
- **Default**: 0

### td_print_eij

- **Type**: Real
- **Description**:
  - <0: don't print $E_{ij}$.
  - \>=0: print the $E_{ij}\ (<\psi_i|H|\psi_j>$) elements which are larger than td_print_eij.
- **Default**: -1

### td_propagator

- **Type**: Integer
- **Description**:
  method of propagator
  - 0: Crank-Nicolson.
  - 1: 4th Taylor expansions of exponential.
  - 2: enforced time-reversal symmetry (ETRS).
- **Default**: 0

### td_vext

- **Type**: Boolean
- **Description**:
  - True: add a laser material interaction (extern laser field).
  - False: no extern laser field.
- **Default**: False

### td_vext_dire

- **Type**: String
- **Description**:
  If `td_vext` is True, the td_vext_dire is a string to set the number of electric fields, like `td_vext_dire 1 2` representing external electric field is added to the x and y axis at the same time. Parameters of electric field can also be written as a string, like `td_gauss_phase 0 1.5707963267948966` representing the Gauss field in the x and y directions has a phase delay of Pi/2. See below for more parameters of electric field.
  - 1: the direction of external light field is along x axis.
  - 2: the direction of external light field is along y axis.
  - 3: the direction of external light field is along z axis.
- **Default**: 1

### td_stype

- **Type**: Integer
- **Description**:
  type of electric field in space domain
  - 0: length gauge.
  - 1: velocity gauge.
- **Default**: 0

### td_ttype

- **Type**: Integer
- **Description**:
  type of electric field in time domain
  - 0: Gaussian type function.
  - 1: Trapezoid function.
  - 2: Trigonometric function.
  - 3: Heaviside function.
  - 4: HHG function.
- **Default**: 0

### td_tstart

- **Type**: Integer
- **Description**:
  number of steps where electric field starts
- **Default**: 1

### td_tend

- **Type**: Integer
- **Description**:
  number of steps where electric field ends
- **Default**: 100

### td_lcut1

- **Type**: Real
- **Description**:
  cut1 of interval in length gauge\
  E = E0 , cut1<x<cut2\
  E = -E0/(cut1+1-cut2) , x<cut1 or cut2<x<1
- **Default**: 0.05

### td_lcut2

- **Type**: Real
- **Description**:
  cut2 of interval in length gauge\
  E = E0 , cut1<x<cut2\
  E = -E0/(cut1+1-cut2) , x<cut1 or cut2<x<1
- **Default**: 0.05

### td_gauss_freq

- **Type**: Real
- **Description**:
  frequency (freq) of Gauss type electric field  (fs^-1)\
  amp\*cos(2pi\*freq(t-t0)+phase)exp(-(t-t0)^2/2sigma^2)
- **Default**: 22.13

### td_gauss_phase

- **Type**: Real
- **Description**:
  phase of Gauss type electric field\
  amp\*(2pi\*freq(t-t0)+phase)exp(-(t-t0)^2/2sigma^2)
- **Default**: 0.0

### td_gauss_sigma

- **Type**: Real
- **Description**:
  sigma of Gauss type electric field  (fs)\
  amp\*cos(2pi\*freq(t-t0)+phase)exp(-(t-t0)^2/2sigma^2)
- **Default**: 30.0

### td_gauss_t0

- **Type**: Real
- **Description**:
  step number of time center (t0) of Gauss type electric field\
  amp\*cos(2pi\*freq(t-t0)+phase)exp(-(t-t0)^2/2sigma^2)
- **Default**: 100

### td_gauss_amp

- **Type**: Real
- **Description**:
  amplitude (amp) of Gauss type electric field  (V/Angstrom)\
  amp\*cos(2pi\*freq(t-t0)+phase)exp(-(t-t0)^2/2sigma^2)
- **Default**: 0.25

### td_trape_freq

- **Type**: Real
- **Description**:
  frequency (freq) of Trapezoid type electric field  (fs^-1)\
  E = amp\*cos(2pi\*freq\*t+phase) t/t1 , t<t1\
  E = amp\*cos(2pi\*freq\*t+phase) , t1<t<t2\
  E = amp\*cos(2pi\*freq\*t+phase) (1-(t-t2)/(t3-t2)) , t2<t<t3\
  E = 0 , t>t3
- **Default**: 1.60

### td_trape_phase

- **Type**: Real
- **Description**:
  phase of Trapezoid type electric field\
  E = amp\*cos(2pi\*freq\*t+phase) t/t1 , t<t1\
  E = amp\*cos(2pi\*freq\*t+phase) , t1<t<t2\
  E = amp\*cos(2pi\*freq\*t+phase) (1-(t-t2)/(t3-t2)) , t2<t<t3\
  E = 0 , t>t3
- **Default**: 0.0

### td_trape_t1

- **Type**: Real
- **Description**:
  step number of time interval 1 (t1) of Trapezoid type electric field\
  E = amp\*cos(2pi\*freq\*t+phase) t/t1 , t<t1\
  E = amp\*cos(2pi\*freq\*t+phase) , t1<t<t2\
  E = amp\*cos(2pi\*freq\*t+phase) (1-(t-t2)/(t3-t2)) , t2<t<t3\
  E = 0 , t>t3
- **Default**: 1875

### td_trape_t2

- **Type**: Real
- **Description**:
  step number of time interval 2 (t2) of Trapezoid type electric field\
  E = amp\*cos(2pi\*freq\*t+phase) t/t1 , t<t1\
  E = amp\*cos(2pi\*freq\*t+phase) , t1<t<t2\
  E = amp\*cos(2pi\*freq\*t+phase) (1-(t-t2)/(t3-t2)) , t2<t<t3\
  E = 0 , t>t3
- **Default**: 5625

### td_trape_t3

- **Type**: Real
- **Description**:
  step number of time interval 3 (t3) of Trapezoid type electric field\
  E = amp\*cos(2pi\*freq\*t+phase) t/t1 , t<t1\
  E = amp\*cos(2pi\*freq\*t+phase) , t1<t<t2\
  E = amp\*cos(2pi\*freq\*t+phase) (1-(t-t2)/(t3-t2)) , t2<t<t3\
  E = 0 , t>t3
- **Default**: 7500

### td_trape_amp

- **Type**: Real
- **Description**:
  amplitude (amp) of Trapezoid type electric field  (V/Angstrom)\
  E = amp\*cos(2pi\*freq\*t+phase) t/t1 , t<t1\
  E = amp\*cos(2pi\*freq\*t+phase) , t1<t<t2\
  E = amp\*cos(2pi\*freq\*t+phase) (1-(t-t2)/(t3-t2)) , t2<t<t3\
  E = 0 , t>t3
- **Default**: 2.74

### td_trigo_freq1

- **Type**: Real
- **Description**:
  frequency 1 (freq1) of Trigonometric type electric field  (fs^-1)\
  amp\*cos(2\*pi\*freq1\*t+phase1)\*sin(2\*pi\*freq2\*t+phase2)^2
- **Default**: 1.164656

### td_trigo_freq2

- **Type**: Real
- **Description**:
  frequency 2 (freq2) of Trigonometric type electric field  (fs^-1)\
  amp\*cos(2\*pi\*freq1\*t+phase1)\*sin(2\*pi\*freq2\*t+phase2)^2
- **Default**: 0.029116

### td_trigo_phase1

- **Type**:Real
- **Description**:
  phase 1 (phase1) of Trigonometric type electric field\
  amp\*cos(2\*pi\*freq1\*t+phase1)\*sin(2\*pi\*freq2\*t+phase2)^2
- **Default**: 0.0

### td_trigo_phase2

- **Type**: Real
- **Description**:
  phase 2 (phase2) of Trigonometric type electric field\
  amp\*cos(2\*pi\*freq1\*t+phase1)\*sin(2\*pi\*freq2\*t+phase2)^2
- **Default**: 0.0

### td_trigo_amp

- **Type**: Real
- **Description**:
  amplitude (amp) of Trigonometric type electric field (V/Angstrom)\
  amp\*cos(2\*pi\*freq1\*t+phase1)\*sin(2\*pi\*freq2\*t+phase2)^2
- **Default**: 2.74

### td_heavi_t0

- **Type**: Real
- **Description**:
  step number of switch time (t0) of Heaviside type electric field\
  E = amp , t<t0\
  E = 0.0 , t>t0
- **Default**: 100

### td_heavi_amp

- **Type**: Real
- **Description**:
  amplitude (amp) of Heaviside type electric field  (V/Angstrom)\
  E = amp , t<t0\
  E = 0.0 , t>t0
- **Default**: 2.74

### out_dipole

- **Type**: Boolean
- **Description**:
  - True: output dipole.
  - False: do not output dipole.
- **Default**: False

### out_current

- **Type**: Boolean
- **Description**:output current in real time TDDFT simulations with the velocity gauge
  - True: output current.
  - False: do not output current.
- **Default**: False

### out_current_k

- **Type**: Boolean
- **Description**:output tddft current for all k points.
  - True: output tddft current for all k points.
  - False: output current in total.
- **Default**: False

### out_efield

- **Type**: Boolean
- **Description**: output TDDFT Efield or not(V/Angstrom)
  - True: output efield.
  - False: do not output efield.
- **Default**: False

### out_vecpot

- **Type**: Boolean
- **Description**: output TDDFT Vector potential or not(a.u.)
  - True: output Vector potential in file "OUT.suffix/At.dat"
  - False: do not output Vector potential.
- **Default**: False

### init_vecpot_file

- **Type**: Boolean
- **Description**: Init vector potential through file or not
  - True: init vector potential from file "At.dat".(a.u.) It consists of four columns, representing istep and vector potential on each direction.
  - False: calculate vector potential by integral of Efield
- **Default**: False

### ocp

- **Type**: Boolean
- **Availability**:
  - For PW and LCAO codes. if set to 1, occupations of bands will be setting of "ocp_set".
  - For TDDFT in LCAO codes. if set to 1, occupations will be constrained since second ionic step.
  - For OFDFT, this feature can't be used.
- **Description**:
- True: fix the occupations of bands.
- False: do not fix the occupations of bands.
- **Default**: False

### ocp_set

- **Type**: String
- **Description**: If ocp is True, the ocp_set is a string to set the number of occupancy, like '1 10 * 1 0 1' representing the 13 band occupancy, 12th band occupancy 0 and the rest 1, the code is parsing this string into an array through a regular expression.
- **Default**: none

[back to top](#full-list-of-input-keywords)

## Variables useful for debugging

### t_in_h

- **Type**: Boolean
- **Description**: Specify whether to include kinetic term in obtaining the Hamiltonian matrix.
  - 0: No.
  - 1: Yes.
- **Default**: 1

### vl_in_h

- **Type**: Boolean
- **Description**: Specify whether to include local pseudopotential term in obtaining the Hamiltonian matrix.
  - 0: No.
  - 1: Yes.
- **Default**: 1

### vnl_in_h

- **Type**: Boolean
- **Description**: Specify whether to include non-local pseudopotential term in obtaining the Hamiltonian matrix.
  - 0: No.
  - 1: Yes.
- **Default**: 1

### vh_in_h

- **Type**: Boolean
- **Description**: Specify whether to include Hartree potential term in obtaining the Hamiltonian matrix.
  - 0: No.
  - 1: Yes.
- **Default**: 1

### vion_in_h

- **Type**: Boolean
- **Description**: Specify whether to include local ionic potential term in obtaining the Hamiltonian matrix.
  - 0: No.
  - 1: Yes.
- **Default**: 1

### test_force

- **Type**: Boolean
- **Description**: Specify whether to output the detailed components in forces.
  - 0: No.
  - 1: Yes.
- **Default**: 0

### test_stress

- **Type**: Boolean
- **Description**: Specify whether to output the detailed components in stress.
  - 0: No.
  - 1: Yes.
- **Default**: 0

- **Type**: Boolean
- **Description**: Specify whether to set the colorful output in terminal.
  - 0: No.
  - 1: Yes.
- **Default**: 0

### test_skip_ewald

- **Type**: Boolean
- **Description**: Specify whether to skip the calculation of the ewald energy.
  - 0: No.
  - 1: Yes.
- **Default**: 0

[back to top](#full-list-of-input-keywords)

## Electronic conductivities

Frequency-dependent electronic conductivities can be calculated with Kubo-Greenwood formula [Phys. Rev. B 83, 235120 (2011)].

Onsager coefficients:

$L_{mn}(\omega)=(-1)^{m+n}\frac{2\pi e^2\hbar^2}{3m_e^2\omega\Omega}$

$\times\sum_{ij\alpha\mathbf{k}}W(\mathbf{k})\left(\frac{\epsilon_{i\mathbf{k}}+\epsilon_{j\mathbf{k}}}{2}-\mu\right)^{m+n-2} \times |\langle\Psi_{i\mathbf{k}}|\nabla_\alpha|\Psi_{j\mathbf{k}}\rangle|^2$

$\times[f(\epsilon_{i\mathbf{k}})-f(\epsilon_{j\mathbf{k}})]\delta(\epsilon_{j\mathbf{k}}-\epsilon_{i\mathbf{k}}-\hbar\omega).$

They can also be computed by $j$-$j$ correlation function.

$L_{mn}=\frac{2e^{m+n-2}}{3\Omega\hbar\omega}\Im[\tilde{C}_{mn}(\omega)]$
Guassian smearing:
$\tilde{C}_{mn}=\int_0^\infty C_{mn}(t)e^{-i\omega t}e^{-\frac{1}{2}s^2t^2}dt$
Lorentzian smearing:
$\tilde{C}_{mn}=\int_0^\infty C_{mn}(t)e^{-i\omega t}e^{-\gamma t}dt$

$C_{mn}(t)=-2\theta(t)\Im\left\{Tr\left[\sqrt{\hat f}\hat{j}_m(1-\hat{f})e^{i\frac{\hat{H}}{\hbar}t}\hat{j}_ne^{-i\frac{\hat{H}}{\hbar}t}\sqrt{\hat f}\right]\right\}$,

where $j_1$ is electric flux and $j_2$ is thermal flux.

Frequency-dependent electric conductivities: $\sigma(\omega)=L_{11}(\omega)$.

Frequency-dependent thermal conductivities: $\kappa(\omega)=\frac{1}{e^2T}\left(L_{22}-\frac{L_{12}^2}{L_{11}}\right)$.

DC electric conductivities: $\sigma = \lim_{\omega\to 0}\sigma(\omega)$.

Thermal conductivities: $\kappa = \lim_{\omega\to 0}\kappa(\omega)$.

### cal_cond

- **Type**: Boolean
- **Availability**: [basis_type](#basis_type) = `pw`
- **Description**: Whether to calculate electronic conductivities.
- **Default**: False

### cond_che_thr

- **Type**: Real
- **Availability**: [esolver_type](#esolver_type) = `sdft`
- **Description**: Control the error of Chebyshev expansions for conductivities.
- **Default**: 1e-8

### cond_dw

- **Type**: Real
- **Availability**: [basis_type](#basis_type) = `pw`
- **Description**: Frequency interval ($\mathrm{d}\omega$) for frequency-dependent conductivities.
- **Default**: 0.1
- **Unit**: eV

### cond_wcut

- **Type**: Real
- **Availability**: [basis_type](#basis_type) = `pw`
- **Description**: Cutoff frequency for frequency-dependent conductivities.
- **Default**: 10.0
- **Unit**: eV

### cond_dt

- **Type**: Real
- **Availability**: [basis_type](#basis_type) = `pw`
- **Description**: Time interval ($\mathrm{d}t$) to integrate Onsager coefficients.
- **Default**: 0.02
- **Unit**: a.u.

### cond_dtbatch

- **Type**: Integer
- **Availability**: [esolver_type](#esolver_type) = `sdft`
- **Description**: exp(iH\*dt\*cond_dtbatch) is expanded with Chebyshev expansion to calculate conductivities. It is faster but costs more memory.
  - If `cond_dtbatch = 0`: Autoset this parameter to make expansion orders larger than 100.
- **Default**: 0

### cond_smear

- **Type**: Integer
- **Description**: Smearing method for conductivities
  - 1: Gaussian smearing
  - 2: Lorentzian smearing
- **Default**: 1

### cond_fwhm

- **Type**: Real
- **Availability**: [basis_type](#basis_type) = `pw`
- **Description**: FWHM for conductivities. For Gaussian smearing, $\mathrm{FWHM}=2\sqrt{2\ln2}s$; for Lorentzian smearing, $\mathrm{FWHM}=2\gamma$.
- **Default**: 0.4
- **Unit**: eV

### cond_nonlocal

- **Type**: Boolean
- **Availability**: [basis_type](#basis_type) = `pw`
- **Description**: Whether to consider nonlocal potential correction when calculating velocity matrix $\bra{\psi_i}\hat{v}\ket{\psi_j}$.
  - True:  $m\hat{v}=\hat{p}+\frac{im}{\hbar}[\hat{V}_{NL},\hat{r}]$.
  - False: $m\hat{v}\approx\hat{p}$.
- **Default**: True

[back to top](#full-list-of-input-keywords)

## Implicit solvation model

These variables are used to control the usage of implicit solvation model. This approach treats the solvent as a continuous medium instead of individual "explicit" solvent molecules, which means that the solute is embedded in an implicit solvent and the average over the solvent degrees of freedom becomes implicit in the properties of the solvent bath.

### imp_sol

- **Type**: Boolean
- **Description**: calculate implicit solvation correction
- **Default**: False

### eb_k

- **Type**: Real
- **Availability**: `imp_sol` is true.
- **Description**: the relative permittivity of the bulk solvent, 80 for water
- **Default**: 80

### tau

- **Type**: Real
- **Description**: The effective surface tension parameter that describes the cavitation, the dispersion, and the repulsion interaction between the solute and the solvent which are not captured by the electrostatic terms
- **Default**: 1.0798e-05
- **Unit**: $Ry/Bohr^{2}$

### sigma_k

- **Type**: Real
- **Description**: the width of the diffuse cavity that is implicitly determined by the electronic structure of the solute
- **Default**: 0.6

### nc_k

- **Type**: Real
- **Description**: the value of the electron density at which the dielectric cavity forms
- **Default**: 0.00037
- **Unit**: $Bohr^{-3}$

[back to top](#full-list-of-input-keywords)

## Quasiatomic Orbital (QO) analysis

These variables are used to control the usage of QO analysis. QO further compress information from LCAO: usually PW basis has dimension in million, LCAO basis has dimension below thousand, and QO basis has dimension below hundred.

### qo_switch

- **Type**: Boolean
- **Description**: whether to let ABACUS output QO analysis required files
- **Default**: 0

### qo_basis

- **Type**: String
- **Description**: specify the type of atomic basis
  - `pswfc`: use the pseudowavefunction in pseudopotential files as atomic basis. To use this option, please make sure in pseudopotential file there is pswfc in it.
  - `hydrogen`: generate hydrogen-like atomic basis (or with Slater screening).
  - `szv`: use the first set of zeta for each angular momentum from numerical atomic orbitals as atomic basis.

  *warning: to use* `pswfc` *, please use norm-conserving pseudopotentials with pseudowavefunctions, SG15 pseudopotentials cannot support this option.*
  *Developer notes: for ABACUS-lcao calculation, it is the most recommend to use `szv` instead of `pswfc` which is originally put forward in work of QO implementation on PW basis. The information loss always happens if `pswfc` or `hydrogen` orbitals are not well tuned, although making kpoints sampling more dense will mitigate this problem, but orbital-adjust parameters are needed to test system-by-system in this case.*
- **Default**: `szv`

### qo_strategy

- **Type**: String \[String...\](optional)
- **Description**: specify the strategy to generate radial orbitals for each atom type. If one parameter is given, will apply to all atom types. If more than one parameters are given but fewer than number of atom type, those unspecified atom type will use default value.

  For `qo_basis hydrogen`
  - `minimal-nodeless`: according to principle quantum number of the highest occupied state, generate only nodeless orbitals, for example Cu, only generate 1s, 2p, 3d and 4f orbitals (for Cu, 4s is occupied, thus $n_{max} = 4$)
  - `minimal-valence`: according to principle quantum number of the highest occupied state, generate only orbitals with highest principle quantum number, for example Cu, only generate 4s, 4p, 4d and 4f orbitals.
  - `full`: similarly according to the maximal principle quantum number, generate all possible orbitals, therefore for Cu, for example, will generate 1s, 2s, 2p, 3s, 3p, 3d, 4s, 4p, 4d, 4f.
  - `energy-full`: will generate hydrogen-like orbitals according to Aufbau principle. For example the Cu (1s2 2s2 2p6 3s2 3p6 3d10 4s1), will generate these orbitals.
  - `energy-valence`: from the highest n (principal quantum number) layer and n-1 layer, generate all occupied and possible ls (angular momentum quantum number) for only once, for example Cu, will generate 4s, 3d and 3p orbitals.

  For `qo_basis pswfc` and `qo_basis szv`
  - `all`: use all possible pseudowavefunctions/numerical atomic orbital (of first zeta) in pseudopotential/numerical atomic orbital file.
  - `s`/`p`/`d`/...: only use s/p/d/f/...-orbital(s).
  - `spd`: use s, p and d orbital(s). Any unordered combination is acceptable.

  *warning: for* `qo_basis hydrogen` *to use* `full`, *generation strategy may cause the space spanned larger than the one spanned by numerical atomic orbitals, in this case, must filter out orbitals in some way*
- **Default**: for `hydrogen`: `energy-valence`, for `pswfc` and `szv`: `all`

### qo_screening_coeff

- **Type**: Real \[Real...\](optional)
- **Description**: rescale the shape of radial orbitals, available for both `qo_basis hydrogen` and `qo_basis pswfc`. cases but has different meaning.

  For `qo_basis pswfc`
  For each atom type, screening factor $e^{-\eta|\mathbf{r}|}$ is multiplied to the pswfc to mimic the behavior of some kind of electron. $\eta$ is the screening coefficient. If only one value is given, then will apply to each atom type. If not enough values are given, will apply default value to rest of atom types. This parameter plays important role in controlling the spread of QO orbitals together with `qo_thr`.

  For `qo_basis hydrogen`
  If any float number is given, will apply Slater screening to all atom types. Slater screening is a classic and empirical method roughly taking many-electron effect into account for obtaining more accurate results when evaluating electron affinity and ionization energy. The Coulomb potential then becomes $V(r) = -\frac{Z-\sigma}{r}$. For example the effective nuclear charge for Cu 3d electrons now reduces from 29 to 7.85, 4s from 29 to 3.70, which means Slater screening will bring about longer tailing effect. If no value is given, will not apply Slater screening.
- **Default**: 0.1
- **Unit**: Bohr^-1

### qo_thr

- **Type**: Real
- **Description**: the convergence threshold determining the cutoff of generated orbital. Lower threshold will yield orbital with larger cutoff radius.
- **Default**: 1.0e-6

## PEXSI

These variables are used to control the usage of PEXSI (Pole Expansion and Selected Inversion) method in calculations.

### pexsi_npole

- **Type**: Integer
- **Description**: the number of poles used in the pole expansion method, should be a even number.
- **Default**: 40

### pexsi_inertia

- **Type**: Boolean
- **Description**: whether inertia counting is used at the very beginning.
- **Default**: True

### pexsi_nmax

- **Type**: Integer
- **Description**: maximum number of PEXSI iterations after each inertia counting procedure.
- **Default**: 80

### pexsi_comm

- **Type**: Boolean
- **Description**: whether to construct PSelInv communication pattern.
- **Default**: True

### pexsi_storage

- **Type**: Boolean
- **Description**: whether to use symmetric storage space used by the Selected Inversion algorithm for symmetric matrices.
- **Default**: True

### pexsi_ordering

- **Type**: Integer
- **Description**: ordering strategy for factorization and selected inversion. 0: Parallel ordering using ParMETIS, 1: Sequential ordering using METIS, 2: Multiple minimum degree ordering
- **Default**: 0

### pexsi_row_ordering

- **Type**: Integer
- **Description**: row permutation strategy for factorization and selected inversion, 0: No row permutation, 1: Make the diagonal entry of the matrix larger than the off-diagonal entries.
- **Default**: 1

### pexsi_nproc

- **Type**: Integer
- **Description**: number of processors for PARMETIS. Only used if pexsi_ordering == 0.
- **Default**: 1

### pexsi_symm

- **Type**: Boolean
- **Description**: whether the matrix is symmetric.
- **Default**: True

### pexsi_trans

- **Type**: Boolean
- **Description**: whether to factorize the transpose of the matrix.
- **Default**: False

### pexsi_method

- **Type**: Integer
- **Description**: the pole expansion method to be used. 1 for Cauchy Contour Integral method, 2 for Moussa optimized method.
- **Default**: 1

### pexsi_nproc_pole

- **Type**: Integer
- **Description**: the point parallelizaion of PEXSI. Recommend two points parallelization.
- **Default**: 1

### pexsi_temp

- **Type**: Real
- **Description**: temperature in Fermi-Dirac distribution, in Ry, should have the same effect as the smearing sigma when smearing method is set to Fermi-Dirac.
- **Default**: 0.015

### pexsi_gap

- **Type**: Real
- **Description**: spectral gap, this can be set to be 0 in most cases.
- **Default**: 0

### pexsi_delta_e

- **Type**: Real
- **Description**: an upper bound for the spectral radius of $S^{-1} H$.
- **Default**: 20

### pexsi_mu_lower

- **Type**: Real
- **Description**: initial guess of lower bound for mu.
- **Default**: -10

### pexsi_mu_upper

- **Type**: Real
- **Description**: initial guess of upper bound for mu.
- **Default**: 10

### pexsi_mu

- **Type**: Real
- **Description**: initial guess for mu (for the solver).
- **Default**: 0

### pexsi_mu_thr

- **Type**: Real
- **Description**: stopping criterion in terms of the chemical potential for the inertia counting procedure.
- **Default**: 0.05

### pexsi_mu_expand

- **Type**: Real
- **Description**: if the chemical potential is not in the initial interval, the interval is expanded by this value.
- **Default**: 0.3

### pexsi_mu_guard

- **Type**: Real
- **Description**: safe guard criterion in terms of the chemical potential to reinvoke the inertia counting procedure.
- **Default**: 0.2

### pexsi_elec_thr

- **Type**: Real
- **Description**: stopping criterion of the PEXSI iteration in terms of the number of electrons compared to numElectronExact.
- **Default**: 0.001

### pexsi_zero_thr

- **Type**: Real
- **Description**: if the absolute value of CCS matrix element is less than this value, it will be considered as zero.
- **Default**: 1e-10

[back to top](#full-list-of-input-keywords)

## Linear Response TDDFT

These parameters are used to solve the excited states using. e.g. LR-TDDFT.

### xc_kernel

- **Type**: String
- **Description**: The exchange-correlation kernel used in the calculation. 
Currently supported: `RPA`, `LDA`, `PBE`, `HSE`, `HF`.
- **Default**: LDA

### lr_init_xc_kernel

- **Type**: String
- **Description**: The method to initalize the xc kernel. 
  - "default": Calculate xc kerenel ($f_\text{xc}$) from the ground-state charge density.
  - "file": Read the xc kernel $f_\text{xc}$ on grid from the provided files. The following words should be the paths of ".cube" files, where the first 1 (*[nspin](#nspin)==1*) or 3 (*[nspin](#nspin)==2*, namely spin-aa, spin-ab and spin-bb) will be read in. The parameter [xc_kernel](#xc_kernel) will be invalid. Now only LDA-type kernel is supproted as the potential will be calculated by directly multiplying the transition density.
  - "from_charge_file": Calculate fxc from the charge density read from the provided files. The following words should be the paths of ".cube" files, where the first [nspin]($nspin) files will be read in. 
- **Default**: "default"

### lr_solver

- **Type**: String
- **Description**: The method to solve the Casida equation $AX=\Omega X$ in LR-TDDFT under Tamm-Dancoff approximation (TDA), where $A_{ai,bj}=(\epsilon_a-\epsilon_i)\delta_{ij}\delta_{ab}+(ai|f_{Hxc}|bj)+\alpha_{EX}(ab|ij)$ is the particle-hole excitation matrix and $X$ is the transition amplitude.
  - `dav`/`dav_subspace`/ `cg`: Construct $AX$ and diagonalize the Hamiltonian matrix iteratively with Davidson/Non-ortho-Davidson/CG algorithm.
  - `lapack`: Construct the full $A$ matrix and directly diagonalize with LAPACK.
  - `spectrum`: Calculate absorption spectrum only without solving Casida equation. The `OUT.${suffix}/` directory should contain the
  files for LR-TDDFT eigenstates and eigenvalues, i.e. `Excitation_Energy.dat` and `Excitation_Amplitude_${processor_rank}.dat`
   output by setting `out_wfc_lr` to true.
- **Default**: dav

### lr_thr

- **Type**: Real
- **Description**: The convergence threshold of iterative diagonalization solver fo LR-TDDFT. It is a pure-math number with the same as [pw_diag_thr](#pw_diag_thr), but since the Casida equation is a one-shot eigenvalue problem, it is also the convergence threshold of LR-TDDFT.
- **Default**: 1e-2

### nocc

- **Type**: Integer
- **Description**: The number of occupied orbitals (up to HOMO) used in the LR-TDDFT calculation. 
  - Note: If the value is illegal ( > [nelec](#nelec)\/2 or <= 0), it will be autoset to [nelec](#nelec)\/2.
- **Default**: [nband](#nband)

### nvirt

- **Type**: Integer
- **Description**: The number of virtual orbitals (staring from LUMO) used in the LR-TDDFT calculation.
- **Default**: 1

### lr_nstates

- **Type**: Integer
- **Description**:  The number of 2-particle states to be solved
- **Default**: 0

### lr_unrestricted
- **Type**: Boolean
- **Description**: Whether to use unrestricted construction for LR-TDDFT (the matrix size will be doubled).
  - True:  Always use unrestricted LR-TDDFT. 
  - False: Use unrestricted LR-TDDFT only when the system is open-shell.
- **Default**: False

### abs_wavelen_range

- **Type**: Real Real
- **Description**: The range of the wavelength for the absorption spectrum calculation.
- **Default**: 0.0 0.0

### out_wfc_lr

- **Type**: Boolean
- **Description**: Whether to output the eigenstates (excitation energy) and eigenvectors (excitation amplitude) of the LR-TDDFT calculation.
The output files are `OUT.${suffix}/Excitation_Energy.dat` and `OUT.${suffix}/Excitation_Amplitude_${processor_rank}.dat`.
- **Default**: False

### abs_broadening
- **Type**: Real
- **Description**: The broadening factor $\eta$ for the absorption spectrum calculation.
- **Default**: 0.01

### ri_hartree_benchmark
- **Type**: String
- **Description**: Whether to use the localized resolution-of-identity (LRI) approximation for the **Hartree** term of kernel in the $A$ matrix of LR-TDDFT for benchmark (with FHI-aims or another ABACUS calculation). Now it only supports molecular systems running with a single processor, and a large enough supercell should be used to make LRI C, V tensors contain only the R=(0 0 0) cell. 
  - `aims`: The `OUT.${suffix}`directory should contain the FHI-aims output files: RI-LVL tensors`Cs_data_0.txt` and `coulomb_mat_0.txt`, and KS eigenstates from FHI-aims: `band_out`and `KS_eigenvectors.out`. The Casida equation will be constructed under FHI-aims' KS eigenpairs.
    - LRI tensor files (`Cs_data_0.txt` and `coulomb_mat_0.txt`)and Kohn-Sham eigenvalues (`bands_out`): run FHI-aims with periodic boundary conditions and with `total_energy_method rpa` and `output librpa`.
    - Kohn-Sham eigenstates under aims NAOs (`KS_eigenvectors.out`): run FHI-aims with `output eigenvectors`.
    - If the number of atomic orbitals of any atom type in FHI-aims is different from that in ABACUS, the `aims_nbasis` should be set.
  - `abacus`: The `OUT.${suffix}`directory should contain the RI-LVL tensors `Cs` and `Vs` (written by setting `out_ri_cv` to 1). The Casida equation will be constructed under ABACUS' KS eigenpairs, with the only difference that the Hartree term is constructed with RI approximation.
  - `none`: Construct the Hartree term by Poisson equation and grid integration as usual.
- **Default**: none

### aims_nbasis
- **Type**: A number(ntype) of Integers
- **Availability**: `ri_hartree_benchmark` = `aims`
- **Description**: Atomic basis set size for each atom type (with the same order as in `STRU`) in FHI-aims.
- **Default**: {} (empty list, where ABACUS use its own basis set size)

## Reduced Density Matrix Functional Theory

ab-initio methods and the xc-functional parameters used in RDMFT.
The physical quantities that RDMFT temporarily expects to output are the kinetic energy, total energy, and 1-RDM of the system in the ground state, etc.

### rdmft

- **Type**: Boolean
- **Description**: Whether to perform rdmft calculation (reduced density matrix funcional theory)
- **Default**: false

### rdmft_power_alpha

- **Type**: Real
- **Description**: The alpha parameter of power-functional(or other exx-type/hybrid functionals) which used in RDMFT, g(occ_number) = occ_number^alpha
- **Default**: 0.656

[back to top](#full-list-of-input-keywords)

#ifndef ATOM_H
#define ATOM_H

#include "atom_pseudo.h"
#include "../module_io/output.h"
class Atom
{
public:

    // constructor and destructor
    Atom();
    ~Atom();

    Atom_pseudo ncpp;
    double mass; // the mass of atom
    std::vector<ModuleBase::Vector3<int>> mbl; // whether the atoms can move or not
	bool flag_empty_element = false; // whether is the empty element for bsse.	Peize Lin add 2021.04.07

    std::vector<int> iw2m; // use iw to find m
    std::vector<int> iw2n; // use iw to find n
    std::vector<int> iw2l; // use iw to find L
	std::vector<int> iw2_ylm;
	std::vector<bool> iw2_new;
    int nw; // number of local orbitals (l,n,m) of this type

    void set_index(void);

    int type; // Index of atom type
    int na; // Number of atoms in this type.

    int nwl; // max L(Angular momentum) (for local basis)
    double Rcut; //pengfei Li 16-2-29
    std::vector<int> l_nchi; // number of chi for each L
    int stapos_wf; // start position of wave functions

    std::string label; // atomic symbol
    std::vector<ModuleBase::Vector3<double>> tau;// Cartesian coordinates of each atom in this type.
    std::vector<ModuleBase::Vector3<double>> dis;// direct displacements of each atom in this type in current step  liuyu modift 2023-03-22
	std::vector<ModuleBase::Vector3<double>> taud;// Direct coordinates of each atom in this type.
    std::vector<ModuleBase::Vector3<double>> vel;// velocities of each atom in this type.
    std::vector<ModuleBase::Vector3<double>> force; // force acting on each atom in this type.
    std::vector<ModuleBase::Vector3<double>> lambda; // Lagrange multiplier for each atom in this type. used in deltaspin
    std::vector<ModuleBase::Vector3<int>> constrain; // constrain for each atom in this type. used in deltaspin
    std::string label_orb; // atomic Element symbol in the orbital file of lcao

	std::vector<double> mag;
	std::vector<double> angle1;//spin angle, added by zhengdy-soc
	std::vector<double> angle2;
    std::vector<ModuleBase::Vector3<double>> m_loc_;
    // Coulomb potential v(r) = z/r
    // It is a local potentail, and has no non-local potential parts.
    bool coulomb_potential = false;
    void print_Atom(std::ofstream &ofs);
    void update_force(ModuleBase::matrix &fcs);
#ifdef __MPI
    void bcast_atom(void);
    void bcast_atom2(void);
#endif

};

#endif //Atomspec


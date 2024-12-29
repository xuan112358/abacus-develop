#include "atom_pair.h"
#include <complex>
#include <cassert>
#include "module_base/blas_connector.h"

namespace hamilt
{

//----------------------------------------------------
// atom pair class
//----------------------------------------------------

// destructor
template <typename T>
AtomPair<T>::~AtomPair()
{
    this->values.clear();
}

template <typename T>
AtomPair<T>::AtomPair(const int& atom_i_, const int& atom_j_, const Parallel_Orbitals* paraV_, T* existed_matrix)
    : atom_i(atom_i_), atom_j(atom_j_), paraV(paraV_)
{
    assert(this->paraV != nullptr);
    this->row_ap = this->paraV->atom_begin_row[atom_i];
    this->col_ap = this->paraV->atom_begin_col[atom_j];
    if (this->row_ap == -1 || this->col_ap == -1)
    {
        throw std::string("Atom-pair not belong this process");
    }
    this->row_size = this->paraV->get_row_size(atom_i);
    this->col_size = this->paraV->get_col_size(atom_j);
    this->R_index.resize(0);
    this->R_index.push_back(ModuleBase::Vector3<int>(0, 0, 0));
    this->current_R = 0;
    if (existed_matrix != nullptr)
    {
        BaseMatrix<T> tmp(row_size, col_size, existed_matrix);
        this->values.push_back(tmp);
    }
    else
    {
        BaseMatrix<T> tmp(row_size, col_size);
        this->values.push_back(tmp);
    }
}

template <typename T>
AtomPair<T>::AtomPair(const int& atom_i_,
                      const int& atom_j_,
                      const int& rx,
                      const int& ry,
                      const int& rz,
                      const Parallel_Orbitals* paraV_,
                      T* existed_matrix)
    : atom_i(atom_i_), atom_j(atom_j_), paraV(paraV_)
{
    assert(this->paraV != nullptr);
    this->row_ap = this->paraV->atom_begin_row[atom_i];
    this->col_ap = this->paraV->atom_begin_col[atom_j];
    if (this->row_ap == -1 || this->col_ap == -1)
    {
        throw std::string("Atom-pair not belong this process");
    }
    this->row_size = this->paraV->get_row_size(atom_i);
    this->col_size = this->paraV->get_col_size(atom_j);
    this->R_index.resize(0);
    this->R_index.push_back(ModuleBase::Vector3<int>(rx, ry, rz));
    this->current_R = 0;
    if (existed_matrix != nullptr)
    {
        BaseMatrix<T> tmp(row_size, col_size, existed_matrix);
        this->values.push_back(tmp);
    }
    else
    {
        BaseMatrix<T> tmp(row_size, col_size);
        this->values.push_back(tmp);
    }
}

template <typename T>
AtomPair<T>::AtomPair(const int& atom_i_,
                      const int& atom_j_,
                      const ModuleBase::Vector3<int> &R_index,
                      const Parallel_Orbitals* paraV_,
                      T* existed_matrix)
    : atom_i(atom_i_), atom_j(atom_j_), paraV(paraV_)
{
    assert(this->paraV != nullptr);
    this->row_ap = this->paraV->atom_begin_row[atom_i];
    this->col_ap = this->paraV->atom_begin_col[atom_j];
    if (this->row_ap == -1 || this->col_ap == -1)
    {
        throw std::string("Atom-pair not belong this process");
    }
    this->row_size = this->paraV->get_row_size(atom_i);
    this->col_size = this->paraV->get_col_size(atom_j);
    this->R_index.resize(0);
    this->R_index.push_back(ModuleBase::Vector3<int>(R_index));
    this->current_R = 0;
    if (existed_matrix != nullptr)
    {
        BaseMatrix<T> tmp(row_size, col_size, existed_matrix);
        this->values.push_back(tmp);
    }
    else
    {
        BaseMatrix<T> tmp(row_size, col_size);
        this->values.push_back(tmp);
    }
}

// direct save whole matrix of atom-pair
template <typename T>
AtomPair<T>::AtomPair(const int& atom_i_,
                      const int& atom_j_,
                      const int* row_atom_begin,
                      const int* col_atom_begin,
                      const int& natom,
                      T* existed_matrix)
    : atom_i(atom_i_), atom_j(atom_j_)
{
    assert(row_atom_begin != nullptr && col_atom_begin != nullptr);
    this->row_ap = row_atom_begin[atom_i];
    this->col_ap = col_atom_begin[atom_j];
    this->row_size = row_atom_begin[atom_i + 1] - row_atom_begin[atom_i];
    this->col_size = col_atom_begin[atom_j + 1] - col_atom_begin[atom_j];
    this->R_index.resize(0);
    this->R_index.push_back(ModuleBase::Vector3<int>(0, 0, 0));
    this->current_R = 0;
    if (existed_matrix != nullptr)
    {
        BaseMatrix<T> tmp(row_size, col_size, existed_matrix);
        this->values.push_back(tmp);
    }
    else
    {
        BaseMatrix<T> tmp(row_size, col_size);
        this->values.push_back(tmp);
    }
}
//
template <typename T>
AtomPair<T>::AtomPair(const int& atom_i_,
                      const int& atom_j_,
                      const int& rx,
                      const int& ry,
                      const int& rz,
                      const int* row_atom_begin,
                      const int* col_atom_begin,
                      const int& natom,
                      T* existed_matrix)
    : atom_i(atom_i_), atom_j(atom_j_)
{
    assert(row_atom_begin != nullptr && col_atom_begin != nullptr);
    this->row_ap = row_atom_begin[atom_i];
    this->col_ap = col_atom_begin[atom_j];
    this->row_size = row_atom_begin[atom_i + 1] - row_atom_begin[atom_i];
    this->col_size = col_atom_begin[atom_j + 1] - col_atom_begin[atom_j];
    this->R_index.resize(0);
    this->R_index.push_back(ModuleBase::Vector3<int>(rx, ry, rz));
    this->current_R = 0;
    if (existed_matrix != nullptr)
    {
        BaseMatrix<T> tmp(row_size, col_size, existed_matrix);
        this->values.push_back(tmp);
    }
    else
    {
        BaseMatrix<T> tmp(row_size, col_size);
        this->values.push_back(tmp);
    }
}

template <typename T>
AtomPair<T>::AtomPair(const int& atom_i_, const int& atom_j_) : atom_i(atom_i_), atom_j(atom_j_)
{
}

// copy constructor
template <typename T>
AtomPair<T>::AtomPair(const AtomPair<T>& other, T* data_pointer)
    : R_index(other.R_index),
      paraV(other.paraV),
      current_R(other.current_R),
      atom_i(other.atom_i),
      atom_j(other.atom_j),
      row_ap(other.row_ap),
      col_ap(other.col_ap),
      row_size(other.row_size),
      col_size(other.col_size)
{
    if(data_pointer == nullptr)
    {
        this->values = other.values;
    }
    else
    {
        this->values.reserve(other.values.size());
        for(int value=0;value<other.values.size();++value)
        {
            hamilt::BaseMatrix<T> tmp(row_size, col_size, data_pointer);
            this->values.push_back(tmp);
            data_pointer += this->get_size();
        }
    }
}

//allocate
template <typename T>
void AtomPair<T>::allocate(T* data_array, bool is_zero)
{
    if(data_array == nullptr)
    {
        for(int value=0;value<this->values.size();++value)
        {
            this->values[value].allocate(nullptr, is_zero);
        }
    }
    else
    {
        for(int value=0;value<this->values.size();++value)
        {
            this->values[value].allocate(data_array, is_zero);
            data_array += this->get_size();
        }
    }
}

// set_zero
template <typename T>
void AtomPair<T>::set_zero()
{
    for(auto& value : values)
    {
        value.set_zero();
    }
}

// The copy assignment operator
template <typename T>
AtomPair<T>& AtomPair<T>::operator=(const AtomPair<T>& other)
{
    if (this != &other)
    {
        R_index = other.R_index;
        values = other.values;
        paraV = other.paraV;
        current_R = other.current_R;
        atom_i = other.atom_i;
        atom_j = other.atom_j;
        row_ap = other.row_ap;
        col_ap = other.col_ap;
        row_size = other.row_size;
        col_size = other.col_size;
    }
    return *this;
}

// move constructor
template <typename T>
AtomPair<T>::AtomPair(AtomPair<T>&& other) noexcept
    : R_index(std::move(other.R_index)),
      values(std::move(other.values)),
      paraV(other.paraV),
      current_R(other.current_R),
      atom_i(other.atom_i),
      atom_j(other.atom_j),
      row_ap(other.row_ap),
      col_ap(other.col_ap),
      row_size(other.row_size),
      col_size(other.col_size)
{
    other.paraV = nullptr;
}

// move assignment operator
template <typename T>
AtomPair<T>& AtomPair<T>::operator=(AtomPair<T>&& other) noexcept
{
    if (this != &other)
    {
        R_index = std::move(other.R_index);
        values = std::move(other.values);
        paraV = other.paraV;
        other.paraV = nullptr;
        current_R = other.current_R;
        atom_i = other.atom_i;
        atom_j = other.atom_j;
        row_ap = other.row_ap;
        col_ap = other.col_ap;
        row_size = other.row_size;
        col_size = other.col_size;
    }
    return *this;
}

template <typename T>
bool AtomPair<T>::operator<(const AtomPair<T>& other) const
{
    if (atom_i < other.atom_i)
    {
        return true;
    }
    else if (atom_i == other.atom_i)
    {
        return atom_j < other.atom_j;
    }
    else
    {
        return false;
    }
}

// get col_size
template <typename T>
int AtomPair<T>::get_col_size() const
{
    return this->col_size;
}

// get row_size
template <typename T>
int AtomPair<T>::get_row_size() const
{
    return this->row_size;
}

// get atom_i
template <typename T>
int AtomPair<T>::get_atom_i() const
{
    return this->atom_i;
}

// get atom_j
template <typename T>
int AtomPair<T>::get_atom_j() const
{
    return this->atom_j;
}

// set size
template <typename T>
void AtomPair<T>::set_size(const int& col_size_in, const int& row_size_in)
{
    this->col_size = col_size_in;
    this->row_size = row_size_in;
    for (int i = 0; i < this->values.size(); i++)
    {
        this->values[i].set_size(row_size_in, col_size_in);
    }
}

// get paraV for check
template <typename T>
const Parallel_Orbitals* AtomPair<T>::get_paraV() const
{
    return this->paraV;
}

// identify
template <typename T>
bool AtomPair<T>::identify(const AtomPair<T>& other) const
{
    if (this->atom_i == other.atom_i && this->atom_j == other.atom_j)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// identify
template <typename T>
bool AtomPair<T>::identify(const int& atom_i_, const int& atom_j_) const
{
    if (this->atom_i == atom_i_ && this->atom_j == atom_j_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// get_HR_values for no-const AtomPair
template <typename T>
BaseMatrix<T>& AtomPair<T>::get_HR_values(int rx_in, int ry_in, int rz_in)
{
    // find existed R index
    const int r_index = this->find_R(rx_in, ry_in, rz_in);
    if (r_index != -1)
    {
        // if found, return this->values[current_R]
        return this->values[r_index];
    }
    // if not found, add a new BaseMatrix for this R index
    R_index.push_back(ModuleBase::Vector3<int>(rx_in, ry_in, rz_in));
    values.push_back(BaseMatrix<T>(this->row_size, this->col_size));
    values.back().allocate(nullptr, true);
    // return the last BaseMatrix reference in values
    return this->values.back();
}

// get_HR_values for const AtomPair
template <typename T>
const BaseMatrix<T>& AtomPair<T>::get_HR_values(int rx_in, int ry_in, int rz_in) const
{
    // if current_R is -1, R index has not been fixed, find existed R index
    const int r_index = this->find_R(rx_in, ry_in, rz_in);
    if (r_index != -1)
    {
        // if found, return this->values[current_R]
        return this->values[r_index];
    }
    // if not found, throw a error message
    else
    {
        throw std::string("AtomPair::get_HR_values: R index not found");
    }
}

// get_HR_values with index
template <typename T>
BaseMatrix<T>& AtomPair<T>::get_HR_values(const int& index) const
{
    if (index >= this->values.size())
    {
        throw std::string("AtomPair::get_HR_values: index out of range");
    }
    return const_cast<BaseMatrix<T>&>(this->values[index]);
}

// find_R
template <typename T>
int AtomPair<T>::find_R(const int& rx_in, const int& ry_in, const int& rz_in) const
{
    for (int i = 0; i < this->R_index.size(); i++)
    {
        if (R_index[i].x == rx_in && R_index[i].y == ry_in && R_index[i].z == rz_in)
        {
            this->current_R = i;
            return i;
        }
    }
    return (-1);
}

// find_R
template <typename T>
int AtomPair<T>::find_R(const ModuleBase::Vector3<int>& R_in) const
{
    for (int i = 0; i < this->R_index.size(); i++)
    {
        if (R_index[i] == R_in)
        {
            this->current_R = i;
            return i;
        }
    }
    return (-1);
}

// find_matrix
template <typename T>
const BaseMatrix<T>* AtomPair<T>::find_matrix(const int& rx_in, const int& ry_in, const int& rz_in) const
{
    const int r_index = this->find_R(rx_in, ry_in, rz_in);
    if(r_index == -1)
    {
        return nullptr;
    }
    else
    {
        return &(this->values[r_index]);
    }
}

// find_matrix
template <typename T>
BaseMatrix<T>* AtomPair<T>::find_matrix(const int& rx_in, const int& ry_in, const int& rz_in)
{
    const int r_index = this->find_R(rx_in, ry_in, rz_in);
    if(r_index == -1)
    {
        return nullptr;
    }
    else
    {
        return &(this->values[r_index]);
    }
}

// find_matrix
template <typename T>
const BaseMatrix<T>* AtomPair<T>::find_matrix(const ModuleBase::Vector3<int>& R_in) const {
    const int r_index = this->find_R(R_in);
    if(r_index == -1)
    {
        return nullptr;
    }
    else
    {
        return &(this->values[r_index]);
    }
}

template <typename T>
BaseMatrix<T>* AtomPair<T>::find_matrix(const ModuleBase::Vector3<int>& R_in){
    const int r_index = this->find_R(R_in);
    if(r_index == -1)
    {
        return nullptr;
    }
    else
    {
        return &(this->values[r_index]);
    }
}

template <typename T>
void AtomPair<T>::convert_add(const BaseMatrix<T>& target, int rx_in, int ry_in, int rz_in)
{
    BaseMatrix<T>& matrix = this->get_HR_values(rx_in, ry_in, rz_in);
    // memory type is 2d-block
    // for 2d-block memory type, the data of input matrix is expected storing in a linear array
    // so we can use pointer to access data, and will be separate to 2d-block in this function
    matrix.add_array(target.get_pointer());
}

// function merge
template <typename T>
void AtomPair<T>::merge(const AtomPair<T>& other, bool skip_R)
{
    if (other.atom_i != atom_i || other.atom_j != atom_j)
    {
        throw std::string("AtomPair::merge: atom pair not match");
    }
    int rx = 0, ry = 0, rz = 0;
    for (int i = 0; i < other.R_index.size(); i ++)
    {
        if (!skip_R)
        {
            rx = other.R_index[i].x;
            ry = other.R_index[i].y;
            rz = other.R_index[i].z;
        }
        const BaseMatrix<T>& matrix_tmp = other.get_HR_values(i);
        //if not found, push_back this BaseMatrix to this->values
        if (this->find_R(rx, ry, rz) == -1)
        {
            this->R_index.push_back(ModuleBase::Vector3<int>(rx, ry, rz));
            this->values.push_back(matrix_tmp);
        }
        //if found but not allocated, skip this BaseMatrix values
        else if (this->values[current_R].get_pointer() == nullptr || matrix_tmp.get_pointer() == nullptr)
        {
            continue;
        }
        // if found and allocated, add data
        else
        {
            this->values[current_R].add_array(matrix_tmp.get_pointer());
        }
    }
}

// merge_to_gamma
template <typename T>
void AtomPair<T>::merge_to_gamma()
{
    // reset R_index to (0, 0, 0)
    this->R_index.clear();
    this->R_index.resize(0);
    this->R_index.push_back(ModuleBase::Vector3<int>(0, 0, 0));
    // merge all values to first BaseMatrix
    BaseMatrix<T> tmp(this->row_size, this->col_size);
    bool empty = true;
    for (int i = 0; i < this->values.size(); i++)
    {
        if(this->values[i].get_pointer() != nullptr)
        {
            if(empty)
            {
                tmp.allocate(nullptr, true);
                empty = false;
            }
            tmp.add_array(this->values[i].get_pointer());
        }
    }
    this->values.clear();
    this->values.push_back(tmp);

    this->current_R = 0;
}

template <typename T>
void AtomPair<T>::add_to_matrix(std::complex<T>* hk,
                                const int ld_hk,
                                const std::complex<T>& kphase,
                                const int hk_type) const
{
    const BaseMatrix<T>& matrix = values[current_R];
    T* hr_tmp = matrix.get_pointer();
    std::complex<T>* hk_tmp = hk;
    T* hk_real_pointer = nullptr;
    T* hk_imag_pointer = nullptr;
    const int ld_hk_2 = ld_hk * 2;
    // row major
    if (hk_type == 0)
    {
        hk_tmp += this->row_ap * ld_hk + this->col_ap;
        for (int mu = 0; mu < this->row_size; mu++)
        {
            hk_real_pointer = (T*)hk_tmp;
            hk_imag_pointer = hk_real_pointer+1;
            BlasConnector::axpy(this->col_size, kphase.real(), hr_tmp, 1, hk_real_pointer, 2);
            BlasConnector::axpy(this->col_size, kphase.imag(), hr_tmp, 1, hk_imag_pointer, 2);
            hk_tmp += ld_hk;
            hr_tmp += this->col_size;
        }
    }
    // column major
    else if (hk_type == 1)
    {
        hk_tmp += this->col_ap * ld_hk + this->row_ap;
        for (int mu = 0; mu < this->row_size; mu++)
        {
            hk_real_pointer = (T*)hk_tmp;
            hk_imag_pointer = hk_real_pointer+1;
            BlasConnector::axpy(this->col_size, kphase.real(), hr_tmp, 1, hk_real_pointer, ld_hk_2);
            BlasConnector::axpy(this->col_size, kphase.imag(), hr_tmp, 1, hk_imag_pointer, ld_hk_2);
            hk_tmp ++;
            hr_tmp += this->col_size;
        }
    }
}

// add_to_matrix
template <typename T>
void AtomPair<T>::add_to_matrix(T* hk, const int ld_hk, const T& kphase, const int hk_type) const
{
    const BaseMatrix<T>& matrix = values[current_R];
    T* hr_tmp = matrix.get_pointer();
    T* hk_tmp = hk;
    // row major
    if (hk_type == 0)
    {
        hk_tmp += this->row_ap * ld_hk + this->col_ap;
        for (int mu = 0; mu < this->row_size; mu++)
        {
            BlasConnector::axpy(this->col_size, kphase, hr_tmp, 1, hk_tmp, 1);
            /*for (int nu = 0; nu < this->col_size; nu++)
            {
                hk_tmp[nu] += matrix.get_value(mu, nu) * kphase;
            }*/
            hk_tmp += ld_hk;
            hr_tmp += this->col_size;
        }
    }
    // column major
    else if (hk_type == 1)
    {
        hk_tmp += this->col_ap * ld_hk + this->row_ap;
        for (int mu = 0; mu < this->row_size; mu++)
        {
            BlasConnector::axpy(this->col_size, kphase, hr_tmp, 1, hk_tmp, ld_hk);
            /*for (int mu = 0; mu < this->row_size; mu++)
            {
                hk_tmp[mu] += matrix.get_value(mu, nu) * kphase;
            }*/
            ++hk_tmp;
            hr_tmp += this->col_size;
        }
    }
}

template <typename T>
void AtomPair<T>::add_from_matrix(const std::complex<T>* hk,
                                const int ld_hk,
                                const std::complex<T>& kphase,
                                const int hk_type)
{
    const BaseMatrix<T>& matrix = values[current_R];
    T* hr_tmp = matrix.get_pointer();
    const std::complex<T>* hk_tmp = hk;
    const T* hk_real_pointer = nullptr;
    const T* hk_imag_pointer = nullptr;
    const int ld_hk_2 = ld_hk * 2;
    // row major
    if (hk_type == 0)
    {
        hk_tmp += this->row_ap * ld_hk + this->col_ap;
        for (int mu = 0; mu < this->row_size; mu++)
        {
            hk_real_pointer = (T*)hk_tmp;
            hk_imag_pointer = hk_real_pointer+1;
            BlasConnector::axpy(this->col_size, kphase.real(), hk_real_pointer, 2, hr_tmp, 1);
            BlasConnector::axpy(this->col_size, kphase.imag(), hk_imag_pointer, 2, hr_tmp, 1);
            hk_tmp += ld_hk;
            hr_tmp += this->col_size;
        }
    }
    // column major
    else if (hk_type == 1)
    {
        hk_tmp += this->col_ap * ld_hk + this->row_ap;
        for (int mu = 0; mu < this->row_size; mu++)
        {
            hk_real_pointer = (T*)hk_tmp;
            hk_imag_pointer = hk_real_pointer+1;
            BlasConnector::axpy(this->col_size, kphase.real(), hk_real_pointer, ld_hk_2, hr_tmp, 1);
            BlasConnector::axpy(this->col_size, kphase.imag(), hk_imag_pointer, ld_hk_2, hr_tmp, 1);
            hk_tmp ++;
            hr_tmp += this->col_size;
        }
    }
}

// add_to_matrix
template <typename T>
void AtomPair<T>::add_from_matrix(const T* hk, const int ld_hk, const T& kphase, const int hk_type)
{
    const BaseMatrix<T>& matrix = values[current_R];
    T* hr_tmp = matrix.get_pointer();
    const T* hk_tmp = hk;
    // row major
    if (hk_type == 0)
    {
        hk_tmp += this->row_ap * ld_hk + this->col_ap;
        for (int mu = 0; mu < this->row_size; mu++)
        {
            BlasConnector::axpy(this->col_size, kphase, hk_tmp, 1, hr_tmp, 1);
            /*for (int nu = 0; nu < this->col_size; nu++)
            {
                hk_tmp[nu] += matrix.get_value(mu, nu) * kphase;
            }*/
            hk_tmp += ld_hk;
            hr_tmp += this->col_size;
        }
    }
    // column major
    else if (hk_type == 1)
    {
        hk_tmp += this->col_ap * ld_hk + this->row_ap;
        for (int mu = 0; mu < this->row_size; mu++)
        {
            BlasConnector::axpy(this->col_size, kphase, hk_tmp, ld_hk, hr_tmp, 1);
            /*for (int mu = 0; mu < this->row_size; mu++)
            {
                hk_tmp[mu] += matrix.get_value(mu, nu) * kphase;
            }*/
            ++hk_tmp;
            hr_tmp += this->col_size;
        }
    }
}

// add_to_array
template <typename T>
void AtomPair<T>::add_to_array(T* array, const T& kphase) const
{
    const BaseMatrix<T>& matrix = values[current_R];
    for (int i = 0; i < this->row_size * this->col_size; i++)
    {
        array[i] += matrix.get_pointer()[i] * kphase;
    }
}

// add_to_array
template <typename T>
void AtomPair<T>::add_to_array(std::complex<T>* array, const std::complex<T>& kphase) const
{
    const BaseMatrix<T>& matrix = values[current_R];
    for (int i = 0; i < this->row_size * this->col_size; i++)
    {
        array[i] += matrix.get_pointer()[i] * kphase;
    }
}

template <typename T>
std::tuple<std::vector<int>, T*> AtomPair<T>::get_matrix_values(int ir) const
{
    if(ir<0) ir = this->current_R;
    return std::tuple<std::vector<int>, T*>({this->row_ap, this->row_size, this->col_ap, this->col_size}, this->values[ir].get_pointer());
}

// interface for get (rx, ry, rz) of index-th R-index in this->R_index, the return should be ModuleBase::Vector3<int>
template <typename T>
ModuleBase::Vector3<int> AtomPair<T>::get_R_index(const int& index) const
{
    if (index >= R_index.size() || index < 0)
    {
        std::cout << "Error: index out of range in s" << std::endl;
        return ModuleBase::Vector3<int>(-1, -1, -1);
    }
    else
    {
        // return the ModuleBase::Vector3<int> of R_index[index]
        return R_index[index];
    }
}

template <typename T>
ModuleBase::Vector3<int> AtomPair<T>::get_R_index() const
{
    // return the ModuleBase::Vector3<int> of R_index[index]
    return R_index[current_R];
}

// get_value
template <typename T>
T& AtomPair<T>::get_value(const int& i) const
{
#ifdef __DEBUG
    assert(i < this->row_size * this->col_size);
    assert(i >= 0);
    assert(current_R < this->values.size());
    assert(current_R >= 0);
#endif
    return this->values[current_R].get_pointer()[i];
}

// get_value
template <typename T>
T& AtomPair<T>::get_value(const int& row, const int& col) const
{
#ifdef __DEBUG
    assert(row < this->row_size && row >= 0);
    assert(col < this->col_size && col >= 0);
    assert(current_R < this->values.size());
    assert(current_R >= 0);
#endif
    return this->values[current_R].get_value(row, col);
}

// get_pointer
template <typename T>
T* AtomPair<T>::get_pointer(int ir) const
{
    if(ir<0)
    { 
        ir = current_R;
    }
#ifdef __DEBUG
    assert(current_R < this->values.size());
    assert(current_R >= 0);
#endif
    return this->values[ir].get_pointer();
}

// get_memory_size
template <typename T>
size_t AtomPair<T>::get_memory_size() const
{
    size_t memory_size = sizeof(*this);
    for (int i = 0; i < this->values.size(); i++)
    {
        memory_size += this->values[i].get_memory_size();
    }
    return memory_size;
}

// T of AtomPair can be double or complex<double>
template class AtomPair<double>;
template class AtomPair<std::complex<double>>;

} // namespace hamilt

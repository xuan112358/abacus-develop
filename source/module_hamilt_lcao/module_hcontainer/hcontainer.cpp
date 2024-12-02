#include "hcontainer.h"

namespace hamilt
{

// class HContainer


// destructor
template <typename T>
HContainer<T>::~HContainer()
{
    if(this->allocated)
    {
        delete[] this->wrapper_pointer;
    }
}

// copy constructor
template <typename T>
HContainer<T>::HContainer(const HContainer<T>& HR_in, T* data_array)
{
    this->sparse_ap = HR_in.sparse_ap;
    this->sparse_ap_index = HR_in.sparse_ap_index;
    this->gamma_only = HR_in.gamma_only;
    this->paraV = HR_in.paraV;
    this->current_R = -1;
    this->wrapper_pointer = data_array;
    this->allocated = false;
    this->atom_pairs = HR_in.atom_pairs;
    // data of HR_in will not be copied, please call add() after this constructor to copy data.
    this->allocate(this->wrapper_pointer, true);
    // tmp terms not copied
}

// move constructor
template <typename T>
HContainer<T>::HContainer(HContainer<T>&& HR_in)
{
    this->atom_pairs = std::move(HR_in.atom_pairs);
    this->sparse_ap = std::move(HR_in.sparse_ap);
    this->sparse_ap_index = std::move(HR_in.sparse_ap_index);
    this->gamma_only = HR_in.gamma_only;
    this->paraV = HR_in.paraV;
    this->current_R = -1;
    // tmp terms not moved
}

// simple constructor
template <typename T>
HContainer<T>::HContainer(int natom)
{
    this->gamma_only = false;
    this->current_R = -1;
    this->sparse_ap.resize(natom);
    this->sparse_ap_index.resize(natom);
}

// use unitcell to initialize atom_pairs
template <typename T>
HContainer<T>::HContainer(const UnitCell& ucell_, const Parallel_Orbitals* paraV)
{
    this->gamma_only = false;
    this->current_R = -1;
    std::vector<int> atom_begin_row(ucell_.nat+1, 0);
    std::vector<int> atom_begin_col(ucell_.nat+1, 0);
    int begin = 0;
    for(int i=0;i<ucell_.nat;++i)
    {
        int it = ucell_.iat2it[i];
        begin += ucell_.atoms[it].nw;
        atom_begin_row[i+1] = begin;
        atom_begin_col[i+1] = begin;
    }
    if(paraV == nullptr)
    {
        // initialize atom_pairs and sparse_ap
        this->atom_pairs.resize(ucell_.nat * ucell_.nat, AtomPair<T>(0, 0));
        this->sparse_ap.resize(ucell_.nat);
        this->sparse_ap_index.resize(ucell_.nat);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < ucell_.nat; i++)
        {
            this->sparse_ap[i].resize(ucell_.nat);
            this->sparse_ap_index[i].resize(ucell_.nat);
            for (int j = 0; j < ucell_.nat; j++)
            {
                //AtomPair<T> atom_ij(i, j, atom_begin_row.data(), atom_begin_col.data(), ucell_.nat);
                this->atom_pairs[i * ucell_.nat + j] = AtomPair<T>(i, j, atom_begin_row.data(), atom_begin_col.data(), ucell_.nat);
                this->sparse_ap[i][j] = j;
                this->sparse_ap_index[i][j] = i * ucell_.nat + j;
            }
        }
    }
    else
    {
        this->paraV = paraV;
        // initialize atom_pairs and sparse_ap
        this->sparse_ap.resize(ucell_.nat);
        this->sparse_ap_index.resize(ucell_.nat);
        for (int i = 0; i < ucell_.nat; i++)
        {
            for (int j = 0; j < ucell_.nat; j++)
            {
                //check if atom_pair(i, j) is empty in this process
                if(paraV->get_row_size(i) <= 0 || paraV->get_col_size(j) <= 0)
                {
                    continue;
                }
                AtomPair<T> atom_ij(i, j, paraV);
                this->insert_pair(atom_ij);
            }
        }
    }
    this->allocate(nullptr, true);
}

//HContainer(const Parallel_Orbitals* paraV, T* data_pointer = nullptr);
template <typename T>
HContainer<T>::HContainer(const Parallel_Orbitals* paraV_in, T* data_pointer, const std::vector<int>* ijr_info)
{
    this->current_R = -1;

    // use HContainer as a wrapper(!nullptr) or container(nullptr)
    this->wrapper_pointer = data_pointer;

#ifdef __DEBUG
    assert(paraV_in != nullptr);
#endif

    // save Parallel_Orbitals pointer
    this->paraV = paraV_in;
    // initialize sparse_ap
    int natom = paraV->atom_begin_row.size();
    this->sparse_ap.resize(natom);
    this->sparse_ap_index.resize(natom);
    if(ijr_info != nullptr)
    {
        this->insert_ijrs(ijr_info);
        // allocate memory
        this->allocate(data_pointer, false);
    }
}

// allocate
template <typename T>
void HContainer<T>::allocate(T* data_array, bool is_zero)
{
    size_t nnr = this->get_nnr();
    if(this->allocated)
    {// delete existed memory of this->wrapper_pointer
        delete[] this->wrapper_pointer;
        this->allocated = false;
    }
    if(data_array == nullptr)
    {
        // use this->wrapper_pointer as data_array
        this->allocated = true;
        this->wrapper_pointer = new T[nnr];
        ModuleBase::GlobalFunc::ZEROS(this->wrapper_pointer, nnr);
        data_array = this->wrapper_pointer;
    }
    else
    {
        // use data_array to replace this->wrapper_pointer
        this->wrapper_pointer = data_array;
        if(is_zero)
        {
            ModuleBase::GlobalFunc::ZEROS(this->wrapper_pointer, nnr);
        }
    }
    for(int it=0;it<this->atom_pairs.size();it++)
    {
        this->atom_pairs[it].allocate(data_array, false);
        // move data_array pointer for the next AtomPair
        data_array += this->atom_pairs[it].get_R_size() * this->atom_pairs[it].get_size();
    }
}

// set_zero
template <typename T>
void HContainer<T>::set_zero()
{
    for(auto& it : this->atom_pairs)
    {
        it.set_zero();
    }
}

template <typename T>
AtomPair<T>* HContainer<T>::find_pair(int atom_i, int atom_j) const
{
    if(atom_i >= this->sparse_ap.size())
    {
        ModuleBase::WARNING_QUIT("HContainer::insert_pair", "atom_i out of range");
    }
    // search atom_i and atom_j in sparse_ap
    auto it = std::lower_bound(this->sparse_ap[atom_i].begin(), this->sparse_ap[atom_i].end(), atom_j);
    if (it != this->sparse_ap[atom_i].end() && *it == atom_j)
    {
        AtomPair<T>* tmp_pointer = const_cast<AtomPair<T>*>(&this->atom_pairs[this->sparse_ap_index[atom_i][it-this->sparse_ap[atom_i].begin()]]);
        return tmp_pointer;
    }
    else
    {
        return nullptr;
    }
}

// find_matrix
template <typename T>
const BaseMatrix<T>* HContainer<T>::find_matrix(int atom_i, int atom_j, int rx, int ry, int rz) const
{
    AtomPair<T>* tmp = this->find_pair(atom_i, atom_j);
    if(tmp == nullptr)
    {
        return nullptr;
    }
    else
    {
        if(this->gamma_only)
        {
            return tmp->find_matrix(0, 0, 0);
        }
        else
        {
            return tmp->find_matrix(rx, ry, rz);
        }
    }
}

template <typename T>
BaseMatrix<T>* HContainer<T>::find_matrix(int atom_i, int atom_j, int rx, int ry, int rz)
{
    AtomPair<T>* tmp = this->find_pair(atom_i, atom_j);
    if(tmp == nullptr)
    {
        return nullptr;
    }
    else
    {
        if(this->gamma_only)
        {
            return tmp->find_matrix(0, 0, 0);
        }
        else
        {
            return tmp->find_matrix(rx, ry, rz);
        }
    }
}

template <typename T>
BaseMatrix<T>* HContainer<T>::find_matrix(int atom_i, int atom_j, const ModuleBase::Vector3<int>& R_index) 
{
    AtomPair<T>* tmp = this->find_pair(atom_i, atom_j);
    if(tmp == nullptr)
    {
        return nullptr;
    }
    else
    {
        if(this->gamma_only)
        {
            return tmp->find_matrix(0, 0, 0);
        }
        else
        {
            return tmp->find_matrix(R_index);
        }
    }
}

template <typename T>
const BaseMatrix<T>* HContainer<T>::find_matrix(int atom_i, int atom_j,const ModuleBase::Vector3<int>& R_index) const
{
    AtomPair<T>* tmp = this->find_pair(atom_i, atom_j);
    if(tmp == nullptr)
    {
        return nullptr;
    }
    else
    {
        if(this->gamma_only)
        {
            return tmp->find_matrix(0, 0, 0);
        }
        else
        {
            return tmp->find_matrix(R_index);
        }
    }
}

template <typename T>
int HContainer<T>::find_matrix_offset(int atom_i, int atom_j, int rx, int ry, int rz) const
{
    const BaseMatrix<T>* tmp = this->find_matrix(atom_i, atom_j, rx, ry, rz);
    if(tmp == nullptr)
    {
        return -1;
    }
    else
    {
        return tmp->get_pointer() - this->wrapper_pointer;
    }
}

template <typename T>
int HContainer<T>::find_matrix_offset(int atom_i, int atom_j, const ModuleBase::Vector3<int>& R_index) const
{
    const BaseMatrix<T>* tmp = this->find_matrix(atom_i, atom_j, R_index);
    if(tmp == nullptr)
    {
        return -1;
    }
    else
    {
        return tmp->get_pointer() - this->wrapper_pointer;
    }
}

// get_atom_pair with atom_ij
template <typename T>
AtomPair<T>& HContainer<T>::get_atom_pair(int atom_i, int atom_j) const
{
    if(atom_i >= this->sparse_ap.size())
    {
        ModuleBase::WARNING_QUIT("HContainer::insert_pair", "atom_i out of range");
    }
    // search atom_i and atom_j in sparse_ap
    auto it = std::lower_bound(this->sparse_ap[atom_i].begin(), this->sparse_ap[atom_i].end(), atom_j);
    if (it != this->sparse_ap[atom_i].end() && *it == atom_j)
    {
        AtomPair<T>* tmp = const_cast<AtomPair<T>*>(&this->atom_pairs[this->sparse_ap_index[atom_i][it-this->sparse_ap[atom_i].begin()]]);
        return *tmp;
    }
    else
    {
        ModuleBase::WARNING_QUIT("HContainer", "atom pair not found in get_atom_pair");
        return const_cast<AtomPair<T>&>(this->atom_pairs[0]);
    }
}

// get_atom_pair with index
template <typename T>
AtomPair<T>& HContainer<T>::get_atom_pair(int index) const
{
#ifdef __DEBUG
    if (this->current_R > -1)
    {
        if (index >= this->tmp_atom_pairs.size() || index < 0)
        {
            std::cout << "Error: index out of range in get_atom_pair" << std::endl;
            exit(1);
        }
    }
    else
    {
        if (index >= this->atom_pairs.size() || index < 0)
        {
            std::cout << "Error: index out of range in get_atom_pair" << std::endl;
            exit(1);
        }
    }
#endif
    if (this->current_R > -1)
    {
        return const_cast<AtomPair<T>&>(*this->tmp_atom_pairs[index]);
    }
    else
    {
        return const_cast<AtomPair<T>&>(this->atom_pairs[index]);
    }
}

// add function
template <typename T>
void HContainer<T>::add(const HContainer<T>& other)
{
    for(int iap=0;iap<other.size_atom_pairs();iap++)
    {
        auto tmp = other.get_atom_pair(iap);
        this->insert_pair(tmp);
    }
}

template <typename T>
bool HContainer<T>::fix_R(int rx_in, int ry_in, int rz_in) const
{
    // clear and reallocate the memory of this->tmp_atom_pairs
    //this->tmp_atom_pairs.clear();
    //this->tmp_atom_pairs.shrink_to_fit();
    //this->tmp_atom_pairs.reserve(this->atom_pairs.size());
    this->tmp_atom_pairs.resize(this->atom_pairs.size());
    int iter = 0;
    // find (rx, ry, rz) in this->atom_pairs[i].R_values
    for (auto it = this->atom_pairs.begin(); it != this->atom_pairs.end(); ++it)
    {
        if (it->find_R(rx_in, ry_in, rz_in) != -1)
        {
            // push bach the pointer of AtomPair to this->tmp_atom_pairs
            const AtomPair<T>* tmp_pointer = &(*it);
            this->tmp_atom_pairs[iter++] = tmp_pointer;
            //this->tmp_atom_pairs.push_back(tmp_pointer);
        }
    }
    if (iter == 0)
    {
        std::cout << "Error: no atom pair found in fix_R" << std::endl;
        this->current_R = -1;
        return false;
    }
    else
    {
        //set current_R
        this->current_R = this->find_R(rx_in, ry_in, rz_in);
        this->tmp_atom_pairs.resize(iter);
        return true;
    }
}

// unfix_R
template <typename T>
void HContainer<T>::unfix_R() const
{
    this->current_R = -1;
    this->tmp_atom_pairs.clear();
    this->tmp_atom_pairs.shrink_to_fit();
}

// fix_gamma
template <typename T>
void HContainer<T>::fix_gamma()
{
    // every AtomPair in this->atom_pairs has the (0, 0, 0) cell index
    // fix every AtomPair in this->atom_pairs to only center cell
    this->gamma_only = true;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int it =0; it< this->atom_pairs.size(); ++it)
    {
        this->atom_pairs[it].merge_to_gamma();
    }
    // in gamma_only case, R_index is not needed, tmp_R_index should be empty
    if (this->current_R != -1)
    {
        this->current_R = -1;
        this->tmp_R_index.clear();
        this->tmp_R_index.shrink_to_fit();
    }
}

// find_R
template <typename T>
int HContainer<T>::find_R(const int& rx_in, const int& ry_in, const int& rz_in) const
{
    // search (rx, ry, rz) in this->tmp_R_index
    if (this->tmp_R_index.empty())
    {
        return -1;
    }
    for (int i = 0; i < this->tmp_R_index.size(); i++)
    {
        if (this->tmp_R_index[i].x == rx_in && this->tmp_R_index[i].y == ry_in
            && this->tmp_R_index[i].z == rz_in)
        {
            return i;
        }
    }
    return -1;
}

template <typename T>
int HContainer<T>::find_R(const ModuleBase::Vector3<int>& R_in) const
{
    // search R_in in this->tmp_R_index
    if (this->tmp_R_index.empty())
    {
        return -1;
    }
    for (int i = 0; i < this->tmp_R_index.size(); i++)
    {
        if (this->tmp_R_index[i] == R_in)
        {
            return i;
        }
    }
    return -1;
}

// size_R_loop, return the number of different cells in this->atom_pairs
template <typename T>
size_t HContainer<T>::size_R_loop() const
{
    // R index is fixed
    if (this->current_R > -1 && this->tmp_R_index.size() > 0)
    {
        return this->tmp_R_index.size();
    }
    /**
     * start a new iteration of loop_R
     * there is different R-index in this->atom_pairs[i].R_values
     * search them one by one and store them in this->tmp_R_index
     */
    this->tmp_R_index.clear();
    for (auto it = this->atom_pairs.begin(); it != this->atom_pairs.end(); ++it)
    {
        /**
         * search (rx, ry, rz) with (it->R_values[i].x, it->R_values[i].y, it->R_values[i].z)
         * if (rx, ry, rz) not found in this->tmp_R_index,
         * insert the (rx, ry, rz) into end of this->tmp_R_index
         * no need to sort this->tmp_R_index, using find_R() to find the (rx, ry, rz) -> int in tmp_R_index
         */
        for (int iR = 0; iR < it->get_R_size(); iR++)
        {
            const ModuleBase::Vector3<int> r_vec = it->get_R_index(iR);
            int it_tmp = this->find_R(r_vec);
            if (it_tmp == -1)
            {
                this->tmp_R_index.push_back(ModuleBase::Vector3<int>(r_vec));
            }
        }
    }
    return this->tmp_R_index.size();
}

template <typename T>
void HContainer<T>::loop_R(const size_t& index, int& rx, int& ry, int& rz) const
{
#ifdef __DEBUG
    if (index >= this->tmp_R_index.size())
    {
        std::cout << "Error: index out of range in loop_R" << std::endl;
        exit(1);
    }
#endif
    // set rx, ry, rz
    rx = this->tmp_R_index[index].x;
    ry = this->tmp_R_index[index].y;
    rz = this->tmp_R_index[index].z;
    return;
}

// get_AP_size
template <typename T>
size_t HContainer<T>::size_atom_pairs() const
{
    // R index is fixed
    if (this->current_R > -1)
    {
        return this->tmp_atom_pairs.size();
    }
    // R index is not fixed
    else
    {
        return this->atom_pairs.size();
    }
}

// data() interface with atom_i and atom_j
template <typename T>
T* HContainer<T>::data(int atom_i, int atom_j) const
{
    AtomPair<T>* atom_ij = this->find_pair(atom_i, atom_j);
    if (atom_ij != nullptr)
    {
        return atom_ij->get_pointer();
    }
    else
    {
        std::cout << "Error: atom pair not found in data()" << std::endl;
        return nullptr;
    }
}

// data() interface with atom_i and atom_j ad R_pointer
template <typename T>
T* HContainer<T>::data(int atom_i, int atom_j, int* R_pointer) const
{
    AtomPair<T>* atom_ij = this->find_pair(atom_i, atom_j);
    if (atom_ij != nullptr)
    {
        return atom_ij->get_HR_values(R_pointer[0], R_pointer[1], R_pointer[2]).get_pointer();
    }
    else
    {
        std::cout << "Error: atom pair not found in data()" << std::endl;
        return nullptr;
    }
}

// insert_pair
template <typename T>
void HContainer<T>::insert_pair(const AtomPair<T>& atom_ij)
{
    int atom_i = atom_ij.get_atom_i();
    if(atom_i >= this->sparse_ap.size())
    {
        ModuleBase::WARNING_QUIT("HContainer::insert_pair", "atom_i out of range");
    }
    int atom_j = atom_ij.get_atom_j();
    // find atom_ij in this->atom_pairs
    // 1. find the index of atom_j in sparse_ap[atom_i]
    auto it = std::lower_bound(this->sparse_ap[atom_i].begin(),
                               this->sparse_ap[atom_i].end(), atom_j);
    if (it != this->sparse_ap[atom_i].end() && *it == atom_j)
    {
        // 2. merge atom_ij
        this->atom_pairs[this->sparse_ap_index[atom_i][it-this->sparse_ap[atom_i].begin()]].merge(atom_ij, this->gamma_only);
    }
    else
    {
        // 3. insert atom_ij
        // check paraV pointer
        if (this->paraV != atom_ij.get_paraV())
        {
            ModuleBase::WARNING_QUIT("HContainer::insert_pair", "atom_ij has different paraV pointer as HContainer");
        }
        else
        { //insert atom_ij, and set paraV pointer for HContainer if atom_ij has paraV pointer
            this->atom_pairs.push_back(atom_ij);
            //if gamma_only case, merge atom_ij to gamma_only
            if(this->gamma_only)
            {
                this->atom_pairs.back().merge_to_gamma();
            }
            // update sparse_ap
            int index = it - this->sparse_ap[atom_i].begin();
            if(it != this->sparse_ap[atom_i].end())
            {
                this->sparse_ap[atom_i].insert(this->sparse_ap[atom_i].begin() + index, atom_j);
                this->sparse_ap_index[atom_i].insert(this->sparse_ap_index[atom_i].begin() + index, this->atom_pairs.size() - 1);
            }
            else
            {
                this->sparse_ap[atom_i].push_back(atom_j);
                this->sparse_ap_index[atom_i].push_back(this->atom_pairs.size() - 1);
            }
        }
    }
}

//operator() is not implemented now, this interface is too expensive to access data
/*template <typename T>
T& HContainer<T>::operator()(int atom_i, int atom_j, int rx_in, int ry_in, int rz_in, int mu, int nu) const
{
        return this->get_atom_pair(atom_i, atom_j).get_HR_values(rx_in, ry_in, rz_in).get_value(mu, nu);
}*/

//get_current_R
template <typename T>
int HContainer<T>::get_current_R() const
{
    return this->current_R;
}

//is_gamma_only
template <typename T>
bool HContainer<T>::is_gamma_only() const
{
    return this->gamma_only;
}

//get_memory_size
template <typename T>
size_t HContainer<T>::get_memory_size() const
{
    size_t memory = sizeof(*this);
    memory += this->atom_pairs.capacity() * sizeof(AtomPair<T>);
    memory += this->sparse_ap.capacity() * sizeof(std::vector<int>);
    memory += this->sparse_ap_index.capacity() * sizeof(std::vector<int>);
    for(int i=0;i<this->atom_pairs.size();++i)
    {
        memory += this->atom_pairs[i].get_memory_size();
    }
    for(int i=0;i<this->sparse_ap.size();++i)
    {
        memory += this->sparse_ap[i].capacity() * sizeof(int);
        memory += this->sparse_ap_index[i].capacity() * sizeof(int);
    }
    memory += this->tmp_atom_pairs.capacity() * sizeof(AtomPair<T>*);
    memory += this->tmp_R_index.capacity() * sizeof(ModuleBase::Vector3<int>);
    if(this->allocated)
    {
        memory += this->get_nnr() * sizeof(T);
    }
    return memory;
}

// synchronize
template <typename T>
void HContainer<T>::shape_synchron( const HContainer<T>& other)
{
    // check paraV pointer
    if (this->paraV != other.paraV)
    {
        ModuleBase::WARNING_QUIT("HContainer::synchronize", "paraV pointer not match");
    }
    // synchronize atom_pairs
    for (int i = 0; i < other.atom_pairs.size(); ++i)
    {
        const int iat1 = other.atom_pairs[i].get_atom_i();
        const int iat2 = other.atom_pairs[i].get_atom_j();
        AtomPair<T>* tmp_pointer = this->find_pair(iat1, iat2);
        if (tmp_pointer == nullptr)
        {
            this->insert_pair(other.atom_pairs[i]);
            // the new AtomPair should be zero
            this->atom_pairs.back().set_zero();
        }
        else
        {
            for(int ir = 0;ir < other.atom_pairs[i].get_R_size();++ir)
            {
                const ModuleBase::Vector3<int> R_vec = other.atom_pairs[i].get_R_index(ir);
                if(tmp_pointer->find_R(R_vec) != -1)
                {
                    // do nothing
                }
                else
                {
                    // insert the new BaseMatrix
                    tmp_pointer->get_HR_values(R_vec.x, R_vec.y, R_vec.z);
                }
            }
        }
    }
}

// get_IJR_info
template <typename T>
std::vector<int> HContainer<T>::get_ijr_info() const
{
    // get number of atom pairs
    std::vector<int> ijr_info;
    ijr_info.push_back(this->atom_pairs.size());
    // loop atom pairs
    for (int i = 0; i < this->atom_pairs.size(); ++i)
    {
        // get atom_i and atom_j
        const int atom_i = this->atom_pairs[i].get_atom_i();
        const int atom_j = this->atom_pairs[i].get_atom_j();
        // push back atom_i, atom_j
        ijr_info.push_back(atom_i);
        ijr_info.push_back(atom_j);
        // get number of R
        const int number_R = this->atom_pairs[i].get_R_size();
        ijr_info.push_back(number_R);
        // loop R
        for (int ir = 0; ir < number_R; ++ir)
        {
            const ModuleBase::Vector3<int> R_vec = this->atom_pairs[i].get_R_index(ir);
            ijr_info.push_back(R_vec.x);
            ijr_info.push_back(R_vec.y);
            ijr_info.push_back(R_vec.z);
        }
    }
    return ijr_info;
}

template<typename T>
void HContainer<T>::insert_ijrs(const std::vector<int>* ijrs)
{
    if (this->paraV == nullptr)
    {
        ModuleBase::WARNING_QUIT("HContainer::insert_ijrs", "paraV pointer can not be nullptr!");
    }
    // get number of atom pairs
    const int number_ap = (*ijrs)[0];
    // loop AtomPairs and unpack values
    const int* ijr_p = ijrs->data() + 1;
    for (int i = 0; i < number_ap; ++i)
    {
        // get atom_i and atom_j
        const int atom_i = *ijr_p++;
        const int atom_j = *ijr_p++;
        // get number of R
        const int number_R = *ijr_p++;
        for (int k = 0; k < number_R; ++k)
        {
            int r_index[3];
            // get R index
            r_index[0] = *ijr_p++;
            r_index[1] = *ijr_p++;
            r_index[2] = *ijr_p++;
            //insert this IJ AtomPair
            AtomPair<T> tmp_ap(atom_i, atom_j, r_index[0], r_index[1], r_index[2], this->paraV);
            this->insert_pair(tmp_ap);
        }
    }
}

template <typename T>
void HContainer<T>::insert_ijrs(const std::vector<int>* ijrs, const UnitCell& ucell, const int npol)
{
    // prepare the row_index and col_index for construct AtomPairs, they are
    // same, name as orb_index
    std::vector<int> orb_index(ucell.nat + 1);
    orb_index[0] = 0;
    for (int i = 1; i < orb_index.size(); i++) {
        int type = ucell.iat2it[i - 1];
        orb_index[i] = orb_index[i - 1] + ucell.atoms[type].nw * npol;
    }

    // get number of atom pairs
    const int number_ap = (*ijrs)[0];
    // loop AtomPairs and unpack values
    const int* ijr_p = ijrs->data() + 1;
    for (int i = 0; i < number_ap; ++i)
    {
        // get atom_i and atom_j
        const int atom_i = *ijr_p++;
        const int atom_j = *ijr_p++;
        // get number of R
        const int number_R = *ijr_p++;
        for (int k = 0; k < number_R; ++k)
        {
            int r_index[3];
            // get R index
            r_index[0] = *ijr_p++;
            r_index[1] = *ijr_p++;
            r_index[2] = *ijr_p++;
            //insert this IJ AtomPair
            AtomPair<T> tmp_ap(atom_i, atom_j, r_index[0], r_index[1], r_index[2], orb_index.data(), orb_index.data(), ucell.nat);
            this->insert_pair(tmp_ap);
        }
    }
}

// T of HContainer can be double or complex<double>
template class HContainer<double>;
template class HContainer<std::complex<double>>;

} // end namespace hamilt
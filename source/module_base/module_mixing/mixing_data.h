#ifndef MIXING_DATA_H_
#define MIXING_DATA_H_
#include <vector>

#include "module_base/module_container/ATen/tensor.h"
namespace Base_Mixing
{

/**
 * @brief data for Mixing class
 *
 */
class Mixing_Data
{
  public:
    Mixing_Data() = default;
    /**
     * @brief Construct a new Mixing_Data object
     *
     * @param ndim store ndim vectors for mixing
     * @param length the length of each vector
     * @param type_size size of type
     *
     */
    Mixing_Data(const int& ndim, const std::size_t& length, const size_t& type_size);

    /**
     * @brief Destroy the Mixing_Data object
     *
     */
    ~Mixing_Data();

    /**
     * @brief resize the data
     *
     * @param ndim store ndim vectors for mixing
     * @param length the length of each vector
     * @param type_size size of type
     *
     */
    void resize(const int& ndim, const std::size_t& length, const size_t& type_size);

    /**
     * @brief push data to the tensor
     *
     */
    template <typename FPTYPE>
    void push(const FPTYPE* data_in)
    {
        this->start = (this->start + 1) % this->ndim_tot;
        this->ndim_use = std::min(this->ndim_use + 1, this->ndim_tot);
        ++this->ndim_history;
        FPTYPE* FP_startdata = static_cast<FPTYPE*>(this->data) + this->start * this->length;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for (std::size_t i = 0; i < length; ++i)
        {
            FP_startdata[i] = data_in[i];
        }
    }

    /**
     * @brief reset mixing
     *
     */
    void reset()
    {
        this->ndim_use = 0;
        this->ndim_history = 0;
        this->start = -1;
    }

    /**
     * @brief get the index of i-th vector
     *
     */
    int index_move(const int& n) const
    {
        return (n + this->start + ndim_tot) % ndim_tot;
    }

  public:
    // Tensor pointer to store the data
    void* data = nullptr;
    // the number of vectors for mixing
    int ndim_tot = 0;
    // the length of each vector
    std::size_t length = 0;
    // the start index for vector: start = this->index_move(0)
    int start = -1;
    // the number of used vectors for mixing
    int ndim_use = 0;
    // the number of history vectors
    int ndim_history = 0;
};

} // namespace Base_Mixing
#endif
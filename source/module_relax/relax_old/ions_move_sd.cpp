#include "ions_move_sd.h"

#include "module_parameter/parameter.h"
#include "ions_move_basic.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include<vector>

using namespace Ions_Move_Basic;

Ions_Move_SD::Ions_Move_SD()
{
    this->energy_saved = 1.0e10;
    this->grad_saved = nullptr;
    this->pos_saved = nullptr;
}
Ions_Move_SD::~Ions_Move_SD()
{
    delete[] grad_saved;
    delete[] pos_saved;
}

void Ions_Move_SD::allocate()
{
    ModuleBase::TITLE("Ions_Move_SD", "allocate");
    assert(dim > 0);
    delete[] grad_saved;
    delete[] pos_saved;
    this->grad_saved = new double[dim];
    this->pos_saved = new double[dim];
    ModuleBase::GlobalFunc::ZEROS(grad_saved, dim);
    ModuleBase::GlobalFunc::ZEROS(pos_saved, dim);
}

void Ions_Move_SD::start(UnitCell& ucell, const ModuleBase::matrix& force, const double& etot_in)
{
    ModuleBase::TITLE("Ions_Move_SD", "start");

    assert(dim > 0);
    assert(grad_saved != nullptr);
    assert(pos_saved != nullptr);

    std::vector<double> pos(dim);
    std::vector<double> grad(dim);
    std::vector<double> move(dim);
    ModuleBase::GlobalFunc::ZEROS(pos.data(), dim);
    ModuleBase::GlobalFunc::ZEROS(grad.data(), dim);
    ModuleBase::GlobalFunc::ZEROS(move.data(), dim);

    // 1: ediff = 0
    // 0: ediff < 0
    bool judgement = false;
    setup_etot(etot_in, judgement);
    setup_gradient(ucell, force, pos.data(), grad.data());

    if (istep == 1 || etot_in <= energy_saved)
    {
        printf("in cheak_converged");
        printf("pos[0]: %f\n", pos[0]);
        energy_saved = etot_in;
        for (int i = 0; i < dim; i++) {
            pos_saved[i] = pos[i];
}
        for (int i = 0; i < dim; i++)
        {
            grad_saved[i] = grad[i];
        }
        // normalize the gradient, in convinience to
        // move atom.
        double norm = dot_func(grad_saved, grad_saved, dim);
        norm = sqrt(norm);
        for (int i = 0; i < dim; i++)
        {
            grad_saved[i] /= norm;
        }
    }

    Ions_Move_Basic::check_converged(ucell, grad.data());
    if (Ions_Move_Basic::converged)
    {
        Ions_Move_Basic::terminate(ucell);
    }
    else
    {
        this->cal_tradius_sd();
        for (int i = 0; i < dim; i++)
        {
            move[i] = -grad_saved[i] * trust_radius;
        }
        move_atoms(ucell, move.data(), pos_saved);
        Ions_Move_Basic::update_iter++;
    }

    return;
}

void Ions_Move_SD::cal_tradius_sd() const
{
    static int accepted_number = 0;

    if (Ions_Move_Basic::istep == 1)
    {
        Ions_Move_Basic::trust_radius = Ions_Move_Basic::relax_bfgs_init;
    }
    else if (Ions_Move_Basic::istep > 1)
    {
        if (Ions_Move_Basic::ediff < 0.0)
        {
            accepted_number++;
            if (accepted_number > 3 && accepted_number % 3 == 1)
            {
                Ions_Move_Basic::trust_radius *= 1.5;
            }
        }
        else if (Ions_Move_Basic::ediff >= 0.0) // == 0 means no accept!
        {
            accepted_number = 0;
            Ions_Move_Basic::trust_radius *= 0.5;
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("Ions_Move_SD::cal_tradius_sd", "istep < 1!");
    }
    if (PARAM.inp.out_level == "ie")
    {
        std::cout << " SD RADIUS (Bohr)     : " << trust_radius << std::endl;
    }
    return;
}

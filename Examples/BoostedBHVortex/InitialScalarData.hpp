/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

#include "ComplexScalarField.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include <cmath>

//! Class which creates a constant scalar field given params for initial
//! matter config
class InitialScalarData
{
  public:
    struct params_t
    {
        double mass;
        double amplitude;
        double omega;
        std::array<double, CH_SPACEDIM> center;
    };

    //! The constructor for the class
    InitialScalarData(params_t a_params, double dx)
        : m_params(a_params), m_dx(dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        ComplexScalarField<>::Vars<data_t> vars;
        VarsTools::assign(vars, 0.);
        const data_t R = coords.get_radius();
        const data_t R2 = R * R;
        const data_t rho2 = simd_max(R2 - coords.z * coords.z, 1e-8);
        const data_t rho = sqrt(rho2);
        const data_t cos_phi = coords.x / rho;
        const data_t sin_phi = coords.y / rho;

        data_t BesselJ = m_params.amplitude * j1(m_params.omega * rho);

        const data_t omega0 = sqrt(m_params.omega * m_params.omega +
                                   m_params.mass * m_params.mass);

        vars.phi_Re = BesselJ * cos_phi;
        vars.phi_Im = BesselJ * sin_phi;
        vars.Pi_Re = BesselJ * m_params.mass * sin_phi;
        vars.Pi_Im = -BesselJ * m_params.mass * cos_phi;

        current_cell.store_vars(vars);
    }

  protected:
    const params_t m_params;
    const double m_dx;
};

#endif /* INITIALSCALARDATA_HPP_ */

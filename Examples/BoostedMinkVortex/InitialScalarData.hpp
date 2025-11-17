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
#include <complex>

#include "ADMFixedBGVars.hpp"

//! Class which creates a constant scalar field given params for initial
//! matter config
template <class background_t> class InitialScalarData
{

    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  protected:
    const background_t m_background; //!< The metric background

  public:
    struct params_t
    {
        double mass;
        double amplitude;
        double velocity;
        double radius;
        double d_to_bh;
        std::array<double, CH_SPACEDIM> center;
    };

    //! The constructor for the class
    InitialScalarData(params_t a_params, double dx, background_t a_background)
        : m_params(a_params), m_dx(dx), m_background(a_background)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        ComplexScalarField<>::Vars<data_t> vars;
        VarsTools::assign(vars, 0.);

        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, coords);

        data_t xx = (coords.x - m_params.d_to_bh);
        double yy = coords.y;
        double zz = coords.z;
        double rho = simd_max(sqrt(xx * xx + yy * yy), 1.e-6);

        const data_t cos_phi = xx / rho;
        const data_t sin_phi = yy / rho;
        const double k = 1.84 / m_params.radius;
        const double omega = sqrt(k * k + m_params.mass * m_params.mass);

        const data_t v = m_params.velocity;
        const data_t A = m_params.amplitude;

        data_t J0 = j0(k * rho);
        data_t J1 = j1(k * rho);
        data_t J2 = jn(2, k * rho);
        data_t dJ1_dr = 0.5 * k * (J0 - J2);
        data_t J1_by_r = J1 / rho;

        vars.phi_Re = A * J1 * cos_phi;
        vars.phi_Im = A * J1 * sin_phi;

        data_t phi_Re = J1 * cos_phi;
        data_t phi_Im = J1 * sin_phi;

        // d/dx components
        Tensor<1, data_t> dphi_Re;
        dphi_Re[0] = (dJ1_dr * cos_phi * cos_phi + J1_by_r * sin_phi * sin_phi);
        dphi_Re[1] = (dJ1_dr * sin_phi * cos_phi - J1_by_r * cos_phi * sin_phi);
        dphi_Re[2] = 0;

        Tensor<1, data_t> dphi_Im;
        dphi_Im[0] = (dJ1_dr * sin_phi * cos_phi - J1_by_r * cos_phi * sin_phi);
        dphi_Im[1] = (dJ1_dr * sin_phi * sin_phi + J1_by_r * cos_phi * cos_phi);
        dphi_Im[2] = 0;

        // β · ∇φ terms
        data_t beta_grad_phiRe = 0;
        data_t beta_grad_phiIm = 0;
        // data_t beta_grad_phiRe =  v * dphi_Re[0];
        // data_t beta_grad_phiIm =  v * dphi_Im[0];

        vars.Pi_Re = A * (omega * phi_Im - beta_grad_phiRe) / metric_vars.lapse;
        vars.Pi_Im =
            A * (-omega * phi_Re - beta_grad_phiIm) / metric_vars.lapse;
        // FOR1(i)
        //{
        //   vars.Pi_Re +=
        //       -A * metric_vars.shift[i] * dphi_Re[i] / metric_vars.lapse;
        //   vars.Pi_Im +=
        //       -A * metric_vars.shift[i] * dphi_Im[i] / metric_vars.lapse;
        // }

        current_cell.store_vars(vars);
    }

  protected:
    const params_t m_params;
    const double m_dx;
};

#endif /* INITIALSCALARDATA_HPP_ */

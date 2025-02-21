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

#include "ADMFixedBGVars.hpp"

//! Class which creates a string + BH params for initial
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
        double *f_profile;
        double *df_profile;
        double pot_eta;
        double pot_lambda;
        double spacing;
        double d_to_bh;
        double velocity;
        std::array<double, CH_SPACEDIM> center; //!< The center of the BH
    };

    //! The constructor for the class
    InitialScalarData(params_t a_params, double a_dx, background_t a_background)
        : m_params(a_params), m_dx(a_dx), m_background(a_background)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

        ComplexScalarField<>::Vars<data_t> vars;
        VarsTools::assign(vars, 0.);
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, coords);

        const double v = m_params.velocity;
        const double v2 = v * v;
        const double boost2 = 1.0 / (1 - v2);
        const double boost = sqrt(boost2);

        double xx = (coords.x - m_params.d_to_bh);
        double yy = coords.y;
        double zz = coords.z;
        double rho = sqrt(xx * xx + yy * yy);
        double sinth = yy / rho;
        double costh = xx / rho;

        double pot_lambda_ref = 2.0;
        double pot_eta_ref = 0.1;

        int indxL = static_cast<int>(floor(rho * sqrt(m_params.pot_lambda / pot_lambda_ref) * m_params.pot_eta / pot_eta_ref / m_params.spacing));
        int indxH = static_cast<int>(ceil(rho  * sqrt(m_params.pot_lambda / pot_lambda_ref) * m_params.pot_eta / pot_eta_ref / m_params.spacing));

        double f_L = *(m_params.f_profile + indxL);
        double f_H = *(m_params.f_profile + indxH);
        double df_L = *(m_params.df_profile + indxL);
        double df_H = *(m_params.df_profile + indxH);

        double fvals = f_L + (rho * sqrt(m_params.pot_lambda / pot_lambda_ref) * m_params.pot_eta / pot_eta_ref / m_params.spacing - indxL) * (f_H - f_L);
        double dfvals = df_L + (rho * sqrt(m_params.pot_lambda / pot_lambda_ref) * m_params.pot_eta / pot_eta_ref/ m_params.spacing - indxL) * (df_H - df_L);
        dfvals *= sqrt(m_params.pot_lambda / pot_lambda_ref) * m_params.pot_eta / pot_eta_ref;

        Tensor<1, data_t> grad_f;
        grad_f[0] = dfvals * xx / rho;
        grad_f[1] = dfvals * yy / rho;
        grad_f[2] = 0.0;
        double phi_Re = sqrt(2.0) * m_params.pot_eta * costh * fvals;
        double phi_Im = sqrt(2.0) * m_params.pot_eta * sinth * fvals;

        Tensor<1, data_t> dphi_Re;
        dphi_Re[0] = yy / rho / rho * phi_Im +
                     grad_f[0] * sqrt(2.0) * m_params.pot_eta * xx / rho;
        dphi_Re[2] = 0;
        dphi_Re[1] = -xx / rho * phi_Im +
                     grad_f[1] * sqrt(2.0) * m_params.pot_eta * xx / rho;

        Tensor<1, data_t> dphi_Im;
        dphi_Im[0] = -yy / rho / rho * phi_Re +
                     grad_f[0] * sqrt(2.0) * m_params.pot_eta * yy / rho;
        dphi_Im[2] = 0;
        dphi_Im[1] = xx / rho / rho * phi_Re +
                     grad_f[1] * sqrt(2.0) * m_params.pot_eta * yy / rho;

        // set the field vars
        vars.phi_Re = phi_Re;
        vars.phi_Im = phi_Im;

        // // data_t gamma_factor = 1.0/sqrt(1.0 - m_params.velocity*m_params.velocity);
        vars.Pi_Re = 0.0;//gamma_factor * m_params.velocity / metric_vars.lapse * dphi_Re[0];
        vars.Pi_Im = 0.0;//gamma_factor * m_params.velocity / metric_vars.lapse * dphi_Im[0];
        FOR1(i)
        {
            vars.Pi_Re += -metric_vars.shift[i] * dphi_Re[i] / metric_vars.lapse;
            vars.Pi_Im += -metric_vars.shift[i] * dphi_Im[i] / metric_vars.lapse;
        }

        current_cell.store_vars(vars);
    }

  protected:
    const params_t m_params;
    const double m_dx;
};

#endif /* INITIALSCALARDATA_HPP_ */
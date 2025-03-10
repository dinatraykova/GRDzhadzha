/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CIRCULATION_HPP_
#define CIRCULATION_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the momentum flux S_i with type matter_t and writes it to the
//! grid, see https://arxiv.org/pdf/2104.13420.pdf for details
template <class matter_t, class background_t> class Circulation
{
    // Use the variable definition in the matter class
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const matter_t m_matter;                        //!< The matter object
    const double m_dx;                              //!< The grid spacing
    const background_t m_background;                //!< The metric background
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center

  public:
    Circulation(matter_t a_matter, background_t a_background, double a_dx,
                std::array<double, CH_SPACEDIM> a_center)
        : m_matter(a_matter), m_deriv(a_dx), m_dx(a_dx),
          m_background(a_background), m_center(a_center)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and derivs
        const auto vars = current_cell.template load_vars<MatterVars>();
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);

        // get the metric vars from the background
        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, coords);

        // some useful quantities
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);
        const emtensor_t<data_t> emtensor = m_matter.compute_emtensor(
            vars, metric_vars, d1, gamma_UU, chris_phys.ULL);
        const data_t R = coords.get_radius();
        data_t rho2 =
            simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
        data_t rho = sqrt(rho2);

        data_t phi2 = (vars.phi_Re * vars.phi_Re + vars.phi_Im * vars.phi_Im);
        Tensor<1, data_t> dl;
        dl[0] = -coords.y / rho;
        dl[1] = coords.x / rho;
        dl[2] = 0;

        Tensor<1, data_t> SiP;
        FOR1(i)
        {
            SiP[i] = vars.phi_Re * d1.phi_Im[i] - vars.phi_Im * d1.phi_Re[i];
        }

        data_t circ = 0.;
        FOR2(i, j) { circ += gamma_UU[i][j] * dl[i] * SiP[j]; }
        circ = circ / phi2;

        current_cell.store_vars(circ, c_circ);
    }
};

#endif /* CIRCULATION_HPP_ */

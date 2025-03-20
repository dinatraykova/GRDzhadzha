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
    const matter_t m_matter;                         //!< The matter object
    const double m_dx;                               //!< The grid spacing
    const background_t m_background;                 //!< The metric background
    const std::array<double, CH_SPACEDIM> m_center1; //!< Circle center 1
    const std::array<double, CH_SPACEDIM> m_center2; //!< Circle center 2
    const std::array<double, CH_SPACEDIM> m_center3; //!< Circle center 3

  public:
    Circulation(matter_t a_matter, background_t a_background, double a_dx,
                std::array<double, CH_SPACEDIM> a_center1,
                std::array<double, CH_SPACEDIM> a_center2,
                std::array<double, CH_SPACEDIM> a_center3)
        : m_matter(a_matter), m_deriv(a_dx), m_dx(a_dx),
          m_background(a_background), m_center1(a_center1),
          m_center2(a_center2), m_center3(a_center3)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and derivs
        const auto vars = current_cell.template load_vars<MatterVars>();
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);

        // define coords for each circle
        Coordinates<data_t> coords1(current_cell, m_dx, m_center1);
        Coordinates<data_t> coords2(current_cell, m_dx, m_center2);
        Coordinates<data_t> coords3(current_cell, m_dx, m_center3);

        // some useful quantities
        using namespace TensorAlgebra;
        data_t rho1 =
            sqrt(simd_max(coords1.x * coords1.x + coords1.y * coords1.y, 1e-8));
        data_t rho2 =
            sqrt(simd_max(coords2.x * coords2.x + coords2.y * coords2.y, 1e-8));
        data_t rho3 =
            sqrt(simd_max(coords3.x * coords3.x + coords3.y * coords3.y, 1e-8));

        data_t phi2 = simd_max(
            vars.phi_Re * vars.phi_Re + vars.phi_Im * vars.phi_Im, 1.e-2);

        Tensor<1, data_t> dl1;
        dl1[0] = -coords1.y / rho1;
        dl1[1] = coords1.x / rho1;
        dl1[2] = 0;

        Tensor<1, data_t> dl2;
        dl2[0] = -coords2.y / rho2;
        dl2[1] = coords2.x / rho2;
        dl2[2] = 0;

        Tensor<1, data_t> dl3;
        dl3[0] = -coords3.y / rho3;
        dl3[1] = coords3.x / rho3;
        dl3[2] = 0;

        Tensor<1, data_t> SiP;
        FOR1(i)
        {
            SiP[i] = vars.phi_Re * d1.phi_Im[i] - vars.phi_Im * d1.phi_Re[i];
        }

        data_t circ1 = 0.;
        FOR2(i, j) { circ1 += delta(i, j) * dl1[i] * SiP[j]; }
        circ1 = circ1 / phi2;

        data_t circ2 = 0.;
        FOR2(i, j) { circ2 += delta(i, j) * dl2[i] * SiP[j]; }
        circ2 = circ2 / phi2;

        data_t circ3 = 0.;
        FOR2(i, j) { circ3 += delta(i, j) * dl3[i] * SiP[j]; }
        circ3 = circ3 / phi2;

        current_cell.store_vars(circ1, c_circ1);
        current_cell.store_vars(circ2, c_circ2);
        current_cell.store_vars(circ3, c_circ3);
    }
};

#endif /* CIRCULATION_HPP_ */

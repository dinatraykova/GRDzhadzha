/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDGRIDSTAGGINGCRITERION_HPP_
#define FIXEDGRIDSTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

#include <iostream>
using namespace std;

// template <class matter_t, class background_t> class
// FixedGridsTaggingCriterion
template <class matter_t> class FixedGridsTaggingCriterion
{
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

  protected:
    const FourthOrderDerivatives m_deriv;
    const matter_t m_matter;
    const double m_dx;
    //    const background_t m_background;
    const double m_L;
    const int m_level;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_d_to_bh;
    const double m_velocity;
    const double m_time;
    const int m_max_vortex_lvl;
    const double m_vortex_refine_threshold;

  public:
    FixedGridsTaggingCriterion(
        const matter_t a_matter, const double a_dx,
        //			       background_t a_background,
        const int a_level, const double a_L,
        const std::array<double, CH_SPACEDIM> a_center, const double a_d_to_bh,
        const double a_velocity, const double a_time,
        const int a_max_vortex_lvl, const double a_vortex_refine_threshold)
        : m_matter(a_matter), m_deriv(a_dx),
          m_dx(a_dx), // m_background(a_background),
          m_L(a_L), m_level(a_level), m_center(a_center), m_d_to_bh(a_d_to_bh),
          m_velocity(a_velocity), m_time(a_time),
          m_max_vortex_lvl(a_max_vortex_lvl),
          m_vortex_refine_threshold(a_vortex_refine_threshold){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<MatterVars>();
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);

        data_t criterion = 0.0;
        // make sure the inner part is regridded around the horizon
        // take L as the length of full grid, so tag inner 1/2
        // of it, which means inner \pm L/4
        double ratio = pow(2.0, -(m_level + 2.0));
        double ratio_string = pow(2.0, -(m_level + 2.0));
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        const data_t max_abs_xy = simd_max(abs(coords.x), abs(coords.y));
        const data_t max_abs_xyz_bh = simd_max(max_abs_xy, abs(coords.z));
        const data_t max_abs_xyz_string = simd_max(
            abs(coords.x - (m_d_to_bh - m_velocity * m_time)), abs(coords.y));

        auto regrid_bh = simd_compare_lt(max_abs_xyz_bh, m_L * ratio);
        auto regrid_string =
            simd_compare_lt(max_abs_xyz_string, m_L * ratio_string);

        data_t gradPhi2 = 0.;
        FOR1(i)
        {
            gradPhi2 +=
                (d1.phi_Re[i] * d1.phi_Re[i] + d1.phi_Im[i] * d1.phi_Im[i]);
        }
        criterion = gradPhi2 / m_vortex_refine_threshold;

        // Refine if inside either window or if the gradient is large enough
        const auto force_refine = (regrid_bh | regrid_string);
        criterion = simd_conditional(force_refine, 100.0, criterion);
        // criterion = simd_conditional(regrid_bh, 100.0, criterion);

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* FIXEDGRIDSTAGGINGCRITERION_HPP_ */

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

class FixedGridsTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_L;
    const int m_level;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_d_to_bh;
    const double m_velocity;
    const double m_time;

  public:
    FixedGridsTaggingCriterion(const double dx, const int a_level,
                               const double a_L,
                               const std::array<double, CH_SPACEDIM> a_center,
                               const double a_d_to_bh, const double a_velocity, const double a_time)
        : m_dx(dx), m_L(a_L), m_level(a_level), m_center(a_center),
          m_d_to_bh(a_d_to_bh), m_velocity(a_velocity), m_time(a_time){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        data_t criterion = 0.0;
        // make sure the inner part is regridded around the horizon
        // take L as the length of full grid, so tag inner 1/2
        // of it, which means inner \pm L/4
        double ratio = pow(2.0, -(m_level + 2.0));
        double ratio_string = pow(2.0, -(m_level + 2.0));
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        const data_t max_abs_xy = simd_max(abs(coords.x), abs(coords.y));
        const data_t max_abs_xyz_bh = simd_max(max_abs_xy, abs(coords.z));
        const data_t max_abs_xyz_string =
            simd_max(abs(coords.x - (m_d_to_bh - m_velocity * m_time)), abs(coords.z));

        auto regrid_bh = simd_compare_lt(max_abs_xyz_bh, m_L * ratio);
        auto regrid_string =
            simd_compare_lt(max_abs_xyz_string, m_L * ratio_string);

/*        Tensor<1, data_t> d1_phi_Re, d1_phi_Im;
        FOR(idir) m_deriv.diff1(d1_phi_Re, current_cell, idir, c_phi_Re);
        FOR(idir) m_deriv.diff1(d1_phi_Im, current_cell, idir, c_phi_Im);

        data_t mod_d1_phi = 0;
        FOR(idir)
        {
            mod_d1_phi += d1_phi_Re[idir] * d1_phi_Re[idir] + d1_phi_Im[idir] * d1_phi_Im[idir];
        }

        criterion = m_dx * (sqrt(mod_d1_phi) / 0.05);
*/
        auto phi_Re = current_cell.load_vars(c_phi_Re);
        auto phi_Im = current_cell.load_vars(c_phi_Im);
	if (m_level < 3){
	criterion = 4e-2/sqrt(phi_Re * phi_Re + phi_Im * phi_Im);}
        criterion = simd_conditional(regrid_bh, 100.0, criterion);
//        criterion = simd_conditional(regrid_string, 100.0, criterion);

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* FIXEDGRIDSTAGGINGCRITERION_HPP_ */
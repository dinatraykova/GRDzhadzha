/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MATTEREVOLUTION_HPP_
#define MATTEREVOLUTION_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//!  Calculates RHS of matter variables only, metric vars assumed analytic
/*!
     The class calculates the RHS evolution for the matter variables.
     It does not assume a specific form of matter or background -
     it is templated over a matter class matter_t, and over a background metric,
     background_t
*/

template <class matter_t, class background_t> class MatterEvolution
{
  public:
    //! Inherit the variable definitions from the Matter vars
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    //  Need d2 of certain matter vars
    template <class data_t>
    using MatterDiff2Vars = typename matter_t::template Diff2Vars<data_t>;

    // This is used for the non evolved ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    //!  Constructor of class MatterEvolution
    MatterEvolution(
        matter_t a_matter, background_t a_background, double sigma, double dx,
        double a_L, std::array<double, CH_SPACEDIM> a_center,
        const InitialScalarData<BoostedMink>::params_t a_initial_params,
        double tau_target, double Nlayer)
        : m_matter(a_matter), m_background(a_background), m_sigma(sigma),
          m_deriv(dx), m_dx(dx), m_L(a_L), m_center(a_center),
          m_initial_params(a_initial_params), m_tau_target(tau_target),
          m_Nlayer(Nlayer)
    {
    }

    //!  The compute member which calculates the RHS at each point in the box
    //!  \sa matter_rhs_equation()
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy matter data from chombo gridpoint into local variable
        const auto matter_vars = current_cell.template load_vars<MatterVars>();

        // compute the background metric vars
        MetricVars<data_t> metric_vars;
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        m_background.compute_metric_background(metric_vars, coords);

        // compute derivs for matter grid vars
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);
        const auto d2 = m_deriv.template diff2<MatterDiff2Vars>(current_cell);
        const auto advec = m_deriv.template advection<MatterVars>(
            current_cell, metric_vars.shift);

        // the RHS
        MatterVars<data_t> matter_rhs;
        VarsTools::assign(matter_rhs, 0.); // All components set to zero

        // add evolution of matter fields and dissipation
        m_matter.matter_rhs(matter_rhs, matter_vars, metric_vars, d1, d2,
                            advec);

        // ---- Damping at x/y face ----
        const double Lx = m_L;
        const double Ly = m_L;
        const double Lz = m_L / 2.;

        const data_t xx = coords.x + m_center[0];
        const double yy = coords.y + m_center[1];

        // distances to faces
        const data_t dxL = simd_max(0.0, xx);
        const data_t dxR = simd_max(0.0, Lx - xx);
        const data_t dyL = simd_max(0.0, yy);
        const data_t dyR = simd_max(0.0, Ly - yy);

        // nearest distance to any face
        const data_t d_near = simd_min(simd_min(dxL, dxR), simd_min(dyL, dyR));

        const double d0 = m_Nlayer * m_dx;      // start of damping
        const double W = m_Nlayer * m_dx;       // thickness of layer
        const double tau_target = m_tau_target; // int_sigma ds across layer

        // smooth step that ramps up ONLY when approaching the boundary
        auto clamp01 = [](data_t x) { return simd_min(1.0, simd_max(0.0, x)); };
        // u = 0 for d_near >= d0 (no damping), u → 1 as d_near -> 0 (toward the
        // face)
        data_t u = clamp01((d0 - d_near) / W);
        auto smooth01 = [](data_t s) { return s * s * (3.0 - 2.0 * s); }; // C1
        const data_t ramp = smooth01(u);

        // keep total optical depth ~ constant when changing W
        const double sigma_peak = tau_target / (W + 1e-14);
        const data_t sigma_iso = sigma_peak * ramp; // isotropic scalar strength

        // shift-aware factors (Minkowski: α=1, β=(v,0))
        auto s = [](double lam) { return std::max(0.0, std::min(1.0, lam)); };
        const double vx = 0.5;
        const double f_xL = s(1.0 + vx); // n=-x
        const double f_xR = s(1.0 - vx); // n=+x
        const double f_y = 1.0;          // β^y=0

        // outgoing characteristics (c=1)
        const data_t wpx_Re_xL = matter_vars.Pi_Re - d1.phi_Re[0];
        const data_t wpx_Re_xR = matter_vars.Pi_Re + d1.phi_Re[0];
        const data_t wpy_Re_yL = matter_vars.Pi_Re - d1.phi_Re[1];
        const data_t wpy_Re_yR = matter_vars.Pi_Re + d1.phi_Re[1];

        const data_t wpx_Im_xL = matter_vars.Pi_Im - d1.phi_Im[0];
        const data_t wpx_Im_xR = matter_vars.Pi_Im + d1.phi_Im[0];
        const data_t wpy_Im_yL = matter_vars.Pi_Im - d1.phi_Im[1];
        const data_t wpy_Im_yR = matter_vars.Pi_Im + d1.phi_Im[1];

        // combine opposite faces (balances x/y), then damp with sigma_iso
        const data_t wdir_Re_x = 0.5 * (f_xL * wpx_Re_xL + f_xR * wpx_Re_xR);
        const data_t wdir_Re_y = 0.5 * (f_y * wpy_Re_yL + f_y * wpy_Re_yR);
        const data_t wdir_Im_x = 0.5 * (f_xL * wpx_Im_xL + f_xR * wpx_Im_xR);
        const data_t wdir_Im_y = 0.5 * (f_y * wpy_Im_yL + f_y * wpy_Im_yR);

        // Π-only damping (isotropic strength)
        matter_rhs.Pi_Re += -sigma_iso * (wdir_Re_x + wdir_Re_y);
        matter_rhs.Pi_Im += -sigma_iso * (wdir_Im_x + wdir_Im_y);

        // ---- End damping layer ----

        m_deriv.add_dissipation(matter_rhs, current_cell, m_sigma);

        // Write the rhs into the output vars for this cell
        current_cell.store_vars(matter_rhs);
    }

  protected:
    const matter_t m_matter;              //!< The matter object
    const background_t m_background;      //!< The metric background
    const FourthOrderDerivatives m_deriv; //!< An object for calculating
                                          //!< derivatives of the vars
    const double m_sigma;                 //!< Sigma for dissipation
    const double m_dx;                    //!< Grid spacing
    const double m_L;                     //!< Grid length
    const std::array<double, CH_SPACEDIM> m_center; //!< Grid center
    const InitialScalarData<BoostedMink>::params_t m_initial_params;
    const double m_tau_target;
    const double m_Nlayer;
};

#endif /* MATTEREVOLUTION_HPP_ */

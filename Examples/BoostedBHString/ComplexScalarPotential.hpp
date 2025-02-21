/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXSCALARPOTENTIAL_HPP_
#define COMPLEXSCALARPOTENTIAL_HPP_

#include "simd.hpp"

class ComplexScalarPotential
{
  protected:
    //    const double m_mu;
    const InitialScalarData<BoostedBH>::params_t m_initial_params;

  public:
    //! The constructor
    ComplexScalarPotential(
        const InitialScalarData<BoostedBH>::params_t a_initial_params)
        : m_initial_params(a_initial_params)
    {
    }

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi_re,
                           data_t &dVdphi_im, const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        const double lambda = m_initial_params.pot_lambda;
        const double eta = m_initial_params.pot_eta;

        V_of_phi = 0.25 * lambda *
                   pow(0.5 * (pow(vars.phi_Re, 2.0) + pow(vars.phi_Im, 2.0)) -
                           pow(eta, 2.0),
                       2.0);

        dVdphi_re = 0.5 * lambda *
                    (0.5 * (pow(vars.phi_Re, 2.0) + pow(vars.phi_Im, 2.0)) -
                     pow(eta, 2.0)) *
                    vars.phi_Re;
        dVdphi_im = 0.5 * lambda *
                    (0.5 * (pow(vars.phi_Re, 2.0) + pow(vars.phi_Im, 2.0)) -
                     pow(eta, 2.0)) *
                    vars.phi_Im;
    }
};

#endif /* COMPLEXSCALARPOTENTIAL_HPP_ */
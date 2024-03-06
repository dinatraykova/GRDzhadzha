/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "FixedBGSimulationParametersBase.hpp"
#include "GRParmParse.hpp"

// Problem specific includes:
#include "BoostedBH.hpp"
#include "InitialProcaData.hpp"

class SimulationParameters : public FixedBGSimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : FixedBGSimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        // Initial PF
        pp.load("proca_mass", proca_mass);
        pp.load("field_amplitude", proca_amplitude);
        pp.load("proca_damping", proca_damping);

        // BH data
        pp.load("bh_mass", bg_params.mass, 1.0);
        pp.load("bh_velocity", bg_params.velocity, 0.0);
        pp.load("bh_center", bg_params.center, center);

        // Volume extraction radii
        pp.load("inner_r", inner_r, extraction_params.extraction_radii[0]);
        pp.load("outer_r", outer_r, extraction_params.extraction_radii[1]);
    }

    void check_params()
    {
        warn_parameter("proca_mass", proca_mass,
                       proca_mass < 0.2 / coarsest_dx / dt_multiplier,
                       "oscillations of proca field do not appear to be "
                       "resolved on coarsest level");
        warn_parameter("inner_r", inner_r,
                       extraction_params.extraction_radii[0] == inner_r,
                       "should be equal to first extraction radius");
        warn_parameter("outer_r", outer_r,
                       extraction_params.extraction_radii[1] == outer_r,
                       "should be equal to second extraction radius");
        warn_parameter("bh_mass", bg_params.mass, bg_params.mass >= 0.0,
                       "should be >= 0.0");
        FOR(idir)
        {
            std::string name = "bh_center[" + std::to_string(idir) + "]";
            warn_parameter(
                name, bg_params.center[idir],
                (bg_params.center[idir] >= 0) &&
                    (bg_params.center[idir] <= (ivN[idir] + 1) * coarsest_dx),
                "should be within the computational domain");
        }
    }

    // Problem specific parameters - the radii for the integrations
    double inner_r, outer_r;
    double proca_mass, proca_damping, proca_amplitude;
    // Collection of parameters necessary for the initial conditions
    // InitialProcaData::params_t initial_params;
    // Collection of parameters necessary for the metric background
    BoostedBH::params_t bg_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

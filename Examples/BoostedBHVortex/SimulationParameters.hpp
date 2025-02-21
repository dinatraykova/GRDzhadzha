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
#include "InitialScalarData.hpp"

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
        // Initial SF
        pp.load("scalar_amplitude", initial_params.amplitude, 0.1);
        pp.load("scalar_mass", initial_params.mass, 0.5);
        pp.load("scalar_omega", initial_params.omega, 0.5);
        pp.load("scalar_center", initial_params.center, center);

        // BH data
        pp.load("bh_mass", bg_params.mass, 1.0);
        pp.load("bh_velocity", bg_params.velocity, 0.0);
        pp.load("bh_center", bg_params.center, center);

        // Circle extraction params
        pp.load("lineout_num_points", lineout_num_points, 36);
        pp.load("r_circle", r_circle, 20.);
        pp.load("circle_center", circle_center, center);
    }

    void check_params()
    {
        warn_parameter("scalar_mass", initial_params.mass,
                       initial_params.mass < 0.2 / coarsest_dx / dt_multiplier,
                       "oscillations of scalar field do not appear to be "
                       "resolved on coarsest level");
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

    // Collection of parameters necessary for the initial conditions
    InitialScalarData::params_t initial_params;
    // Collection of parameters necessary for the metric background
    BoostedBH::params_t bg_params;
    // Problem specific parameters
    int lineout_num_points;
    double r_circle;
    std::array<double, CH_SPACEDIM> circle_center;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

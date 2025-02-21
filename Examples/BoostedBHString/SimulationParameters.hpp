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
#include <fstream>

class SimulationParameters : public FixedBGSimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : FixedBGSimulationParametersBase(pp)
    {
        pout() << "A" << endl;

        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        // Initial SF
        pp.load("pot_eta", initial_params.pot_eta);
        pp.load("pot_lambda", initial_params.pot_lambda);
        pp.load("d_to_bh", initial_params.d_to_bh);

        // BH data
        pp.load("bh_mass", bg_params.mass, 1.0);
        pp.load("bh_velocity", bg_params.velocity, 0.0);
        pp.load("bh_center", bg_params.center, center);
        initial_params.center = bg_params.center;
        initial_params.velocity = bg_params.velocity;

        // Reading data
        int lines;
        pp.load("lines", lines);
        pp.load("spacing", initial_params.spacing);
        double f_profile[lines];
        double df_profile[lines];

        double tmp_data;
        ifstream read_file("f_profile.csv");

        for (int i = 0; i < lines; ++i)
        {
            read_file >> tmp_data;
            f_profile[i] = tmp_data;
        }

        for (int i = 0; i < lines - 1; ++i)
        {
            if (i == 0)
            {
                df_profile[i] = 0.;
            }
            else
            {
                df_profile[i] = pow(2. * initial_params.spacing, -1.0) *
                                (f_profile[i + 1] - f_profile[i - 1]);
            }
        }
        initial_params.f_profile = f_profile;
        initial_params.df_profile = df_profile;
    }

    void check_params()
    {
        // warn_parameter("scalar_mass", initial_params.mass,
        //                initial_params.mass < 0.2 / coarsest_dx / dt_multiplier,
        //                "oscillations of scalar field do not appear to be "
        //                "resolved on coarsest level");
        // warn_parameter("bh_mass", bg_params.mass, bg_params.mass >= 0.0,
        //                "should be >= 0.0");
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
    InitialScalarData<BoostedBH>::params_t initial_params;
    // Collection of parameters necessary for the metric background
    BoostedBH::params_t bg_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_chi,
    c_circ1,
    c_circ2,
    c_circ3,
    c_rhoEnergy,
    c_rhoLinMom,
    c_rhoLinMomY,
    c_sourceLinMom,
    c_sourceLinMomY,
    c_fluxEnergy,
    c_fluxEnergyY,
    c_fluxLinMom,
    c_fluxLinMomY,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi",         "circ1",      "circ2",        "circ3",         "rhoEnergy",
    "rhoLinMom",   "rhoLinMomY", "sourceLinMom", "sourceLinMomY", "fluxEnergy",
    "fluxEnergyY", "fluxLinMom", "fluxLinMomY"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */

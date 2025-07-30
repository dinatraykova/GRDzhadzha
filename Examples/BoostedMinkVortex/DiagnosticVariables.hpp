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
    c_circ4,
    c_rhoEnergy,
    c_rhoParticle,
    c_rhoDensity,
    c_rhoLinMom,
    c_rhoLinMomY,
    c_sourceLinMom,
    c_sourceLinMomY,
    c_fluxLinMom,
    c_fluxLinMomY,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi",          "circ1",         "circ2",      "circ3",      "circ4",
    "rhoEnergy",    "rhoParticle",   "rhoDensity", "rhoLinMom",  "rhoLinMomY",
    "sourceLinMom", "sourceLinMomY", "fluxLinMom", "fluxLinMomY"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */

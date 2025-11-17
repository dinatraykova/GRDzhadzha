/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "BoostedBHScalarLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"

// For RHS update
#include "BoostedMink.hpp"
#include "MatterEvolution.hpp"

// For tag cells
// #include "FixedGridsTaggingCriterion.hpp"
#include "FixedGridsTaggingCriterionNew.hpp"

// Problem specific includes
#include "Circulation.hpp"
#include "ComplexScalarField.hpp"
#include "ComplexScalarPotential.hpp"
#include "CustomExtraction.hpp"
// #include "EnergyConservation.hpp"
#include "Densities.hpp"
#include "ExcisionDiagnostics.hpp"
#include "ExcisionEvolution.hpp"
#include "FluxExtraction.hpp"
#include "InitialScalarData.hpp"
#include "LinearMomConservation.hpp"
#include "LinearMomConservationY.hpp"

// Initial data for field and metric variables
void BoostedBHScalarLevel::initialData()
{
    CH_TIME("BoostedBHScalarLevel::initialData");
    if (m_verbosity)
        pout() << "BoostedBHScalarLevel::initialData " << m_level << endl;

    // First set everything to zero, then set the value of the conformal factor
    // This is just for the diagnostics
    SetValue set_zero(0.0);
    BoostedMink boosted_bh(m_p.bg_params, m_dx); // just calculates chi
    auto compute_pack = make_compute_pack(set_zero, boosted_bh);
    BoxLoops::loop(compute_pack, m_state_diagnostics, m_state_diagnostics,
                   SKIP_GHOST_CELLS);

    // Now set the actual evolution variables
    InitialScalarData<BoostedMink> initial_sf(m_p.initial_params, m_dx,
                                              boosted_bh);
    BoxLoops::loop(initial_sf, m_state_new, m_state_new, FILL_GHOST_CELLS,
                   disable_simd());

    // excise evolution vars within horizon, turn off simd vectorisation
    // BoxLoops::loop(ExcisionEvolution<ScalarFieldWithPotential, BoostedMink>(
    //									    m_dx,
    //m_p.bg_params.center,
    // boosted_bh), 		   m_state_new, m_state_new, SKIP_GHOST_CELLS,
    // disable_simd());
}

void BoostedBHScalarLevel::specificPostTimeStep()
{
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd());

    // At any level, but after the timestep on the minimum extraction level
    int min_level = 0;
    if (m_p.activate_extraction == 1)
        min_level = m_p.extraction_params.min_extraction_level();

    bool calculate_diagnostics = at_level_timestep_multiple(min_level);
    if (calculate_diagnostics)
    {
        fillAllGhosts();
        ComplexScalarPotential potential(m_p.initial_params);
        ScalarFieldWithPotential scalar_field(potential);
        BoostedMink boosted_bh(m_p.bg_params, m_dx);
        Densities<ScalarFieldWithPotential, BoostedMink> densities(
            scalar_field, boosted_bh, m_dx, m_p.center);
        int direction = 0; // we want the x direction for the momentum
        LinearMomConservation<ScalarFieldWithPotential, BoostedMink>
            linear_momenta(scalar_field, boosted_bh, direction, m_dx,
                           m_p.center);
        LinearMomConservationY<ScalarFieldWithPotential, BoostedMink>
            linear_momenta_y(scalar_field, boosted_bh, 1, m_dx, m_p.center);
        Circulation<ScalarFieldWithPotential, BoostedMink> circulation(
            scalar_field, boosted_bh, m_dx, m_p.circle1_center,
            m_p.circle2_center, m_p.circle3_center, m_p.circle4_center);
        BoxLoops::loop(make_compute_pack(densities, linear_momenta,
                                         linear_momenta_y, circulation),
                       m_state_new, m_state_diagnostics, SKIP_GHOST_CELLS);

        // excise within/outside specified radii, no simd
        if (m_p.activate_excision == 1)
        {
            BoxLoops::loop(
                ExcisionDiagnostics<ScalarFieldWithPotential, BoostedMink>(
                    m_dx, m_p.center, boosted_bh, m_p.inner_r, m_p.outer_r),
                m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
                disable_simd());
        }
    }
}

// Things to do in RHS update, at each RK4 step
void BoostedBHScalarLevel::specificEvalRHS(GRLevelData &a_soln,
                                           GRLevelData &a_rhs,
                                           const double a_time)
{
    // Calculate right hand side with matter_t = ScalarField
    // and background_t = BoostedBH
    ComplexScalarPotential potential(m_p.initial_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoostedMink boosted_bh(m_p.bg_params, m_dx);
    MatterEvolution<ScalarFieldWithPotential, BoostedMink> my_evolution(
        scalar_field, boosted_bh, m_p.sigma, m_dx, m_p.L, m_p.center,
        m_p.initial_params, m_p.tau_target, m_p.Nlayer);
    // MatterEvolution<ScalarFieldWithPotential, BoostedMink> my_evolution(
    //    scalar_field, boosted_bh, m_p.sigma, m_dx, m_p.center);
    a_soln.exchange(); // MPI halos
    // boundary_conditions.fill_exc(a_soln); // your existing BC hook
    BoxLoops::loop(my_evolution, a_soln, a_rhs, SKIP_GHOST_CELLS);
}

// Note that for the fixed grids this only happens on the initial timestep
void BoostedBHScalarLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state)
{
    ComplexScalarPotential potential(m_p.initial_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoostedMink boosted_bh(m_p.bg_params, m_dx);

    FixedGridsTaggingCriterion<ScalarFieldWithPotential> my_tagging(
        scalar_field, m_dx, m_level, m_p.L, m_p.center,
        m_p.initial_params.d_to_bh, m_p.bg_params.velocity, m_time,
        m_p.max_vortex_lvl, m_p.vortex_refine_threshold);
    BoxLoops::loop(my_tagging, current_state, tagging_criterion);
}
/*void BoostedBHScalarLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state)
{
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.L, m_p.center),
                   current_state, tagging_criterion);
                   }*/

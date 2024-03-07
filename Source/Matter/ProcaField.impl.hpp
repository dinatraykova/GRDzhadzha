/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(PROCAFIELD_HPP_)
#error "This file should only be included through ProcaField.hpp"
#endif

#ifndef PROCAFIELD_IMPL_HPP_
#define PROCAFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> ProcaField::compute_emtensor(
    const vars_t<data_t> &vars, const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<2, data_t> &gamma_UU,
    const Tensor<3, data_t> &chris_phys_ULL) const
{
    emtensor_t<data_t> out;

    // Some useful quantities
    const double msquared = pow(m_proca_mass, 2.0);

    // D_i A_j
    Tensor<2, data_t> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] += -chris_phys_ULL[k][i][j] * vars.Avec[k]; }
    }

    // D_i A_j - D_j A_i (NB exterior derivative, so christoffel symbols cancel)
    Tensor<2, data_t> diff_DA;
    FOR2(i, j) { diff_DA[i][j] = d1.Avec[j][i] - d1.Avec[i][j]; }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i, j)
    {
        out.Sij[i][j] = msquared * (vars.Avec[i] * vars.Avec[j] +
                                    0.5 * metric_vars.gamma[i][j] * vars.Avec0 *
                                        vars.Avec0);

        FOR2(k, l)
        {
            out.Sij[i][j] += gamma_UU[k][l] * (diff_DA[i][k] * diff_DA[j][l]) -
                             metric_vars.gamma[i][k] * metric_vars.gamma[j][l] *
                                 vars.Evec[k] * vars.Evec[l] +
                             0.5 * metric_vars.gamma[k][l] *
                                 metric_vars.gamma[i][j] * vars.Evec[k] *
                                 vars.Evec[l] -
                             0.5 * gamma_UU[k][l] * metric_vars.gamma[i][j] *
                                 msquared * vars.Avec[k] * vars.Avec[l];
            FOR2(m, n)
            {
                out.Sij[i][j] += -0.5 * metric_vars.gamma[i][j] *
                                 gamma_UU[k][m] * gamma_UU[l][n] * DA[k][l] *
                                 diff_DA[m][n];
            }
        }
    }

    // S = Tr_S_ij
    out.S = 0.0;
    FOR2(i, j) { out.S += out.Sij[i][j] * gamma_UU[i][j]; }

    // S_i (note lower index) = n^a T_a0
    FOR1(i)
    {
        out.Si[i] = msquared * vars.Avec0 * vars.Avec[i];

        FOR1(j) { out.Si[i] += vars.Evec[j] * diff_DA[i][j]; }
    }

    // rho = n^a n^b T_ab
    out.rho = 0.5 * msquared * (vars.Avec0 * vars.Avec0);
    FOR2(i, j)
    {
        out.rho +=
            0.5 * metric_vars.gamma[i][j] * vars.Evec[i] * vars.Evec[j] +
            0.5 * gamma_UU[i][j] * msquared * vars.Avec[i] * vars.Avec[j];

        FOR2(k, l)
        {
            // This is 0.5 * B_i B^i
            out.rho += 0.5 * gamma_UU[i][k] * gamma_UU[j][l] * DA[k][l] *
                       diff_DA[i][j];
        }
    }

    return out;
}

// Adds VF evolution to the RHS
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void ProcaField::matter_rhs(rhs_vars_t<data_t> &total_rhs,
                            const vars_t<data_t> &vars,
                            const MetricVars<data_t> &metric_vars,
                            const vars_t<Tensor<1, data_t>> &d1,
                            const diff2_vars_t<Tensor<2, data_t>> &d2,
                            const vars_t<data_t> &advec) const
{
    // calculate full spatial christoffel symbols
    using namespace TensorAlgebra;
    const auto gamma_UU = compute_inverse(metric_vars.gamma);
    const auto chris_phys = compute_christoffel(metric_vars.d1_gamma, gamma_UU);

    // evolution equations for vector fields A_0, A_i (note indices down) and
    // the conjugate momentum E^i (index up)
    total_rhs.Avec0 = metric_vars.lapse * metric_vars.K * vars.Avec0 +
                      advec.Avec0 - metric_vars.lapse * vars.Zvec;

    FOR2(i, j)
    {
        total_rhs.Avec0 +=
            -gamma_UU[i][j] * (vars.Avec[i] * metric_vars.d1_lapse[j] +
                               metric_vars.lapse * d1.Avec[i][j]);

        FOR1(k)
        {
            total_rhs.Avec0 += gamma_UU[i][j] * metric_vars.lapse *
                               chris_phys.ULL[k][i][j] * vars.Avec[k];
        }
    }

    FOR1(i)
    {
        total_rhs.Avec[i] = -metric_vars.lapse * d1.Avec0[i] -
                            vars.Avec0 * metric_vars.d1_lapse[i] +
                            advec.Avec[i];
        FOR1(j)
        {
            total_rhs.Avec[i] +=
                -metric_vars.lapse * metric_vars.gamma[i][j] * vars.Evec[j] +
                vars.Avec[j] * metric_vars.d1_shift[j][i];
        }
    }

    // variable for term (D_i A_j - D_j A_i)
    // NB Christoffel symbols cancel and take care with
    // indices - the second index is the derivative index
    Tensor<2, data_t> diff_DA;
    FOR2(i, j) { diff_DA[i][j] = d1.Avec[j][i] - d1.Avec[i][j]; }

    // NB This is for E^i with indices up
    FOR1(i)
    {
        total_rhs.Evec[i] =
            metric_vars.lapse * metric_vars.K * vars.Evec[i] + advec.Evec[i];
        FOR1(j)
        {
            total_rhs.Evec[i] +=
                gamma_UU[i][j] * (metric_vars.lapse * d1.Zvec[j] +
                                  metric_vars.lapse * pow(m_proca_mass, 2.0) *
                                      vars.Avec[j]) -
                vars.Evec[j] * metric_vars.d1_shift[i][j];
        }

        FOR3(j, k, l)
        {
            total_rhs.Evec[i] +=
                gamma_UU[k][j] * gamma_UU[i][l] *
                (metric_vars.d1_lapse[k] * diff_DA[l][j] +
                 metric_vars.lapse * (d2.Avec[j][l][k] - d2.Avec[l][j][k]));

            FOR1(m)
            {
                total_rhs.Evec[i] += -gamma_UU[k][j] * gamma_UU[i][l] *
                                     metric_vars.lapse *
                                     (chris_phys.ULL[m][k][l] * diff_DA[m][j] +
                                      chris_phys.ULL[m][k][j] * diff_DA[l][m]);
            }
        }
    }

    // evolution equations for vector field Avec and its conjugate momentum Evec
    total_rhs.Zvec = metric_vars.lapse * (pow(m_proca_mass, 2.0) * vars.Avec0 -
                                          m_proca_damping * vars.Zvec) +
                     advec.Zvec;

    FOR1(i)
    {
        total_rhs.Zvec += metric_vars.lapse * d1.Evec[i][i];
        FOR1(j)
        {
            total_rhs.Zvec +=
                metric_vars.lapse * chris_phys.ULL[i][i][j] * vars.Evec[j];
        }
    }
}

#endif /* PROCAFIELD_IMPL_HPP_ */

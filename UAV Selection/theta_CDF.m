function out = theta_CDF(gamma_val, gamma_hat_AU ,gamma_hat_UB, m_AU, m_UB)

out = 1 - (gammainc(m_AU, gamma_val./gamma_hat_AU).*gammainc(m_UB, gamma_val./gamma_hat_UB))./(gamma(m_AU).*gamma(m_UB));

end


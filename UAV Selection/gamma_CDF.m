function out = gamma_CDF(gam, gamma_hat_AU ,gamma_hat_UB, m_AU, m_UB)

    K = @(t, p) besselk(t-p+1, 2*sqrt((gam.^2)./(gamma_hat_AU.*gamma_hat_UB)));

    x = 2*exp(-gam*(1./gamma_hat_AU + 1./gamma_hat_UB));
    xx = 0;

    for p = 0: m_AU-1
        for t = 0: p + m_UB-1

            v = ((gam.^(p+m_UB))*(gamma_hat_UB.^((t-p-2*m_UB+1)/2)))./...
                ((gamma_hat_AU.^((p+t+1)/2)).*factorial(p).*gamma(m_UB));
            xx = xx + nchoosek(p+m_UB-1, t) .* v .* K(t, p);
%             disp(K(t, p))

        end
    end


    out = 1 - x .* xx;
    % out = xx;
    
end


function [c, c_heur, convex_factor] = getC(M_list, h, m, L, sigma_z, step_length)
    %Convergence constants
    K = length(h);
    n = length(M_list);
    c = zeros(n, 1);
    c_heur = zeros(n, 1);
    for i = 1:n
        M = M_list(i);
        [p, eta] = power_control(h, sigma_z, M);
        c(i) = 1-2*step_length(i)/(K*sqrt(eta))*m*L/(m+L)*sum(sqrt(p).*abs(h));
        c_heur(i) = 1-2*step_length(i)/(K*sqrt(eta))*0.5*sum(sqrt(p).*abs(h));
    end
    %Also gonna try with a convex proof instead of strongly convex
    convex_factor = zeros(n, 1);
    for i = 1:n
        M = M_list(i);
        [p, eta] = power_control(h, sigma_z, M);
        convex_factor(i) = K*sqrt(eta)/(2*step_length(i)*sum(sqrt(p).*abs(h)));
    end
end
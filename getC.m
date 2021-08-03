function [c] = getC(M_list, h, m, L, sigma_z, step_length)
    %Convergence constants
    K = length(h);
    n = length(M_list);
    c = zeros(n, 1);
    for i = 1:n
        M = M_list(i);
        [p, eta] = power_control(h, sigma_z, M);
        c(i) = 1-2*step_length/(K*sqrt(eta))*m*L/(m+L)*sum(sqrt(p).*abs(h));
    end
end
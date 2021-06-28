function [c] = optimalSchedule(budget, h, sigma_z, step_length, m, L)
    %Convergence constants
    num_devices = length(h);
    c = zeros(budget/10, 1);
    for i = 1:budget/10
        num_tx = i;
        [p, eta] = power_control(h, sigma_z, num_tx);
        c(i) = 1-2*step_length/(num_devices*sqrt(eta))*m*L/(m+L)*sum(sqrt(p).*abs(h));
    end
end
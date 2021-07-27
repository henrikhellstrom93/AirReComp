function [c, schedule] = optimalSchedule(budget, h, sigma_z, step_length, m, L, sigma_g, beta_g_init, true_beta)
    %Convergence constants
    num_devices = length(h);
    n = budget/10;
    d = length(beta_g_init);
    c = zeros(n, 1);
    eta = zeros(n, 1);
    for i = 1:n
        num_tx = i;
        [p, eta_i] = power_control(h, sigma_z, num_tx);
        eta(i) = eta_i;
        c(i) = 1-2*step_length/(num_devices*sqrt(eta_i))*m*L/(m+L)*sum(sqrt(p).*abs(h));
    end
    
    r_zero = norm(true_beta-beta_g_init)^2;
    
    %Optimization problem
    counter = (1:n)';

    cvx_begin quiet
        variable schedule(n)
        objective = L/2*r_zero;
        for i = 1:n
            objective = objective*c(i)^schedule(i);
            parenthesis = sigma_g'*sigma_g*sum(p.*abs(h).^2)+d*sigma_z^2/(i*num_devices);
            objective = objective + L*step_length^2/(2*num_devices*eta(i))*parenthesis/(1-c(i));
        end
        minimize( objective )
        subject to
            schedule'*counter <= budget
            schedule >= zeros(n,1)
    cvx_end
    
    schedule = round(schedule);
end
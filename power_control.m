%ASSUMES P_max and h have been sorted by quality indicator
function [p, eta] = power_control(h, sigma_z, M)
    n = length(h);
    
    P_max = ones(n, 1);
    
    %Find eta by taking the smallest eta_tilde
    eta_tilde = zeros(n,1);
    for k = 1:n
        sum1 = 0;
        sum2 = 0;
        for i = 1:k
            sum1 = sum1 + P_max(i)*abs(h(i))^2;
            sum2 = sum2 + sqrt(P_max(i))*abs(h(i));
        end
        eta_tilde(k) = ((sigma_z^2/M+sum1)/sum2)^2;
    end
    eta = min(eta_tilde);
    
    %Find p by channel inversion or max power
    p = zeros(n,1);
    for i = 1:n
        if P_max(i) > eta/abs(h(i))^2
            %Possible to invert channel
            p(i) = eta/abs(h(i))^2;
        else
            %Maximum power used
            p(i) = P_max(i);
        end
    end
end
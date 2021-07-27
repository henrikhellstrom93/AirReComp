%The schedule is a vector of length max_tx where entry i is how many
%communication rounds will have i uplink transmissions.
function [sigma_g] = calculateSigmaG(grad_len, num_rounds, X, Y, beta_g_init, step_length)
    sigma_gradient = zeros(grad_len, 1);
    sigma_g_samples = 0;
    beta_g = beta_g_init;
    data_size = size(Y);
    num_devices = data_size(2);
    for a = 1:num_rounds        
        %Least-squares loss gradient
        grad_l = getGradient(X, Y, beta_g);

        %Transmit gradients over MAC
        rcv_grad = sum(grad_l, 2)/4;
        
        % Model update
        beta_g = beta_g - step_length*rcv_grad;
        
        %For upper bound analysis
        for i = 1:num_devices
            sigma_g_samples = sigma_g_samples + 1;
            for j = 1:grad_len
                sigma_gradient(j) = sigma_gradient(j) + (rcv_grad(j)-grad_l(j,i))^2;
            end
        end
    end
    sigma_g = sigma_gradient/sigma_g_samples;
end
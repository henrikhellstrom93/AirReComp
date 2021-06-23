function [rcv_vec] = otaComputation(grad_l, sigma_z, num_tx)
    dimensions = size(grad_l);
    num_devices = dimensions(2);
    % Generate channel
    h = generateH(num_devices);
    %Normalize
    [grad_l_n, mu, sigma] = normalize(grad_l);

    grad_g = zeros(2, 1);
    %Generate power control
    [p, eta] = power_control(h, sigma_z, num_tx);
    for i = 1:num_tx
        for j = 1:2
            %Generate AWGN
            z = normrnd(0,sigma_z,1,1);
            %Apply fading and noise over MAC
            grad_g(j) = grad_g(j) + real(grad_l_n(j,:)*(abs(h).*sqrt(p)) + z)/(num_devices*sqrt(eta));
        end
    end
    grad_g = grad_g/num_tx;
    
    %Denormalize
    rcv_vec = denormalize(grad_g, mu, sigma);
end
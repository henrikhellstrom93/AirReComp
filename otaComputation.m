function [rcv_vec] = otaComputation(grad_l, sigma_z, num_tx, h)
    dimensions = size(grad_l);
    num_devices = dimensions(2);
    %Fix vector norm to support power control
    [grad_l_n, norms] = setVectorNorms(grad_l, sqrt(2));

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
    rcv_vec = grad_g/num_tx;
    
    %Lengthen to perserve original norm
    rcv_vec = rcv_vec*mean(norms)/sqrt(2);
end
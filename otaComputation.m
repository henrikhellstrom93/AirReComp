% Returns (d+1) x 1 vector
function [rcv_vec] = otaComputation(grad_l, sigma_z, num_tx, h)
    dims = size(grad_l);
    d = dims(1)-1;
    K = dims(2);
    new_norm = sqrt(K);
    %Fix vector norm to support power control
    [grad_l_n, norms] = setVectorNorms(grad_l, new_norm);

    grad_g = zeros(d+1, 1);
    %Generate power control
    [p, eta] = power_control(h, sigma_z, num_tx);
    for i = 1:num_tx
        for j = 1:d+1
            %Generate AWGN
            z = normrnd(0,sigma_z,1,1);
            %Apply fading and noise over MAC
            grad_g(j) = grad_g(j) + real(grad_l_n(j,:)*(abs(h).*sqrt(p)) + z)/(K*sqrt(eta));
        end
    end
    rcv_vec = grad_g/num_tx;
    
    %Change length to perserve original norm
    rcv_vec = rcv_vec*mean(norms)/new_norm;
end
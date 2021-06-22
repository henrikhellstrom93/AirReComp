function [grad_l] = getGradient(X, Y, beta_g)
    num_devices = min(size(X));
    grad_l = zeros(2, num_devices);
    grad_l(1,:) = 2*sum((beta_g(1)*X+beta_g(2)-Y).*X);
    grad_l(2,:) = 2*sum((beta_g(1)*X+beta_g(2)-Y));
end
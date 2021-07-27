% Returns (d+1) x K matrix consisting of K gradients of dimension (d+1)
function [grad_l] = getGradient(X_struct, y_struct, beta_g)
    dims = size(X_struct);
    D_k = dims(1);
    d = dims(2)-1;
    K = dims(3);
    grad_l = zeros(d+1, K);
    for k = 1:K
        X = X_struct(:,:,k);
        y = y_struct(:,k);
        grad_l(:,k) = 2*X'*(X*beta_g-y)/D_k;
    end
end
%X_struct is a 3-dimensional array where the first two elements correspond
%to a D_k x (d+1) matrix and the third element is a selector of the matrices
%y_struct is a D_k x K matrix
function [X_struct, y_struct] = splitDataset(X, y, K)
    %split dataset evenly between the devices
    dims = size(X);
    D = dims(1);
    d = dims(2)-1;
    D_k = floor(D/K);
    X_struct = zeros(D_k, d+1);
    y_struct = zeros(D_k, K);
    for k = 1:K
        X_struct(:,:,k) = X(1+(k-1)*D_k:k*D_k,:);
        y_struct(:,k) = y(1+(k-1)*D_k:k*D_k);
    end
end
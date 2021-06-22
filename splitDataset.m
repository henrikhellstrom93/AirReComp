function [X Y] = generateDataset(x, y, num_devices)
    %split dataset evenly between the devices
    n = length(x);
    n_k = n/num_devices;
    X = zeros(n_k, num_devices);
    Y = zeros(n_k, num_devices);
    for k = 1:num_devices
        X(:,k) = x((k-1)*n_k+1:k*n_k);
        Y(:,k) = y((k-1)*n_k+1:k*n_k);
    end
end
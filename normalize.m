function [vecs_n, mu, sigma] = normalize(vecs)
    dims = size(vecs);
    num_devices = dims(2);
    mu = zeros(num_devices, 1);
    sigma = zeros(num_devices, 1);
    for k = 1:num_devices
        mu(k) = mean(vecs(:,k));
        vecs(:,k) = vecs(:,k) - mu(k);
        sigma(k) = std(vecs(:,k));
        vecs(:,k) = vecs(:,k)/sigma(k);
    end
    vecs_n = vecs;
end
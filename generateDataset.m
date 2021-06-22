function [x y] = generateDataset(beta_1, beta_2, n, sigma_w)
    %Generate dataset
    x = 10*rand(n, 1);
    y  = beta_1*x+beta_2+sigma_w*randn(n,1);
end
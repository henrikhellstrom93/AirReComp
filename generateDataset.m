% X is a D x (d+1) matrix consisting of D datapoints of dimension d plus a column of ones.
% y is a D x 1 matrix consisting of D datapoints of dimension 1.
function [X, y] = generateDataset(beta, D, sigma_w)
    d = length(beta)-1;
    X = 10*rand(D, d);
    X = [ ones(D, 1) X ]; % Size is now D x (d+1)
    y  = X*beta + sigma_w*randn(D,1);
end
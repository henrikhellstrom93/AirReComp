clc
clear
% Parameters
num_rounds = 10000;
num_devices = 2;
num_samples = 1000;
sigma_w = 1; % Dataset noise
sigma_z = 1; % Channel noise
[x, y] = generateDataset(2, 1, num_samples, sigma_w);
[X, Y] = splitDataset(x, y, num_devices);

% Strong convexity and Lipschitz smoothness
m = num_samples+sum(x.^2)-sqrt(num_samples^2-2*num_samples*sum(x.^2)+4*sum(x).^2+sum(x.^2)^2);
L = num_samples+sum(x.^2)+sqrt(num_samples^2-2*num_samples*sum(x.^2)+4*sum(x).^2+sum(x.^2)^2);
step_length = 1/L;

% Initialize
beta_g = 10*rand(2,1);
loss = zeros(num_rounds, 1);

for r = 1:num_rounds
    % Least-squares loss on global dataset
    loss(r) = sum((y-x*beta_g(1)-beta_g(2)).^2);
    grad_l = getGradient(X, Y, beta_g);
    
    % Generate channel and noise
    h = generateH(num_devices);
    %Normalize
    [grad_l_n, mu, sigma] = normalize(grad_l);
    
    %Over-the-air computation
    grad_g = zeros(2, 1);
    %Generate power control
    [p, eta] = power_control(h, sigma_z, 1);
    for j = 1:2
        %Generate AWGN
        z = normrnd(0,sigma_z/2,1,1) + 1j*normrnd(0,sigma_z/2,1,1);
        %Apply fading and noise over MAC
        grad_g(j) = real(grad_l_n(j,:)*(abs(h).*sqrt(p)) + z)/(num_devices*sqrt(eta));
    end
    
    %Denormalize
    grad_l_dn = denormalize(grad_g, mu, sigma);
    
    % Model update
    beta_g = beta_g - step_length*grad_l_dn;
end

start = 10;
figure;
plot(start:num_rounds, loss(start:end))

y_hat = x*beta_g(1)+beta_g(2);
figure;
plot(x, y, ".")
hold on;
plot(x, y_hat, 'LineWidth', 2);

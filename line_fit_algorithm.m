clc
clear
% Parameters
budget = 100;
K = 2; % Number of devices
D = 1000; % Number of samples
sigma_w = 3; % Dataset noise
true_beta = [2;1]; % beta(1) = constant term
d = length(true_beta)-1; % input data dimension
[X, y] = generateDataset(true_beta, D, sigma_w);
[X_struct, y_struct] = splitDataset(X, y, K);

% Strong convexity and Lipschitz smoothness
% m = num_samples+sum(x.^2)-sqrt(num_samples^2-2*num_samples*sum(x.^2)+4*sum(x)^2+sum(x.^2)^2);
% L = num_samples+sum(x.^2)+sqrt(num_samples^2-2*num_samples*sum(x.^2)+4*sum(x)^2+sum(x.^2)^2);
% 
% % Normalized loss
% m = m/num_samples;
% L = L/num_samples;

% Initialize
% beta_g_init = 10*rand(2,1);
beta_g_init = zeros(d+1,1);

% Generate channel
h = generateH(K);

%% Everything above this line can be generated once and kept for multiple experiments
clc
sigma_z = 1; % Channel noise
num_tests = 1;
mean_err_tests = 0;
mean_loss = zeros(budget, 1);
step_length = 0.04;
schedule_type = "constant";
c_2 = 0;
T = 0;
if schedule_type == "constant"
    c_1 = 1;
elseif schedule_type == "lin_increase"
    c_1 = 0.005;
elseif schedule_type == "lin_decrease"
    c_1 = 10;
    c_2 = 0.5;
elseif schedule_type == "cutoff"
    c_1 = 4;
    T = 20;
    c_2 = 1;
end
retransmission_schedule = setTransmissionSchedule(schedule_type, budget, c_1, c_2, T);
calculate_bound = false;
if calculate_bound == true
    sigma_gradient = calculateSigmaG(length(true_beta), budget, X, Y, beta_g_init, step_length);
    [c, opt_schedule] = optimalSchedule(budget, h, sigma_z, step_length, m, L, sigma_gradient, beta_g_init, true_beta);
end

for t = 1:num_tests
    t
    %Initialization
    remaining_budget = budget;
    loss = zeros(budget, 1);
    beta_g = beta_g_init;
    mean_err = 0;
    num_tx = 1;

    r = 1;
    while remaining_budget > 0
        %Select number of uplink transmissions
        num_tx = retransmission_schedule(r);
        
        %Least-squares loss on global dataset
        loss(r) = (X*beta_g-y)'*(X*beta_g-y)/D;
        grad_l = getGradient(X_struct, y_struct, beta_g);

        %Transmit gradients over MAC
        rcv_grad = otaComputation(grad_l, sigma_z, num_tx, h);
        
        % Model update
        beta_g = beta_g - step_length*rcv_grad;
        
        %End condition
        remaining_budget = remaining_budget - num_tx;
        r = r + 1;
    end

    %Fill remaining losses with straight line
    loss(r:end) = loss(r-1)*ones(budget-r+1,1);
    mean_loss = mean_loss + loss;
end
mean_loss = mean_loss/num_tests;
loss = mean_loss;
loss(end)

%Plotting
plotting = true;
optimal_beta = (X'*X)\X'*y;
optimal_loss = (X*true_beta-y)'*(X*true_beta-y)/D;
start = 1;
finish = budget;
if plotting == true
    figure;
    plot(start:finish, loss(start:finish))
    hold on;
    plot(start:finish, optimal_loss*ones(finish-start+1, 1))
    legend("Fitted", "True line")
end
if schedule_type == "constant"
    filename = append("data/", "num_tx=", int2str(num_tx), "sigma_z=", int2str(sigma_z));
else
    filename = append("data/", append("schedule=", schedule_type));
end
save(filename, 'loss')

%%
y_hat = x*beta_g(1)+beta_g(2);
y_true = x*true_beta(1)+true_beta(2);
figure;
plot(x, y, ".")
hold on;
plot(x, y_hat, 'LineWidth', 2);
% hold on;
% plot(x, y_true, 'LineWidth', 2);
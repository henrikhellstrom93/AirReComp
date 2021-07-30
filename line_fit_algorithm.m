clc
clear
% Parameters
budget = 10000;
K = 10; % Number of devices
D = 1000; % Number of samples
sigma_w = 3; % Dataset noise
true_beta = [2;1]; % beta(1) = constant term
d = length(true_beta)-1; % input data dimension
[X, y] = generateDataset(true_beta, D, sigma_w);
[X_struct, y_struct] = splitDataset(X, y, K);

% Strong convexity and Lipschitz smoothness
if d == 1
    m = D+sum(X(:,2).^2)-sqrt(D^2-2*D*sum(X(:,2).^2)+4*sum(X(:,2))^2+sum(X(:,2).^2)^2);
    L = D+sum(X(:,2).^2)+sqrt(D^2-2*D*sum(X(:,2).^2)+4*sum(X(:,2))^2+sum(X(:,2).^2)^2);

    % Normalized loss
    m = m/D;
    L = L/D;
end

% Initialize
% beta_g_init = 10*rand(2,1);
beta_g_init = zeros(d+1,1);

% Generate channel
h = generateH(K);

%% Everything above this line can be generated once and kept for multiple experiments
clc
sigma_z = 1; % Channel noise
M = 1;
[table1, table2, table3] = power_control_lookup(h, 0, M, K, m, L);
num_tests = 1;
mean_err_tests = 0;
mean_loss = zeros(budget, 1);
grad_history = zeros(ceil(budget/M), 1);
% step_length = 0.02;
step_length = 1/L;
step_reduce = 0.9;
schedule_type = "constant";
c_2 = 0;
T = 0;
if schedule_type == "constant"
    c_1 = M;
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
        grad_history(r) = rcv_grad(1);
        
        %End condition
        remaining_budget = remaining_budget - num_tx;
        r = r + 1;
        if mod(r, 300) == 0
            step_length = step_length*step_reduce;
        end
    end

    %Fill remaining losses with straight line
    loss(r:end) = loss(r-1)*ones(budget-r+1,1);
    mean_loss = mean_loss + loss;
end
mean_loss = mean_loss/num_tests;
loss = mean_loss;
loss(end);

%Plotting
plotting = true;
optimal_beta = (X'*X)\X'*y;
optimal_loss = (X*optimal_beta-y)'*(X*optimal_beta-y)/D;
start = 10;
finish = budget;
if plotting == true
    figure;
    plot(start:finish, loss(start:finish))
    hold on;
    plot(start:finish, optimal_loss*ones(finish-start+1, 1))
    legend("OTA FL fit", "Optimal fit")
end
if schedule_type == "constant"
    filename = append("data/", "num_tx=", int2str(num_tx), "sigma_z=", int2str(sigma_z), "budget=", int2str(budget));
else
    filename = append("data/", append("schedule=", schedule_type));
end
save(filename, 'loss', 'retransmission_schedule')

%Printing
fit_loss = mean(loss(budget/M-budget/M*0.05:budget/M));
disp([ 'Loss gap: ', num2str(fit_loss-optimal_loss) ])
disp([ 'rcv_grad size: ', num2str( abs( rcv_grad(1) ) ) ])


%%
y_hat = x*beta_g(1)+beta_g(2);
y_true = x*true_beta(1)+true_beta(2);
figure;
plot(x, y, ".")
hold on;
plot(x, y_hat, 'LineWidth', 2);
% hold on;
% plot(x, y_true, 'LineWidth', 2);
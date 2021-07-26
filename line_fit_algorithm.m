clc
clear
% Parameters
budget = 1000;
num_devices = 4;
num_samples = 1000;
num_tests = 10;
sigma_w = 1; % Dataset noise
sigma_z = 5; % Channel noise
true_beta = [2;1];
[x, y] = generateDataset(true_beta(1), true_beta(2), num_samples, sigma_w);
[X, Y] = splitDataset(x, y, num_devices);

% Strong convexity and Lipschitz smoothness
m = num_samples+sum(x.^2)-sqrt(num_samples^2-2*num_samples*sum(x.^2)+4*sum(x)^2+sum(x.^2)^2);
L = num_samples+sum(x.^2)+sqrt(num_samples^2-2*num_samples*sum(x.^2)+4*sum(x)^2+sum(x.^2)^2);

% Initialize
% beta_g_init = 10*rand(2,1);
beta_g_init = [0; 0];

% Generate channel
h = generateH(num_devices);

%% Everything above this line can be generated once and kept for multiple experiments
clc
mean_err_tests = 0;
mean_loss = zeros(budget, 1);
schedule_type = "constant";
step_length = 1/(4*L);
max_tx = 1;
calculate_bound = true;
if calculate_bound == true
    sigma_gradient = calculateSigmaG(length(true_beta), budget, X, Y, beta_g_init, step_length);
end

[c, opt_schedule] = optimalSchedule(budget, h, sigma_z, step_length, m, L, sigma_gradient, beta_g_init, true_beta);

for t = 1:num_tests
    t
    %Initialization
    remaining_budget = budget;
    loss = zeros(budget, 1);
    beta_g = beta_g_init;
    mean_err = 0;
    retransmission_schedule = setTransmissionSchedule(schedule_type, max_tx, budget);
    num_tx = 1;

    r = 1;
    while remaining_budget > 0
        %Select number of uplink transmissions
        if retransmission_schedule(num_tx) > 0
            retransmission_schedule(num_tx) = retransmission_schedule(num_tx) - 1;
        else
            while retransmission_schedule(num_tx) <= 0
                num_tx = num_tx + 1;
            end
        end
        
        %Least-squares loss on global dataset
        loss(r) = sum((y-x*beta_g(1)-beta_g(2)).^2);
        grad_l = getGradient(X, Y, beta_g);

        %Transmit gradients over MAC
        rcv_grad = otaComputation(grad_l, sigma_z, num_tx, h);
        mean_err = mean_err + abs(mean(grad_l, 2)-rcv_grad);
        % Model update
        beta_g = beta_g - step_length*rcv_grad;
        
        %End condition
        remaining_budget = remaining_budget - num_tx;
        r = r + 1;
        
        %Update step size
        %step_length = step_length*5^(-1/budget);
    end
    mean_err_tests = mean_err_tests + mean_err/r;

    %Fill remaining losses with straight line
    loss(r:end) = loss(r-1)*ones(budget-r+1,1);
    mean_loss = mean_loss + loss;
end
mean_err_tests = mean_err_tests/num_tests
mean_loss = mean_loss/num_tests;
loss = mean_loss;

%Plotting
true_loss = sum((y-x*true_beta(1)-true_beta(2)).^2);
start = 10;
finish = budget;
figure;
plot(start:finish, loss(start:finish))
hold on;
plot(start:finish, true_loss*ones(finish-start+1, 1))
legend("Fitted", "True line")
if schedule_type == "constant"
    filename = append("data/", "num_tx=", int2str(num_tx), "sigma_z=", int2str(sigma_z))
elseif schedule_type == "simple"
    filename = append("data/", "schedule=simple")
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
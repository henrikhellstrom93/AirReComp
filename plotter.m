clear
clc

num_tx_list = [1 2 4 8];
num_tests = length(num_tx_list);
budget = 1000;
sigma_z = 10;
path = "data/";

losses = zeros(budget, num_tests);
for i = 1:num_tests
    num_tx = num_tx_list(i);
    filename = append(path, "num_tx=", int2str(num_tx), "sigma_z=", int2str(sigma_z));
    a = load(filename);
    losses(:,i) = a.loss;
end

start = 1;
finish = budget;

for i = 1:num_tests
    plot(start:finish, losses(start:finish, i))
    hold on;
end
xlabel("Communication round")
ylabel("Least-squares loss")
legend("M=1", "M=2", "M=4", "M=8")
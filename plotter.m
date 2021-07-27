clear
clc

num_tx_list = [1 2 4 8 20 50];
schedule = "constant";
num_tests = length(num_tx_list);
if schedule == "simple"
    num_tests = num_tests + 1;
end
budget = 500;
sigma_z = 100;
path = "data/";

losses = zeros(budget, num_tests);
for i = 1:num_tests
    filename = "";
    if i == num_tests && schedule == "simple"
        filename = append(path, "schedule=simple")
    else
        num_tx = num_tx_list(i);
        filename = append(path, "num_tx=", int2str(num_tx), "sigma_z=", int2str(sigma_z));
    end
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
legend("M=1", "M=2", "M=4", "M=8", "schedule")

losses(500,:)
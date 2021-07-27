clear
clc

num_tx_list = [1 2 4 8];
schedules = ["constant", "constant", "constant", "constant", "lin_increase", "lin_decrease", "cutoff"];
num_tests = length(schedules);
budget = 400;
sigma_z = 20;
path = "data/";

losses = zeros(budget, num_tests);
for i = 1:num_tests
    filename = "";
    if schedules(i) ~= "constant"
        filename = append(path, append("schedule=", schedules(i)))
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
legend("M=1", "M=2", "M=4", "M=8", "lin\_inc", "lin\_dec", "cutoff")

losses(budget,:)
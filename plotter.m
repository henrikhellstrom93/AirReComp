clear
clc

num_tx_list = [1 2 4 8];
% num_tx_list = [1 2];
%schedules = ["constant", "constant", "constant", "constant", "lin_increase", "lin_decrease", "cutoff"];
schedules = ["constant", "constant", "constant", "constant", "cutoff"]
% schedules = ["constant", "constant"]
num_curves = length(schedules);
budget = 1000;
sigma_z = floor(sqrt(20));
alpha = 1;
padding = false; %Pads the data to show them in the same timescale
path = "data/";

losses = zeros(budget, num_curves);
for i = 1:num_curves
    filename = "";
    if schedules(i) ~= "constant"
        filename = append(path, append("schedule=", schedules(i)))
    else
        num_tx = num_tx_list(i);
        filename = append(path, "num_tx=", int2str(num_tx), "sigma_z=", int2str(sigma_z), "budget=", int2str(budget), "alpha=", num2str(alpha));
    end
    a = load(filename);
    losses(:,i) = a.loss;
    schedule = a.retransmission_schedule;
    padded_loss = zeros(length(a.loss), 1);
    if padding == true
        schedule_index = 1;
        padding_index = 1;
        while padding_index <= budget
            if schedule(schedule_index) == 0
                schedule_index = schedule_index + 1;
                if schedule_index > budget
                    break
                end
            elseif schedule(schedule_index) == 1
                padded_loss(padding_index) = a.loss(schedule_index);
                schedule_index = schedule_index + 1;
                padding_index = padding_index + 1;
            else
                padded_loss(padding_index) = a.loss(schedule_index);
                padding_index = padding_index + 1;
                num_repetitions = schedule(schedule_index)-1;
                for j = 1:num_repetitions
                    padded_loss(padding_index) = a.loss(schedule_index);
                    padding_index = padding_index+1;
                end
                schedule_index = schedule_index + 1;
            end
        end
        losses(:,i) = padded_loss;
    end
end

start = 1;
finish = budget/(alpha+1);

for i = 1:num_curves
    plot(start:finish, losses(start:finish, i))
    hold on;
end
xlabel("Communication round")
ylabel("Least-squares loss")
legend("M=1", "M=2", "M=4", "M=8", "cutoff", "M=32")
ylim([9, 100])

losses(budget,:)
%The schedule is a vector of length budget, where entry i is the number of
%uplink transmissions at communication round i.
function [schedule] = setTransmissionSchedule(type, budget, c_1, c_2, c_3, T_1, T_2)
    schedule = zeros(budget,1);
    rem_budget = budget;
    i = 1;
    while rem_budget > 0
        if type == "constant"
            if rem_budget < ceil(c_1)
                break
            end
            schedule(i) = ceil(c_1);
        elseif type == "lin_increase"
            schedule(i) = ceil(c_1*i);
        elseif type == "lin_decrease"
            if ceil(c_1 - c_2*i) > 0
                schedule(i) = ceil(c_1 - c_2*i);
            else
                schedule(i) = 1;
            end
        elseif type == "cutoff"
            if i <= T_1
                schedule(i) = c_1;
            elseif i <= T_2
                schedule(i) = c_2;
            else
                schedule(i) = c_3;
            end
        end
        rem_budget = rem_budget - schedule(i);
        i = i + 1;
    end
end
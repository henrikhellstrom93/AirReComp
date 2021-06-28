function [schedule] = setTransmissionSchedule(type, max_tx, budget)
    schedule = zeros(max_tx,1);
    if type == "simple"
        for i = 1:max_tx
            schedule(i) = floor(budget/max_tx/i);
        end
    end
    if type == "constant"
        schedule(max_tx) = budget;
    end
end
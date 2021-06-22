function [h] = generateH(num_devices)
    h = normrnd(0,1/2,num_devices,1) + 1j*normrnd(0,1/2,num_devices,1);
    [quality_indicator, sort_index] = sort(abs(h).^2);
    h_sort = zeros(num_devices,1);
    for i = 1:num_devices
        h_sort(i) = h(sort_index(i));
    end
    h = h_sort;
end
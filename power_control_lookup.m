% Returns a lookup table for the power control, given the noise std
function [table, table2, table3] = power_control_lookup(h, range, M, K, m, L)
    table = zeros(length(range), 1);
    table2 = zeros(length(range), 1);
    for i = 1:length(range)
        [p, eta] = power_control(h, range(i), M);
        table(i) = K*(m+L)*sqrt(eta)/(2*m*L*abs(h)'*sqrt(p));
        table2(i) = 2*sqrt(eta)*abs(h)'*sqrt(p)/((m+L)*abs(h).^2'*p);
        table3(i) = sqrt(eta);
    end
end


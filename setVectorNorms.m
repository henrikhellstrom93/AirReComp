function [ret_vecs, norms] = setVectorNorms(vecs, new_norm)
    dims = size(vecs);
    num_devices = dims(2);
    ret_vecs = zeros(2, num_devices);
    norms = zeros(num_devices, 1);
    for k = 1:num_devices
        norms(k) = norm(vecs(:,k));
        ret_vecs(:,k) = vecs(:,k)/norms(k)*new_norm;
    end
end
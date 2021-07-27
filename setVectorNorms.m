% ret_vecs is a (d+1) x K matrix consisting of K vectors of length (d+1)
% norms is a K x 1 vector consisting of the norms of the vecs input vectors
function [ret_vecs, norms] = setVectorNorms(vecs, new_norm)
    dims = size(vecs);
    d = dims(1)-1;
    K = dims(2);
    ret_vecs = zeros(d+1, K);
    norms = zeros(K, 1);
    for k = 1:K
        norms(k) = norm(vecs(:,k));
        ret_vecs(:,k) = vecs(:,k)/norms(k)*new_norm;
    end
end
function [vec_dn] = denormalize(vec_n, mu, sigma)
    vec_dn = vec_n*mean(sigma);
    vec_dn = vec_dn + mean(mu);
end
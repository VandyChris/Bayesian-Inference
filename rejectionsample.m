function smpl = rejectionsample( f, fmax, bound, n )
% f: the destination function for sample, f(N*D matrix) = N*1 matrix
% f: the maximum value of f
% bound: constains the bound of each variable, which should be a matrix of
% [x1_lower_bound, x1_upper_bound; x2_lower_bound; ...]
% n: number of sampling points you want to include in smpl
%
% smpl: n * D matrix, where D is the dimension
% here the instrumental distribution is simply the uniform distribution
% within the given bounds
% reference: https://en.wikipedia.org/wiki/Rejection_sampling
%%

n_candidate = 500;
D = size(bound, 1);

    function y = sample_candidate()
        y = zeros(n_candidate, D);
        for i = 1:D
            y(:, i) = unifrnd(bound(i, 1), bound(i, 2), [n_candidate, 1]);
        end
    end
smpl = zeros(n+n_candidate, D);
n_smpl = 0; % current non-zero samples in smpl
while n_smpl < n
    smpl_pool = sample_candidate();
    p = f(smpl_pool)./fmax;
    U = unifrnd(0, 1, [n_candidate, 1]);
    n_accept = sum(U<p);
    ratio = n_accept / n_candidate;
    disp(['Accpt ratio: ', num2str(ratio)]);
    smpl(n_smpl+1:n_smpl+n_accept, :) = smpl_pool(U<p, :);
    n_smpl = n_smpl + n_accept;
end

smpl = smpl(1:n, :);


end
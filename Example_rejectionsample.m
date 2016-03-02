clc;
clear;
close all;
%% 1-D normal
mu = 5;
sigma = 0.5;
f = @(x) 1/sigma/sqrt(2*pi)*exp(-(x-mu).^2/2/sigma^2);
bound = [mu-4*sigma, mu+4*sigma];
fmax = f(mu);
n = 3000;
tic
smpl = rejectionsample( f, fmax, bound, n );
toc

figure(1)
[fi, xi] = ksdensity(smpl);
x = linspace(bound(1), bound(2), 500);
y = normpdf(x, mu, sigma);
plot(x, y, 'red', 'linewidth', 2, 'linestyle', '--')
hold on
plot(xi, fi, 'linewidth', 2);
legend({'True distribution', 'Rejection sampling'})
xlabel('QoI')
ylabel('PDF')
title('1-D problem')

%% multi normal
mu = [2, 3, 4, 5];
corr = [1, 0.174, 0.451, 0.082;
        0.174, 1, -0.80, 0.059;
        0.451, -0.80, 1, -0.004;
        0.082, 0.059, -0.004, 1.0];
std = [2,1,1,1];
sigma = zeros(4, 4);
for i = 1:4;
    for j = 1:4;
        sigma(i, j) = std(i) * std(j) * corr(i, j);
    end
end

k = size(mu, 2);

f =@(x) mvnpdf(x, mu, sigma);
bound = [mu(1) - 4*std(1), mu(1) + 4*std(1);
         mu(2) - 4*std(2), mu(2) + 4*std(2);
         mu(3) - 4*std(3), mu(3) + 4*std(3);
         mu(4) - 4*std(4), mu(4) + 4*std(4)];
fmax = f(mu);
n = 3000;
tic
smpl = rejectionsample( f, fmax, bound, n );
toc

for i = 1:k
    figure(i+1)
    [fi, xi] = ksdensity(smpl(:, i));
    x = linspace(bound(i, 1), bound(i, 2), 500);
    y = normpdf(x, mu(i), std(i));
    plot(x, y, 'red', 'linestyle', '--', 'linewidth', 2);
    hold on
    plot(xi, fi, 'linewidth', 2);
    
    legend({'True distribution', 'Rejection sampling'})
    xlabel(strcat('QoI', num2str(i)))
    ylabel('PDF')
    title(strcat('4-D problem, QoI', num2str(i)))
end

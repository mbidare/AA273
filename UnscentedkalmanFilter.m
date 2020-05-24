function [mu,sigma]= UnscentedkalmanFilter(f, g, y, u, n,lambda, Q, R)

num = length(y(1,:));
dim = n;

mu = zeros(dim,num);
sigma = zeros(dim,dim,num);

weights = generate_weights(n, lambda);
for i = 2:num
    
    prediction = predict(f,mu(:,i-1), u(:,i-1), sigma(:,:,i-1), n, lambda);
    [mu(:,i),sigma(:,:,i)] = update(g, prediction, y(i), n, lambda, weights, Q, R);
    
end
end

function sigpoints = computeSigmaPoints(mu,sigma, n, lambda)

half_sigma = sqrtm((lambda+n)*sigma);

sigpoints = zeros(3,2*n+1);
sigpoints(:,1) = mu;
for j=2:n+1
    sigpoints(:,j) = mu + half_sigma(:,j-1);
    sigpoints(:,n+j) = mu - half_sigma(:,j-1);
end
end

function weights = generate_weights(n, lambda)
% generate weights for finding sigma points
weights = ones(1,2*n+1);
weights(1) = lambda/(n+lambda);
weights(2:end) = 1/(2*(lambda+n));
end

function prediction = predict(f,mu, u, sigma, n, lambda)

sigpoints = computeSigmaPoints(mu,sigma, n, lambda);

prediction = zeros(3,2*n+1);
for j=1:2*n+1
    prediction(:,j) = f(sigpoints(:,j), u(:,j));
end
end

function [mu,sigma] = update(g, prediction, y, n, lambda, weights, Q, R)

% mean %
mu_bar = prediction*weights';

% cov %
cov_bar = Q;
for j=1:2*n+1
    temp = (prediction(:,j)-mu_bar);
    cov_bar = cov_bar + weights(j)*(temp*temp');
end

updated_sigpoints = computeSigmaPoints(mu_bar,cov_bar, n, lambda);

sig_measure = zeros(1,2*n+1);
for j=1:2*n+1
    sig_measure(j) = g(updated_sigpoints(:,j));
end

y_mean = sig_measure*weights';

% cov Y %
cov_y = R;
for j=1:2*n+1
    temp = (sig_measure(j)-y_mean);
    cov_y = cov_y + weights(j)*(temp*temp');
end

% cov XY %
cov_xy = zeros(3,1);
for j=1:2*n+1
    temp = (updated_sigpoints(:,j) - mu_bar)*(sig_measure(j)-y_mean)';
    cov_xy = cov_xy + weights(j)*temp;
end

% update %
mu = mu_bar + cov_xy*inv(cov_y)*(y - y_mean);
sigma = cov_bar - cov_xy*inv(cov_y)*(cov_xy');
end

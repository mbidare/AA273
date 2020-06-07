function [mu,sigma]= UnscentedkalmanFilter(y, u, C, Q, R,dt)
 
lambda = 2;
 
num = length(y(1,:));
dim = 6;
 
mu = zeros(dim,num);
sigma = zeros(dim,dim,num);
sigma(:,:,1) = .1*eye(dim);
 
for i = 2:num
    A = zeros(6,6); A(1,1) = 1; A(3,3) = 1; A(5,5) = 1; A(5,6) = dt; A(6,6)=1;
    B = zeros(6,3); B(1,1) = dt*cos(u(2)); B(1,2) = -dt*u(1)*sin(u(2));
    B(3,1) = dt*sin(u(2)); B(3,2) = dt*u(1)*cos(u(2)); B(6,3) = dt;
    [pred_mean, pred_cov] =  predict(A, B, mu(:,i-1), u(:,i-1), sigma(:,:,i-1), dim, lambda, Q, dt);
    [mu(:,i),sigma(:,:,i)] = update (C, pred_mean, pred_cov, y(:,i), dim, lambda, Q, R);
    
end
end
 
function sigpoints = computeSigmaPoints(mu,sigma, n, lambda)
% This function computes the sigma points used by the unscented kalman
% filter
half_sigma = sqrtm((lambda+n)*sigma);
 
sigpoints = zeros(n,2*n+1);
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
 
function [pred_mean, pred_cov] = predict(A, B, mu, u, sigma, n, lambda, Q, dt)
sigpoints = computeSigmaPoints(mu,sigma, n, lambda);
prediction = zeros(n,2*n+1);
for j=1:2*n+1
    prediction(:,j) = act_dyn(sigpoints(:,j), u, j, dt);
end
 
weights = generate_weights(n, lambda);
 
% mean %
pred_mean = prediction*weights';
 
% cov %
pred_cov = Q;
for j=1:2*n+1
    temp = (prediction(:,j)-pred_mean);
    pred_cov = pred_cov + weights(j)*(temp*temp');
end
 
 
end
 
function [mu,sigma] = update(C, pred_mean, pred_cov, y, n, lambda, Q, R)
% This function updates the unscented kalman prediction of the robot's
% position
 
weights = generate_weights(n, lambda);
 
updated_sigpoints = computeSigmaPoints(pred_mean,pred_cov, n, lambda);
 
sig_measure = zeros(n,2*n+1);
for j=1:2*n+1
    sig_measure(:,j) = meas_model(C,updated_sigpoints(:,j));
end
 
y_mean = sig_measure*weights';
 
% cov Y %
cov_y = R;
for j=1:2*n+1
    temp = (sig_measure(:,j)-y_mean);
    cov_y = cov_y + weights(j)*(temp*temp');
end
 
% cov XY %
cov_xy = zeros(length(pred_mean),length(y_mean));
for j=1:2*n+1
    temp = (updated_sigpoints(:,j) - pred_mean)*(sig_measure(:,j)-y_mean)';
    cov_xy = cov_xy + weights(j)*temp;
end
 
% update %
mu = pred_mean + cov_xy*inv(cov_y)*(y - y_mean);
sigma = pred_cov - cov_xy*inv(cov_y)*(cov_xy');
end


function [mu,sigma] = ExtendedkalmanFilter(y,u,A,B,C,g,Q)
%Extended kalman filter
%y, u are fat vectors
%A, C are aproprate jacobians
%act_dyn is non-linear dynamics
%g is non-linear estimator
num = length(y);
dim = length(A,1);
mu = zeros(dim,num);
sigma = zeros(dim,dim,num);
mu(:,1) = y(:,1);

for j = 2:num
    %predict
    muPred = act_dyn(mu(:,j-1),u(j-1),dt, A, B);
    sigPred = A*sigma(:,:,j-1)*A' + Q;
    %update
    mu(:,j) = muPred+sigPred*C'*(C*sigPred*C'+R)^(-1)*(y(:,j)-g(muPred,u(:,j)));
    sigma(:,:,j) = sigPred - sigPred*C'*(C*sigPred*C'+R)^(-1)*C*sigPred;
end

end
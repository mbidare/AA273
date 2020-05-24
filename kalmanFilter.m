function [mu,sigma] = kalmanFilter(y,u,A,B,C,Q,R,Cmean)
%Plain old kalman filter
%y, u are fat vectors
num = length(y);
dim = length(A,1);
mu = zeros(dim,num);
sigma = zeros(dim,dim,num);
mu(:,1) = y(:,1);

for j = 2:num
    %predict
    muPred = lin_dyn(mu(:,j-1),u(j-1),dt, A, B);
    sigPred = A*sigma(:,:,j-1)*A' + Q;
    %update
    mu(:,j) = muPred+sigPred*C'*(C*sigPred*C'+R)^(-1)*(y(:,j)-(C*muPred-Cmean));
    sigma(:,:,j) = sigPred - sigPred*C'*(C*sigPred*C'+R)^(-1)*C*sigPred;
end

end


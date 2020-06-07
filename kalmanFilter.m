function [mu,sigma] = kalmanFilter(y,u,C,Q,R,Cmean,dt)
%Plain old kalman filter
%y, u are fat vectors
num = length(y);
dim = 6;
mu = zeros(dim,num);
sigma = zeros(dim,dim,num);
mu(:,1) = y(:,1);

for j = 2:num
    %predict
    A = zeros(6,6); A(1,1) = 1; A(3,3) = 1; A(5,5) = 1; A(5,6) = dt; A(6,6)=1;
    B = zeros(6,3); B(1,1) = dt*cos(u(2)); B(1,2) = -dt*u(1)*sin(u(2));
    B(3,1) = dt*sin(u(2)); B(3,2) = dt*u(1)*cos(u(2)); B(6,3) = dt/mass;
    muPred = lin_dyn(mu(:,j-1),u(:,j-1),dt, A, B);
    sigPred = A*sigma(:,:,j-1)*A' + Q;
    %update
    mu(:,j) = muPred+sigPred*C'*(C*sigPred*C'+R)^(-1)*(y(:,j)-(C*muPred+Cmean));
    sigma(:,:,j) = sigPred - sigPred*C'*(C*sigPred*C'+R)^(-1)*C*sigPred;
end

end


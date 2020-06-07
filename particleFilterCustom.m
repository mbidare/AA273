function [muPF] = particleFilterCustom(y,u,numParticles,dim,Q,R,dt,C)
num = length(y);
particles = zeros(dim,numParticles,num);
muPF = zeros(dim,num);
muPF(:,1) = mean(particles(:,:,1),2);

wb = waitbar(0,'Particle Filter (This can be slow depending on settings)');
for i = 2:num
    xPred = zeros(dim,numParticles);
    wp = mvnrnd(zeros(dim,1),Q,numParticles)';
    for j = 1:numParticles
       xPred(:,j) = act_dyn(particles(:,j,i-1),u(:,i-1),i-1,dt) + wp(:,j);
    end
   wBar = zeros(1,numParticles);
   dimR = size(R,1);
   eta = 1/sqrt((2*pi)^dimR*det(R));
   vp = mvnrnd([0;0;0;0;0;0],R,numParticles)';
   for j = 1:numParticles
       gVal = meas_model(C,particles(:,j,i-1)) ; %+ vp(:,j); %g(particles(:,j,i-1)) + vp(:,j);
       wBar(j) = eta*exp(-0.5*(y(:,i)-gVal)'*R^(-1)*(y(:,i)-gVal));
   end
   w = wBar./(sum(wBar));
   wSum = w;
   for j = 2:numParticles
       wSum(j) = sum(w(1:j));
   end
   for j = 1:numParticles
       pointer = rand;
       [~,indx] = find(wSum>pointer,1);
       particles(:,j,i) = xPred(:,j);
   end
   muPF(:,i) = mean(particles(:,:,i),2);
   waitbar(i/num,wb,'Particle Filter(This can be slow depending on settings)');
end
close(wb)
end
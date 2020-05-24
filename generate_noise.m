function generate_noise(Qgauss,Qbi1,Qbi2,mu,numSteps,s)
    gaussNoise = gaussian(zeros(6,1),Qgauss,numSteps);
    uniformNoise = zeros(6,numSteps);
    uniformNoise(2,:) = unifrom(-0.75,0.75,[1,numSteps]);
    uniformNoise(4,:) = unifrom(-0.75,0.75,[1,numSteps]);
    uniformNoise(6,:) = unifrom(-0.1,0.1,[1,numSteps]);
    expNoise = zeros(6,numSteps);
    expNoise(2,:) = exponteialNoise(0.002,[1,numSteps]);
    expNoise(4,:) = exponteialNoise(0.002,[1,numSteps]);
    expNoise(6,:) = exponteialNoise(0.0004,[1,numSteps]);
    brnNoise = zeros(6,numSteps);
    brnNoise(2,:) = 0.0005*coloredNoise('brown',numSteps);
    brnNoise(4,:) = 0.0005*coloredNoise('brown',numSteps);
    brnNoise(6,:) = 0.0001*coloredNoise('brown',numSteps);
    cauchNoise = zeros(6,numSteps);
    cauchNoise(2,:) = 0.25*cauchy(10,numSteps);
    cauchNoise(4,:) = 0.25*cauchy(10,numSteps);
    cauchNoise(6,:) = 0.05*cauchy(10,numSteps);
    biModalNoise = bimodal(mu,Qbi1,-1*mu,Qbi2,numSteps);
    if s
        save('noise')
    end
end

function [noise] = gaussian(mean,cov,num)
noise = mvnrnd(mean,cov,num)';
end

function [noise] = unifrom(min,max,num)
noise = min + (max-min)*rand(num);
end

function [noise] = exponteialNoise(mean,num)
noise = exprnd(mean,num);
end

function [noise] = coloredNoise(col,num)
brn =  dsp.ColoredNoise('Color',col,'SamplesPerFrame',num);
noise = brn();
noise = noise(randperm(length(noise)));
end

function [noise] = cauchy(dof,num)
noise = trnd(dof,num,1);
end

function [noise] = bimodal(mu1,cov1,mu2,cov2,num)
a = mvnrnd(mu1,cov2,num)';
b = mvnrnd(mu2,cov1,num)';
noise = [a;b];
noise = noise(randperm(length(noise)));
end
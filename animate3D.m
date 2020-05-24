function [] = animate3D(state,pauseTime,sigma,ellipseIndx)
%posititions are state(1,:) state(3,:), state(5,:)
num = size(state,2);
figure
plot3(state(1,1:4),state(3,1:4),state(5,1:4),'k');
grid on
hold on
drawnow
pause(pauseTime);
xlabel('X(m)')
ylabel('Y(m)')
zlabel('Z(m)')
for i = 5:5:num-5
    plot3(state(1,i:i+5),state(3,i:i+5),state(5,i:i+5),'k');
    drawnow
    pause(pauseTime);
    if(any(ismember([i:i+5],ellipseIndx)))
        ell = sigma(:,:,i);
        C = ell(1:2:5,1:2:5);
        mu = [state(1,i);state(3,i);state(5,i)];
        P = 0.95;
        s = -2*log(1-P);
        [eigvec,eigval] = eig(s*C);
        [X,Y,Z] = ellipsoid(0,0,0,1,1,1);
        XYZ = [X(:),Y(:),Z(:)]*sqrt(eigval)*eigvec';
        X(:) = (XYZ(:,1)+mu(1));
        Y(:) = (XYZ(:,2)+mu(2));
        Z(:) = (XYZ(:,3)+mu(3));
        h=surf(X,Y,Z,'FaceAlpha',0.5);
        h.EdgeColor = 'none';
    end
end
    
end


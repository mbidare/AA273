function plotting(mu,x,disp_name)
    figure
    plot3(x(1,:),x(3,:),x(5,:))
    hold on
    plot3(mu(1,:),mu(3,:),mu(5,:))
    legend('Actual Trajectory','Prediced Trajectory')
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title(disp_name)
end
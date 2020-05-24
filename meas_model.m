% Measurement Model used to calculate y in main and in EKF (can be non-lin)
function y = meas_model(H,x)
    y = H*x;
end
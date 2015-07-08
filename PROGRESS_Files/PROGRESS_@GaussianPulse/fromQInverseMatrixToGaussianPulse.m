function [ gaussianPulseArray ] = fromQInverseMatrixToGaussianPulse( QInverseMatrix )
%fromQMatrixToGaussianPulse Compute the gaussian pulsed beam from its complex q
% matrices        
% The code is also vectorized. The gaussian pulsed beam array could be computed
% from multiple inputs.

% QMatrix:  Size 2X2XN for N gaussian pulsed beams

gaussianPulseArray = GaussianPulse(QInverseMatrix);
end


addpath('../../FastDMF')
addpath('../../FastDMF/data')
addpath('../../jidt/demos/octave')  
javaaddpath('../../jidt/infodynamics.jar');

% Calculate BOLD data from dmf model
%regions = 46:5:90;
%regions = 1:5:90;
regions = [44,45,46,47];
example

% Read Data
data = b(regions,:)';

TC = totalCorrelation(data)

DTC = dualTotalCorrelation(data)

O = TC - DTC
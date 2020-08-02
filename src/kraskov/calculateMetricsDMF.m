addpath('../../FastDMF')
addpath('../../FastDMF/data')
addpath('../../jidt/demos/octave')  
javaaddpath('../../jidt/infodynamics.jar');

% Calculate BOLD data from dmf model
%regions = 46:5:90;
%regions = 1:5:90;
%regions = [44,45,46,47];

dmn_ids = [12 79 16 75 18 73 33 58 34 57];
example

% Read Data
data = b(dmn_ids,:)';

[TC, TC_list] = totalCorrelation(data);

[DTC, DTC_list] = dualTotalCorrelation(data);

O = TC - DTC;

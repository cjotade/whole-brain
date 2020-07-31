addpath('../../jidt/demos/octave')  
javaaddpath('../../jidt/infodynamics.jar');

% Read Data
datafile = '../../data/muestras_1e5/gaussian_9/samples50.txt';
data = load(datafile);

[TC, TC_list] = totalCorrelation(data);

[DTC, DTC_list] = dualTotalCorrelation(data);

O = TC - DTC;
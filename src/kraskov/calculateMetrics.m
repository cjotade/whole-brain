addpath('../../jidt/demos/octave')  
javaaddpath('../../jidt/infodynamics.jar');

% Read Data
datafile = '../../data/muestras/gaussian_7/samples1.txt';
data = load(datafile);

TC = totalCorrelation(data)

DTC = dualTotalCorrelation(data)

O = TC - DTC
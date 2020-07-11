addpath('../../jidt/demos/octave')  
javaaddpath('../../jidt/infodynamics.jar');

% Read Data
datafile = '../../data/muestras/10000/samples9.txt';
data = load(datafile);

TC = totalCorrelation(data)

DTC = dualTotalCorrelation(data)

O = TC - DTC
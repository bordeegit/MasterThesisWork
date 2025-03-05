clear 
close all

set(groot,'DefaultAxesFontSize', 15);
set(groot,'DefaultLineLineWidth', 1.5);
set(groot,'DefaultTextInterpreter', 'Latex');
set(groot,'DefaultAxesTickLabelInterpreter', 'Latex');  
set(groot,'DefaultLegendInterpreter', 'Latex');

flightData = 'FlightData\Kitemill_90S.mat';
processingScriptPath = which("Kitemill_TL");
windEstimationLS(flightData, processingScriptPath)

flightData = 'FlightData\Standard_LinY.mat'; 
processingScriptPath = which('SoftKite_TL');
windEstimationLS(flightData, processingScriptPath)

addpath(genpath(pwd));
flightData = which('Log_12m2_fix.mat'); 
processingScriptPath = which('RealFlight_TL');
windEstimationLS(flightData, processingScriptPath)



% Kitemill_TL
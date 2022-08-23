clear all
close all


%% physical array parameters
configFBGarray.pulseLength = 220; 
configFBGarray.compensatorLength = 220; %%length in ns (standard choice is == pulseLength)
configFBGarray.channelDelayLeadIn = 10; %% 10 ns delay for lead in fiber.
configFBGarray.numFBG = 25;
configFBGarray.fbgReflectivity = 0.01; %%% reflectivity of 1 is perfect relection, 0 is no reflection.
configFBGarray.channelDelay = 220*ones(1,configFBGarray.numFBG-1); %% channel length in ns (25m corresponds to 220ns)
configFBGarray.samplingRateADC = 100e6; %%% current sampling rate is 100MHz
configFBGarray.samplesPerOpticalWavePeriod = 10;  %%% this is a prameter of the simulator. Determines how many samples per sampling period of the ADC.
configFBGarray.A = 0.9994*ones(1,configFBGarray.numFBG-1); %%% corresponds to 0.2 dB loss per km

%% dependant parameters (i.e. depend on previous parameters). Don't need to be set, they are automatically calculated
configFBGarray.samplingPeriodADC = 1/configFBGarray.samplingRateADC*1e9;
configFBGarray.Ts = configFBGarray.samplingPeriodADC/configFBGarray.samplesPerOpticalWavePeriod;
configFBGarray.tp = 0:configFBGarray.Ts:configFBGarray.pulseLength-configFBGarray.Ts;
configFBGarray.f = 1/configFBGarray.samplingPeriodADC;

%% modulation parameters
configFBGarray.fm = 1.25/(configFBGarray.pulseLength);  %% modulation frequency is such that there is 1.25 periods per 220ns
configFBGarray.m = 2*pi*configFBGarray.fm*configFBGarray.tp;

configFBGarray

%% array independant parameters
%theta = zeros(N,configFBGarray.numFBG-1);
%psi = zeros(N,2*(configFBGarray.numFBG-1));
numTimeFrames2sim = configFBGarray.numFBG;
    
N = 10000; %number of pulses launched


cd = 220:20:520;
bins = 0.025:0.025:3;
hist0 = zeros(length(cd),length(bins),configFBGarray.numFBG-1);
hist120 = zeros(length(cd),length(bins),configFBGarray.numFBG-1);
hist240 = zeros(length(cd),length(bins),configFBGarray.numFBG-1);
robustSTD = zeros(length(cd), configFBGarray.numFBG-1);
for ii = 1:length(cd)

    theta = 0.5*pi*rand(N,configFBGarray.numFBG-1);
    psi = 2*pi*rand(N,2*(configFBGarray.numFBG-1))-pi;

    configFBGarray.channelDelay(15) = cd(ii);

    %% start simulator
    opticalWave = simFBGArrayStepwise2(psi,theta,numTimeFrames2sim, configFBGarray);
    [adc0, adc120, adc240] = simOpticalSensors(opticalWave, configFBGarray);
    [theta0, theta120, theta240] = decodeADC(adc0, adc120, adc240, configFBGarray);
    [hist0(ii,:,:), hist120(ii,:,:), hist240(ii,:,:), robustSTD(ii,:)] = histEstimationError(theta, theta0, theta120, theta240, bins);

end

Yspacing = 0.5;
figure
hold on
for ii = 2:length(cd)
    plot(robustSTD(1,1:22) - (ii-1)*Yspacing,'b', 'lineWidth', 3);
    plot(robustSTD(ii,1:22) - (ii-1)*Yspacing, 'r' )
    
end

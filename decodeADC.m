function [theta0, theta120, theta240] = decodeADC(adc0, adc120, adc240, config)

samplesPerChannel = config.pulseLength/config.samplingPeriodADC;
samplesLeadIn = config.channelDelayLeadIn/config.samplingPeriodADC;


%% points spacing should not be changed. The 5 points must cover 1 period of the modualtion frequency. Offset can be tuned however. 
points = round(samplesPerChannel/10:samplesPerChannel/5:(samplesPerChannel-samplesPerChannel/10));
numTimeFrames = floor(size(adc0,2)/samplesPerChannel);
numPulses = size(adc0,1);

theta0 = zeros(numPulses,numTimeFrames);
theta120 = zeros(numPulses,numTimeFrames);
theta240 = zeros(numPulses,numTimeFrames);

channelStartPoints = [2*samplesLeadIn, cumsum(config.channelDelay)/config.samplingPeriodADC+2*samplesLeadIn];
channelStartPoints = [channelStartPoints, channelStartPoints(end) + (1:numTimeFrames - length(config.channelDelay)-1)*samplesPerChannel];
ii = length(channelStartPoints);
while channelStartPoints(ii)+samplesPerChannel>size(adc0,2)
    ii = ii - 1;
end
channelStartPoints(ii+1:end) = [];
numTimeFrames = length(channelStartPoints);

offset = zeros(1,numTimeFrames);
for ii = 1:length(config.channelDelay)
    if config.channelDelay(ii) < 220 
        offset(ii) = -1;
    end
end

figure
subplot(311)
hold on
plot(adc0(1,:))
plot(channelStartPoints, adc0(1,channelStartPoints), 'r*');
for ii = 1:numTimeFrames
    plot(channelStartPoints(ii)+points, adc0(1,channelStartPoints(ii)+points+offset(ii)), 'k*');
end
subplot(312)
hold on
plot(adc120(1,:))
plot(channelStartPoints, adc120(1,channelStartPoints), 'r*');
for ii = 1:numTimeFrames
    plot(channelStartPoints(ii)+points, adc120(1,channelStartPoints(ii)+points+offset(ii)), 'k*');
end
subplot(313)
hold on
plot(adc240(1,:))
plot(channelStartPoints, adc240(1,channelStartPoints), 'r*');
for ii = 1:numTimeFrames
    plot(channelStartPoints(ii)+points, adc240(1,channelStartPoints(ii)+points+offset(ii)), 'k*');
end





for ii = 1:numPulses
    for jj = 1:numTimeFrames
        theta0(ii,jj) = estimatePhase5pnt2(adc0(ii,channelStartPoints(jj)+points+offset(jj)));
        theta120(ii,jj) = estimatePhase5pnt2(adc120(ii,channelStartPoints(jj)+points+offset(jj)));
        theta240(ii,jj) = estimatePhase5pnt2(adc240(ii,channelStartPoints(jj)+points+offset(jj)));
    end
end


xxx = (config.samplesPerOpticalWavePeriod:config.samplesPerOpticalWavePeriod:config.channelDelay(1)) - 9;
for ii = 1:numTimeFrames
    theta0(:,ii) = theta0(:,ii) + 2*pi*config.fm*(xxx(points(1)+offset(ii))+0.5);
    theta120(:,ii) = theta120(:,ii) + 2*pi*config.fm*(xxx(points(1)+offset(ii))+0.5);
    theta240(:,ii) = theta240(:,ii) + 2*pi*config.fm*(xxx(points(1)+offset(ii))+0.5);
end

% - 2*pi*config.fm*(points(1)-1)


% + 2*pi*config.fm*config.samplingPeriodADC*points(1)



% xxx = (config.samplesPerOpticalWavePeriod:config.samplesPerOpticalWavePeriod:config.channelDelay(1)) - 9;
% figure
% hold on
% plot(xxx,adc0(ii,channelStartPoints(jj)+(1:22)))
% plot(xxx(points),adc0(ii,channelStartPoints(jj)+points),'*')
% A = max(adc0(ii,channelStartPoints(jj)+(1:22))) - min(adc0(ii,channelStartPoints(jj)+(1:22)));
% O = (max(adc0(ii,channelStartPoints(jj)+(1:22))) + min(adc0(ii,channelStartPoints(jj)+(1:22))))/2;
% plot(xxx(points(1):points(end)),A/2*cos(2*pi*config.fm*((0:points(end)-points(1))*config.samplingPeriodADC) - theta0(ii,jj))+O)
% %plot(xxx,A/2*cos(2*pi*config.fm*(xxx-1) + theta0(ii,jj) - 2*pi*config.fm*(points(1)-1)) + O )
% plot(xxx,A/2*cos(2*pi*config.fm*((0:21)*config.samplingPeriodADC) - theta0(ii,jj) - 2*pi*config.fm*xxx(points(1)))+O)




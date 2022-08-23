clear all
close all

N = 22;
t = 0:0.0625:0.0625*(N-1);

mod_fun = 2*pi*1*[t, zeros(1,2*N)] + 2*pi*2*[zeros(1,N), t, zeros(1,N)] + 2*pi/3*[zeros(1,N), ones(1,N), zeros(1,N)] + 2*pi*3*[zeros(1,2*N), t] + 4*pi/3*[zeros(1,2*N), ones(1,N)];
mod_fun = mod_fun';

ph = 0:0.0001:4*pi;
ph = ph';
%ph = (pi-0.01)*ones(10000,1);

Niter = 150;
meanEstimationError = zeros(Niter,1);
meanEstimationErrorCur = zeros(Niter,1);
meanEstimationErrorFreq = zeros(Niter,1);

for iterSim = 1:Niter
    measurementNoise = (iterSim-1)/100;
    v = randn(3*N,length(ph));
    
    Aa = zeros(length(ph),1);
    Ab = zeros(length(ph),1);
    Ac = zeros(length(ph),1);
    Acur = zeros(length(ph),1);
    
    for ii = 1:length(ph)
        m = sin(mod_fun + ph(ii)) + measurementNoise*v(:,ii);
        
        Ma = fft(m(1:1+15));
        Mb = fft(m(23:23+15));
        Mc = fft(m(45:45+15));
        Aa(ii) = angle(Ma(2));
        Ab(ii) = angle(Mb(3));
        Ac(ii) = angle(Mc(4));
        
        
        S1 = m(1);
        S2 = m(5);
        S3 = m(9);
        S4 = m(13);
        S5 = m(17);
        CR = 2*S3 - S5 - S1;
        SR = 2*(S2 - S4);
        Acur(ii) = -atan2(CR,SR);

    end
    
    Afreq = Aa + pi/2;
    
    Aa = mod(Aa + pi/2, 2*pi);
    Ab = mod(Ab - 2*pi/3 + pi/2, 2*pi);
    Ac = mod(Ac - 4*pi/3 + pi/2, 2*pi);
    Amat = wrapMed([Aa, Ab, Ac]);
    Amed = median(Amat, 2);
    
    dAa = (Aa-Amed);
    dAb = (Ab-Amed);
    dAc = (Ac-Amed);
    dAmat = wrap(wrapMed([dAa, dAb, dAc]));
    maxDist2med = max([dAa, dAb, dAc], [], 2);
    poi = find(maxDist2med>2);
    meanA = mean(Amat, 2);
    meanA(poi) = Amed(poi);
    
    uA = unwrap(wrap(meanA));
    uAcur = unwrap(Acur);
    uAfreq = unwrap(Afreq);
    figure
    hold on
    plot(Aa)
    plot(Ab)
    plot(Ac)
    
    meanEstimationError(iterSim) = mean(abs(uA-ph));
    meanEstimationErrorCur(iterSim) = mean(abs(uAcur-ph));
    meanEstimationErrorFreq(iterSim) = mean(abs(uAfreq-ph));
    
%     if meanEstimationError(iterSim) > 5
%         fprintf( 'here')
%     end
    
    makePlot = 0;
    if makePlot == 1
        uAa = unwrap(Aa);
        uAb = unwrap(Ab);
        uAc = unwrap(Ac);
        
        figure
        hold on
        plot(uAa)
        plot(uAb)
        plot(uAc)
        plot(uA)
    end
    
end

idx1 = 38;
idx2 = 71;
idx3 = 112;
figure
hold on
measurementNoise = 0:0.01:Niter/100-0.01;
plot(measurementNoise(1:idx1), meanEstimationErrorCur(1:idx1))
plot(measurementNoise(1:idx2), meanEstimationErrorFreq(1:idx2))
plot(measurementNoise(1:idx3), meanEstimationError(1:idx3))
xlabel('Standard Deviation of Measurement Noise')
ylabel('Mean Estimation Error')
title('Mean Estimation Error' );
legend('5 point', 'DFT', 'Robust DFT' ); 


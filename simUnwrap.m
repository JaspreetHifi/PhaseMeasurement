close all
clear all

%t = 0:0.01:1.25;
%t = linspace(0,1.25,22);
%t2 = linspace(0,2.5,44);
%t3 = linspace(0,3.75,66);
t = 0:0.0625:0.0625*(22-1);
t2 = 0:0.0625:0.0625*(44-1);
t3 = 0:0.0625:0.0625*(66-1);



f1 = 1;
f2 = 2;
ph0 = 0:0.001:4*pi;
%ph0 = pi*ones(1,10000)-0.01;

%phSearch = 0:0.01:2*pi-0.01;
%phSearch2 = 0:0.01:2*pi-0.01;

%ddHist = zeros(1,length(ph)-1);
mod_fun = 2*pi*f1*[t2(1:22), zeros(1,22)] + 2*pi*f2*[zeros(1,22), t2(23:44)] + pi/2*[zeros(1,22), ones(1,22)];
mod_fun2 = 2*pi*0.5*[t2(1:22), zeros(1,22)] + 2*pi*1*[zeros(1,22), t2(23:44)] + pi/4*[zeros(1,22), ones(1,22)];
mod_fun3 = 2*pi*0.5*[t3(1:22), zeros(1,44)] + 2*pi*1*[zeros(1,22), t, zeros(1,22)] ...
    + 0*pi/4*[zeros(1,22), ones(1,22), zeros(1,22)] + 2*pi*2*[zeros(1,44), t] ...
    + pi/4*[zeros(1,44), ones(1,22)];

hist_erratic = zeros(1,length(ph0));
hist_erraticA5a = zeros(1,length(ph0));
hist_erraticA5b = zeros(1,length(ph0));


for iter_sim = 1:90
    fprintf('Iteration %d/%d. ', iter_sim, 1000);
    
    %proc_noise = (iter_sim-1)/10;
    proc_noise = 0;
    meas_noise = (iter_sim-1)/100;
    
    ph = ph0 + proc_noise*randn(1,length(ph0));
    
    for ii = 1:length(ph)
        v = meas_noise*randn(1,length(t));
        v2 = meas_noise*randn(1,length(t2));
        v3 = meas_noise*randn(1,length(t3));
        
        m = sin(2*pi*f1*t + ph(ii) ) + v;
        %        m2 = sin(2*pi*f1*t.^2 + 2*pi*f1*t + ph(ii) ) + v;
        m3 = sin(mod_fun + ph(ii)) + v2;
        m4 = sin(mod_fun2 + ph(ii)) + v2;
        m5 = sin(mod_fun3 + ph(ii)) + v3;
        %m3b = sin(2*pi*f2*t + pi + ph(ii)) + v;
        %m3 = [m, m3b];
        %         objective = zeros(length(phSearch),1);
        %         objective2 = zeros(length(phSearch),1);
        %         for jj = 1:length(phSearch)
        %             bf = sin(2*pi*f1*t + phSearch(jj));
        %             objective(jj) = sum((bf-m).^2);
        %         end
        %         for jj = 1:length(phSearch2)
        %             bf2 = sin(2*pi*f1*t.^2 + 2*pi*f1*t + phSearch2(jj));
        %             objective2(jj) = sum((bf2-m2).^2);
        %         end
        %         [~, minInd] = min(objective);
        %         A2(ii) = phSearch(minInd);
        %
        %         [~, minInd] = min(objective2);
        %         A4(ii) = phSearch2(minInd);
        tic
        M = fft(m(1:16));
        A2(ii) = angle(M(2));
        timer1(ii) = toc;
        M5b = fft(m5(23:23+15));
        M5c = fft(m5(45:45+15));
        A5a(ii) = angle(M5b(2));
        A5b(ii) = angle(M5c(3));
        
        %plot(t,m)
        S1 = m(1);
        S2 = m(5);
        S3 = m(9);
        S4 = m(13);
        S5 = m(17);
        CR = 2*S3 - S5 - S1;
        SR = 2*(S2 - S4);
        A(ii) = -atan2(CR,SR);
        
%         if (ph(ii) - unwrap(A(ii))) < 5
%            fprintf('here');
%         end
    end
    fprintf('Mean excution time: %1.4e\n', mean(timer1) )
    A2 = A2 + pi/2;
    A5a = A5a + pi/2;
    A5b = A5b + pi/4;
    
    if mean(ph-unwrap(A2)) < -6
        A2 = A2 - 2*pi;
    elseif mean(ph-unwrap(A2)) > 6
        A2 = A2 + 2*pi;
    end
    
    %     if mean(ph-unwrap(A4)) < -6
    %         A4 = A4 - 2*pi;
    %     elseif mean(ph-unwrap(A4)) > 6
    %         A4 = A4 + 2*pi;
    %     end
    makePlot = 0;
    if makePlot == 1
        figure
        
        subplot(3,1,1)
        hold on
        plot(ph)
        plot(A)
        plot(A2)
        %plot(A2)
        %plot(A4)
        axis('tight')
        title( 'Estimated Phase' );
        xlabel('Sample')
        ylabel('Phase')
        legend( 'actual', 'mthd 1', 'mthd 2', 'mthd 3', 'mthd 4' );
        
        subplot(3,1,2)
        hold on
        plot(ph)
        plot(unwrap(A))
        plot(unwrap(A2))
        %plot(unwrap(A2))
        %plot((A4))
        axis('tight')
        title( 'Estimated Phase (Unwrapped)' );
        xlabel('Sample')
        ylabel('Phase')
        legend( 'actual', 'mthd 1', 'mthd 2', 'mthd 3', 'mthd 4' );
        
        subplot(3,1,3)
        hold on
        plot(ph-ph0)
        plot(ph-unwrap(A))
        plot(ph-unwrap(A2))
        %plot(ph-unwrap(A2))
        %plot(ph-(A4))
        axis('tight')
        title( 'Difference Between Estimated and Actual Phase' );
        xlabel('Sample')
        ylabel('Diff (rad)')
        legend( 'actual', 'mthd 1', 'mthd 2', 'mthd 3', 'mthd 4' );
    end
    meanError(iter_sim) = mean(abs(ph0-unwrap(A)));
    meanError2(iter_sim) = mean(abs(ph0-unwrap(A2))); 
    
    
    
    tabulateErraticStrain = 0;
    if tabulateErraticStrain == 1
        Aunwrap = unwrap(A);
        A5aunwrap = unwrap(A5a);
        A5bunwrap = unwrap(A5b);
        
        ii = 1;
        while ii <= length(A)
            %         if Aunwrap(ii)-ph0(ii) > 2*pi
            %             Aunwrap(ii:end) = Aunwrap(ii:end) - 2*pi;
            %             hist_erratic(ii) = hist_erratic(ii) + 1;
            %         elseif Aunwrap(ii) - ph0(ii) < -2*pi
            %             Aunwrap(ii:end) = Aunwrap(ii:end) + 2*pi;
            %             hist_erratic(ii) = hist_erratic(ii) + 1;
            %         end
            
            
            if A5aunwrap(ii)-ph0(ii) > 2*pi
                A5aunwrap(ii:end) = A5aunwrap(ii:end) - 2*pi;
                hist_erraticA5a(ii) = hist_erraticA5a(ii) + 1;
            elseif A5aunwrap(ii) - ph0(ii) < -2*pi
                A5aunwrap(ii:end) = A5aunwrap(ii:end) + 2*pi;
                hist_erraticA5a(ii) = hist_erraticA5a(ii) + 1;
            end
            
            if A5bunwrap(ii)-ph0(ii) > 2*pi
                A5bunwrap(ii:end) = A5bunwrap(ii:end) - 2*pi;
                hist_erraticA5b(ii) = hist_erraticA5b(ii) + 1;
            elseif A5bunwrap(ii) - ph0(ii) < -2*pi
                A5bunwrap(ii:end) = A5bunwrap(ii:end) + 2*pi;
                hist_erraticA5b(ii) = hist_erraticA5b(ii) + 1;
            end
            
            ii = ii + 1;
        end
    end
    %figure
%     medWin = 10;
%     Aunwrapped = unwrap(A);
%     for ii = 1 + medWin:length(A)-medWin
%         Amed(ii-medWin) = median(Aunwrapped(ii-medWin:ii+medWin));
%     end
%     temp1 = find(abs(diff(Amed))> 2);
%     hist_erratic(temp1) = hist_erratic(temp1) + 1;
% 
%     A5aunwrapped = unwrap(A5a);
%     for ii = 1 + medWin:length(A)-medWin
%         Amed(ii-medWin) = median(A5aunwrapped(ii-medWin:ii+medWin));
%     end
%     temp1 = find(abs(diff(Amed))> 2);
%     hist_erraticA5a(temp1) = hist_erraticA5a(temp1) + 1;
%     
%     A5bunwrapped = unwrap(A);
%     for ii = 1 + medWin:length(A)-medWin
%         Amed(ii-medWin) = median(A5bunwrapped(ii-medWin:ii+medWin));
%     end
%     temp1 = find(abs(diff(Amed))> 2);
%     hist_erraticA5b(temp1) = hist_erraticA5b(temp1) + 1;
    

    
    %num_erratic1(iter_sim) = length(temp1);
    %num_erratic2(iter_sim) = length(find(abs(diff(unwrap(A2)) - diff(ph))> 5));
    %num_erratic3(iter_sim) = length(find(abs(diff(unwrap(A2)) - diff(ph))> 5));
    %num_erratic4(iter_sim) = length(find(abs(diff((A4)) - diff(ph))> 5));
    
    %mean_err1(iter_sim) = mean(abs(ph-unwrap(A)));
    %mean_err2(iter_sim) = mean(abs(ph-unwrap(A2)));
    %mean_err3(iter_sim) = mean(abs(ph-unwrap(A2)));
    %mean_err4(iter_sim) = mean(abs(ph-(A4)));
    
end

% binSize = 100;
% histAmed = zeros(1,floor(length(hist_erratic)/binSize));
% for ii = 1:floor(length(hist_erratic)/binSize)
%     histAmed(ii) = sum(hist_erratic((ii-1)*binSize+1:ii*binSize));
% end

binSize = 100;
histAmedA5a = zeros(1,floor(length(hist_erraticA5a)/binSize));
for ii = 1:floor(length(hist_erraticA5a)/binSize)
    histAmedA5a(ii) = sum(hist_erraticA5a((ii-1)*binSize+1:ii*binSize));
end

binSize = 100;
histAmedA5b = zeros(1,floor(length(hist_erraticA5b)/binSize));
for ii = 1:floor(length(hist_erraticA5b)/binSize)
    histAmedA5b(ii) = sum(hist_erraticA5b((ii-1)*binSize+1:ii*binSize));
end

figure
plot(ph(1:binSize:end-binSize),histAmed)

figure
plot(ph(1:binSize:end-binSize),histAmedA5a)
figure
plot(ph(1:binSize:end-binSize),histAmedA5b)
% figure
% hold on
% plot(num_erratic1)
% plot(num_erratic2)
% %plot(num_erratic3)
% %plot(num_erratic4)
%
% figure
% hold on
% plot(mean_err1)
% plot(mean_err2)
% %plot(mean_err3)
% %plot(mean_err4)





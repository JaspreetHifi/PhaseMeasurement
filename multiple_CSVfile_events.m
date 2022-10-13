%%
clear all

clc
ch = 130;
fileList = dir( 'D:\Hifi\events\crossPlotData\Event*.csv');
% fileList = dir( 'D:\Hifi\events\CPd\July201.xlsx');

numFiles = length(fileList);

STCdate = [];
STCch = [];

pigdate = [];
pigch = [];

PSDdate = [];
PSDch = [];

Acousticdate = [];
Acousticch = [];

Velocitydate = [];
Velocitych = [];

for ii = 1:numFiles
    fprintf( 'Loading file %d/%d: %s.\n', ii, numFiles, fileList(ii).name );
    [~, ~, raw] = xlsread(['D:\Hifi\events\crossPlotData\' fileList(ii).name]);
    for ii=2:size(raw,1)
        if length(raw{ii,2}) == 4
            if raw{ii,2}(2:4) == 'STC'
                STCdate = [STCdate;datenum([raw(ii,3)])];
                STCch = [STCch;raw(ii,5)];
            elseif raw{ii,2}(2:4) == 'Pig'
                pigdate = [pigdate;datenum([raw(ii,3)])];
                pigch = [pigch;raw(ii,5)];
            end
        else
            if raw{ii,2}(2:5) == 'Acou'
                Acousticdate = [Acousticdate;datenum([raw(ii,3)])];
                Acousticch = [Acousticch;raw(ii,5)];
            elseif raw{ii,2}(2:5) == 'Leak'
                if length(raw{ii,2}) == 5
                    PSDdate = [PSDdate;datenum([raw(ii,3)])];
                    PSDch = [PSDch;raw(ii,5)];
                elseif length(raw{ii,2}) == 11
                    Velocitydate = [Velocitydate;datenum([raw(ii,3)])];
                    Velocitych = [Velocitych;raw(ii,5)];
                end
            end
            
        end
    end
    
    
end

m=1;
for ii=1:length(STCdate)
    if ~ischar(cell2mat(STCch(ii)))
        STCdate1(m) = STCdate(ii);STCch1(m) = cell2mat(STCch(ii));m=m+1;
    else
        chGroup = str2num(cell2mat(STCch(ii)));
        STCdate1(m:m+length(chGroup)-1) = STCdate(ii);STCch1(m:m+length(chGroup)-1) = chGroup;m = m+length(chGroup);
    end
end
chGroup=[];

m=1;
for ii=1:length(pigdate)
    if ~ischar(cell2mat(pigch(ii)))
        pigdate1(m) = pigdate(ii);pigch1(m) = cell2mat(pigch(ii));m=m+1;
    else
        chGroup = str2num(cell2mat(pigch(ii)));
        pigdate1(m:m+length(chGroup)-1) = pigdate(ii);pigch1(m:m+length(chGroup)-1) = chGroup;m = m+length(chGroup);
    end
end
chGroup=[];

m=1;
for ii=1:length(PSDdate)
    if ~ischar(cell2mat(PSDch(ii)))
        PSDdate1(m) = PSDdate(ii);PSDch1(m) = cell2mat(PSDch(ii));m=m+1;
    else
        chGroup = str2num(cell2mat(PSDch(ii)));
        PSDdate1(m:m+length(chGroup)-1) = PSDdate(ii);PSDch1(m:m+length(chGroup)-1) = chGroup;m = m+length(chGroup);
    end
end
chGroup=[];

m=1;
for ii=1:length(Acousticdate)
    if ~ischar(cell2mat(Acousticch(ii)))
        Acousticdate1(m) = Acousticdate(ii);Acousticch1(m) = cell2mat(Acousticch(ii));m=m+1;
    else
        chGroup = str2num(cell2mat(Acousticch(ii)));
        Acousticdate1(m:m+length(chGroup)-1) = Acousticdate(ii);Acousticch1(m:m+length(chGroup)-1) = chGroup;m = m+length(chGroup);
    end
end
chGroup=[];

m=1;
for ii=1:length(Velocitydate)
    if ~ischar(cell2mat(Velocitych(ii)))
        Velocitydate1(m) = Velocitydate(ii);Velocitych1(m) = cell2mat(Velocitych(ii));m=m+1;
    else
        chGroup = str2num(cell2mat(Velocitych(ii)));
        Velocitydate1(m:m+length(chGroup)-1) = Velocitydate(ii);Velocitych1(m:m+length(chGroup)-1) = chGroup;m = m+length(chGroup);
    end
end
chGroup=[];
%%
figure,
plot(STCdate1,STCch1,'*r');hold on
% plot(PSDdate1,PSDch1,'*k'); hold on
% plot(pigdate1,pigch1,'*c'); hold on
% plot(Acousticdate1,Acousticch1,'*','color',[0, 0.5, 0]); hold on
% plot(Acousticdate1,Acousticch1,'*g'); hold on
% plot(Velocitydate1,Velocitych1,'*','color',[0.9100    0.4100    0.1700]);
axis tight
% legend('STC events','PSD events','Acoustic events','Velocity events')
% legend('PSD events','Acoustic events')
%  xlim([datenum('15-Aug-2022 00:00:00'),datenum('19-Aug-2022 09:00:00')])
%   title('Pipestone events - June 12 - July 5');
xData1 = linspace(STCdate1(1),STCdate1(end),100);
%   xData1 = linspace(datenum('15-Aug-2022 00:00:00'),datenum('19-Aug-2022 09:00:00'),50);
ax = gca;ax.XTick = xData1;ax.XTickLabel = datestr(xData1);ax.XTickLabelRotation = 270;
% ylim([1,ch]);
xlabel('Date');ylabel('Channel');grid on; box on

%%

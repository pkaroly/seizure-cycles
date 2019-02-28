clear
clc
close all

% Patients
filePath = 'STMatlab/';
allPts = dir([filePath 'Patient*']);
Npt = length(allPts);
id = cell(1,Npt);

resWin = [6 12 18 24 ...        % go up in 6 hour blocks to one day
    48:24:30*24 ...                 % go up in one day blocks to one month
    840:7*24:52*7*24 ...           % then go up in one week blocks to one year
    78*7*24 104*7*24];              % 18 months, 2 years

circR = nan(Npt, length(resWin));
circRtest = circR;
circRtest2 = circR;
gender = nan(Npt,1);
nPtSz = nan(Npt,1);

dayPeak = zeros(1,24);
weekPeak = zeros(1,7);

daySz = 0;
totalSz = 0;

circDaySz = 0;
circTotalSz = 0;

vmTest = 0;
vmCount = 0;

propCluster = nan(Npt,1);

avDur = nan(Npt,1);
for ind = 1:Npt
    
    load([filePath allPts(ind).name]);
    id{ind} = curPt;
    
    SzTimes = 24 * (szTime - szTime(1));  % seizure times in hours
    [SzTimes, I] = sort(SzTimes);
    szHour = szHour(I);
    szTime = szTime(I);
    szDur = szDur(I);
    
    % remove seizures at 1am
    invalid = szHour == 1;
    SzTimes(invalid) = [];
    szTime(invalid) = [];
    szDur(invalid) = [];
    
    % remove seizure times before 2007
    szYear = datevec(szTime);
    invalid = szYear(:,1) < 2007 & szYear(:,1) > 2017;
    SzTimes(invalid) = [];
    szTime(invalid) = [];
    
    % remove seizures in same hour
    [~, valid] = unique(floor(SzTimes));
    SzTimes = SzTimes(valid);
    szTime = szTime(valid);
    
    if length(SzTimes) < 100
        continue;
    end
    
    % check time
    SzDate = datevec(szTime);
    Hour = SzDate(:,4);
    daySz = daySz + sum(Hour > 9 & Hour < 21);
    totalSz = totalSz + length(Hour);
    
    % check max res
    maxRes = max(SzTimes) / 10;
    resWinCur = resWin;
    resWinCur(resWin > maxRes) = [];
    
    % save pt gender
    gender(ind) = ptGender;
    nPtSz(ind) = length(SzTimes);
    avDur(ind) = (max(szTime) - min(szTime)) / 365;
    
    % windows in hr
    count = 1;
    for n = resWinCur
        
        % remove clusters
        if n > 48
            leadTime = 24;
        else
            leadTime = 1;
        end
        ISI = [leadTime ; diff(SzTimes)];
        SzTimeLead = SzTimes(ISI > leadTime);
        
        % actual variance
        SzPhase = 2 * pi * mod(SzTimes, n) / n; 
%         SzPhase = 2 * pi * mod(SzTimeLead, n) / n;
        circR(ind, count) = circ_r(SzPhase);
        circRtest(ind, count) = circ_rtest(SzPhase);
        circRtest2(ind, count) = circ_otest(SzPhase);
        
%         p = watson(SzPhase');
%         vmCount = vmCount + 1;
%         if p < 0.05
%             vmTest = vmTest + 1;
%         end
        
        if circRtest(ind, count) < 0.05 && n == 24
            SzDate = datevec(szTime);
            Hour = SzDate(:,4);
            circDaySz = circDaySz + sum(Hour > 9 & Hour < 21);
            circTotalSz = circTotalSz + length(Hour);
            [~,I] = max(histcounts(Hour,0:24));
            dayPeak(I) = dayPeak(I) + 1;
        end
        
        if circRtest(ind, count) < 0.05 && n == 24*7
            Day = weekday(szTime);
            [~,I] = max(histcounts(Day,1:8));
            weekPeak(I) = weekPeak(I) + 1;
        end
        
        count = count + 1;
    end
    
end

%% plots
font = 'arial';
fsize = 8;

count = 1;
x = resWin;
close all
figure;

% yyaxis left
invalid = sum(isnan(circR),2) == size(circR,2);
circR(invalid,:) = [];
gender(invalid) = [];
nPtSz(invalid) = [];
avDur(invalid) = [];
circRtest(invalid,:) = [];
circRtest2(invalid,:) = [];
id =  categorical(id);
id(invalid) = [];

[~,ind] = sort(circR(:,10),'descend');

imagesc(circR(:,:));
C = brewermap(255,'PuBu');
colormap(C);

%%
resWin = [6 12 18 24 ...        % go up in 6 hour blocks to one day
    48:24:28*24 ...                 % go up in one day blocks to one month
    840:7*24:52*7*24 ...           % then go up in one week blocks to one year
    78*7*24 104*7*24];              % 18 months, 2 years
xvals = [4, 10, 17, 31, 39, 51];
xlab = {'day', 'week', '2wk', 'month', '3 month', '6 month'};
set(gca,'box','off','tickdir','out', ...
    'ycolor','k', 'xlim',[1 40],...
    'xtick',xvals, 'xticklabels',xlab,...
    'ytick',[500 1000], ...
    'fontsize',fsize,'fontname','arial');
ylabel('Patients','fontsize',fsize);
xlabel('Time Scale','fontsize',fsize);

hold on;
yyaxis right
rMu = nanmean(circR);
rSEM = 1.96 * nanstd(circR) ./ sqrt(sum(~isnan(circR)));
plot(rMu,'linewidth',1.5,'color','k');
plot(rMu + rSEM,'--','linewidth',1,'color','k');
plot(rMu - rSEM,'--','linewidth',1,'color','k');
ylabel('R-value','fontsize',fsize);
% plot(nanmean(pVals),'linewidth',1.5,'color','w');
% ylabel('Mean','fontsize',fsize);

colorbar('color','k','fontsize',fsize,'ytick',[0.01 0.5 0.9],'yticklabel',[0 0.5 1]);
set(gca,'box','off','ytick',[]);
set(gca,'xcolor','k','ycolor','k','fontsize',fsize);
% set(gcf,'color',[63 63 63]/255,'inverthardcopy','off');

set(gcf,'paperunits','centimeters','paperposition',[0 0 15 10]);
% print(gcf,'ST_RSig','-dpng','-r300');
print(gcf,'ST_RVals','-dpng','-r300');

%%
% vals
simThresh = 0.2644; % from simulation
alpha = 0.05 / 80;
a0 = sum(sum(circRtest2 < alpha,2) > 0);
a1 = sum(sum(circRtest2 < alpha,2) > 1);
c = sum(sum(circRtest2(:,24:end) < alpha,2) > 0);
b = sum(sum(circRtest2(:,10) < alpha,2) > 0);
a = sum(sum(circRtest2(:,4) < alpha,2) > 0);

% over 3 weeks
c2 = sum(sum(circRtest(:,24:39) < 0.0033,2) > 0);
% one week
b2 = sum(sum(circRtest(:,10) < 0.05,2) > 0);
% one day
a2 = sum(sum(circRtest(:,4) < 0.05,2) > 0);

fprintf('any: %d, any+1: %d, day: %d, week: %d, long: %d\n',a0,a1,a,b,c);
fprintf('any: %.0f, any+1: %.0f, day: %.0f, week: %.0f, long: %.0f\n',100*[a0,a1,a,b,c]/sum(~invalid));

%%
nM = sum(gender == 1);
nF = sum(gender == 0);
pM = sum(gender(circR(:,4) > simThresh) == 1) / nM;
pW = sum(gender(circR(:,4) > simThresh) == 0) / nF;

p = a / sum(~invalid);
r = (pM - pW) /  sqrt((p / nM) + (p / nF));

NactualPt = Npt - sum(isnan(avDur));
N = sum(~isnan(gender));
mDay = sum(gender(circR(:,4) > simThresh) == 1) / a;
ciMDay = CI(mDay,a);
fDay = sum(gender(circR(:,4) > simThresh) == 0) / a;
ciFDay = CI(fDay,a);
% disp(mDay - fDay)
mWeek = sum(gender(circR(:,10) > simThresh) == 1) / b;
ciMWeek = CI(mWeek,b);
fWeek = sum(gender(circR(:,10) > simThresh) == 0) / b;
ciFWeek = CI(fWeek, b);
% disp(mWeek - fWeek)
[ind, ~] = find(circR(:,24:end) > simThresh);
mMonth = sum(gender(unique(ind)) == 1) / c;
ciMMonth = CI(mMonth,c);
fMonth = sum(gender(unique(ind)) == 0) / c;
ciFMonth = CI(fMonth,c);
% disp(mMonth - fMonth)

%%
close all
B = bar([mDay, fDay ; mWeek, fWeek ; mMonth, fMonth]);
hold on
line([0.85 0.85], [mDay+ciMDay mDay-ciMDay],'color','k','linewidth',1.5);
line([1.15 1.15], [fDay+ciFDay fDay-ciFDay],'color','k','linewidth',1.5);
line([1.85 1.85], [mWeek+ciMWeek mWeek-ciMWeek],'color','k','linewidth',1.5);
line([2.15 2.15], [fWeek+ciFWeek fWeek-ciFWeek],'color','k','linewidth',1.5);
line([2.85 2.85], [mMonth+ciMMonth mMonth-ciMMonth],'color','k','linewidth',1.5);
line([3.15 3.15], [fMonth+ciFMonth fMonth-ciFMonth],'color','k','linewidth',1.5);

B(1).FaceColor = [0 0 0];
B(2).FaceColor = [0.5 0.5 0.5];
set(gca,'box','off',...
    'ytick',[0 0.6],'xticklabel',{'Daily', 'Weekly', '>3 Weeks'},...
    'fontsize',fsize,'fontname',font);
ylabel('Gender ratio','fontsize',fsize,'fontname',font);
legend({'male', 'female'},'location','northoutside');
set(gcf,'paperunits','centimeters','paperposition',[0 0 8 6]);
print(gcf,'MaleFemale','-dtiff','-r300');

%%
close all
subplot(211)
bar(dayPeak, 'facecolor','k');
set(gca,'box','off',...
    'ytick',[0 120],...
    'xtick',[0 12 24],'xticklabel',{'midnight','midday','midnight'},...
'fontsize',fsize,'fontname',font)
xlabel('Hour','fontsize',fsize,'fontname',font)
title('A. Peak for Circadian Patients','fontsize',fsize,'fontname',font,'fontweight','bold');

N = sum(dayPeak);
CIs = N * CI(dayPeak ./ N, N);
hold on;
line(repmat(1:24,2,1), [dayPeak - CIs ; dayPeak + CIs],'color','k','linewidth',1.5);

subplot(212);
bar(weekPeak,'facecolor','k');
set(gca,'box','off',...
    'ytick',[0 60],...
    'xticklabel',{'S','M','T','W','T','F','S'},...
    'fontsize',fsize,'fontname',font)
xlabel('Day','fontsize',fsize,'fontname',font)
title('B. Peak for Circaseptan Patients','fontsize',fsize,'fontname',font,'fontweight','bold');
ylabel('N patients','fontsize',fsize,'fontname',font);

N = sum(weekPeak);
CIs = N * CI(weekPeak ./ N, N);
hold on;
line(repmat(1:7,2,1), [weekPeak - CIs ; weekPeak + CIs],'color','k','linewidth',1.5);

set(gcf,'paperunits','centimeters','paperposition',[0 0 8 8]);
print(gcf,'Peaks','-dtiff','-r300');
clc
clear
close all

col1 = [204 0 0]/255;
col2 = [255 128 0]/255;
col3 = [0 102 204]/255;

load('CycleIDs');
load('Profiles');
load('StudyIDs');

id = categorical(Profiles.Unlinked_ID);
types = cell(1,10);

count = 1;
allDat = nan(1,10);
dayDat = allDat;
wkDat = allDat;
longDat = allDat;

keepInd = ismember(id,studyId);

for cat = 2:16
    temp = categorical(Profiles.(cat));
    
    catName = char(temp(1));
    catName = strrep(catName,'_',' ');
    
    if strcmp(catName,'Stroke') || ...
        strcmp(catName,'Heart Attack') || ...
        strcmp(catName,'Alzheimers')
    
        continue;
    end
    
    temp = temp(keepInd);
    idTemp = id(keepInd);
    tempID = idTemp(~isundefined(temp));
    
    
    if strcmp(catName,'Lack Of Oxygen During Birth')  || ...
            strcmp(catName,'Maternal Drug Or Alcohol Abuse')  || ...
            strcmp(catName, 'Brain Malformations')  || ...
            strcmp(catName,'Brain Injury During Fetal Development')
        types{count} = 'Foetal';
        allDat(count) = nansum([allDat(count),sum(~isundefined(temp)) - 1]);
        dayDat(count) = nansum([dayDat(count),sum(ismember(daily,tempID))]);
        wkDat(count) = nansum([wkDat(count),sum(ismember(weekly,tempID))]);
        longDat(count) = nansum([longDat(count),sum(ismember(long,tempID))]);
    else
        dayDat(count) = sum(ismember(daily,tempID));
        wkDat(count) = sum(ismember(weekly,tempID));
        longDat(count) = sum(ismember(long,tempID));
        types{count} = catName;
        allDat(count) = sum(~isundefined(temp)) - 1;
        count = count + 1;
    end
    
    %     idByCond{cat-1} = id(~isundefined(temp));
end
%%
XLab = types;
XLab{4} = 'Hematomas';
XLab{9} = 'Electrolyte';

A = getCustomAxesPos(3,1,0,0.05);

YLim = [-0.1 0.1];
YLab = [-10 10];

a = allDat / sum(allDat);
CI1 = CI(a,sum(allDat),0.05/12/4);
CI2 = CI(a,sum(allDat),0.01/12/4);

axes(A(1))
b = dayDat / sum(dayDat);
bar(b-a,'facecolor',col1,'facealpha',0.6);
set(gca,'box','off','xtick',[],'xlim',[0 11],...
    'ylim',YLim,'ytick',YLim,'yticklabel',[],...
    'fontname','arial','fontsize',8);
Pos = get(gca,'position');
set(gca,'position',[Pos(1),0.75,Pos(3),0.2]);
text(6,0.9*YLim(2), 'Daily Cycles','fontname','arial','fontsize',8,'fontweight','bold');

sig = abs(a-b) > CI1 & abs(a-b) < CI2 ;
text(find(sig),repmat(0.8*YLim(1),1,sum(sig)),'*');
sig =  abs(a-b) > CI2 ;
text(find(sig)-0.2,repmat(0.8*YLim(1),1,sum(sig)),'**');

axes(A(2))
c = wkDat / sum(wkDat);
% abs(a-c) > CIs
bar(c-a,'facecolor',col2,'facealpha',0.6);
set(gca,'box','off','xtick',[],'xlim',[0 11],...
    'ylim',YLim,'ytick',YLim,'yticklabel',[],...
    'fontname','arial','fontsize',8);
Pos = get(gca,'position');
set(gca,'position',[Pos(1),0.5,Pos(3),0.2]);
text(6,0.9*YLim(2), 'Weekly Cycles','fontname','arial','fontsize',8,'fontweight','bold');

sig = abs(a-c) > CI1 & abs(a-c) < CI2 ;
text(find(sig),repmat(0.8*YLim(1),1,sum(sig)),'*');
sig =  abs(a-c) > CI2 ;
text(find(sig)-0.2,repmat(0.8*YLim(1),1,sum(sig)),'**');

axes(A(3))
d = longDat / sum(longDat);
% abs(a-d) > CIs
bar(d-a,'facecolor',col3,'facealpha',0.6);
set(gca,'box','off','ylim',YLim,'ytick',YLim,'yticklabel',YLab,...
    'xlim',[0 11],'xtick',1:length(XLab),'xticklabel',XLab,...
    'fontname','arial','fontsize',8);
Pos = get(gca,'position');
set(gca,'position',[Pos(1),0.25,Pos(3),0.2]);
xtickangle(60);
text(6,0.9*YLim(2), 'Long Cycles','fontname','arial','fontsize',8,'fontweight','bold');

sig = abs(a-d) > CI1 & abs(a-d) < CI2 ;
text(find(sig),repmat(0.8*YLim(1),1,sum(sig)),'*');
sig =  abs(a-d) > CI2 ;
text(find(sig)-0.2,repmat(0.8*YLim(1),1,sum(sig)),'**');

ylabel('Difference (%)');
% xlabel('Seizure type');

set(gcf,'paperunits','centimeters','paperposition',[0 0 8 10]);
print(gcf,'-dpng','-r300','SyndromeProp');


%%
figure;
bar(allDat/sum(allDat),'facecolor','k','facealpha',0.6);
% histogram(typeLong,'normalization','pdf','facecolor',col3);
set(gca,'box','off','ylim',[0 0.5],'ytick',[0.2 0.4],'yticklabel',[20 40],...
    'xlim',[0 11],'xtick',1:length(allDat),'xticklabel',XLab,...
    'fontname','arial','fontsize',8);
xtickangle(60);
ylabel('Proportion of People (%)');

set(gcf,'paperunits','centimeters','paperposition',[0 0 8 8]);
print(gcf,'-dpng','-r300','AllSyndromeProp');
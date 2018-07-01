clc
clear
close all

col1 = [204 0 0]/255;
col2 = [255 128 0]/255;
col3 = [0 102 204]/255;

load('CycleIDs');
load('SeizureTypes');
load('StudyIDs');

remove = ismember(type,'type') | isundefined(type) ...
    | ismember(type,'Gelastic') | ismember(type,'Unknown') ...
    | ismember(type,'Other');
type(remove) = [];
id(remove) = [];

type = removecats(type);

typeAll = type(ismember(id,studyId));
typeDay = type(ismember(id,daily));
typeWk = type(ismember(id,weekly));
typeLong = type(ismember(id,long));

A = getCustomAxesPos(3,1,0,0.05);

YLim = [-0.1 0.2];
YLab = [-10 20];
XLab = categories(type);

XLab{7} = 'Infantile Spasms';
XLab{10} = 'Secondary Gen.';

% axes(A(1))
a = histcounts(typeAll,'normalization','probability');
CI1 = CI(a,length(typeAll),0.05/12/4);
CI2 = CI(a,length(typeAll),0.01/12/4);

axes(A(1))
b = histcounts(typeDay,'normalization','probability');
bar(b-a,'facecolor',col1,'facealpha',0.6);
% histogram(typeDay,'normalization','pdf','facecolor',col1);
set(gca,'box','off','xtick',[],'ylim',YLim,'ytick',YLim,'yticklabel',[],...
    'fontname','arial','fontsize',8);
Pos = get(gca,'position');
set(gca,'position',[Pos(1),0.75,Pos(3),0.2]);
text(6,0.9*YLim(2), 'Daily Cycles','fontname','arial','fontsize',8,'fontweight','bold');

sig = abs(a-b) > CI1 & abs(a-b) < CI2 ;
text(find(sig),repmat(0.8*YLim(1),1,sum(sig)),'*');
sig =  abs(a-b) > CI2 ;
text(find(sig)-0.2,repmat(0.8*YLim(1),1,sum(sig)),'**');

axes(A(2))
c = histcounts(typeWk,'normalization','probability');
% abs(a-c) > CIs
% histogram(typeWk,'normalization','pdf','facecolor',col2);
bar(c-a,'facecolor',col2,'facealpha',0.6);
set(gca,'box','off','xtick',[],'ylim',YLim,'ytick',YLim,'yticklabel',[],...
    'fontname','arial','fontsize',8);
Pos = get(gca,'position');
set(gca,'position',[Pos(1),0.5,Pos(3),0.2]);
text(6,0.9*YLim(2), 'Weekly Cycles','fontname','arial','fontsize',8,'fontweight','bold');

sig = abs(a-c) > CI1 & abs(a-c) < CI2 ;
text(find(sig),repmat(0.8*YLim(1),1,sum(sig)),'*');
sig =  abs(a-c) > CI2 ;
text(find(sig)-0.2,repmat(0.8*YLim(1),1,sum(sig)),'**');

axes(A(3))
d = histcounts(typeDay,'normalization','probability');
% abs(a-d) > CIs
bar(d-a,'facecolor',col3,'facealpha',0.6);
% histogram(typeLong,'normalization','pdf','facecolor',col3);
set(gca,'box','off','ylim',YLim,'ytick',YLim,'yticklabel',YLab,...
    'xtick',1:length(XLab),'xticklabel',XLab,...
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
print(gcf,'-dpng','-r300','TypeProp');
%%

figure;
histogram(typeAll,'normalization','probability','facecolor','k');
set(gca,'box','off','fontname','arial','fontsize',8,...
    'ytick',[0 0.10 0.20],'yticklabel',[0 10 20]);
ylabel('Proportion of Seizures (%)');

set(gcf,'paperunits','centimeters','paperposition',[0 0 8 8]);
print(gcf,'-dpng','-r300','AllTypeProp');

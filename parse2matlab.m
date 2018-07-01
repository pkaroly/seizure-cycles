close all
clc
clear

load('pt_seizure_dates');
load('pt_gender');
mkdir('STMatlab')

% remove headers
Date_Time(1) = [];
Unlinked_ID(1) = [];
length_hr(1) = [];
length_min(1) = [];
length_sec(1) = [];

% remove zero durations
invalid = 3600*length_hr + 60*length_min + length_sec <= 0;
Date_Time(invalid) = [];
Unlinked_ID(invalid) = [];
length_hr(invalid) = [];
length_min(invalid) = [];
length_sec(invalid) = [];

% remove long durations
invalid = 3600*length_hr + 60*length_min + length_sec > 3600;
Date_Time(invalid) = [];
Unlinked_ID(invalid) = [];
length_hr(invalid) = [];
length_min(invalid) = [];
length_sec(invalid) = [];

patients = categories(Unlinked_ID);
Np = length(patients);

Nsz = zeros(1,Np);

count = 0;
for iPt = 1:Np
    
    % get the times for this patient
    curPt = patients{iPt};
    save_name = sprintf('STMatlab/Patient_%04d.mat',count);
    
    Nsz(iPt) = sum(Unlinked_ID == curPt);
    if sum(Unlinked_ID == curPt) < 50
        continue;
    end
    
    dates = datetime(Date_Time(Unlinked_ID == curPt));
    temp = datevec(dates);
    hours = length_hr(Unlinked_ID == curPt);
    mins = length_min(Unlinked_ID == curPt);
    secs = length_sec(Unlinked_ID == curPt);
    
    szHour = temp(:,4);
    szTime = datenum(dates);
    szDur = 3600*hours + 60*mins + secs;
    
    % get gender (1 = male, 0 = female, nan = NA)
    ptGender = gender(strcmp(id,curPt));
    
    % save times for this patient
    save(save_name,'szTime','szHour','szDur','curPt', 'ptGender');
      
    count = count + 1;
end
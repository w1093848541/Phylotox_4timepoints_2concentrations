Path = 'C:\Users\10938\Desktop\Jacobson Lab\2021.10-\4 time points';
File = dir(fullfile(Path,'*tsv'));
FileNames = {File.name};
PeakList = readtable('NAME.txt');
Mass_List = table2array(PeakList(1:end,1));
n = length(Mass_List);
FileNumber = 1;
Filename = FileNames(FileNumber);
Fname = cell2mat(Filename);
fid = fopen(Fname);
GetColumn = readtable(Fname,'FileType','text');
NumC = size(GetColumn,2);
result = zeros(n+1,NumC-7);
DataSet = textscan(fid,repmat('%s',1,NumC));
fclose(fid);
stop = length(DataSet{1,1,1});
t = 1;
for t = 1:n
    result(t+1,1)= Mass_List (t);
end
s = 1;
% 
% for s = 1:n
%     mass = Mass_List(s);
%     i = 7;
%     for i = 7:stop
%         data_mass1 = DataSet{1,1}{i,1};
%         data_mass = str2num(data_mass1);
%          diff = mass - data_mass;
%         ppm = abs(diff/mass);
%         if ppm <= 0.000002
%             m = 2;
%             for m = 2:NumC-7
%                 result(s+1,m) = str2num(DataSet{1,m+7}{i,1});
%             end
%             end
%         end
% end 
result1 = num2cell(result);
s = 2;
for s = 2:NumC-7
    result1{1,s} = DataSet{1,s}{4,1};
end
result2 = result1;
result3 = readmatrix('PQN_predictedresult.csv');
QC=0;
E4=0;
E3=0;
E2=0;
E1=0;
C4=0;
C3=0;
C2=0;
C1=0;
C0=0;
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'QC'
        QC = QC +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E4'
        E4 = E4 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E4+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E3'
        E3 = E3 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E4+E3+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E2'
        E2 = E2 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E4+E3+E2+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E1'
        E1 = E1 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E4+E3+E2+E1+1} = result1 {m,s};
        end
    end
end     
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'C4'
        C4 = C4 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E4+E3+E2+E1+C4+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'C3'
        C3 = C3 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E4+E3+E2+E1+C4+C3+1} = result1 {m,s};
        end
    end
end  
  s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'C2'
        C2 = C2 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E4+E3+E2+E1+C4+C3+C2+1} = result1 {m,s};
        end
    end
end    
  s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'C1'
        C1 = C1 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E4+E3+E2+E1+C4+C3+C2+C1+1} = result1 {m,s};
        end
    end
end 
  s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'C0'
        C0 = C0 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E4+E3+E2+E1+C4+C3+C2+C1+C0+1} = result1 {m,s};
        end
    end
end 



Col = size (result3,2);
Row = size (result3,1);
result3 = num2cell(result3);
ColResult2 = 2;
for ColResult2 = 2:Col
    RowResult2 = 2;
    for RowResult2 = 2: Row+1
        result2 {RowResult2,ColResult2} = result3{RowResult2-1,ColResult2};
    end
end


resultforANOVA = cell2mat(result3.');

t = 1;

n = length(result3);
s = size(result3,2);

AnovaTPpass = 1;
AnovaECpass = 1;
AnovaIApass = 1;


for t = 1:n
    QCt = zeros(QC,1);
    E4t = zeros(E4,1);
    E3t = zeros(E3,1);
    E2t = zeros(E2,1);
    E1t = zeros(E1,1);
    C4t = zeros(C4,1);
    C3t = zeros(C3,1);
    C2t = zeros(C2,1);
    C1t = zeros(C1,1);
    C0t = zeros(C0,1);
    filt = 1;
    for filt = 1:QC
        QCt(filt,1) = resultforANOVA(filt,t);
    end
    filt = 1+QC;
    for filt = 1+QC:QC+E4
        E4t(filt,1) = resultforANOVA(filt,t);
    end
    filt = 1 + QC+E4;
    for filt = 1+QC+E4:QC+E4+E3
        E3t(filt,1) = resultforANOVA(filt,t);
    end
    filt = 1 +QC+ E4 + E3;
    for filt = 1+QC+E4+E3:QC+E4+E3+E2
        E2t(filt,1) = resultforANOVA(filt,t);
    end 
    filt = 1 + QC + E4 + E3 + E2;
    for filt = 1+E4+E3+E2+QC:QC+E4+E3+E2+E1
        E1t(filt,1) = resultforANOVA(filt,t);
    end
    filt = 1 + QC + E4 + E3 + E2 + E1;
    for filt = 1+E4+E3+E2+E1+QC:QC+C4+E4+E3+E2+E1
        C4t(filt,1) = resultforANOVA(filt,t);
    end
    filt =1+QC+E4+E3+E2+E1+C4;
    for filt = 1+QC+E4+E3+E2+E1+C4:QC+E4+E3+E2+E1+C4+C3
        C3t(filt,1) = resultforANOVA(filt,t);
    end
    filt = 1+QC+E4+E3+E2+E1+C4+C3;
    for filt = 1+QC+E4+E3+E2+E1+C4+C3: 1+QC+E4+E3+E2+E1+C4+C3+C2
        C2t(filt,1) = resultforANOVA(filt,t);
    end
    filt = 1+QC+E4+E3+E2+E1+C4+C3+C2;
    for filt = 1+QC+E4+E3+E2+E1+C4+C3+C2: 1+QC+E4+E3+E2+E1+C4+C3+C2+C1
        C1t(filt,1) = resultforANOVA(filt,t);
    end
    filt = 1+QC+E4+E3+E2+E1+C4+C3+C2+C1;
    for filt = 1+QC+E4+E3+E2+E1+C4+C3+C2+C1: 1+QC+E4+E3+E2+E1+C4+C3+C2+C1+C0
        C0t(filt,1) = resultforANOVA(filt,t);
    end
    
    
    
    E4t(E4t==0) = [];
    E3t(E3t==0) = [];
    E2t(E2t==0) = [];
    E1t(E1t==0) = [];
    C4t(C4t==0) = [];
    C3t(C3t==0) = [];
    C2t(C2t==0) = [];
    C1t(C1t==0) = [];
    
    A_Data = [E4t;E3t;E2t;E1t;C4t;C3t;C2t;C1t];
    ValidE4 = length(E4t);
    ValidE3 = length(E3t);
    ValidE2 = length(E2t);
    ValidE1 = length(E1t);
    ValidC4 = length(C4t);
    ValidC3 = length(C3t);
    ValidC2 = length(C2t);
    ValidC1 = length(C1t);
    ValidE = ValidE4 +ValidE3 + ValidE2 + ValidE1;
    ValidC = ValidC4 +ValidC3 + ValidC2 + ValidC1;
    
    AnovaTP = zeros(ValidE+ValidC,1);
    AnovaEC = zeros(ValidE+ValidC,1);
    
    for TP = 1:ValidE4
        AnovaTP(TP,1) = 4;
    end
    for TP = ValidE4+1:ValidE4+ValidE3
        AnovaTP(TP,1) = 3;
    end
    for TP = ValidE4+ValidE3+1:ValidE4+ValidE3+ValidE2
        AnovaTP(TP,1) = 2;
    end
    for TP = ValidE4+ValidE3+ValidE2+1:ValidE4+ValidE3+ValidE2+ValidE1
        AnovaTP(TP,1)= 1;
    end
    for TP = ValidE4+ValidE3+ValidE2+ValidE1+1:ValidE4+ValidE3+ValidE2+ValidE1+ValidC4
        AnovaTP(TP,1) = 4;
    end
    for TP = ValidE4+ValidE3+ValidE2+ValidE1+ValidC4+1:ValidE4+ValidE3+ValidE2+ValidE1+ValidC4+ValidC3
        AnovaTP(TP,1) = 3;
    end
    for TP = ValidE4+ValidE3+ValidE2+ValidE1+ValidC4+ValidC3+1:ValidE4+ValidE3+ValidE2+ValidE1+ValidC4+ValidC3+ValidC2
        AnovaTP(TP,1) = 2;
    end
    for TP = ValidE4+ValidE3+ValidE2+ValidE1+ValidC4+ValidC3+ValidC2+1:ValidE4+ValidE3+ValidE2+ValidE1+ValidC4+ValidC3+ValidC2+ValidC1
        AnovaTP(TP,1) = 1;
    end
    
    
    for EC = 1:ValidE
        AnovaEC(EC,1) = 1;
    end
    for EC = 1+ValidE:ValidE+ValidC
        AnovaEC(EC,1) = 2;
    end
    
     AnovaP = anovan(A_Data,{AnovaTP AnovaEC},'model',2,'varnames',{'Time Point','Exposed vs. Control'},'display','off');
    
    if AnovaP(1,1) <= 0.05
        ts = 1;
        for ts = 1:s
            AnovaTPresult(AnovaTPpass,ts) = cell2mat(result3(t,ts));
            AnovaTPresult(AnovaTPpass,s+1) = AnovaP(1,1);
        end
        AnovaTPpass = AnovaTPpass+1;
    end
    
    if AnovaP(2,1) <= 0.05
        ts = 1;
        for ts = 1:s
            AnovaECresult(AnovaECpass,ts) = cell2mat(result3(t,ts));
            AnovaECresult(AnovaECpass,s+1) = AnovaP(2,1);
        end
        AnovaECpass = AnovaECpass+1;
    end
       
    if AnovaP(3,1) <= 0.05
        ts = 1;
        for ts = 1:s
            AnovaIAresult(AnovaIApass,ts) = cell2mat(result3(t,ts));
            AnovaIAresult(AnovaIApass,s+1) = AnovaP(3,1);
        end
        AnovaIApass = AnovaIApass+1;
    end
end

nWrite = ['ANOVA_TP.xlsx'];

xlswrite(nWrite,AnovaTPresult);

nWrite = ['ANOVA_EC.xlsx'];

xlswrite(nWrite,AnovaECresult);

nWrite = ['ANOVA_IA.xlsx'];

xlswrite(nWrite,AnovaIAresult);

Anova_all = [AnovaTPresult;AnovaECresult;AnovaIAresult];
p_position = size(Anova_all,2);
Anova_all(:,p_position) = [];

nWrite = ['ANOVA_all.xlsx'];
xlswrite(nWrite,Anova_all);

result4 = readmatrix('glog_predictedresult.csv');



Col = size (result4,2);
Row = size (result4,1);
result4 = num2cell(result4);



Means_SDs = zeros(n+1,21);
t = 1;
for t = 1:n
    Means_SDs(t,1) = result4{t,1};
end


for s = 1:n
    QCs = zeros(QC,1);
    m = 1;
    for m = 1:QC
        QCs(m,1) = result4{s,m+1};
    end
    QCs = QCs(~isnan(QCs));
    average = mean(QCs);
    Means_SDs(s,2) = average;
    S = std(QCs);
    Means_SDs(s,3) = S;
end

for s = 1:n
    E4s = zeros(E4,1);
    m = 1;
    for m = 1:E4
        E4s(m,1) = result4{s,m+QC+1};
    end
    E4s = E4s(~isnan(E4s));
    average = mean(E4s);
    Means_SDs(s,4) = average;
    S = std(E4s);
    Means_SDs(s,5) = S;
end


for s = 1:n
    E3s = zeros(E3,1);
    m = 1;
    for m = 1:E3
        E3s(m,1) = result4{s,m+QC+E4+1};
    end
    E3s = E3s(~isnan(E3s));
    average = mean(E3s);
    Means_SDs(s,6) = average;
    S = std(E3s);
    Means_SDs(s,7) = S;
end


for s = 1:n
    E2s = zeros(E2,1);
    m = 1;
    for m = 1:E2
        E2s(m,1) = result4{s,m+QC+E4+E3+1};
    end
    E2s = E2s(~isnan(E2s));
    average = mean(E2s);
    Means_SDs(s,8) = average;
    S = std(E2s);
    Means_SDs(s,9) = S;
end


for s =1:n
    E1s = zeros(E1,1);
    m = 1;
    for m = 1:E1
        E1s(m,1) = result4{s,m+QC+E4+E3+E2+1};
    end
    E1s = E1s(~isnan(E1s));
    average = mean(E1s);
    Means_SDs(s,10) = average;
    S = std(E1s);
    Means_SDs(s,11) = S;
end


for s = 1:n
    C4s = zeros(C4,1);
    m = 1;
    for m = 1:C4
        C4s(m,1) = result4{s,m+QC+E4+E3+E2+E1+1};
    end
    C4s = C4s(~isnan(C4s));
    average = mean(C4s);
    Means_SDs(s,12) = average;
    S = std(C4s);
    Means_SDs(s,13) = S;
end


for s = 1:n
    C3s = zeros(C3,1);
    m = 1;
    for m = 1:C3
        C3s(m,1) = result4{s,m+QC+E4+E3+E2+E1+C4+1};
    end
    C3s = C3s(~isnan(C3s));
    average = mean(C3s);
    Means_SDs(s,14) = average;
    S = std(C3s);
    Means_SDs(s,15) = S;
end


for s = 1:n
    C2s = zeros(C2,1);
    m = 1;
    for m = 1:C2
        C2s(m,1) = result4{s,m+QC+E4+E3+E2+E1+C4+C3+1};
    end
    C2s = C2s(~isnan(C2s));
    average = mean(C2s);
    Means_SDs(s,16) = average;
    S = std(C2s);
    Means_SDs(s,17) = S;
end


for s = 1:n
    C1s = zeros(C1,1);
    m = 1;
    for m = 1:C1
        C1s(m,1) = result4{s,m+QC+E4+E3+E2+E1+C4+C3+C2+1};
    end
    C1s = C1s(~isnan(C1s));
    average = mean(C1s);
    Means_SDs(s,18) = average;
    S = std(C1s);
    Means_SDs(s,19) = S;
end


for s = 1:n
    C0s = zeros(C0,1);
    m = 1;
    for m = 1:C0
        C0s(m,1) = result4{s,m+QC+E4+E3+E2+E1+C4+C3+C2+C1+1};
    end
    C0s = C0s(~isnan(C0s));
    average = mean(C0s);
    Means_SDs(s,20) = average;
    S = std(C0s);
    Means_SDs(s,21) = S;
end

s = 1;
for s = 1:n
    
    m = 2;
    for m = 2:21
        if isnan(Means_SDs (s,m))
            Means_SDs (s,m) = 0;
    end
    end
end


Means_SDs(all(Means_SDs==0,2),:) = [];
Means_SDs(:,all(Means_SDs==0,1))= [];

nWrite = ['means&SD_afterGLOG.xlsx'];

xlswrite(nWrite,Means_SDs);
       

 result5 = xlsread('glog_predictedresult.csv');


SigPeak_Info = xlsread('ANOVA_TP.xlsx','Sheet1');%Change the significant peaks document here
SigPeak_Number = size(SigPeak_Info,1);
Info = zeros(SigPeak_Number,21);
s = 1;

 for s = 1:SigPeak_Number
     mass = SigPeak_Info(s,1);
     m = 1;
     for m = 1:n
     if (Means_SDs(m,1)-mass) <= 0.000002
         Info(s,:) = Means_SDs(m,:); 

        end
     end 
 end

 for s = 1:SigPeak_Number
      mass = SigPeak_Info(s,1);
     m = 1;
     for m = 1:n
     if (result5(m,1)-mass) <= 0.000002
         SigpeaksinformationTP(s,:) = result5(m,:); 

        end
     end 
 end
 
 nWrite = ['sigpeaks_TP.xlsx'];

xlswrite(nWrite,SigpeaksinformationTP);

 
 
 
 
%Heat maps 
SampleType = size(Info,2);

NumSigPeak = size(Info,1);
yvalues = ["QCs","Exposed TP4","Exposed TP3","Exposed TP2","Exposed TP1","Control TP4","Control TP3","Control TP2","Control TP1","TP0"];
xvalues = cell(1,NumSigPeak);
t =1 ;
for t = 1:NumSigPeak
xvalues{1,t} = num2str(Info(t,1));
end
data = zeros(NumSigPeak,20);
t = 1;
for t = 1:NumSigPeak
m = 1;
for m = 1:20
data(t,m) = Info(t,m+1);
end
end
t= 1;
for t = 1:NumSigPeak
m = 1;
for m = 1:20
data(t,m) = Info(t,m+1);
end
end
Sample_Type = yvalues;
Metabolite_MZ = xvalues;
realdata = zeros(NumSigPeak,10);
t = 1;
for t = 1:10
m = 1;
for m = 1:NumSigPeak
realdata(m,t) = data(m,2*t-1);
end
end
realdata = realdata';

%below is the part for putting endpoints in order
Value_Order = realdata';
realdata = realdata';
t =1 ;
for t = 1:NumSigPeak
m = 1;
row_intensity = zeros (1,10);
for m = 1:10
row_intensity(1,m) = realdata(t,m);
end
[Xsorted,Xidx] = sort(row_intensity);
[Xsorted,Xidx2] = sort(Xidx);
Value_Order (t,:) = Xidx2(1,:);
end
Value_Order = Value_Order';
Zscore_final = zeros(NumSigPeak,10);

t =1 ;
for t = 1:NumSigPeak
m = 1;
Zscore_intensity = zeros(1,10);
for m = 1:10
Zscore_intensity(1,m)  = realdata(t,m);
end
Zscore_rows = zscore(Zscore_intensity);
Zscore_final(t,:) = Zscore_rows(1,:);
end
ZZZ = Zscore_final';

%clustergram
cgo = clustergram (ZZZ);
Metabolites = Info(:,1);
Metabolites = Metabolites';
set(cgo,'ColumnLabels',Metabolites);
ytest = cellstr(yvalues);
set(cgo,'RowLabels',ytest);
set(cgo,'Linkage','complete','Dendrogram',3)
set(cgo,'Colormap',autumn);


SigPeak_Info = xlsread('ANOVA_EC.xlsx','Sheet1');%Change the significant peaks document here
SigPeak_Number = size(SigPeak_Info,1);
Info = zeros(SigPeak_Number,21);
s = 1;

 for s = 1:SigPeak_Number
     mass = SigPeak_Info(s,1);
     m = 1;
     for m = 1:n
     if (Means_SDs(m,1)-mass) <= 0.000002
         Info(s,:) = Means_SDs(m,:); 

        end
     end 
 end
 
 for s = 1:SigPeak_Number
      mass = SigPeak_Info(s,1);
     m = 1;
     for m = 1:n
     if (result5(m,1)-mass) <= 0.000002
         SigpeaksinformationEC(s,:) = result5(m,:); 

        end
     end 
 end
 
 nWrite = ['sigpeaks_EC.xlsx'];

xlswrite(nWrite,SigpeaksinformationEC);

 
 
 
 
%Heat maps 
SampleType = size(Info,2);

NumSigPeak = size(Info,1);
yvalues = ["QCs","Exposed TP4","Exposed TP3","Exposed TP2","Exposed TP1","Control TP4","Control TP3","Control TP2","Control TP1","TP0"];
xvalues = cell(1,NumSigPeak);
t =1 ;
for t = 1:NumSigPeak
xvalues{1,t} = num2str(Info(t,1));
end
data = zeros(NumSigPeak,20);
t = 1;
for t = 1:NumSigPeak
m = 1;
for m = 1:20
data(t,m) = Info(t,m+1);
end
end
t= 1;
for t = 1:NumSigPeak
m = 1;
for m = 1:20
data(t,m) = Info(t,m+1);
end
end
Sample_Type = yvalues;
Metabolite_MZ = xvalues;
realdata = zeros(NumSigPeak,10);
t = 1;
for t = 1:10
m = 1;
for m = 1:NumSigPeak
realdata(m,t) = data(m,2*t-1);
end
end
realdata = realdata';

%below is the part for putting endpoints in order
Value_Order = realdata';
realdata = realdata';
t =1 ;
for t = 1:NumSigPeak
m = 1;
row_intensity = zeros (1,10);
for m = 1:10
row_intensity(1,m) = realdata(t,m);
end
[Xsorted,Xidx] = sort(row_intensity);
[Xsorted,Xidx2] = sort(Xidx);
Value_Order (t,:) = Xidx2(1,:);
end
Value_Order = Value_Order';
Zscore_final = zeros(NumSigPeak,10);

t =1 ;
for t = 1:NumSigPeak
m = 1;
Zscore_intensity = zeros(1,10);
for m = 1:10
Zscore_intensity(1,m)  = realdata(t,m);
end
Zscore_rows = zscore(Zscore_intensity);
Zscore_final(t,:) = Zscore_rows(1,:);
end
ZZZ = Zscore_final';

%clustergram
cgo = clustergram (ZZZ);
Metabolites = Info(:,1);
Metabolites = Metabolites';
set(cgo,'ColumnLabels',Metabolites);
ytest = cellstr(yvalues);
set(cgo,'RowLabels',ytest);
set(cgo,'Linkage','complete','Dendrogram',3)
set(cgo,'Colormap',autumn);


SigPeak_Info = xlsread('ANOVA_IA.xlsx','Sheet1');%Change the significant peaks document here
SigPeak_Number = size(SigPeak_Info,1);
Info = zeros(SigPeak_Number,21);
s = 1;

 for s = 1:SigPeak_Number
     mass = SigPeak_Info(s,1);
     m = 1;
     for m = 1:n
     if (Means_SDs(m,1)-mass) <= 0.000002
         Info(s,:) = Means_SDs(m,:); 

        end
     end 
 end
 
 for s = 1:SigPeak_Number
      mass = SigPeak_Info(s,1);
     m = 1;
     for m = 1:n
     if (result5(m,1)-mass) <= 0.000002
         SigpeaksinformationIA(s,:) = result5(m,:); 

        end
     end 
 end
 
 nWrite = ['sigpeaks_IA.xlsx'];

xlswrite(nWrite,SigpeaksinformationIA);

 
 
 
 
%Heat maps 
SampleType = size(Info,2);

NumSigPeak = size(Info,1);
yvalues = ["QCs","Exposed TP4","Exposed TP3","Exposed TP2","Exposed TP1","Control TP4","Control TP3","Control TP2","Control TP1","TP0"];
xvalues = cell(1,NumSigPeak);
t =1 ;
for t = 1:NumSigPeak
xvalues{1,t} = num2str(Info(t,1));
end
data = zeros(NumSigPeak,20);
t = 1;
for t = 1:NumSigPeak
m = 1;
for m = 1:20
data(t,m) = Info(t,m+1);
end
end
t= 1;
for t = 1:NumSigPeak
m = 1;
for m = 1:20
data(t,m) = Info(t,m+1);
end
end
Sample_Type = yvalues;
Metabolite_MZ = xvalues;
realdata = zeros(NumSigPeak,10);
t = 1;
for t = 1:10
m = 1;
for m = 1:NumSigPeak
realdata(m,t) = data(m,2*t-1);
end
end
realdata = realdata';

%below is the part for putting endpoints in order
Value_Order = realdata';
realdata = realdata';
t =1 ;
for t = 1:NumSigPeak
m = 1;
row_intensity = zeros (1,10);
for m = 1:10
row_intensity(1,m) = realdata(t,m);
end
[Xsorted,Xidx] = sort(row_intensity);
[Xsorted,Xidx2] = sort(Xidx);
Value_Order (t,:) = Xidx2(1,:);
end
Value_Order = Value_Order';
Zscore_final = zeros(NumSigPeak,10);

t =1 ;
for t = 1:NumSigPeak
m = 1;
Zscore_intensity = zeros(1,10);
for m = 1:10
Zscore_intensity(1,m)  = realdata(t,m);
end
Zscore_rows = zscore(Zscore_intensity);
Zscore_final(t,:) = Zscore_rows(1,:);
end
ZZZ = Zscore_final';

%clustergram
cgo = clustergram (ZZZ);
Metabolites = Info(:,1);
Metabolites = Metabolites';
set(cgo,'ColumnLabels',Metabolites);
ytest = cellstr(yvalues);
set(cgo,'RowLabels',ytest);
set(cgo,'Linkage','complete','Dendrogram',3)
set(cgo,'Colormap',autumn);
    
    
    
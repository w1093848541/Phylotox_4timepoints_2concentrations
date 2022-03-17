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

for s = 1:n
    mass = Mass_List(s);
    i = 7;
    for i = 7:stop
        data_mass1 = DataSet{1,1}{i,1};
        data_mass = str2num(data_mass1);
         diff = mass - data_mass;
        ppm = abs(diff/mass);
        if ppm <= 0.000002
            m = 2;
            for m = 2:NumC-7
                result(s+1,m) = str2num(DataSet{1,m+7}{i,1});
            end
            end
        end
end
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
% Fname(end-3:end) = [];
% nWrite = [Fname '_all.xlsx'];
% xlswrite(nWrite,result2);



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

% Below is the t-test 


resultfort = cell2mat(result3.');
t = 1;



n = length(result3);
s = size(result3,2);
ttestallT4 = 1;
ttestallT3 = 1;
ttestallT2 = 1;
ttestallT1 = 1;

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
        QCt(filt,1) = resultfort(filt,t);
    end
    filt = 1+QC;
    for filt = 1+QC:QC+E4
        E4t(filt,1) = resultfort(filt,t);
    end
    filt = 1 + QC+E4;
    for filt = 1+QC+E4:QC+E4+E3
        E3t(filt,1) = resultfort(filt,t);
    end
    filt = 1 +QC+ E4 + E3;
    for filt = 1+QC+E4+E3:QC+E4+E3+E2
        E2t(filt,1) = resultfort(filt,t);
    end 
    filt = 1 + QC + E4 + E3 + E2;
    for filt = 1+E4+E3+E2+QC:QC+E4+E3+E2+E1
        E1t(filt,1) = resultfort(filt,t);
    end
    filt = 1 + QC + E4 + E3 + E2 + E1;
    for filt = 1+E4+E3+E2+E1+QC:QC+C4+E4+E3+E2+E1
        C4t(filt,1) = resultfort(filt,t);
    end
    filt =1+QC+E4+E3+E2+E1+C4;
    for filt = 1+QC+E4+E3+E2+E1+C4:QC+E4+E3+E2+E1+C4+C3
        C3t(filt,1) = resultfort(filt,t);
    end
    filt = 1+QC+E4+E3+E2+E1+C4+C3;
    for filt = 1+QC+E4+E3+E2+E1+C4+C3: 1+QC+E4+E3+E2+E1+C4+C3+C2
        C2t(filt,1) = resultfort(filt,t);
    end
    filt = 1+QC+E4+E3+E2+E1+C4+C3+C2;
    for filt = 1+QC+E4+E3+E2+E1+C4+C3+C2: 1+QC+E4+E3+E2+E1+C4+C3+C2+C1
        C1t(filt,1) = resultfort(filt,t);
    end
    filt = 1+QC+E4+E3+E2+E1+C4+C3+C2+C1;
    for filt = 1+QC+E4+E3+E2+E1+C4+C3+C2+C1: 1+QC+E4+E3+E2+E1+C4+C3+C2+C1+C0
        C0t(filt,1) = resultfort(filt,t);
    end
    
    
    
    E4t(E4t==0) = [];
    E3t(E3t==0) = [];
    E2t(E2t==0) = [];
    E1t(E1t==0) = [];
    C4t(C4t==0) = [];
    C3t(C3t==0) = [];
    C2t(C2t==0) = [];
    C1t(C1t==0) = [];
    
    [h,p] = ttest2(E4t,C4t,'Vartype','unequal');
    meanE4t = mean(E4t(:));
    meanC4t = mean(C4t(:));
    log2FC = log2(meanE4t/meanC4t);
       ts = 1;
       for ts = 1:s
       ttest2resultT4(ttestallT4,ts) = cell2mat(result3(t,ts));
       ttest2resultT4(ttestallT4,s+1) = p;
       ttest2resultT4(ttestallT4,s+2) = log2FC;
       end
       ttestallT4 = ttestallT4 +1;
      
    
    
    
    [h,p] = ttest2(E3t,C3t,'Vartype','unequal');
    meanE3t = mean(E3t(:));
    meanC3t = mean(C3t(:));
    log2FC = log2(meanE3t/meanC3t);
        ts = 1;
        for ts = 1:s
       ttest2resultT3(ttestallT3,ts) = cell2mat(result3(t,ts));
       ttest2resultT3(ttestallT3,s+1) = p;
       ttest2resultT3(ttestallT3,s+2) = log2FC;
        end
       ttestallT3 = ttestallT3 +1;
        
    
   
    
    [h,p] = ttest2(E2t,C2t,'Vartype','unequal');
    meanE2t = mean(E2t(:));
    meanC2t = mean(C2t(:));
    log2FC = log2(meanE2t/meanC2t);
        ts = 1;
        for ts = 1:s
       ttest2resultT2(ttestallT2,ts) = cell2mat(result3(t,ts));
       ttest2resultT2(ttestallT2,s+1) = p;
        ttest2resultT2(ttestallT2,s+2) = log2FC;
        end
       ttestallT2 = ttestallT2 +1;
       
    
    
    [h,p] = ttest2(E1t,C1t,'Vartype','unequal');
    meanE1t = mean(E1t(:));
    meanC1t = mean(C1t(:));
    log2FC = log2(meanE1t/meanC1t);
        ts = 1;
        for ts = 1:s
       ttest2resultT1(ttestallT1,ts) = cell2mat(result3(t,ts));
       ttest2resultT1(ttestallT1,s+1) = p;
        ttest2resultT1(ttestallT1,s+2) = log2FC;
        end
       ttestallT1 = ttestallT1 +1;
       
    
       
    
    
 
       
    
end    

%Output the t-test p-values and log2FC%

nWrite = ['ttest2resultT4_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT4);


nWrite = ['ttest2resultT3_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT3);


nWrite = ['ttest2resultT2_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT2);

nWrite = ['ttest2resultT1_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT1);


%Opearte the FDR test and output the results with q-values, saved for
%volcano plots


Pvaluesall  =  xlsread('ttest2resultT4_ALL.xlsx','Sheet1');
Pvalues = Pvaluesall(:,Col+1);
log2FC = Pvaluesall(:,Col+2);
ids = Pvaluesall(:,1);
fdr = mafdr(Pvalues);
[fdr,q,priori,R2] = mafdr(Pvalues,'Method','polynomial');

qvalue_result4 = [ids,fdr,q,log2FC];

priori_R24 = [priori,R2];

nWrite = ['qvalues_allT4.xlsx'];

xlswrite(nWrite,qvalue_result4);

nWrite = ['priori_R2T4.xlsx'];

xlswrite(nWrite,priori_R24);


Pvaluesall  =  xlsread('ttest2resultT3_ALL.xlsx','Sheet1');
Pvalues = Pvaluesall(:,Col+1);
log2FC = Pvaluesall(:,Col+2);
ids = Pvaluesall(:,1);
fdr = mafdr(Pvalues);
[fdr,q,priori,R2] = mafdr(Pvalues,'Method','polynomial');

qvalue_result3 = [ids,fdr,q,log2FC];

priori_R23 = [priori,R2];

nWrite = ['qvalues_allT3.xlsx'];

xlswrite(nWrite,qvalue_result3);
nWrite = ['priori_R2T3.xlsx'];

xlswrite(nWrite,priori_R23);

Pvaluesall  =  xlsread('ttest2resultT2_ALL.xlsx','Sheet1');
Pvalues = Pvaluesall(:,Col+1);
log2FC = Pvaluesall(:,Col+2);
ids = Pvaluesall(:,1);
fdr = mafdr(Pvalues);
[fdr,q,priori,R2] = mafdr(Pvalues,'Method','polynomial');

qvalue_result2 = [ids,fdr,q,log2FC];

priori_R22 = [priori,R2];

nWrite = ['qvalues_allT2.xlsx'];

xlswrite(nWrite,qvalue_result2);
nWrite = ['priori_R2T2.xlsx'];

xlswrite(nWrite,priori_R22);

Pvaluesall  =  xlsread('ttest2resultT1_ALL.xlsx','Sheet1');
Pvalues = Pvaluesall(:,Col+1);
log2FC = Pvaluesall(:,Col+2);
ids = Pvaluesall(:,1);
fdr = mafdr(Pvalues);
[fdr,q,priori,R2] = mafdr(Pvalues,'Method','polynomial');

qvalue_result1 = [ids,fdr,q,log2FC];

priori_R21 = [priori,R2];

nWrite = ['qvalues_allT1.xlsx'];

xlswrite(nWrite,qvalue_result1);


nWrite = ['priori_R2T1.xlsx'];

xlswrite(nWrite,priori_R21);


%Collect the metabolites that passed the FDR test

FDRpass = 1;

for npass = 1:n
    if qvalue_result1(npass,3) <= 0.05 || qvalue_result2(npass,3) <= 0.05 ||qvalue_result3(npass,3) <= 0.05 ||qvalue_result4(npass,3) <= 0.05
       FDRpass_all(FDRpass,:) = ttest2resultT1(FDRpass,:);
      
       FDRpass = FDRpass + 1;
    end
end
p_position = size(FDRpass_all,2);
p_position = p_position -1;
FDRpass_all(:,p_position) = [];
FDRpass_all(:,p_position) = [];


% Below is the output for each time point


% FDRpass1 = 1;
% 
% for npass = 1:n
%     if qvalue_result1(npass,3) <= 0.05 
%        FDRpass_1(FDRpass1,:) = ttest2resultT1(FDRpass1,:);
%       
%        FDRpass1 = FDRpass1 + 1;
%     end
% end
% 
% FDRpass2 = 1;
% 
% for npass = 1:n
%     if qvalue_result2(npass,3) <= 0.05
%        FDRpass_2(FDRpass2,:) = ttest2resultT2(FDRpass2,:);
%       
%        FDRpass2 = FDRpass2 + 1;
%     end
% end
% 
% 
% FDRpass3 = 1;
% 
% for npass = 1:n
%     if qvalue_result3(npass,3) <= 0.05
%        FDRpass_3(FDRpass3,:) = ttest2resultT3(FDRpass3,:);
%       
%        FDRpass3 = FDRpass3 + 1;
%     end
% end
% 
% FDRpass4 = 1;
% 
% for npass = 1:n
%     if qvalue_result4(npass,3) <= 0.05
%        FDRpass_4(FDRpass4,:) = ttest2resultT4(FDRpass4,:);
%       
%        FDRpass4 = FDRpass4 + 1;
%     end
% end
% FDRpass_1(:,p_position) = [];
% FDRpass_1(:,p_position) = [];
% 
% FDRpass_2(:,p_position) = [];
% FDRpass_2(:,p_position) = [];
% 
% FDRpass_3(:,p_position) = [];
% FDRpass_3(:,p_position) = [];
% 
% FDRpass_4(:,p_position) = [];
% FDRpass_4(:,p_position) = [];











% nWrite = ['qpass_T1.xlsx'];
% 
% xlswrite(nWrite,FDRpass_1);
% 
% nWrite = ['qpass_T2.xlsx'];
% 
% xlswrite(nWrite,FDRpass_2);
% 
% nWrite = ['qpass_T3.xlsx'];
% 
% xlswrite(nWrite,FDRpass_3);
% 
% nWrite = ['qpass_T4.xlsx'];
% 
% xlswrite(nWrite,FDRpass_4);

nWrite = ['qpass_ALL.xlsx'];

xlswrite(nWrite,FDRpass_all);



%Now, calculate the means and SD of gloged data


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




%Combine the significant peaks information with gloged data

SigPeak_Info = xlsread('qpass_ALL.xlsx','Sheet1');%Change the significant peaks document here
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
 result5 = xlsread('glog_predictedresult.csv');
 for s = 1:SigPeak_Number
      mass = SigPeak_Info(s,1);
     m = 1;
     for m = 1:n
     if (result5(m,1)-mass) <= 0.000002
         Sigpeaksinformation(s,:) = result5(m,:); 

        end
     end 
 end
 
 nWrite = ['sigpeaks_afterGLOG.xlsx'];

xlswrite(nWrite,Sigpeaksinformation);

 
 
 
 
%Heat maps 
SampleType = size(Info,2);

NumSigPeak = size(Info,1);
yvalues = ["Quality Controls","Exposed TP4","Exposed TP3","Exposed TP2","Exposed TP1","Control TP4","Control TP3","Control TP2","Control TP1","Time Point 0"];
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

metabolitesIDs = Sigpeaksinformation(:,1);
metabolitesIDs = metabolitesIDs';
ZZZs = [metabolitesIDs;ZZZ];

nWrite = ['Zscored_sigdata.xlsx'];

xlswrite(nWrite,ZZZs);


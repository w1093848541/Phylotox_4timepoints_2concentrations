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

%Eliminating metabolites that have more than 50% of missing values in QC
%samples
s = 2;
NumToDelete = 1;
ToDelete= [];
for s = 2:n
    testQ = 2;
    for testQ = 1:QC+1
        Test_QC(1,testQ) = result2{s,testQ};
    end
    numberOfZeros = sum(Test_QC == 0);
    if numberOfZeros > QC/2
       ToDelete (NumToDelete,1) = s;
       NumToDelete = NumToDelete + 1 ;
    end
end
ToDelete = flip(ToDelete);

NumToDelete = 1;
for NumToDelete = 1:length(ToDelete)
    Deleting = ToDelete(NumToDelete,1);
    result2(Deleting,:) = [];
end

Fname(end-3:end) = [];
nWrite = [Fname '_all.xlsx'];
xlswrite(nWrite,result2);
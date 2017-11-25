%% Use this code first to build a histogram containing donor only molecules.
%  In this step, set the LEAKAGE, donor_blank and acceptor_blank values as 0. Observe the position of the donor only peak.

function batchFRET()

close all;
fclose('all');

%read data
pth=input('Directory [default=C:\\User\\tir data\\yyyy\\New Folder]  ');
if isempty(pth)
   	pth='C:\User\tir data\yyyy\New Folder';
end
cd(pth);
fList=dir;
nf=size(fList,1);

fretE = [];
LEAKAGE=0.13;
donor_blank=0;
acceptor_blank=0;
for n = 3:nf
    s=fList(n).name;
    if fList(n).isdir || ~strcmp(s(end-5:end), 'traces')
        continue;
    end
    fid=fopen(s,'r');

    %first line of binary file specifies length of trace
    len=fread(fid,1,'int32');
    disp('The len of the time traces is: ')
    disp(len);

    %number of traces
    Ntraces=fread(fid,1,'int16');
    disp('The number of traces is: ')
    disp(Ntraces/2);

    %raw is a linear array, looks like so
    raw=fread(fid,Ntraces*len,'int16');
    disp('Done reading data.');
    fclose(fid);

    %convert into traces
    index=(1:Ntraces*len);
    Data=zeros(Ntraces,len);
    Data(index)=raw(index);
    
    donor=zeros(Ntraces/2,len);
    acceptor=zeros(Ntraces/2,len);
    total=zeros(Ntraces/2,1);
    
    startFrame=12;
    for i=1:(Ntraces/2)
        donor(i,:)=Data(i*2-1,:);
        acceptor(i,:)=Data(i*2,:);
        tempD=sum(donor(i,(startFrame:startFrame+7)),2);
        tempA=sum(acceptor(i,(startFrame:startFrame+7)),2);
        total(i)=(tempD+tempA)/8.-donor_blank-acceptor_blank;
    end

    hdl = figure;
    hist(total,0:50:2000);
    grid on;
    zoom on;

    fcutoff1=input('low cutoff intensity: ','s');
    cutoff1=str2num(fcutoff1);
    if isempty(cutoff1)
        cutoff1=300;
    end
    fcutoff2=input('high cutoff intensity: ','s');
    cutoff2=str2num(fcutoff2);
    if isempty(cutoff2)
        cutoff2=1000;
    end

    for i = 1:(Ntraces/2)
        if total(i) < cutoff1 || total(i) > cutoff2
            continue;
        end
        tempD=sum(donor(i,(startFrame:startFrame+7)),2)-8*donor_blank;
        tempA=sum(acceptor(i,(startFrame:startFrame+7)),2)-8*acceptor_blank;
        fretE=[fretE;(tempA-LEAKAGE*tempD)/(tempD+tempA-LEAKAGE*tempD)];
    end
    close(hdl);
end

hist(fretE, -0.0875:0.025:1.0125);
[counts,centers]=hist(fretE, -0.0875:0.025:1.0125);
centers=centers';
counts=counts';

save('FRET histogram.dat','counts','-ascii');

fclose('all');


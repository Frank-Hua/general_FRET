%Program allows to skim through traces
%Use 's' to save individual traces as .dat files.
%Use 'b' to go to previous, 'g' to specific

function s_tr()
close all;
fclose('all');

%read data
pth=input('directory [default=C:\\User\\tir data\\yyyy\\New Folder]  ');
if isempty(pth)
    disp('error');
end
cd(pth);
save_file=pth;

fname=input('index # of filename [default=1]  ');
if isempty(fname)
    fname=1;
end
fname=num2str(fname);
disp(['hel' fname '.traces']);
%disp(['film' fname '.traces']);

timeunit=input('time unit [default=0.03 sec]  ');
if isempty(timeunit)
    timeunit=0.03;
end

%select the folder to which the files need to be saved
newfolder = [fname ' selected traces'];
mkdir(newfolder);

fid=fopen(['hel' fname '.traces'],'r');
%fid=fopen(['film' fname '.traces'],'r');

%first line of binary file specifies length of trace
len=fread(fid,1,'int32');
disp('The len of the time traces is: ')
disp(len);

%number of traces
Ntraces=fread(fid,1,'int16');
disp('The number of traces is: ')
disp(Ntraces/2);

%raw is a linear array, looks like it
raw=fread(fid,Ntraces*len,'int16');
disp('Done reading data.');
fclose(fid);


%convert into traces
index=(1:Ntraces*len);
Data=zeros(Ntraces,len);
donor=zeros(Ntraces/2,len);
acceptor=zeros(Ntraces/2,len);
total=zeros(Ntraces/2,1);
Data(index)=raw(index);
for i=1:(Ntraces/2),
    donor(i,:)=Data(i*2-1,:);
    acceptor(i,:)=Data(i*2,:);
    tempD=sum(donor(i,(11:20)),2);
    tempA=sum(acceptor(i,(11:20)),2);
    total(i)=(tempA+tempD)/10.; 
end

time=(0:(len-1))*timeunit;

figure;
hist(total,80);
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
    cutoff2=2500;
end
    
hdl=figure;
i=0;
n=0;
N_mol=Ntraces/2;

while ((Ntraces/2)-i) > 0 ,
    i = i+1 ;
    
    if total(i) < cutoff1 || total(i) > cutoff2
        n = n+1;
        continue;
    end
    
    %trace window
    figure(hdl);
    subplot(2,1,1);
    plot(time,donor(i,:),'g', time,acceptor(i,:),'r', time,acceptor(i,:)+donor(i,:)+200,'k');
    title(['  Molecule ' num2str(i) ' of ' num2str(N_mol)]);
    axis tight;
    temp=axis;
    temp(3)=0;
    %temp(4)=temp(4)*1.1;
    if temp(4) < 800
        temp(4)=800;
    end
    axis(temp);
    grid on;
    zoom on;
    
    subplot(2,1,2);
    fretE=acceptor(i,:)./(acceptor(i,:)+donor(i,:));
    for m=1:len,
        if acceptor(i,m)+donor(i,m)==0
            fretE(m)=-0.5;
        end
    end % This is to avoid undefined fretE interfering with future analysis
    plot(time,fretE,'b');
    axis tight;
    temp=axis;
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp);
    grid on;
    zoom on;

    answer=input('press b-back,g-go,s-save  ','s');
    if answer=='b'
        i=i-2;
    end

    if answer=='g'
        mol= input('which molecule do you choose:  ');
        i= mol-1;
    end

    %to save individual traces
    if answer=='s'
        fname1=[save_file '\' newfolder  '\hel' fname ' tr' num2str(i) '.dat'];
        %fname1=[save_file '\' newfolder  '\film' fname ' tr' num2str(i) '.dat'];
        output=[time' donor(i,:)' acceptor(i,:)'];
        save(fname1,'output','-ascii') ;
        disp(n);
    end

end

disp('number of molecules selected:');
disp(N_mol-n);

close all;
fclose('all');







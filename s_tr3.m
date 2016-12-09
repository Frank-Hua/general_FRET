%Program allows to skim through traces
%'Enter' to click and measure dwelltimes.
%Use 'b' to go to previous, 'g' to specific

function s_tr3()
close all;
fclose('all');

%read data
pth=input('Directory [default=C:\\User\\tir data\\yyyy\\New Folder]  ');
if isempty(pth)
   	pth='C:\User\tir data\yyyy\New Folder';
end
cd(pth);
A=dir;
[nf,dum]=size(A);

donor=[];
acceptor=[];
Ntraces=0;
clicktime=[];

for i=1:nf,
    if A(i).isdir == 0
        s=A(i).name;
        if strcmp(s(end-2:end), 'dat')
            disp(s);
            Data=dlmread(s);
            donor=[donor,Data(:,2)];
            acceptor=[acceptor,Data(:,3)];
            Ntraces=Ntraces+1;
        end
    end
end

donor=donor';
acceptor=acceptor';

len=size(Data,1);
time=Data(:,1);

donor_blank=input('Donor channel baseline [default=0]  ');
acceptor_blank=input('Acceptor channel baseline [default=0]  ');
if isempty(donor_blank)
    donor_blank=0;
end
if isempty(acceptor_blank)
    acceptor_blank=0;
end

leakage=0.12;

donor=donor-donor_blank;
acceptor=acceptor-acceptor_blank;

i=0;
hdl=figure;

while (Ntraces-i) > 0 ,
    i = i+1 ;
    %trace window
    figure(hdl);
    subplot(2,10,[1 9]);
    plot(time,donor(i,:),'g', time,acceptor(i,:)-leakage*donor(i,:),'r', time,donor(i,:)+acceptor(i,:)-leakage*donor(i,:)+200,'k');
    title(['  Molecule ' num2str(i) ' of ' num2str(Ntraces)]);
    axis tight;
    temp=axis;
    temp(3)=-temp(4)*0.05;
    temp(4)=temp(4)*1.1;
    axis(temp);
    grid on;
    zoom on;
    
    subplot(2,10,[11 19]);
    fretE=(acceptor(i,:)-leakage*donor(i,:))./(donor(i,:)+acceptor(i,:)-leakage*donor(i,:));
    for m=1:len,
        if acceptor(i,m)+donor(i,m)==0
            fretE(m)=-0.5;
        end
    end % This is to avoid undefined fretE interfering with future analysis
    %fretE=smooth(fretE,5,'moving');
    plot(time,fretE,'b');
    axis tight;
    temp=axis;
    temp(3)=-0.1;
    temp(4)=1.1;
    axis(temp);
    grid on;
    zoom on;

    answer=input('press b-back,g-go,enter-dwelltime  ','s');
    if isempty(answer)
        answer='f';
    end
    disp(answer);
    
    if answer=='b'
        i=i-2;
    end

    if answer=='g'
        mol= input('which molecule do you choose:  ');
        i= mol-1;
    end
    
    if answer=='f'
        [X,~]=ginput;
        result=X;
        templength=length(result);
        if templength >= 2 && mod(templength,2) == 0
            temp=zeros(templength/2,1);
            temp=result(2:2:end)-result(1:2:end);
            temp=temp';
            clicktime=[clicktime temp];
        end
    end
end

save('dwelltimes.dat','clicktime','-ascii');

input('enter-to continue ','s');
close all;
fclose('all');







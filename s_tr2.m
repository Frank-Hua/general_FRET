%Program allows to skim through traces
%Use 'c' to correct traces, 'Enter' to analyze traces.
%Use 'b' to go to previous, 'g' to specific

function s_tr2()
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
timeunit=Data(2,1);
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
Cor=[];
countsT=[];

donor=donor-donor_blank;
acceptor=acceptor-acceptor_blank;

newfolder = ['HaMMy traces'];
mkdir(newfolder);
cd([pth '\' newfolder]);

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

    answer=input('press b-back,g-go,c-correct,enter-FRET histogram  ','s');
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
    
    if answer=='c'
        [~,Y]=ginput(4);
        correction=zeros(3,1);
        correction(1)=Y(1);
        correction(2)=Y(2);
        correction(3)=(Y(4)-Y(2))/(Y(3)-Y(1));
        Cor=[Cor,correction];
    end
    
    if answer=='f'
        [p,~]=ginput(2);
        p=round(p/timeunit);
        X(1)=p(1);
        X(2)=p(2);
        output=[time(X(1):X(2)) donor(i,(X(1):X(2)))' acceptor(i,(X(1):X(2)))'];
        save(['tr' num2str(i) '.dat'],'output','-ascii');
        
        xbins=(-0.0875:0.025:1.0125);
        [counts,centers] = hist(fretE(X(1):X(2)),xbins);
        counts=counts/(X(2)-X(1)+1);
        subplot(2,10,20);
        plot(counts,centers,'bo');
        axis tight;
        temp=axis;
        temp(3) = -0.1;
        temp(4) = 1.1;
        axis(temp);
        %axis off;
        countsT=[countsT,counts'];
        
        Ave_donor=mean(donor(i,X(1):X(2)));
        Ave_acceptor=mean(acceptor(i,X(1):X(2)));
        delta_donor=donor(i,X(1):X(2))-Ave_donor;
        delta_acceptor=acceptor(i,X(1):X(2))-Ave_acceptor;
        Ave_donor_cross_acceptor=zeros(1,X(2)-X(1)+1);
        for tau=0:(X(2)-X(1))
            normalize_number=0;
            sum_tau_1=0;
            for t=1:(X(2)-X(1)+1-tau)
                sum_tau_1=sum_tau_1+(delta_donor(t)*delta_acceptor(t+tau));
                normalize_number=normalize_number+1;
            end
            Ave_donor_cross_acceptor(tau+1)=sum_tau_1/normalize_number/(Ave_donor+Ave_acceptor);
        end
        subplot(2,10,10)
        t=0:(X(2)-X(1));
        plot(t,Ave_donor_cross_acceptor,'-b');
        title('Cross Correlation');
        axis tight;
        temp=axis;
        temp(1)=-1;
        temp(2)=50;
        temp(4)=2;
        axis(temp);
        grid on;

        input('enter-to continue ','s');
        clf(hdl,'reset');
    end
    
end

cd(pth);
newfolder2 = ['analysis'];
mkdir(newfolder2);
cd([pth '\' newfolder2]);

elevel=zeros(Ntraces,1);
total=elevel;

for i=1:Ntraces
    tempD=sum(donor(i,(21:30))-donor_blank,2);
    tempA=sum(acceptor(i,(21:30))-acceptor_blank,2);
    total(i)=(tempA+tempD)/10; %frank
    elevel(i)=(tempA-leakage*tempD)/(tempD+tempA-leakage*tempD);
end

%{
hdl2=figure;
figure(hdl2);
subplot(2,1,1);
hist(total,10);
zoom on;
%}
save('intensityResult.dat','total','-ascii');
%{
title(s);
subplot(2,1,2);
xbins=(0.05:0.1:0.95);
hist(elevel,xbins);
zoom on;
%}
save('FRETResult.dat','elevel','-ascii');
save('FRETResult_2.dat','countsT','-ascii');

disp(Cor);
cd(pth);

input('enter-to continue ','s');
close all;
fclose('all');







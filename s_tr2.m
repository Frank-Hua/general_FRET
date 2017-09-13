%Use this code to further analyze the selected traces by s_tr.m
%Use 'Enter' to analyze traces. Blank subtraction and donor leakage correction are applied.
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
[nf,~]=size(A);

donor=[];
acceptor=[];
Ntraces=0;

for i=1:nf
    if A(i).isdir == 0
        s=A(i).name;
        if strcmp(s(end-3:end), '.dat')
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

LEAKAGE=input('input donor leakage correction value: ');
if isempty(LEAKAGE)
    LEAKAGE=0;
end

for i=1:Ntraces
    tD(i)=sum(donor(i,(end-18:end-11)),2)/8.;
    tA(i)=sum(acceptor(i,(end-18:end-11)),2)/8.;
end

figure;
plot(tD,tA,'x');
axis square;
temp=axis;
temp(1)=0;
temp(3)=0;
temp(2)=max(temp(2),temp(4));
temp(4)=temp(2);
axis(temp);
grid on;
zoom on;
[donor_blank,acceptor_blank]=ginput(1);
disp(donor_blank);
disp(acceptor_blank);

donor=donor-donor_blank;
acceptor=acceptor-acceptor_blank;
acceptor=acceptor-LEAKAGE*donor;

newfolder = ['HaMMy traces'];
mkdir(newfolder);
cd([pth '\' newfolder]);

hdl=figure;
i=0;
countsT=[];

while (Ntraces-i) > 0
    i = i+1;
    
    %trace window
    figure(hdl);
    subplot(2,10,[1 9]);
    plot(time,donor(i,:),'g', time,acceptor(i,:),'r', time,donor(i,:)+acceptor(i,:)+400,'k');
    title(['  Molecule ' num2str(i) ' of ' num2str(Ntraces)]);
    axis tight;
    temp=axis;
    temp(3)=-temp(4)*0.05;
    temp(4)=temp(4)*1.1;
    if temp(4) < 800
        temp(4)=800;
    end
    axis(temp);
    grid on;
    zoom on;
    
    subplot(2,10,[11 19]);
    fretE = medfilt1(acceptor(i,:),3)./(medfilt1(donor(i,:),3)+medfilt1(acceptor(i,:),3));
    for m=1:len
        if acceptor(i,m)+donor(i,m)<=0
            fretE(m)=-0.2;
        end
    end % This is to avoid undefined fretE interfering with future analysis
    fretE(fretE>1.1)=1.1;
    fretE(fretE<-0.2)=-0.2;
    plot(time,fretE,'b');
    axis tight;
    temp=axis;
    temp(3)=-0.2;
    temp(4)=1.1;
    axis(temp);
    grid on;
    zoom on;

    answer=input('press b-back,g-go,enter-FRET histogram  ','s');
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
        [X,~]=ginput(2);
        X=round(X/timeunit);
        output=[time(X(1):X(2)) donor(i,(X(1):X(2)))' acceptor(i,(X(1):X(2)))'];
        save(['tr' num2str(i) '.dat'],'output','-ascii');
        
        xbins=(-0.0875:0.025:1.0125);
%         xbins=(-0.025:0.05:1.025);
        [counts,centers] = hist(fretE(X(1):X(2)),xbins);
        counts=counts/(X(2)-X(1)+1);
        subplot(2,10,20);
        plot(counts,centers,'bo');
        title('FRET Histogram');
        temp=axis;
        temp(3) = -0.1;
        temp(4) = 1.1;
        axis(temp);
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
    tempD=sum(donor(i,(12:19)),2);
    tempA=sum(acceptor(i,(12:19)),2);
    total(i)=(tempD+tempA)/8; %frank
    elevel(i)=(tempA)/(tempD+tempA);
end

save('intensityResult.dat','total','-ascii');
save('FRETResult.dat','elevel','-ascii');
save('FRETResult_tr.dat','countsT','-ascii');

cd(pth);

close all;
fclose('all');







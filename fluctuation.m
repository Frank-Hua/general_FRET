function fluctuation()

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
leakage=0.116;

hdl=figure;
fluc=zeros(1,nf);
for i=1:nf,
    if A(i).isdir == 0
        s=A(i).name;
        if strcmp(s(end-2:end), 'dat')
            disp(s);
            Data=dlmread(s);
            time=Data(:,1);
            donor=Data(:,2);
            acceptor=Data(:,3);
            len=size(Data,1);
            
            FRET=(acceptor-leakage*donor)./(donor+acceptor-leakage*donor);
            for m=1:len,
                if acceptor(m)+donor(m)-leakage*donor(m)==0
                    FRET(m)=-0.5;
                end
            end % This is to avoid undefined fretE interfering with future analysis

            figure(hdl);
            subplot(2,10,[1 9]);
            plot(time,donor,'g', time,acceptor-leakage*donor,'r', time,donor+acceptor-leakage*donor+200,'k');
            subplot(2,10,[11 19]);
            plot_results(time,FRET,'b');
            
            %input('enter to next trace');
            fluc(i)=std(FRET);
        end
    end
end

close all;
fclose('all');

figure;
plot(fluc(3:end));

disp(mean(fluc(3:end)));
disp(std(fluc(3:end)));
disp(std(fluc(3:end))/sqrt(nf-2));

end

function plot_results(a,b,c)

%fretE=smooth(fretE,5,'moving');
plot(a,b,c);
axis tight;
temp=axis;
temp(3)=-0.1;
temp(4)=1.1;
axis(temp);
grid on;
zoom on;

end

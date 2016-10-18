function batchb()

close all;
fclose('all');

%read data
pth=input('Directory [default=C:\\User\\tir data\\yyyy\\New Folder]  ');
if isempty(pth)
   	pth='C:\User\tir data\yyyy\New Folder';
end
cd(pth);
disp(pth);
A=dir;
[nf,dum]=size(A);

leakage=input('Donor leakage correction [default=0.116]  ');
if isempty(leakage)
    leakage=0.116;
end

countsT=[];

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
            
            xbins=(-0.0875:0.025:1.0125);
            [counts,centers] = hist(FRET,xbins);
            counts=counts/len;
            countsT=[countsT,counts'];
        end
    end
end

save('FRETResult_2.dat','countsT','-ascii');

end
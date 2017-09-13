%% stack HaMMy traces for ebFRET

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
Ntraces=[];

for i=3:nf
    if A(i).isdir == 0
        s=A(i).name;
        if strcmp(s(end-2:end), 'dat')
            disp(s);
            Data=dlmread(s);
            donor=[donor',Data(:,2)']';
            acceptor=[acceptor',Data(:,3)']';
            if isempty(max(Ntraces))
                Ntraces=[Ntraces',ones(1,size(Data,1))]';
            else
                Ntraces=[Ntraces',(max(Ntraces)+1)*ones(1,size(Data,1))]';
            end
        end
    end
end

output=[Ntraces donor acceptor];
save('ebfret.dat','output','-ascii');


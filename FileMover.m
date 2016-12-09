% random pick files

NFolders=2;
NFiles=50;
path='D:\Hua\00projects\p14-Rep-X project with Sarah\08maximum percentage unwound\2016-04-28\20 mM MgCl2';
cd(path);

A=dir;
[nf,dum]=size(A);

if (nf-2)>=NFiles*NFolders
    for i=1:NFolders
        mkdir(['group' num2str(i)]);
        j=0;
        while j<NFiles
            n=3+floor((nf-2)*rand);
            if exist(A(n).name,'file')
                movefile(A(n).name,[path '\group' num2str(i)]);
                j=j+1;
            end
        end
    end
end
      
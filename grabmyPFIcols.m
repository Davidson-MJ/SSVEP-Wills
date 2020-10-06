%grabmyPFIcols

reds = (cbrewer('seq', 'Greys', 10));
drkred = reds(3,:);
ltred = reds(10,:);


blues =(cbrewer('seq', 'Blues', 10));
drkbl= blues(3,:);
ltbl = blues(10,:);

mycols=[];

%15 hz, pfi, 
mycols(1,1).c = drkbl;
%15 hz, pmd, 
mycols(1,2).c = ltbl;

%15 hz, pfi, 
mycols(2,1).c = drkred;
%15 hz, pmd, 
mycols(2,2).c = ltred;
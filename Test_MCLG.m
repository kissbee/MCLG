function  Test_MCLG(dataname, X, Y )
% Input£º     X£ºdata matrix
%             Y£ºtruth label
% Output:     result£ºclustering result

k = 6;

beta=100;
mu=0.001;
alpha=0.001;
lambda=0.001;
gamma=0.005;
theta=1;%fixed


row1=string(zeros(1,13));
row1(1)='beta';
row1(2)='mu';
row1(3)='alpha';
row1(4)='lambda';
row1(5)='gamma';
row1(6)='theta';
row1(7)='ACC';
row1(8)='NMI';
row1(9)='Purity';
row1(10)='Fscore';
row1(11)='Precision';
row1(12)='Recall';
row1(13)='AR';
row1(14)='Entropy';

resultfile = strcat('result/',dataname,'.csv'); % record result(ACC NMI Purity...)
objfile = strcat('result/',dataname,'Obj.csv'); % record objvalue


f1 = fopen(resultfile, 'w+', 'n', 'utf8');      % create table title
fprintf(f1,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',row1(1,1),row1(1,2),row1(1,3),row1(1,4),row1(1,5),row1(1,6),row1(1,7),row1(1,8),row1(1,9),row1(1,10),row1(1,11),row1(1,12),row1(1,13),row1(1,14));
fclose(f1);

f2 = fopen(objfile, 'w+', 'n', 'utf8');         % create table title
fprintf(f2, '%s,%s,%s,%s,%s,%s,%s\n', row1(1,1),row1(1,2),row1(1,3),row1(1,4),row1(1,5),row1(1,6),'obj');
fclose(f2);

f3 = fopen(resultfile, 'a', 'n', 'utf8');       % add running results to file

%
V = size(X,2);      % view number
c = max(Y);         % class number
n = length(Y);      % object number

X_bar = [];
for iv = 1 : V
    X_bar = [X_bar,X{iv}];
end

for iv = 1:V
    fea = X{iv};
    fea = mapminmax(fea,0,1);  % input n * dv,output n * dv
    X{iv} = fea';  % dv * n
end

X_bar = double(mapminmax(X_bar,0,1));%n*d

% Gv
for iv = 1:V
    S = constructS_PNG(X{iv}, k, 1); %% set k,  dv * n
    S1= (S+S')/2;
    
    Gv{iv}=S1; %
    
    D1 = diag(sum(S1));
    Lv{iv} = D1-S1;
    [Fviv, ~, ~]=eig1(Lv{iv}, c, 0);
    Fv{iv} = Fviv;%
end
%disp(size(X_bar))
S_bar = constructS_PNG(X_bar', k, 1); %% set k,  dv * n
S_bar = (S_bar+S_bar')/2;


tic;
Z = MCLG(dataname,X,Y,Gv,S_bar,Fv,beta,mu,alpha,lambda,gamma,theta);%
%Z = Z';
% Repeat kmeans 20 times.
res=newMetrics_Kmeans(Z,Y,c);

fprintf('\tbeta:%.4f,mu:%.4f,alpha:%.4f,lambda:%.4f,gamma:%.4f,theta:%.4f\n', beta,mu,alpha,lambda,gamma,theta);
disp(res);

% write mean data to file
row = [beta,mu,alpha,lambda,gamma,theta,res(1,1),res(1,2),res(1,3),res(1,4),res(1,5),res(1,6),res(1,7),res(1,8)];
fprintf(f3, '%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%f,%f,%f,%f,%f,%f,%f,%f\n', row(1),row(2),row(3),row(4),row(5),row(6),row(7),row(8),row(9),row(10),row(11),row(12),row(13),row(14));

% write std data to file
fprintf(f3, '%s,%s,%s,%s,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f\n',' ',' ',' ',' ',' ',' ',res(2,1),res(2,2),res(2,3),res(2,4),res(2,5),res(2,6),res(2,7),res(2,8));

toc;

fclose(f3);

end


function [res]= newMetrics_Kmeans(F,Y,c)
% F : label matrix
% Y : true label
% c : class number
stream = RandStream.getGlobalStream;
reset(stream);

MAXiter = 100; % Maximum number of iterations for KMeans.
REPlic = 10;   % Number of replications for KMeans.
for idx = 1:20
    pre_Y = kmeans(F,c,'maxiter',MAXiter,'replicates',REPlic,'emptyaction','singleton');
    Res(idx,:) = Clustering8Measure(Y, pre_Y);
end
res = [mean(Res);std(Res)];
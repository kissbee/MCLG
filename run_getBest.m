clc
%1:100Leaves.mat 2:20newsgroups.mat 3:3Sources.mat 4:BBC.mat 5:BBCsport.mat 
%6:BBCsport_3view.mat 7:COIL20.mat 8:Caltech101-20.mat 9:Caltech101-7.mat 
%10:CiteSeer.mat 11:Cora.mat 12:HW.mat 13:HW2sources.mat 14:HandWritten_5view.mat 
%15:Hdigit.mat 16:LandUse_21.mat 17:MNIST.mat 18:MSRC.mat 19:MSRC_v1.mat 
%20:Mfeat.mat 21:Movie.mat 22:NGs.mat 23:NUS.mat 24:ORL.mat 25:ORL_mtv.mat 
%26:Prokaryotic.mat 27:Reuters_1200.mat 28:Scene.mat 29:Scene_15.mat 
%30:Synthetic.mat 31:TC_DS.mat 32:TM_DS.mat 33:ThreeRing.mat 34:ToyData.mat 
%35:TwoMoon.mat 36:Uci_digit.mat 37:Umist_mtv.mat 38:Webkb.mat 39:Wiki_textimage.mat 
%40:WikipediaArticles.mat 41:Yale.mat 42:YaleB.mat 43:Youtube.mat 

% Caltech101-7 ����p=7 8����
data=csvread('result/Prokaryotic.csv',1, 0);
p = 7; %��������Ŀ��һ

% ��������+1
pc = 7;


% lambdaset9 = 10.^(-3:1:3);
% tauset9 = (100:50:500);
% disp(lambdaset9);

%% ��ѯ��һ��ACC/NMI���һ
m = max(data(:,p));                 % ��ȡACC�����е����ֵ
%disp(m);
n=find(data(:,p)==m);               % �������ֵ������
disp(n)                            % ���ֵ��ͬ��������
%disp(data(n,:));
if length(n)>1                      % ������ͬ���ACC��ֵ��ȡ��Сstd��Ӧ�ļ���
    std_n = n+1;                    % �õ���ͬACC����Ӧ��std������
    %disp(std_n);
    tempmin=100;
    for j=1:length(std_n)
        if data(j,p)<=tempmin
            tempmin=data(j,p);
            finalN=std_n(j)-1;
        end
    end
end
%disp(finalN);
if length(n)==1
    finalN = n;
end

%disp(finalN)
fprintf('%d\n', finalN);
for i=1:length(data(finalN,:))
    if(i<pc)
        fprintf('%.6f ', data(finalN,i));
    end
    if(i==pc)
        fprintf('\n');
    end
    if(i>=pc)
        fprintf('%.2f',roundn(data(finalN,i),-4)*100);
        %fprintf('%.4f', data(finalN,i));
        fprintf('(%.2f) ', roundn(data(finalN+1,i)*100,-4));
        %fprintf('(%.4f) ', data(finalN+1,i));
    end
end
fprintf('\n '); % ����

%% ��ѯ�ڶ���ACC
D = sort(data(:,p),'descend');      % ��ȡACC�����У�������

for i=1:length(data(:,1))              % ��õڶ���Ԫ��ֵ
    secendD = D(i);
    if secendD~=m
        break;
    end
end

nn=find(data(:,p)==secendD);    % ���ҵڶ���ֵ������
if length(nn)>1                 % ������ͬ�ڶ���ACC��ֵ��ȡ��һ��n����
    finalN = nn(1);
end
if length(nn)==1
    finalN = nn;
end

fprintf('\n%d\n', finalN);
for i=1:length(data(finalN,:))
    if(i<pc)
        fprintf('%.6f ', data(finalN,i));
    end
    if(i==pc)
        fprintf('\n');
    end
    if(i>=pc)
        fprintf('%.2f', data(finalN,i)*100);
        fprintf('(%.2f) ', data(finalN+1,i)*100);
    end
end
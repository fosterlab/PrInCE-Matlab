%% SVM: protein interactions, generate P-R curve
% To generate all variables, run ROC_PCPSILAC.m and pause after section 5
%
% Training set = inverse_self & Int_matrix
% TP labels = TP_Matrix

% Get data
I = inverse_self & Int_matrix;
labels = TP_Matrix(I(:));
y = labels(:);
y(y>0) = 1;
y(y~=1) = -1;
data1 = Dist.R2(I(:)); % 1 - R^2
data2 = Dist.Euc(I(:));
data3 = Dist.Center(I(:));
data4 = Dist.dtw(I(:));
X = [data1 data2 data3 data4];
%X = data1;
Nd = size(X,2);

% Soft whiten data
eps = 2e-16;
wmu = zeros(Nd,1);
wstd = zeros(Nd,1);
for ii = 1:Nd
  wmu(ii) = mean(X(:,ii));
  wstd(ii) = std(X(:,ii));
  X(:,ii) = (X(:,ii) - wmu(ii)) / 2 / (wstd(ii) + eps);
end

% Make training and testing data
nn = 1000;
if 1 % balance training data
  I1 = find(y==1);
  I1 = I1(randsample(length(I1),nn/2));
  I0 = find(y==-1);
  I0 = I0(randsample(length(I0),nn/2));
  Itrain = ismember(1:length(y),[I1;I0]);
  Ipred = ~Itrain;
else
  Itrain = zeros(size(y));
  Itrain(randsample(length(Itrain),nn)) = 1;
  Itrain = Itrain==1;
  Ipred = ~Itrain;
end
Xtr = X(Itrain,:);
ytr = y(Itrain);
Xnew = X(Ipred,:);
ynew = y(Ipred);

% % Make custom Cost matrix for SVM
% fpw = 10^1; % false positive penalization weight.
% fnw = 10^1; % false negative penalization weight.
% Cost = [0 fpw; fnw 0];

% Feature selection
% f2consider = find(IndFeat(Xtr,ytr) > 2);
% [fsin,history] = sequentialfs(@SVM_class_fun,Xtr(:,f2consider),ytr,'cv',5);
% Xtr2 = Xtr(:,f2consider(fsin));
% Xnew2 = Xnew(:,f2consider(fsin));
% f2consider(fsin)

% Fit svm model
svm = fitcsvm(Xtr,ytr);
[~,score_svm] = predict(svm,Xnew);
[precision_svm,recall_svm,~,auc_svm] = perfcurve(ynew,score_svm(:,2),1,'xcrit','prec');

% Fit naive bayes model
nab = fitcnb(Xtr,ytr);
[~,score_nb] = predict(nab,Xnew);
[precision_nb,recall_nb,~,auc_nb] = perfcurve(ynew,score_nb(:,2),1,'xcrit','prec');

disp([num2str(auc_svm) ', ' num2str(auc_nb)])

%figure
hold on
plot(recall_svm,precision_svm,'r')
plot(recall_nb,precision_nb,'g')
xlabel('Recall','fontsize',13)
ylabel('Precision','fontsize',13)
hl = legend('SVM','Naive Bayes');
set(hl,'fontsize',13)
title('Comparing classifiers SVM and NB','fontsize',13)



%% Prove that the classifiers are working
% Include fake informative variables

clear precnb{ii} recnb{ii} aucnb(ii) precsvm{ii} recsvm{ii} aucsvm(ii)

% Get data
I = inverse_self & Int_matrix;
labels = TP_Matrix(I(:));
y = labels(:);
y(y>0) = 1;
y(y~=1) = -1;
data1 = Dist.R2(I(:)); % 1 - R^2
data2 = y + rand(size(y))*10; % weak fake data
data3 = y + rand(size(y))*5; % medium fake data
data4 = y + rand(size(y))*1; % strong fake data
data5 = rand(size(y))*6; % awful fake data
clear data
% data{1} = data1;
% data{2} = [data1 data2];
% data{3} = [data1 data3];
% data{4} = [data1 data4];
% data{5} = [data1            data5];
% data{6} = [data1 data2      data5];
% data{7} = [data1 data3      data5];
% data{8} = [data1 data4      data5];
data{1} = data1;
data{2} = [data1 data5];
data{3} = [data1 data5*1000];
data{4} = [data1 data5*1000];
data{5} = [data1 data5*1000];
data{6} = [data1 data5*1000];

for ii=1:length(data)
  disp(num2str(ii))
  X = data{ii};
  Nd = size(X,2);
  
  if ii == 4
    % Soft whiten data
    eps = 2e-16;
    wmu = zeros(Nd,1);
    wstd = zeros(Nd,1);
    for jj = 1:Nd
      wmu(jj) = mean(X(:,jj));
      wstd(jj) = std(X(:,jj));
      X(:,jj) = (X(:,jj) - wmu(jj)) / 2 / (wstd(jj) + eps);
    end
  end
  
  if ii == 5 || ii ==6
    % Hard whiten data
    eps = 2e-16;
    mm = zeros(Nd,1);
    mn = zeros(Nd,1);
    mx = zeros(Nd,1);
    for jj = 1:Nd
      mm(jj) = mean(X(:,jj));
      mn(jj) = min(X(:,jj));
      mx(jj) = max(X(:,jj));
      X(:,jj) = (X(:,jj) - mm(jj)) / (mx(jj) - mn(jj));
    end
  end
  
  % Make training and testing data
  nn = 2000;
  if 1 % balance training data
    I1 = find(y==1);
    I1 = I1(randsample(length(I1),nn/2));
    I0 = find(y==-1);
    I0 = I0(randsample(length(I0),nn/2));
    Itrain = ismember(1:length(y),[I1;I0]);
    Ipred = ~Itrain;
  else
    Itrain = zeros(size(y));
    Itrain(randsample(length(Itrain),nn)) = 1;
    Itrain = Itrain==1;
    Ipred = ~Itrain;
  end
  Xtr = X(Itrain,:);
  ytr = y(Itrain);
  Xnew = X(Ipred,:);
  ynew = y(Ipred);
  
  % weight features by F-ratio
  if ii==6
    ww = IndFeat(Xtr,ytr);
    for jj = 1:Nd
      X(:,jj) = X(:,jj) * ww(jj);
    end
  end
  
  % Fit svm model
  svm = fitcsvm(Xtr,ytr);
  [label_svm,score_svm] = predict(svm,Xnew);
  [precsvm{ii},recsvm{ii},~,aucsvm(ii)] = perfcurve(ynew,score_svm(:,2),1,'xcrit','prec');
  
  % Fit naive bayes model
  nab = fitcnb(Xtr,ytr);
  [~,score_nb] = predict(nab,Xnew);
  [precnb{ii},recnb{ii},~,aucnb(ii)] = perfcurve(ynew,score_nb(:,2),1,'xcrit','prec');
  
end

figure
subplot(1,2,1),hold on
plot(recsvm{1},precsvm{1},'k')
plot(recsvm{2},precsvm{2},'m')
plot(recsvm{3},precsvm{3},'r')
plot(recsvm{4},precsvm{4},'color',[.5 .5 .5])
plot(recsvm{5},precsvm{5},'b')
plot(recsvm{6},precsvm{6},'color',[.5 .5 1])
legend('R2','R2+bad','R2+1000*bad','R2+1000*bad+sw','R2+1000*bad+hw','R2+1000*bad+hw','location','northeast')
title('SVM')
subplot(1,2,2),hold on
plot(recnb{1},precnb{1},'k')
plot(recnb{2},precnb{2},'m')
plot(recnb{3},precnb{3},'r')
plot(recnb{4},precnb{4},'color',[.5 .5 .5])
plot(recnb{5},precnb{5},'b')
plot(recnb{6},precnb{6},'color',[.5 .5 1])
title('Naive Bayes')
set(gcf,'units','normalized','position',[.1 .1 .8 .6])



%% Try out DTW matrix
% Make Dist.DTW using a subset of proteins.
% Goal: 1000 interactions, 1000 non-interactions, run time <300s

clear nint
for oi = 1:size(TP_Matrix,1)-201
pri = (1:180) + oi;
I = inverse_self(pri,pri) & Int_matrix(pri,pri);
tpm = TP_Matrix(pri,pri);
nint(oi) = sum(tpm(:));
end
[~,offs] = max(nint);
pri = (1:180) + offs;
I = inverse_self(pri,pri) & Int_matrix(pri,pri);
tpm = TP_Matrix(pri,pri);
Nint = sum(tpm(:));
runhat = 0.007 * numel(I);
disp(['Number of interactions: ' num2str(Nint) ', Estimated run time: ' num2str(runhat) 's'])
labels = tpm(I(:));
y = labels(:);
y(y>0) = 1;
y(y~=1) = -1;


Dist2.Euc = Dist.Euc(pri,pri);
Dist2.Center = Dist.Center(pri,pri);
Dist2.R2 = Dist.R2(pri,pri); % one minus R squared
Dist2.dtw = ones(length(pri),length(pri));
for ii = 1:length(pri)
  if mod(ii,10)==0;disp(num2str(ii));end
  i1 = pri(ii);
  for jj = 1:length(pri)
    j1 = pri(jj);
    Dist2.dtw(ii,jj) = dtw_old(Chromatograms(i1,:),Chromatograms(j1,:));
  end
end



data1 = Dist2.R2(I(:)); % 1 - R^2
data2 = Dist2.Euc(I(:));
data3 = Dist2.Center(I(:));
data4 = Dist2.dtw(I(:));
X = [data1 data2 data3 data4];
Nd = size(X,2);

if 1
  % Soft whiten data
  eps = 2e-16;
  wmu = zeros(Nd,1);
  wstd = zeros(Nd,1);
  for ii = 1:Nd
    wmu(ii) = mean(X(:,ii));
    wstd(ii) = std(X(:,ii));
    X(:,ii) = (X(:,ii) - wmu(ii)) / 2 / (wstd(ii) + eps);
  end
else
  % Hard whiten data
  eps = 2e-16;
  mm = zeros(Nd,1);
  mn = zeros(Nd,1);
  mx = zeros(Nd,1);
  for jj = 1:Nd
    mm(jj) = mean(X(:,jj));
    mn(jj) = min(X(:,jj));
    mx(jj) = max(X(:,jj));
    X(:,jj) = (X(:,jj) - mm(jj)) / (mx(jj) - mn(jj));
  end
end

% Make training and testing data
% balance training data
nn = 1000;
I1 = find(y==1);
I1 = I1(randsample(length(I1),nn/2));
I0 = find(y==-1);
I0 = I0(randsample(length(I0),nn/2));
Itrain = ismember(1:length(y),[I1;I0]);
Ipred = ~Itrain;
Xtr = X(Itrain,:);
ytr = y(Itrain);
Xnew = X(Ipred,:);
ynew = y(Ipred);
% 
% % weight features by F-ratio
% ww = IndFeat(X,y);
% for jj = 1:Nd
%   X(:,jj) = X(:,jj) * ww(jj);
% end

% Fit models on just [R2]
% Fit svm model
svm = fitcsvm(Xtr(:,1),ytr);
[~,score_svm] = predict(svm,Xnew(:,1));
[precision_svm1,recall_svm1,~,pav1] = perfcurve(ynew,score_svm(:,2),1,'xcrit','prec');
[fpr_svm1,tpr_svm1,~,av1] = perfcurve(ynew,score_svm(:,2),1);
% Fit naive bayes model
nab = fitcnb(Xtr(:,1),ytr);
[~,score_nb] = predict(nab,Xnew(:,1));
[precision_nb1,recall_nb1,~,pab1] = perfcurve(ynew,score_nb(:,2),1,'xcrit','prec');
[fpr_nb1,tpr_nb1,~,ab1] = perfcurve(ynew,score_nb(:,2),1);

% Fit models on [R2 Euc C DTW]
% Fit svm model
svm = fitcsvm(Xtr,ytr);
[~,score_svm] = predict(svm,Xnew);
[precision_svm2,recall_svm2,~,pav2] = perfcurve(ynew,score_svm(:,2),1,'xcrit','prec');
[fpr_svm2,tpr_svm2,~,av2] = perfcurve(ynew,score_svm(:,2),1);
% Fit naive bayes model
nab = fitcnb(Xtr,ytr);
[~,score_nb] = predict(nab,Xnew);
[precision_nb2,recall_nb2,~,pab2] = perfcurve(ynew,score_nb(:,2),1,'xcrit','prec');
[fpr_nb2,tpr_nb2,~,ab2] = perfcurve(ynew,score_nb(:,2),1);


disp(['ROC AUC: ' num2str(av2-av1) ', ' num2str(ab2-ab1)])
disp(['PR AUC: ' num2str(pav2-pav1) ', ' num2str(pab2-pab1)])


figure
hold on
plot(recall_svm1,precision_svm1,'r')
plot(recall_nb1,precision_nb1,'b')
plot(recall_svm2,precision_svm2,'m','linewidth',2)
plot(recall_nb2,precision_nb2,'c','linewidth',2)
xlabel('Recall','fontsize',13)
ylabel('Precision','fontsize',13)
hl = legend('SVM, R2','Naive Bayes, R2','SVM, all','Naive Bayes, all');
set(hl,'fontsize',13)
title('Comparing classifiers SVM and NB','fontsize',13)
set(gcf,'units','normalized','position',[.05 .2 .4 .5])

figure
hold on
plot(fpr_svm1,tpr_svm1,'r')
plot(fpr_nb1,tpr_nb1,'b')
plot(fpr_svm2,tpr_svm2,'m','linewidth',2)
plot(fpr_nb2,tpr_nb2,'c','linewidth',2)
xlabel('FPR','fontsize',13)
ylabel('TPR','fontsize',13)
hl = legend('SVM, R2','Naive Bayes, R2','SVM, all','Naive Bayes, all','location','southeast');
set(hl,'fontsize',13)
title('Comparing classifiers SVM and NB','fontsize',13)
set(gcf,'units','normalized','position',[.55 .2 .4 .5])

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
X = [data1 data2 data3];
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


%% scratchpad for updating ROC_PCPSILAC.m

x = score_nb(:,2);
dR1 = log10(linspace(min(10.^x),max(10.^x)*.9,100));
dR2 = log10(linspace(max(10.^x)*.9,max(10.^x)*.99,300));
dR3 = log10(linspace(max(10.^x)*.99,max(10.^x),800));
dR = [dR1 dR2 dR3];

rec = nan(size(dR));
prec = nan(size(dR));
for di = 1:length(dR)
  D = dR(di);
  TP = sum(x>D & ynew==1);
  FP = sum(x>D & ynew==-1);
  TN = sum(x<D & ynew==-1);
  FN = sum(x<D & ynew==1);
  
  rec(di) = TP/(TP+FN);
  prec(di) = TP/(TP+FP);
end

figure
plot(rec,prec)

%%

x = score_svm(:,2);
%x = 1-Xnew(:,1);
dR = zeros(1,2000);
dx1 = (max(x) - min(x))/length(dR);
mindx = dx1/100;
maxdx = dx1*2;
rec = zeros(1,length(dR));
prec = zeros(1,length(dR));
di = 1;
prec0 = 0;
dR0 = 0;
while prec0<0.95 && dR0<max(x)
  
  % build next D
  if di>2
    m = diff(prec([di-2 di-1])) / diff(dR([di-2 di-1]));
    dx = max(mindx, 1/length(dR)/m);
    dx = min(maxdx, dx);
  end
  
  if di==1
    dR(di) = min(x);
  elseif di==2
    dR(di) = dR(di) + dx1;
  else
    dR(di) = dR(di-1) + dx;
  end
  
  
  %D = dR(di);
  TP = sum(x>dR(di) & ynew==1);
  FP = sum(x>dR(di) & ynew==-1);
  TN = sum(x<dR(di) & ynew==-1);
  FN = sum(x<dR(di) & ynew==1);
  
  rec(di) = TP/(TP+FN);
  prec(di) = TP/(TP+FP);
  prec0 = prec(di);
  dR0 = dR(di);
  di = di+1
end
prec = prec(1:di-1);
rec = rec(1:di-1);
dR = dR(1:di-1);

% sort dR ascending
[~,I] = sort(dR,'ascend');
dR = dR(I);
prec = prec(I);
rec = rec(I);

figure
plot(rec,prec)



%% How do I average svm models?

I = inverse_self & Int_matrix;
labels = TP_Matrix(I(:));
y = labels(:);
y(y>0) = 1;
y(y~=1) = -1;
data1 = Dist.R2(I(:)); % 1 - R^2
data2 = Dist.Euc(I(:));
data3 = Dist.Center(I(:));
data4 = Dist.dtw(I(:));
X = [data1 data2 data3];
%X = data1;
Nd = size(X,2);

maxIter = 10;
scoreMatrix = nan(size(X,1),maxIter);

for iter = 1:maxIter
  % % Make training and testing data
  nn = 10000; % length of training data
  % balance training data
  I1 = find(y==1);
  I1 = I1(randsample(length(I1),nn/2));
  I0 = find(y==-1);
  I0 = I0(randsample(length(I0),nn/2));
  Itrain = ismember(1:length(y),[I1;I0]);
  Ipred = ~Itrain;
  Xtr = X(Itrain,:);
  ytr = y(Itrain);
  Xnew = X(Ipred,:);
  
  % Fit svm model
  svm = fitcsvm(Xtr,ytr,'Standardize',true);
  [~,score] = predict(svm,Xnew);
  
  % Combine scores
  scoreMatrix(Ipred,iter) = score(:,2);
end


%% Compare ensemble averaging and voting

disp('Calculating scoreMatrix...')
[scoreMatrix, feats] = scorenb(Dist,possibleInts,TP_Matrix);
scoreMatrix = scoreMatrix(possibleInts(:),:);

desiredPrecision = [0.4 0.6 0.8 0.95];
class = TP_Matrix(possibleInts(:));
nn = 25;
Tol = 0.01; % get within 1% of precision
maxIter = 20; % zoom in 20 times at most

% single model, negative control
score_neg = scoreMatrix(:,2);

% vote, most extreme wins
%score_vote = scoreMatrix;

% average 1
score_avg1 = mean(scoreMatrix(:,1:ceil(size(scoreMatrix,2)/2)),2);
% average 2
score_avg2 = mean(scoreMatrix,2);

disp('Calculating precision/recall...')
scores = {scoreMatrix score_avg1 score_avg2 scoreMatrix(:,1)};
precSave = cell(size(scores));
recSave = cell(size(scores));
scoreSave = cell(size(scores));
for si = 1:length(scores)
  si
  score = scores{si};
  ds = linspace(min(score(:)),max(score(:)),nn); % start off with coarse search

  xcutoff = nan(size(desiredPrecision));
  calcprec = zeros(size(xcutoff));
  calcrec = zeros(size(xcutoff));
  scoreRange = nan(nn*maxIter*length(desiredPrecision),1);  % Save these for
  precRange = nan(size(scoreRange));                        % plotting
  recRange = nan(size(scoreRange));                         % later.
  for di = 1:length(desiredPrecision)
    calcTol = 10^10;
    iter = 0;
    deltaPrec = 1;
    prec0 = zeros(nn,1);
    % Stop zooming in when one of three things happens:
    % i) you get close enough to the desired precision (within calcTol)
    % ii) you've been through maxIter iterations
    % iii) zooming in stops being useful (precision changes by less than deltaPrec b/w iterations)
    while calcTol>Tol && iter<maxIter && deltaPrec>1e-3
      iter=iter+1;
      rec = nan(nn,1);
      prec = nan(nn,1);
      for dd =1:length(ds)
        if si ==1 % vote, i.e. have 2D score
          pos = sum(score>ds(dd),2) > sum(score<ds(dd),2);
          TP = sum(pos & class==1);
          FP = sum(pos & class==0);
          FN = sum(~pos & class==1);
        else % average, i.e. 1D score
          TP = sum(score>ds(dd) & class==1);
          FP = sum(score>ds(dd) & class==0);
          FN = sum(score<ds(dd) & class==1);
        end
        prec(dd) = TP/(TP+FP);
        rec(dd) = TP/(TP+FN);
      end
      deltaPrec = nanmean(abs(prec - prec0));
      
      % Save vectors for plotting
      i1 = find(isnan(scoreRange),1,'first');
      I = i1 : i1+nn-1;
      scoreRange(I) = ds;
      precRange(I) = prec;
      recRange(I) = rec;
      
      % Calculate how close to desiredPrecision(di) you got
      [calcTol,I] = min(abs(prec - desiredPrecision(di)));
      
      % Zoom in on region of interest
      i1 = find(prec>desiredPrecision(di));
      if isempty(i1);
        mx = max(score(:));
      else
        mx = ds(i1(1));
      end
      i2 = find(prec<desiredPrecision(di));
      if isempty(i2);
        mn = min(score(:));
      else
        mn = ds(i2(end));
      end
      ds = linspace(mn,mx,nn);
      prec0 = prec;
    end
    xcutoff(di) = ds(I);
    calcprec(di) = prec(I);
    calcrec(di) = rec(I);
  end
  
  [scoreRange,I] = sort(scoreRange);
  precRange = precRange(I);
  recRange = recRange(I);
  
  precSave{si} = precRange;
  recSave{si} = recRange;
  scoreSave{si} = scoreRange;
end


figure,hold on
plot(recSave{1},precSave{1},'r')
plot(recSave{2},precSave{2},'g')
plot(recSave{3},precSave{3},'c')
plot(recSave{4},precSave{4},'k')
legend('vote','avg1','avg2','neg')
title(['N = ' num2str(size(scoreMatrix,2)) ' models'])
set(gca,'xtick',0:.01:0.15)
axis([0 0.15 0.5 1])
grid on

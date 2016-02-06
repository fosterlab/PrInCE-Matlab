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



%% Try out the "new" feature

% Get data
I = possibleInts(:);
y = TP_Matrix(I);
y(y>0) = 1;
y(y~=1) = -1;

X2{1} = Dist.R2(I);
X2{2} = [Dist.R2(I) Dist.Euc(I) Dist.Center(I) Dist.Ngauss(I) Dist.CoApex(I) Dist.AUC(I)];
X2{3} = [Dist.R2(I) Dist.Euc(I) Dist.Center(I) Dist.Ngauss(I) Dist.CoApex(I) Dist.AUC(I) Dist.RealOverlap(I)];
X2{4} = [Dist.R2(I) Dist.Euc(I) Dist.Center(I) Dist.Ngauss(I) Dist.CoApex(I) Dist.AUC(I) Dist.R2raw(I)];
X2{5} = [Dist.R2(I) Dist.Euc(I) Dist.Center(I) Dist.Ngauss(I) Dist.CoApex(I) Dist.AUC(I) Dist.R2raw(I) Dist.RealOverlap(I)];

clear scores
for jj = 1:length(X2)
  jj
  X = X2{jj};
  Nd = size(X,2);
  
  % Soft whiten data
  eps = 2e-16;
  wmu = zeros(Nd,1);
  wstd = zeros(Nd,1);
  for ii = 1:Nd
    wmu(ii) = nanmean(X(:,ii));
    wstd(ii) = nanstd(X(:,ii));
    X(:,ii) = (X(:,ii) - wmu(ii)) / 2 / (wstd(ii) + eps);
  end
  
  
  Nmodel = 1;
  score = nan(size(X,1),Nmodel);
  feats = nan(Nmodel,size(X,2));
  for iter = 1:Nmodel
    
    % Make training and testing data
    nn = 1000; % length of training data
    % balance training data
    I1 = find(y==1 & sum(isnan(X),2)==0);
    I1 = I1(randsample(length(I1),nn/2));
    I0 = find(y==-1 & sum(isnan(X),2)==0);
    I0 = I0(randsample(length(I0),nn/2));
    Itrain = ismember(1:length(y),[I1;I0]);
    Ipred = ~Itrain;
    Xtr = X(Itrain,:);
    ytr = y(Itrain);
    Xnew = X(Ipred,:);
    
    % Feature selection
    feats(iter,:) = IndFeat(Xtr,ytr);
    f2consider = find(feats(iter,:) > 2);
    
    % Fit Naive Bayes model
    nab = fitcnb(Xtr(:,f2consider),ytr);
    [~,scoretmp] = predict(nab,Xnew(:,f2consider));
    %nab = fitcnb(Xtr,ytr);
    %[~,scoretmp] = predict(nab,Xnew);
    
    scores{jj}(Ipred,iter) = scoretmp(:,2);
  end
end



disp('Calculating precision/recall...')
precSave = cell(size(scores));
recSave = cell(size(scores));
scoreSave = cell(size(scores));
desiredPrecision = [.5 .7 .95];
maxIter = 7;
nn = 25;
for si = 1:length(scores)
  si
  score = scores{si};
  
  xcutoff = nan(size(desiredPrecision));
  calcprec = zeros(size(xcutoff));
  calcrec = zeros(size(xcutoff));
  scoreRange = nan(nn*maxIter*length(desiredPrecision),1);  % Save these for
  precRange = nan(size(scoreRange));                        % plotting
  recRange = nan(size(scoreRange));                         % later.
  fprRange = nan(size(scoreRange));                        % plotting
  tprRange = nan(size(scoreRange));
  for di = 1:length(desiredPrecision)
    ds = linspace(min(score(:)),max(score(:)),nn); % start off with coarse search
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
      fpr = nan(nn,1);
      tpr = nan(nn,1);
      for dd =1:length(ds)
        pos = score>ds(dd);
        TP = sum(pos & y==1);
        FP = sum(pos & y==-1);
        FN = sum(~pos & y==1);
        TN = sum(~pos & y==-1);
        prec(dd) = TP/(TP+FP);
        rec(dd) = TP/(TP+FN);
        fpr(dd) = FP/(FP+TN);
        tpr(dd) = TP/(TP+FN);
      end
      deltaPrec = nanmean(abs(prec - prec0));
      
      % Save vectors for plotting
      i1 = find(isnan(scoreRange),1,'first');
      I = i1 : i1+nn-1;
      scoreRange(I) = ds;
      precRange(I) = prec;
      recRange(I) = rec;
       fprRange(I) = fpr;
      tprRange(I) = tpr;
      
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
  tprRange = tprRange(I);
  fprRange = fprRange(I);
  
  precSave{si} = precRange;
  recSave{si} = recRange;
  scoreSave{si} = scoreRange;
  fprSave{si} = fprRange;
  tprSave{si} = tprRange;
end

colordef black
figure(1),hold on
plot(recSave{1},precSave{1},'--w')
plot(recSave{2},precSave{2},'--r')
plot(recSave{3},precSave{3},'--g')
plot(recSave{4},precSave{4},'--c')
plot(recSave{5},precSave{5},'--m')
legend('R2','All','All+Overlap','All+Raw','All+Overlap+Raw')
title(['N = ' num2str(size(scoreMatrix,2)) ' models'])
axis([0 1 0 1])
grid on

figure(2),hold on
plot(fprSave{1},tprSave{1},'--w')
plot(fprSave{2},tprSave{2},'--r')
plot(fprSave{3},tprSave{3},'--g')
plot(fprSave{4},tprSave{4},'--c')
plot(fprSave{5},tprSave{5},'--m')
plot([0 1],[0 1],'--w')
legend('R2','All','All+Overlap','All+Raw','All+Overlap+Raw')
title(['N = ' num2str(size(scoreMatrix,2)) ' models'])
axis([0 1 0 1])
grid on
colordef white


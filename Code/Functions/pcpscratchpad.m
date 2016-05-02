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



%% Compare importdata for different input files

% Alignment
% _Combined_OutputGaus_rep
f1 = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/Data/GaussBuild/MvsL_Combined_OutputGaus_rep1.csv';
f2 = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Data/GaussBuild/MvsL_Combined_OutputGaus_rep1.csv';
f3 = '/Users/Mercy/Academics/Foster/Jenny_PCPSILAC/PCPSILAC_Analysis/Data/GaussBuild/MvsL_Combined_OutputGaus_rep1.csv';
% _Summary_Gausians_for_individual_proteins_rep
f4 = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/Data/GaussBuild/MvsL_Summary_Gausians_for_individual_proteins_rep1.csv';
f5 = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Data/GaussBuild/MvsL_Summary_Gausians_for_individual_proteins_rep1.csv';
f6 = '/Users/Mercy/Academics/Foster/Jenny_PCPSILAC/PCPSILAC_Analysis/Data/GaussBuild/MvsL_Summary_Gausians_for_individual_proteins_rep1.csv';

d1 = importdata(f1);
d2 = importdata(f2);
d3 = importdata(f3);
d4 = importdata(f4);
d5 = importdata(f5);
d6 = importdata(f6);



%% Import and explore Craig's data

fn = '/Users/Mercy/Downloads/proteinGroups.txt';

% % importdata
% data1 = importdata(fn);
%
%
% % textscan
% fid = fopen(fn);
% C = textscan(fid, '%s','delimiter', '\n');
% fclose(fid);
% data2 = cell(size(C{1}));
% for ii = 1:size(data1,1)
%   data1{ii} = strsplit(C{1}{ii},'\t');
% end
% fclose(fid);
%
%
% % readtable
% T = readtable(fn);

% fgetl
fid = fopen(fn);
data = cell(10000,200);
cc = 0;
line = fgetl(fid);
while line~=-1
  cc = cc+1;
  tabs = find(ismember(line,'	'));
  data{cc,1} = line(1:tabs(1)-1);
  for ii = 1:length(tabs)-1
    data{cc,ii+1} = line(tabs(ii)+1: tabs(ii+1)-1);
  end
  data{cc,ii+1} = line(tabs(end):end);
  line = fgetl(fid);
end
data = data(1:cc,:);
fclose(fid);
N = size(data,1)-1;

% find columns with
%   H/L
%   M/L
tmp = data(1,1:140);
h1 = find(ismember(tmp,{'Ratio H/L'}));
h2 = find(ismember(tmp,{'Ratio H/L 1'}));
h3 = find(ismember(tmp,{'Ratio H/L 2'}));
m1 = find(ismember(tmp,{'Ratio M/L'}));
m2 = find(ismember(tmp,{'Ratio M/L 1'}));
m3 = find(ismember(tmp,{'Ratio M/L 2'}));

mvsl = zeros(size(data,1),3);
hvsl = zeros(size(data,1),3);
for ii = 2:N
  mvsl(ii,1) = str2num(data{ii,m1});
  mvsl(ii,2) = str2num(data{ii,m2});
  mvsl(ii,3) = str2num(data{ii,m3});
  hvsl(ii,1) = str2num(data{ii,h1});
  hvsl(ii,2) = str2num(data{ii,h2});
  hvsl(ii,3) = str2num(data{ii,h3});
end


pp1 = zeros(N-1,1);
pp2 = zeros(N-1,1);
rat = zeros(N-1,1);
for ii = 1:N-1
  pp1(ii) = ttest3(log10(mvsl(ii+1,:)),log10(hvsl(ii+1,:)));
  pp2(ii) = ttest3(mvsl(ii+1,:),hvsl(ii+1,:));
  rat(ii) = mean(mvsl(ii,:))/mean(hvsl(ii,:));
end


% Volcano plot
figure,hold on
scatter(log2(rat),-log10(pp1))
I = abs(log2(rat))>1 & pp1<.01;
scatter(log2(rat(I)),-log10(pp1(I)),'r','filled')
xlabel('Log2 fold change, M/H')
ylabel('-log10 ttest pvalue')



%% Solve Alignment bug
% Summary_gausian_infomration is not doing it's job. It should be saying how many Gaussians were fit
% for each protein name. That's it.

% Alignment

ci = 1;
rep = 1;
ng = Summary_gausian_infomration{ci,rep}.data(:,1); % # of gaussians fit
listofnames = Gaus_import{ci,align_rep}.textdata(:,1);
Isingle = find(ng==1);

for ii = 1:length(Isingle)
  i = Isingle(ii) + 1;
  protName = Summary_gausian_infomration{ci,rep}.textdata{i,2};
  
  ng2(ii) = sum(ismember(listofnames,protName));
end


figure,hold on
plot(ng2 - ng(Isingle)')



%% Illustrate the problem with using Precision for unbalanced data

x = linspace(0,3,101);
sig = 0.2;
mu1 = 1.5;
mu2 = 1.7;
Nrange = [10^3 10^4 10^6];
tts = {'Balanced' 'Unbalanced' 'Really unbalanced!'};

for jj = 1:3
  % case 2: unbalanced data
  N1 = Nrange(jj);
  N2 = 1000;
  data1 = normrnd(mu1,sig,[N1 1]);
  data2 = normrnd(mu2,sig,[N2 1]);
  h1 = hist(data1,x); % neg
  h2 = hist(data2,x); % pos
  figure
  subplot(2,1,1),hold on
  p0 = patch([0 x x(end)],[0 h1 0],'r');
  p1 = patch([0 x x(end)],[0 h2 0],'g');
  p0.FaceAlpha = 0.4;
  p1.FaceAlpha = 0.4;
  grid on
  Prec = nan(size(x));
  Rec = nan(size(x));
  FPR = nan(size(x));
  TPR = nan(size(x));
  for ii = 2:length(x)
    TP = sum(data2>x(ii));
    FP = sum(data1>x(ii));
    TN = sum(data1<x(ii));
    FN = sum(data2<x(ii));
    
    FPR(ii) = FP/(FP+TN);
    TPR(ii) = TP/(TP+FN);
    
    Prec(ii) = TP/(TP+FP);
    Rec(ii) = TP/(TP+FN);
  end
  title(tts{jj})
  if jj == 1;ax = axis;end
  axis(ax)
  ylabel('Counts')
  subplot(2,2,3),hold on
  plot(Rec,Prec)
  axis([0 1 0 1])
  xlabel('Recall')
  ylabel('Precision')
  title('P-R curve')
  subplot(2,2,4),hold on
  plot(FPR,TPR)
  axis([0 1 0 1])
  xlabel('FPR')
  ylabel('TPR')
  title('ROC curve')
end


%% Figure out why Nick achieved high precision with the tissue data, but I don't

% This takes a while!!

y1 = TP_Matrix==1 & possibleInts;
y0 = TP_Matrix==0 & possibleInts;

cols = {'b' 'c' 'm' 'k'};

% 1. Is it due to the way I calculate score??
% Do I get low precision with Dist.Euc? Dist.R2? scoreMatrix?
figure,hold on%,subplot(1,2,1),hold on,subplot(1,2,2), hold on
hold on
for ii = 1:3
  ii
  if ii==1
    tmp = Dist.Euc;
    x = linspace(0, log10(max(tmp(:))),101).^10;
  elseif ii ==2
    tmp = Dist.R2;
    x = linspace(0, max(tmp(:)),101);
  elseif ii ==3
    tmp = 1 - scoreMatrix_nb;
    x = linspace(0, 1,201).^10;
  elseif ii ==4
    tmp = 1 - scoreMatrix_svm;
    x = linspace(min(tmp), max(tmp),201);
  end
  FPR = nan(size(x));
  TPR = nan(size(x));
  Rec = nan(size(x));
  Prec = nan(size(x));
  for xi = 1:length(x)
    xx = x(xi);
    
    TP = sum(y1(:) & tmp(:)<xx);
    FP = sum(y0(:) & tmp(:)<xx);
    TN = sum(y0(:) & tmp(:)>xx);
    FN = sum(y1(:) & tmp(:)>xx);
    
    FPR(xi) = FP/(FP+TN);
    TPR(xi) = TP/(TP+FN);
    
    Prec(xi) = TP/(TP+FP);
    Rec(xi) = TP/(TP+FN);
  end
  stairs(Rec,Prec,'color',cols{ii})
  pause(.0001)
end
%subplot(1,2,1)
xlabel('Recall','fontsize',12)
ylabel('Precision','fontsize',12),set(gca,'fontsize',12)
title('PR Curve - Training on nonbalanced data','fontsize',14)
legend('Feature 1','Feature 2','Naive Bayes','SVM','location','northeast')


I = randsample(find(~isnan(scoreMatrix_svm) & ~isnan(scoreMatrix_nb)),length(1:500:length(scoreMatrix_svm)));
rocout = roc([Dist.R2(I) TP_Matrix(I)==1 & possibleInts(I)]);
rocout2 = roc([1-scoreMatrix_nb(I) TP_Matrix(I)==1 & possibleInts(I)]);
rocout3 = roc([1-scoreMatrix_svm(I) TP_Matrix(I)==1 & possibleInts(I)]);
rocout4 = roc([Dist.Euc(I) TP_Matrix(I)==1 & possibleInts(I)]);

figure,hold on
plot(rocout4.xr,rocout4.yr,'b')
plot(rocout.xr,rocout.yr,'c')
plot(rocout2.xr,rocout2.yr,'m')
plot(rocout3.xr,rocout3.yr,'k')
xlabel('FPR','fontsize',12)
ylabel('TPR','fontsize',12)
set(gca,'fontsize',12)
title('ROC Curve','fontsize',14)
legend('Feature 1','Feature 2','Naive Bayes','SVM','location','southeast')
plot([0 1],[0 1],':r')




%% Understand enrichment analysis
% Complex vs rest-of-the-proteins.
% Parameters:
%   N = total proteins, e.g. 1000
%   Nc = complex size, e.g. 10
%   Nl = how many times the label occurs
%   Nlc = number of labels in this complex
%
% contingency table:
%                label?
%             Y        N
%             _________________
% complex? Y |Nlc      Nc-Nlc
%          N |Nl-Nlc   N-Nc-Nlc

% CONCLUSIONS:
% 1. The minimum size for a complex is 3.
% 2. Medium size complexes are the best, i.e. 3<Nc<50.
% 3. Labels need to occur at least 3 times.
% 4. As long as labels occur >3 times and complexes aren't big, total number of labels doesnt' matter.
% 5. For rare labels (Nl<5), complexes can't be to big (Nc<75).

% hard-coded parameters
N = 1000;

% explored parameters
Nc_range = [2:5 10 15 20 25 35 50 75 100];
Nl_range = 1:2:21;
Nlc_range = 1:2:21;

figure
for ii = 1:length(Nc_range)
  subplot(3,4,ii)
  Nc = Nc_range(ii);
  pval = ones(length(Nl_range), length(Nlc_range)) * 1.5;
  pthreshold = ones(length(Nl_range), length(Nlc_range));
  for jj = 1:length(Nl_range)
    Nl = Nl_range(jj);
    for kk = 1:length(Nlc_range)
      Nlc = Nlc_range(kk);
      if Nlc>Nc || Nlc>Nl; continue; end
      
      T = [Nlc Nc-Nlc; Nl-Nlc N-Nc-Nlc];
      [~,pval(jj,kk)] = fishertest(T);
      
      Ncomp = 75 * 500 / Nl;
      pthreshold(jj,kk) = 0.05 / Ncomp;
      
    end
  end
  
  pval(pval>pthreshold & pval<=1) = 1;
  
  imagesc(Nlc_range,Nl_range,pval),caxis([0 1.5])
  axis xy
  xlabel('labels in complex')
  ylabel('total labels')
  title(num2str(Nc))
end



%% Understand enrichment analysis 2
% Question: For a given complex size and number of label occurrences, in order to be "enriched",
% what is the minimum fraction of protein members that need to have the label?

%Nc_range = [2:5 10 15 20 25 35 50 75 100];
%Nl_range = 1:2:21;
Nc_range = 1:2:100;
Nl_range = 1:25;

pth = 1e-7;

for N = [100 500 1000 5000 10000 100000];
  
  fraction = nan(length(Nc_range),length(Nl_range));
  for ii = 1:length(Nc_range)
    Nc = Nc_range(ii);
    for jj = 1:length(Nl_range)
      Nl = Nl_range(jj);
      
      Nlc = 0;
      pval = 1;
      while pval>pth
        Nlc = Nlc+1;
        if Nlc>Nc || Nlc>Nl || Nc+Nlc>N
          fraction(ii,jj) = 1.5;
          break
        end
        T = [Nlc Nc-Nlc; Nl-Nlc N-Nc-Nlc];
        [~,pval] = fishertest(T);
        fraction(ii,jj) = Nlc/Nc;
      end
      
    end
  end
  
  figure
  imagesc(fraction)
  caxis([0 1.5])
  axis xy
  set(gca,'xtick',1:length(Nl_range),'xticklabel',num2cell(Nl_range),...
    'ytick',1:length(Nc_range),'yticklabel',num2cell(Nc_range))
  colorbar
  xlabel('How many times does the label occur?')
  ylabel('Size of complex')
  title(['What fraction of the complex needs to have the label? N=' num2str(N)])
  
end


%% Play around with geometric accuracy
% CONCLUSIONs:
% - GA is a good measure.
% - It's not penalized by novel predictions, either whole complexes or complex members.
% - Sn is the ratio reference-proteins-in-the-predicted : total-reference-proteins.
% - PPV is the ratio predicted-proteins-in-the-reference : total-predicted-proteins.


clear ref1 pred1 ref2 pred2 ref3 pred3

% what if predicted complexes have NOVEL members?
clear a b c d
a{1} = 1:3;
b{1} = 1:4;
c{1} = 1:100;
ga1a = geomacc(a,a); % no novel members
ga1b = geomacc(b,a); % 1 novel
ga1c = geomacc(c,a); % lots novel
[ga1a ga1b ga1c]
disp('--> G.A. is UNAFFECTED for novel members in predicted complexes.')

% what if references have an unmatched complex?
clear a b c d
a{1} = 1:3;
b{1} = 1:3; b{2} = 4:6;
c{1} = 1:3; c{2} = 4:6; c{3} = 7:100;
ga2a = geomacc(a,a); % no unmatched complexes
ga2b = geomacc(a,b); % 1 unmatched
ga2c = geomacc(a,c); % lots unmatched
[ga2a ga2b ga2c]
disp('--> G.A. is PENALIZED by unmatched reference complexes.')

% what if predicted have an unmatched complex?
clear a b c d
a{1} = 1:3;
b{1} = 1:3; b{2} = 4:6;
c{1} = 1:3; c{2} = 4:6; c{3} = 7; c{4} = 8; c{5} = 9:10; c{6} = 11:15;
ga3a = geomacc(a,a); % no unmatched complexes
ga3b = geomacc(b,a); % 1 unmatched
ga3c = geomacc(c,a); % lots unmatched
[ga3a ga3b ga3c]
disp('--> G.A. is UNAFFECTED by unmatched predicted complexes.')

% what if predicted&complex size increases?
clear a b c d
a{1} = 1:3;
b{1} = 1:4;
c{1} = 1:100;
ga4a = geomacc(a,a); % size=3
ga4b = geomacc(b,b); % size=4
ga4c = geomacc(c,c); % size=100
[ga4a ga4b ga4c]
disp('--> G.A. is UNAFFECTED by larger complex sizes.')


%% Look at chromatograms for the proteins Nick selected for microscopy
%
% IQGAP1 - cdc42
% STK4/MST1 - DLP1
% WD(l)rp47 - PTS1r
%
% P46940 - P60953  ?>  47% precision, 11093 / 14226 interactions
% Q13043 - O00429  ?>  Q13043 not found
% O94967 - P50542  ?>  both found, but not interacting

rawdata = cell(2,1);
txt_val = cell(2,1);
for ii = 1:2
  [rawdata{ii},txt_val{ii}] = xlsread(user.MQfiles{ii});
  
  % Remove first column as the replicate
  replicate = rawdata{ii}(:,1);
  rawdata{ii} = rawdata{ii}(:,2:end);
  
  
  
  %   I1a = find(ismember(txt_val{ii}(:,1),'P46940')) - 1;
  %   I1b = find(ismember(txt_val{ii}(:,1),'P60953')) - 1;
  %   I2a = find(ismember(txt_val{ii}(:,1),'Q13043')) - 1;
  %   I2b = find(ismember(txt_val{ii}(:,1),'O00429')) - 1;
  %   I3a = find(ismember(txt_val{ii}(:,1),'O94967')) - 1;
  %   I3b = find(ismember(txt_val{ii}(:,1),'P50542')) - 1;
  I1a = strmatch('P46940', txt_val{ii}(:,1)) - 1;
  I1b = strmatch('P60953', txt_val{ii}(:,1)) - 1;
  I2a = strmatch('Q13043', txt_val{ii}(:,1)) - 1;
  I2b = strmatch('O00429', txt_val{ii}(:,1)) - 1;
  I3a = strmatch('O94967', txt_val{ii}(:,1)) - 1;
  I3b = strmatch('P50542', txt_val{ii}(:,1)) - 1;
  
  figure
  subplot(2,2,1),hold on
  plot(rawdata{ii}(I1a,:)','g')
  plot(rawdata{ii}(I1b,:)','r')
  subplot(2,2,2),hold on
  plot(rawdata{ii}(I2a,:)','g')
  plot(rawdata{ii}(I2b,:)','r')
  subplot(2,2,3),hold on
  plot(rawdata{ii}(I3a,:)','g')
  plot(rawdata{ii}(I3b,:)','r')
  
end



%% Why does precision drop when combining interaction lists?
% Theory: TP are common between lists (signal), FP are unique (noise).
%
% Run testROC_tissue.m to make the data


key = rand(1,100);

% Make list A
Ia = find(dataOld.data(:,5)==1);
listA = zeros(length(Ia),3);
for ii = 1:length(Ia)
  interaction = dataOld.text{Ia(ii),1};
  listA(ii,1) = sum(double(interaction) .* key(1:length(interaction)));
  listA(ii,2) = dataOld.data(Ia(ii),1);
  listA(ii,3) = dataOld.data(Ia(ii),2);
end
listA = unique(listA,'rows');
TPa = sum(listA(:,2)==1 & listA(:,3)==1);
FPa = sum(listA(:,2)==1 & listA(:,3)==0);

% Make list B
Ib = find(dataOld.data(:,6)==1);
listB = zeros(length(Ib),3);
for ii = 1:length(Ib)
  interaction = dataOld.text{Ib(ii),1};
  listB(ii,1) = sum(double(interaction) .* key(1:length(interaction)));
  listB(ii,2) = dataOld.data(Ib(ii),1);
  listB(ii,3) = dataOld.data(Ib(ii),2);
end
listB = unique(listB,'rows');
TPb = sum(listB(:,2)==1 & listB(:,3)==1);
FPb = sum(listB(:,2)==1 & listB(:,3)==0);

% Make unique(union(A,B))
listAB = unique([listA; listB],'rows');
TPab = sum(listAB(:,2)==1 & listAB(:,3)==1);
FPab = sum(listAB(:,2)==1 & listAB(:,3)==0);

% Calculate the ratio of FP:TP in A, B and combined list
ratioA = FPa / TPa;
ratioB = FPb / TPb;
ratioAB = FPab / TPab;
[ratioA ratioB ratioAB]

% Calculate how what fraction of TP and FP were removed
fracTPremoved = (TPa + TPb - TPab) / (TPa + TPb);
fracFPremoved = (FPa + FPb - FPab) / (FPa + FPb);
[fracTPremoved fracFPremoved]



%% Understand why some Tissue replicates have no interactions
% Compare scoreMatrix between replicates
% Run Initialize in ROC_PCPSILAC

prct = [50 75 90 95 99];

rep_prctile = zeros(6,length(prct)+1);
for ii = 1:(number_of_replicates*number_of_channels)
  sf = [datadir2 'score_rep' num2str(ii) '.mat'];
  load(sf)
  
  clear TP_Matrix possibleInts Protein inverse_self Chromatograms Dist
  
  for jj = 1:length(prct)
    rep_prctile(ii,jj) = prctile(scoreMatrix(:),prct(jj));
  end
  rep_prctile(ii,jj+1) = max(scoreMatrix(:));
end



%% Play with Tukey. Does it work for Ali's 30-group, 2-treatment, 3-replicate data?

data = rand(30,2,3);
data(1,1,:) = data(1,1,:) + 1;
data(2,1,:) = data(2,1,:) + 1;
data(3,1,:) = data(3,1,:) + 1;

cc = 0;
X = zeros(30*2*3,1);
G = zeros(size(X));
for ii = 1:30 % group
  for jj = 1:2 % treatment
    for kk = 1:3 % replicate
      cc = cc+1;
      X(cc) = data(ii,jj,kk);
      G(cc) = jj;
    end
  end
end

[~,~,stats] = anovan(X,G,'display','off');


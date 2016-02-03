function [score, feats] = scorenb(Dist,possibleInts,TP_Matrix)
% Use an SVM classifier to predict protein interactions.

% Get data
I = possibleInts(:);
labels = TP_Matrix(:);
y = labels(:);
y(y>0) = 1;
y(y~=1) = -1;
%data1 = Dist.R2(:); % 1 - R^2
%data2 = Dist.Euc(:);
%data3 = Dist.Center(:);
X = [Dist.R2(:) Dist.Euc(:) Dist.Center(:) Dist.Ngauss(:) Dist.CoApex(:) Dist.AUC(:)];
%X = [Dist.R2(:) Dist.Euc(:) Dist.CoApex(:) Dist.Ngauss(:) Dist.AUC(:)];
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


Nmodel = 15;
score = nan(size(X,1),Nmodel);
feats = nan(Nmodel,size(X,2));
for iter = 1:Nmodel
  
  % Make training and testing data
  nn = 1000; % length of training data
  % balance training data
  I1 = find(y==1 & I);
  I1 = I1(randsample(length(I1),nn/2));
  I0 = find(y==-1);
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
  
  score(Ipred,iter) = scoretmp(:,2);
end

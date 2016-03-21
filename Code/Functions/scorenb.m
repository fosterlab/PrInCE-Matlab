function [score, feats] = scorenb(Dist,possibleInts,TP_Matrix)
% Use a Naive Bayes classifier to predict protein interactions.

% Get data
I = possibleInts(:);
labels = TP_Matrix(:);
y = labels(:);
y(y>0) = 1;
y(y~=1) = -1;
X = [Dist.R2(:) Dist.Euc(:) Dist.Center(:) Dist.Ngauss(:) Dist.CoApex(:) Dist.AUC(:)];% Dist.RawOverlap(:) Dist.R2raw(:)];
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

Nlabel1 = sum(y==1);
trainingLength = min([1000 round(Nlabel1 * 0.8)]);

Nmodel = 15;
score = nan(size(X,1),Nmodel);
feats = nan(Nmodel,size(X,2));
for iter = 1:Nmodel
  % Make training and testing data
  
  % balance training data
%   I1 = find(y==1 & I);
%   I1 = I1(randsample(length(I1),round(trainingLength/2)));
%   I0 = find(y==-1);
%   I0 = I0(randsample(length(I0),round(trainingLength/2)));
%   Itrain = ismember(1:length(y),[I1;I0]);
  
  % non-balanced training data
  Iall = find(I);
  Iall = Iall(randsample(length(Iall),trainingLength));
  Itrain = ismember(1:length(y),Iall);
  sum(y(Iall)==1)
  
  Ipred = ~Itrain;
  Xtr = X(Itrain,:);
  ytr = y(Itrain);
  Xnew = X(Ipred,:);
  
  % Feature selection
  feats(iter,:) = IndFeat(Xtr,ytr);
  if sum(feats(iter,:)<=2)==size(feats,2)
    feats(iter,:) = 3;
  end
  f2consider = find(feats(iter,:) > 2);
  
  % Fit Naive Bayes model
  nab = fitcnb(Xtr(:,f2consider),ytr);
  [~,scoretmp] = predict(nab,Xnew(:,f2consider));
  
  score(Ipred,iter) = scoretmp(:,2);
end


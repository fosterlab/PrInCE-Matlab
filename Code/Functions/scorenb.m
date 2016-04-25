function [score, feats] = scorenb(Dist,possibleInts,TP_Matrix)
% Use a Naive Bayes classifier to predict protein interactions.

% Get data
I = possibleInts(:);
labels = TP_Matrix(:);
y = labels(:);
y(y>0) = 1;
y(y~=1) = -1;
fn = fieldnames(Dist);
X = nan(length(Dist.(fn{1})(:)),length(fn));
for ii = 1:length(fn);
  X(:,ii) = Dist.(fn{ii})(:);
end
Nd = size(X,2);

% Set NaN's to maxiumum value
for ii = 1:Nd
  nanelements = isnan(X(:,ii));
  X(nanelements,ii) = max(X(~nanelements,ii));
end

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
  go = 1;
  train_iter = 0;
  iterMax = 100;
  while go
    train_iter = train_iter+1;
    Iall = find(I);
    Iall = Iall(randsample(length(Iall),trainingLength));
    Itrain = ismember(1:length(y),Iall);
    
    Ipred = ~Itrain;
    Xtr = X(Itrain,:);
    ytr = y(Itrain);
    Xnew = X(Ipred,:);
    
    % Ensure there are at least 3 data points in each class
    go1 = sum(ytr>0)<3 | sum(ytr<0)<3;
    
    % Ensure no class-variable pair has zero variance
    testvar = nan(2,Nd);
    for jj = 1:Nd
      testvar(1,jj) = var(Xtr(ytr>0,jj));
      testvar(2,jj) = var(Xtr(ytr<0,jj));
    end
    go2 = sum(testvar(:)==0) > 0;
    
    go = (go1 | go2) & train_iter<=iterMax;
  end
  
  
  % Feature selection
  feats(iter,:) = IndFeat(Xtr,ytr);
  if sum(feats(iter,:)<=2)==size(feats,2)
    feats(iter,:) = 3;
  end
  f2consider = find(feats(iter,:) > 2)
    
  % Fit Naive Bayes model
  nab = fitcnb(Xtr(:,f2consider),ytr);
  [~,scoretmp] = predict(nab,Xnew(:,f2consider));
  
  score(Ipred,iter) = scoretmp(:,2);
end


function [score, feats] = scoresvm(Dist,possibleInts,TP_Matrix)
% Use an SVM classifier to predict protein interactions.

% Get data
I = possibleInts(:);
y = TP_Matrix(:);
y(y>0) = 1;
y(y~=1) = -1;
fn = fieldnames(Dist);
X = nan(length(Dist.(fn{1})(:)),length(fn));
for ii = 1:length(fn);
  X(:,ii) = Dist.(fn{ii})(:);
end
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
trainingLength = min([8000 round(Nlabel1 * 0.8)]);

iterMax = 15;
score = nan(size(X,1),iterMax);
feats = nan(iterMax,size(X,2));
for iter = 1:iterMax
  % Make training and testing data
  
  % balance training data
  %   I1 = find(y==1 & I);
  %   I1 = I1(randsample(length(I1),nn/2));
  %   I0 = find(y==-1);
  %   I0 = I0(randsample(length(I0),nn/2));
  %   Itrain = ismember(1:length(y),[I1;I0]);
  
  % non-balanced training data
  Iall = find(I);
  Iall = Iall(randsample(length(Iall),trainingLength));
  Itrain = ismember(1:length(y),Iall);
  
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
  
  % Fit svm model
  svm = fitcsvm(Xtr(:,f2consider),ytr);
  [~,scoretmp] = predict(svm,Xnew(:,f2consider));
  
  score(Ipred,iter) = scoretmp(:,2);
end


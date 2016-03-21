function score = scoresvm(Dist,possibleInts,TP_Matrix)
% Use an SVM classifier to predict protein interactions.

% Get data
I = possibleInts(:);
labels = TP_Matrix(:);
y = labels(:);
y(y>0) = 1;
y(y~=1) = -1;
data1 = Dist.R2(:); % 1 - R^2
data2 = Dist.Euc(:);
data3 = Dist.Center(:);
X = [data1 data2 data3];
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


iterMax = 2;
score = nan(size(X,1),iterMax);
for iter = 1:iterMax
  iter
  % Make training and testing data
  nn = 1000; % length of training data
  
  % balance training data
%   I1 = find(y==1 & I);
%   I1 = I1(randsample(length(I1),nn/2));
%   I0 = find(y==-1);
%   I0 = I0(randsample(length(I0),nn/2));
%   Itrain = ismember(1:length(y),[I1;I0]);
  
  % non-balanced training data
  Iall = find(I);
  Iall = Iall(randsample(length(Iall),nn));
  Itrain = ismember(1:length(y),Iall);
  
  Ipred = ~Itrain;
  Xtr = X(Itrain,:);
  ytr = y(Itrain);
  Xnew = X(Ipred,:);
  
  % Fit svm model
  svm = fitcsvm(Xtr,ytr,'Standardize',true);
  [~,scoretmp] = predict(svm,Xnew);
  
  score(Ipred,iter) = scoretmp(:,2);
end


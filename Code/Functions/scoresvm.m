function score = scoresvm(Dist,possibleInts,TP_Matrix)
% Use an SVM classifier to predict protein interactions.

% Get data
I = possibleInts;
labels = TP_Matrix(I(:));
y = labels(:);
y(y>0) = 1;
y(y~=1) = -1;
data1 = Dist.R2(I(:)); % 1 - R^2
data2 = Dist.Euc(I(:));
data3 = Dist.Center(I(:));
%data4 = Dist.dtw(I(:));
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


iterMax = 10;
score = nan(size(X,1),iterMax);
for iter = 1:iterMax

  % % Make training and testing data
  nn = 1000; % length of training data
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
  %CVSVMModel = crossval(svm);
  %FirstModel = CVSVMModel.Trained{1};
  %svm = fitcsvm(X,y,'Standardize',true);
  [~,scoretmp] = predict(svm,Xnew);
  %[precision_svm,recall_svm,~,auc_svm] = perfcurve(ynew,score_svm(:,2),1,'xcrit','prec');
  
  score(Ipred,iter) = scoretmp(:,2);
end
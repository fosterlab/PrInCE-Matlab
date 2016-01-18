


%% Calculate PR curve for just R^2
I = inverse_self & Int_matrix;
data = Dist.R2(I(:));
labels = TP_Matrix(I(:));
y = labels(:) > 0;

RRange = linspace(0, 1,201);
Recall = zeros(size(RRange));
Precision = zeros(size(RRange));
TPR = zeros(size(RRange));
FPR = zeros(size(RRange));
for di = 1:length(RRange)
  ypred2 = data<RRange(di);
  
  TP = sum(ypred2==1 & y==1);
  FP = sum(ypred2==1 & y==0);
  TN = sum(ypred2==0 & y==0);
  FN = sum(ypred2==0 & y==1);
  Recall(di) = TP/(TP+FN);
  Precision(di) = TP/(TP+FP);
  TPR(di) = TP/(TP+FN);
  FPR(di) = FP/(FP+TN);
end

figure
subplot(3,1,1),hold on
x = linspace(0,1,101);
h1 = hist(data(y==1,1),x);
h0 = hist(data(y==0,1),x);
plot(x,h1,'g')
plot(x,h0,'k')
subplot(3,1,2),
plot(Recall,Precision)
axis([0 1 0 1])
subplot(3,1,3),
plot(FPR,TPR)
axis([0 1 0 1])
set(gcf,'units','normalized','position',[.1 .1 .3 .9])



%% SVM: protein interactions, generate P-R curve
% To generate all variables, run ROC_PCPSILAC.m and pause after section 5 
%
% Training set = inverse_self & Int_matrix
% TP labels = TP_Matrix

% Get data
I = inverse_self & Int_matrix;
data1 = Dist.R2(I(:)); % 1 - R^2
data2 = Dist.Euc(I(:));
data3 = Dist.Center(I(:));
labels = TP_Matrix(I(:));
y = labels(:);% > 0;
y(y>0) = 1;
y(y~=1) = -1;
X = [data1 data2 data3];
Nd = size(X,2);
%X = data1;
%X(y==1,:) = repmat([0 0 0],sum(y==1),1); % cheat
%X(y==0,:) = X(y==0,:) + repmat([1 1 1],sum(y==0),1)*0.38; % cheat2

% Make training and testing data
if 1
  Itrain = zeros(size(y));
  Itrain(randsample(length(Itrain),1400)) = 1;
  Itrain = Itrain==1;
else % equalize 1s and 0s labels
  nn = 700;
  I1 = find(y==1);
  I1 = I1(randsample(length(I1),nn));
  I0 = find(y==-1);
  I0 = I0(randsample(length(I0),nn));
  Itrain = ismember(1:length(y),[I1;I0]);
  Ipred = ~Itrain;
end
Xtr = X(Itrain,:);
ytr = y(Itrain);
Xnew = X(Ipred,:);
ynew = y(Ipred);

% Whiten data
eps = 2e-16;
wmu = zeros(Nd,1);
wstd = zeros(Nd,1);
for ii = 1:Nd
  wmu(ii) = mean(Xtr(:,ii));
  wstd(ii) = std(Xtr(:,ii));
  X(:,ii) = (X(:,ii) - wmu(ii)) / 2 / (wstd(ii) + eps);
end

% Fit the model
%model = fitcsvm(Xtr,ytr);
%[label,score] = predict(model,Xnew);
svm = svmtrain(Xtr,ytr,'method','smo','kernel_function','rbf');
shift = svm.ScaleData.shift;
scale = svm.ScaleData.scaleFactor;
Xnew2 = bsxfun(@plus,Xnew,shift);
Xnew2 = bsxfun(@times,Xnew2,scale);
sv = svm.SupportVectors;
alphaHat = svm.Alpha;
bias = svm.Bias;
kfun = svm.KernelFunction;
kfunargs = svm.KernelFunctionArgs;
f = kfun(sv,Xnew2,kfunargs{:})'*alphaHat(:) + bias;
f = -f; % flip the sign to get the score for the +1 class

% Predict
label = svmclassify(svm,Xnew);
[precision_svm,recall_svm]=perfcurve(ynew,f,1,'xcrit','prec');

% Evaluate
TP = sum(label==1 & ynew==1);
FP = sum(label==1 & ynew==-1);
TN = sum(label==-1 & ynew==-1);
FN = sum(label==-1 & ynew==1);
Recall = TP/(TP+FN);
Precision = TP/(TP+FP);

[TP FP TN FN]
[Recall Precision]




%% Naive Bayes: protein interactions, generate P-R curve
% To generate all variables, run ROC_PCPSILAC.m and pause after section 5 

I = inverse_self & Int_matrix;
data1 = Dist.R2(I(:)); % 1 - R^2
data2 = Dist.Euc(I(:));
data3 = Dist.Center(I(:));
labels = TP_Matrix(I(:));
y = labels(:) > 0;
X = [data1 data2 data3];

% Train on Ntr points
% Ntr = length(y) - 499;
% Itr = zeros(size(y));
% Itr(randsample(length(Itr),Ntr)) = 1;
% Itr = Itr==1;
nn = 10000;
I1 = find(y==1);
I1 = I1(randsample(length(I1),nn));
I0 = find(y==0);
I0 = I0(randsample(length(I0),nn*5));
Itr = ismember(1:length(y),[I1;I0]);
Xtr = X(Itr,:);
ytr = y(Itr);

Itest = find(~Itr);
Itest = Itest(randsample(length(Itest),500));
Xtest = X(Itest,:);
ytest = y(Itest);

ypred = zeros(size(ytest));
Dprob = zeros(size(ytest));
for ii = 1:length(ytest)
  x = Xtest(ii,:); % test data

  prob = zeros(2,1);
  for ki = 1:2
    % calculate Bayesian probability 
    % p(Ck | x) = p(C_k) * p(x_1|C_k) * p(x_2|C_k) * ... * p(x_n|C_k)    
    Ck = ki - 1;
    prob(ki) = sum(ytr == Ck) / length(ytr); % prior
    
    % gaussian likelihood parameters
    xmu = mean(Xtr(ytr == Ck,:),1);
    xstd = std(Xtr(ytr == Ck,:),1);
    % beta likelihood parameters
    alpha = ((1-xmu)./xstd.^2 - 1./xmu).*xmu.^2;
    beta = alpha.*(1./xmu - 1);
    
    likelihood = zeros(size(x));
    for xi=1:size(X,2)
      if xi==1 % X(:,1) looks beta-distributed
        likelihood(xi) = pdf('beta',x(xi),alpha(xi),beta(xi));%exp(-(x(xi) - xmu(xi))^2 / xstd(xi)^2); % p(x_i|C_k)
      else % X(:,2:3) look normally-distributed
        likelihood(xi) = pdf('normal',x(xi),xmu(xi),xstd(xi));%exp(-(x(xi) - xmu(xi))^2 / xstd(xi)^2); % p(x_i|C_k)
      end
      prob(ki) = prob(ki)*likelihood(xi);
    end
    %[x xmu likelihood sum(ytr == Ck) / length(ytr)] 
  end
  % what does the classifier say?
  [~,I] = max(prob);
  ypred(ii) = I-1;
  Dprob(ii) = diff(prob);
  
  if mod(ii,10)==0
  disp([num2str(ii) ',    ' num2str(x(1))  ',    ' num2str(ypred(ii)) ',    ' num2str(ytest(ii))])
  end
  
  %[prob' ytest(ii) ypred(ii)]
  %pause
end

TP = sum(ypred==1 & ytest==1);
FP = sum(ypred==1 & ytest==0);
TN = sum(ypred==0 & ytest==0);
FN = sum(ypred==0 & ytest==1);
Recall = TP/(TP+FN);
Precision = TP/(TP+FP);
TPR = TP/(TP+FN);
FPR = FP/(FP+TN);

[TP FP TN FN]
[Recall Precision]


%% Calculate PR curve for Naive Bayes
%Dprob2 = -log(-Dprob); % what?!? why does this make sense!
Dprob2 = Dprob;
DprobRange = linspace(min(Dprob2),max(Dprob2),201);
Recall = zeros(size(DprobRange));
Precision = zeros(size(DprobRange));
TPR = zeros(size(DprobRange));
FPR = zeros(size(DprobRange));
for di = 1:length(DprobRange)
  ypred2 = Dprob2>DprobRange(di);
  
  TP = sum(ypred2==1 & ytest==1);
  FP = sum(ypred2==1 & ytest==0);
  TN = sum(ypred2==0 & ytest==0);
  FN = sum(ypred2==0 & ytest==1);
  Recall(di) = TP/(TP+FN);
  Precision(di) = TP/(TP+FP);
  TPR(di) = TP/(TP+FN);
  FPR(di) = FP/(FP+TN);
end

figure
subplot(3,1,1),hold on
x = linspace(0,1,101);
h1 = hist(X(y==1,1),x);
h0 = hist(X(y==0,1),x);
subplot(3,1,2),
plot(Recall,Precision)
axis([0 1 0 1])
subplot(3,1,3),
plot(FPR,TPR)
axis([0 1 0 1])
set(gcf,'units','normalized','position',[.1 .1 .3 .8])


% % Leave-one-out cross validation
% predClass = zeros(N,1);
% realClass = zeros(size(predClass));
% for vi = 1:N
%   
%   Itr = [1:vi-1 vi+1:N];    % training indices
%   Xtr = X(Itr,:);           % training data, X
%   classtr = class(Itr);     % training data, class
%   x = X(vi,:);              % test data
%   
%   % calculate Bayesian probability p(Ck | x) = p(C_k) * p(x_1|C_k) * p(x_2|C_k) * ... * p(x_n|C_k)
%   prob = zeros(length(class_cat),1);
%   for ki=1:length(class_cat)
%     Ck = class_cat(ki);
%     prob(ki) = sum(classtr == Ck) / length(classtr); % prior
%     
%     % likelihood parameters
%     xmu = mean(Xtr(classtr == Ck,:),1);
%     xstd = std(Xtr(classtr == Ck,:),1);
%     
%     likelihood = zeros(length(x));
%     for xi=1:D
%       likelihood(xi) = exp(-(x(xi) - xmu)^2 / 2 / xstd^2); % p(x_i|C_k)
%       prob(ki) = prob(ki)*likelihood(xi);
%     end
%   end
%   
%   % what does the classifier say?
%   [~,I] = max(prob);
%   predClass(vi) = class_cat(I);
%   
%   % what's the real class?
%   realClass(vi) = class(vi);
% end


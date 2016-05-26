function [threshold, varoutput] = calcPPIthreshold(score, class, desiredPrecision)

%CALCPPITHRESHOLD calculates the threshold needed to
%    determine interactions at the desired precision.
%
%    Calculates precision as a function of a range of
%    score values, and detects when precision crosses
%    the desired level. "Zooms in" on that score
%    value to obtain finer esimate of the required
%    score. In the case of multiple scores leading to
%    the same precision, pick the largest score.
%
%    Created by R. Greg Stacey, March 2016.

nn = 25;
Tol = 0.001; % get within 0.1% precision
maxIter = 20; % zoom in 20 times at most


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
    % ensemble VOTING
    TP = sum(score>ds(dd) & class==1);
    FP = sum(score>ds(dd) & class==0);
    FN = sum(score<ds(dd) & class==1);
    TN = sum(score<ds(dd) & class==0);
    prec(dd) = TP/(TP+FP);
    rec(dd) = TP/(TP+FN);
    fpr(dd) = FP/(FP+TN);
    tpr(dd) = TP/(TP+FN);
  end
  deltaPrec = nanmean(abs(prec - prec0));
  
  % Save vectors for plotting
  I = (iter-1)*nn+1 : iter*nn;
  scoreRange(I) = ds;
  precRange(I) = prec;
  recRange(I) = rec;
  tprRange(I) = tpr;
  fprRange(I) = fpr;
  
  %   % Zoom in on region of interest
  %   i1 = find(prec>desiredPrecision);
  %   if isempty(i1);
  %     mx = max(score(:));
  %   else
  %     mx = ds(i1(1));
  %   end
  %   i2 = find(prec<desiredPrecision);
  %   if isempty(i2);
  %     mn = min(score(:));
  %   else
  %     mn = ds(i2(end));
  %   end
  
  % Check if desiredPrecision is reached exactly
  tmp = prec - desiredPrecision;
  if sum(tmp==0)>0
    I = find(tmp==0);
    dsi = I(end);
    break
  end
  
  % If not, zoom in on the region of interest
  zeroCross = find(sign(tmp(1:end-1)) .* sign(tmp(2:end)) == -1);
  if isempty(zeroCross) % prec never crosses desiredPrecision
    i1 = find(prec>desiredPrecision);
    i2 = find(prec<desiredPrecision);
    if isempty(i1) % prec is always less than desiredPrecision
      mx = max(score(:));
      mn = ds(i2(end));
      dsI = [i2(end) nn];
    elseif isempty(i2) % prec is always greater than desiredPrecision
      mn = min(score(:));
      mx = ds(i1(1));
      dsI = [1 i1(1)];
    else
      error('calcPPIthreshold: error in algorithm')
    end
  else
    mn = ds(zeroCross(end));
    mx = ds(zeroCross(end)+1);
    dsI = [zeroCross(end) zeroCross(end)+1];
  end
  
  % Calculate how close to desiredPrecision(di) you got
  [calcTol,I] = min(abs(prec(dsI) - desiredPrecision));
  dsi = dsI(I);
  
  ds = linspace(mn,mx,nn);
  prec0 = prec;
end
threshold = ds(dsi);
calcprec = prec(dsi);
calcrec = rec(dsi);

[scoreRange,I] = sort(scoreRange);
precRange = precRange(I);
recRange = recRange(I);
tprRange = tprRange(I);
fprRange = fprRange(I);


varoutput.calcprec = calcprec;
varoutput.calcrec = calcrec;
varoutput.scoreRange = scoreRange;
varoutput.precRange = precRange;
varoutput.recRange = recRange;
varoutput.tprRange = tprRange;
varoutput.fprRange = fprRange;
varoutput.Ninteractions = recRange * sum(class==1);

 

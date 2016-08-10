function [threshold, varoutput] = calcPPIthreshold(score, class, desiredPrecision)

%CALCPPITHRESHOLD calculates the score threshold that 
%    results in a desired precision level.
%
%    Calculates precision as a function of score,
%    detects when when precision crosses the desired
%    level. "Zooms in" on that score value to obtain 
%    finer esimate of the required score.
%
%    score, an Nx1 vector is a measure of interaction 
%    confidence and ranges between 0 and 1. Score = 1 
%    denotes highest confidence in the interaction, 
%    score = 0 denotes lowest confidence. Score is 
%    often the output of a binary classifier.
%   
%    class is binary Nx1 vector. class = 1 denotes a
%    known interaction. class = 0 denotes a known
%    non-interaction.
%
%    desiredPrecision is the desired precision level,
%    expressed as a fraction between 0 and 1.
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
  
  % Ensure that prec = desiredPrecision counts as a "zero-crossing"
  tmp = prec - desiredPrecision;
  if sum(tmp==0)>0
    tmp = tmp + 10e-6;
  end
  
  % Zoom in on the region of interest
  % Either i) top end, ii) bottom end, iii) first zero-crossing.
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
    mn = ds(zeroCross(1));
    mx = ds(zeroCross(1)+1);
    dsI = [zeroCross(1) zeroCross(end)+1];
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

 

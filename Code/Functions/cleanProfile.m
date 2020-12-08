function cleandata = cleanProfile(rawdata,Nminconsecutive)
%CLEANPROFILE Pre-process co-fractionation profile
%   cleandata=CLEANPROFILE(rawdata) performs several pre-processing steps
%   on the co-fractionation profile rawdata. The Nx1 vector rawdata is a
%   co-fractionation profile for a single protein. cleandata is a (N+10)x1
%   vector, corresponding to the co-fractionation profile with values
%   padded on both sides. The pre-processing steps are:
%   
%   1. Impute single missing values with the average of neighboring values.
%   2. Add 5 nans to the start and end of rawdata.
%   3. Ensures at least 5 consecutive real values.
%   4. Replace nans with near-zero real values.
%   5. Smooth with a linear boxcar filter.
%
%   See also GAUSSBUILD.



% 0. Initialize

if nargin<1
  error('No input detected.')
end

[n1,n2] = size(rawdata);
if n1~=1 && n2~=1 || (n1==1 && n2==1)
  error('Input must be a 1-dimensional vector, i.e. a single co-fractionation profile.')
end

% transform into row vector if necessary
if n1>n2
  rawdata = rawdata';
end

tmpchrom = rawdata; % dummy variable


%% 1. Replace single missing values with mean of neighbours
steps = 1:5;

if ismember(1,steps)
  tmp = [1 0 1]; % look for a single nan flanked by real values, i.e. [1 0 1]
  tmp2 = strfind(~isnan(tmpchrom),tmp); % find where [1 0 1] occurs in the chromatogram
  for ii = 1:length(tmp2)
    I = tmp2(ii) + 1;
    tmpchrom(I) = mean(tmpchrom([I-1 I+1]));
  end
end


% 2. Add 5 nans to either side of the chromatogram
if ismember(2,steps)
  tmpchrom = [nan(1,5)*0.01 tmpchrom nan(1,5)*0.04];
end


% 3. Consecutive numbers, if less then 5 consecutive number removes chromogram and replace with 0.05
if ismember(3,steps)
  a = ones(1,Nminconsecutive); % the pattern to look for: [1 1 1]
  a2 = strfind(~isnan(tmpchrom),a); % find the pattern
  tmpchrom2 = zeros(size(tmpchrom)); % dummy variable
  for ii = 1:length(a2)
    I = a2(ii) : a2(ii) + Nminconsecutive-1;
    tmpchrom2(I) = tmpchrom(I);
  end
  tmpchrom = tmpchrom2;
end


% 4. Replace nans with near-zero noise
if ismember(4,steps)
  I = find(isnan(tmpchrom));
  tmpchrom(I) = rand(size(I))*0.01;
end


% 5. Smooth with a boxcar filter
%if ismember(5,steps)
%  tmpchrom = smooth(tmpchrom,4);
%end

cleandata = tmpchrom;



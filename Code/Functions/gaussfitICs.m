function ICs = gaussfitICs(cleanchrom)

% Generates initial conditions (ICs) for Gaussian fitting.
%
% Input: nx1 cleaned chromatogram, where n is the number of fractions.
% Output: 1x5 cell, where each element is a vector containing the ICs for a single Gaussian model.
%
% Adapted from Nichollas Scott's Gaus_build_24_1.m.
% Made by Greg Stacey on Nov 25 2015.



%% 0. Initialize

if nargin~=1
  error('Input must be a single variable, i.e. one chromatogram.')
end

[n1,n2] = size(cleanchrom);
if n1~=1 && n2~=1 || (n1==1 && n2==1)
  error('Input must be a 1-dimensional vector, i.e. a single chromatogram.')
end

% transform into row vector if necessary
if n1>n2
  cleanchrom = cleanchrom';
end



%% Find initial conditions for up to 5 Gaussians.

% Here I'm keeping the same syntax as Nick's original code. This could likely be changed to make
% things more readable.
% N5: Peak location (location)
% N6: Peak height

mH = max(cleanchrom); % maximum value of the chromatogram

% find the top 5 peaks of cleanchrom
[N6, N8] = findpeaks(cleanchrom,'NPeaks', 5);
[~,I] = sort(N6, 'descend');
N6 = N6(I);
N8 = N8(I);

% If there were 5 peaks found, use those. If <5 peaks, fill in with random guesses.
N5 = zeros(1,5);
N5(1:length(N6)) = N8;
for ii = length(N6)+1:5
  N5(ii) = random('Normal',N8(1),2);
  N6(ii) = mH;
end

% Store the peak height and locations in the right format
ICs{1} = [N6(1) N5(1) 2];
ICs{2} = [N6(1) N5(1) 2 N6(2) N5(2) 2];
ICs{3} = [N6(1) N5(1) 2 N6(2) N5(2) 2 N6(3) N5(3) 2];
ICs{4} = [N6(1) N5(1) 2 N6(2) N5(2) 2 N6(3) N5(3) 2 N6(4) N5(4) 2];
ICs{5} = [N6(1) N5(1) 2 N6(2) N5(2) 2 N6(3) N5(3) 2 N6(4) N5(4) 2 N6(5) N5(5) 2];


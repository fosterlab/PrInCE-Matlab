function M2 = mclpcp(M, params)

%MCLPCP Prunes protein-protein interaction matrices.
%    M2 = mclpcp(M,P) returns a pruned version of the interaction
%    matrix M according to the MCL algorithm developed by Stign van
%    Dongen. M is an NxN matrix where a non-zero entry in element
%    (i,j) indicates an interaction between protein_i and protein_j. P
%    is a structure that contains parameters for the MCL algorithm.
%
%    The parameters in P are
%       P.p - Exponent power used in the expansion.
%       P.minrep - Before MCL, ignore all interactions with a score
%                  less than this.
%       P.minval - Minimum value below which interactions are pruned.
%
% Adapted from software by Stijn van Dongen obtained here http://micans.org/mcl/
% References:
% 1. Stijn van Dongen, Graph Clustering by Flow Simulation, PhD thesis, University of Utrecht, May
%    2000. ( http://www.library.uu.nl/digiarchief/dip/diss/1895620/inhoud.htm )
% 2. Enright A.J., Van Dongen S., Ouzounis C.A., An efficient algorithm for large-scale detection
%    of protein families, Nucleic Acids Research 30(7):1575-1584 (2002).


if ~isfield(params,'p')
  p = 2;
else
  p = params.p;
end
if ~isfield(params,'minrep')
  minrep = 1;
else
  minrep = params.minrep;
end
if ~isfield(params,'minval')
  minval = 0.000001;
else
  minval = params.minval;
end

Nprot = size(M,1);

% Find and store non-zero entries of M
idx = find(M~=0);
vals = M(idx);
Icorum = M>99;

% Reduce M to binary
M = single(M>=minrep);

% Normalize M
total_unique_proteins_interacting_with = sum(M,2);
M2 = zeros(Nprot,Nprot);
for jj = 1:Nprot
  if total_unique_proteins_interacting_with(jj)>0
    M2(:,jj) = M(jj,:)./total_unique_proteins_interacting_with(jj);
  end
end
clear M

% Initialize
energy = 1;
delta_energy = 1;

%Repeat till the variation is less then 10% OR 10 interations have been
%evaluated. Note the limit on interations is due convergence always being
%outcome from MCL
iteration = 0;
while delta_energy>0.1 && iteration<10
  iteration=iteration+1;
  
  %Set emax
  emax = energy;
  
  m2 = M2 .^ p;       % inflation
  m2(Icorum) = 1;             % MY HACK! Corum interactions are always 1
  I = m2 < minval;            % pruning
  m2(I) = 0;
  dinv = diag(1./sum(m2));    % normalisation
  m2 = m2 * dinv;
  m2(isnan(m2))=0;            % Remove NaN with zeros
  
  % calculate residual energy
  maxs = max(m2);
  sqsums = sum(m2 .^ 2);
  energy = max(maxs - sqsums);
  delta_energy = abs((energy - emax) ./ emax);
  
  %set m2 to mTemp
  M2 = m2;
end

% Return M(M2>0)
idx2 = find(M2~=0);
I = idx(ismember(idx,idx2));
M2(I) = vals(ismember(idx,idx2));

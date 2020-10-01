function mr = matchingratio(predComplex, refComplex)

%MATCHINGRATIO Calculates the matching ratio between two
%    groups of complexes (predicted and reference).
%
% Ref:
% Detecting overlapping protein complexes in protein-protein 
% interaction networks Nature Methods 9, 471?472 (2012) 
% doi:10.1038/nmeth.1938 

Na = length(refComplex);
Nb = length(predComplex);

overlapMatrix = zeros(Na,Nb);

for ii = 1:Na
  for jj = 1:Nb
        
    overlap = length(intersect(refComplex{ii},predComplex{jj}))^2;
    overlapMatrix(ii,jj) = overlap / length(refComplex{ii}) / length(predComplex{jj});
    
  end
end


% Pick the optimal edges b/w predicted and reference complexes
% traverse the smallest dimension
if size(overlapMatrix,1)<size(overlapMatrix,2)
  overlapMatrix = overlapMatrix';
end
sortMatrix = nan(size(overlapMatrix));
for ii = 1:size(overlapMatrix,2)
  [~,sortMatrix(:,ii)] = sort(overlapMatrix(:,ii),'descend');
end
already_picked = zeros(size(overlapMatrix,1),1);
edges = zeros(size(overlapMatrix,2),1);
for ii = 1:size(overlapMatrix,1)
  tmp = zeros(size(overlapMatrix,2),1);
  for jj = 1:size(overlapMatrix,2)
    tmp(jj) = overlapMatrix(sortMatrix(ii,jj),jj);
  end
  [~,Iorder] = sort(tmp,'descend');
  
  for jj = 1:size(overlapMatrix,2)
    Ipred = Iorder(jj);
    Iref = sortMatrix(ii,Ipred);
    
    if edges(Ipred)~=0
      continue
    end
    
    if already_picked(Iref) == 0
      edges(Ipred) = Iref;
      already_picked(Iref) = 1;
    end
  end
end


% Calculate Matching Ratio
mr = 0;
for ii = 1:min(Na,Nb)
  if edges(ii)>0
    mr = mr + overlapMatrix(edges(ii), ii);
  end
end
mr = mr/Na;


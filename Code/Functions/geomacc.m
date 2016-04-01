function ga = geomacc(predComplex, refComplex)

%GEOMACC Calculates the geometric accuracy between two
%    groups of complexes (predicted and reference).
%
% Ref:
% Evaluation of clustering algorithms for protein-
% protein interaction networks. BMC Bioinf. 7, 488 
% (2006) doi: 10.1186/1471-2105-7-488

Na = length(refComplex);
Nb = length(predComplex);

T = zeros(Na,Nb);

Ni = nan(size(refComplex));
if size(Ni,1)==1 && size(Ni,2)>1
  Ni = Ni';
end
for ii = 1:Na
  Ni(ii) = length(refComplex{ii});
  for jj = 1:Nb
    T(ii,jj) = length(intersect(refComplex{ii},predComplex{jj}))^2;
  end
end


% Calculate sensitivity
sn = T ./ repmat(Ni,1,Nb);
Sn_ref = max(sn,[],2);
Sn = sum(Sn_ref .* Ni) / sum(Ni);


% Calculate PPV
ppv = T ./ repmat(sum(T),Na,1);
PPV_pred = max(ppv,[],1);
PPV = nansum(PPV_pred .* sum(T)) / sum(T(:));


% Calculate geometric accuracy, as defined in Brohee and van Helden (2006)
ga = sqrt(Sn * PPV);

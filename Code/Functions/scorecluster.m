function score = scorecluster(Dist)


% if necessary, transform data into list format
if isstruct(Dist)
  fn = fieldnames(Dist);
  Nfields = length(fn);
  data = zeros(length(Dist.(fn{1})(:)),length(fn));
  for ii = 1:Nfields
    data(:,ii) = Dist.(fn{ii})(:);
  end
else
  Nfields = size(Dist,2);
  % if necessary, impute missing data
  for ii = 1:Nfields
    Ibad = isnan(Dist(:,ii));
    Igood = find(~Ibad);
    imp = randsample(Igood,sum(Ibad),1);
    Dist(Ibad,ii) = Dist(imp,ii);
  end
  data = Dist;
end

% soft-whiten and center data
for ii = 1:Nfields
  data(:,ii) = data(:,ii) ./ nanstd(data(:,ii)); % soft-whiten
  data(:,ii) = data(:,ii) - nanmean(data(:,ii)); % center
end

% pca transform
[~, pcascore] = pca(data);
data = pcascore;

% cluster data with with k-means
[idx,C] = kmeans(data,2);

% which cluster is the "most similar"?
aa(1) = sum(nanmean(data(idx==1,:)));
aa(2) = sum(nanmean(data(idx==2,:)));
[~,I] = sort(aa,'descend');
C = C(I,:);
C1 = C(1,:);
C2 = C(2,:);

% project all data onto line between 2 centers
ap = data - repmat(C1,length(data),1);
ab = repmat(C2 - C1,length(data),1);
proj = repmat(C1,length(data),1) + repmat(dot(ap',ab')'./dot(ab',ab')',1,Nfields) .* ab;
score = zeros(length(data),1);
a3 = norm(ab(1,:));
for ii = 1:length(proj)
  a1 = norm(proj(ii,:) - C1);
  a2 = norm(proj(ii,:) - C2);
  s = 1;
  if a2>a3 && a1<a2
    s = -1;
  end
  score(ii) = norm(proj(ii,:) - C1) ./ a3 * s;
end

%score1 = sqrt(sum(abs(data - repmat(C1,size(data,1),1)).^2,2));
%score = score.*score1;

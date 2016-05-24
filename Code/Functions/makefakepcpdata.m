
%% Initialize

datadir = [user.maindir 'Input/']; % where data files live
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist([user.maindir 'Output/'], 'dir'); mkdir([user.maindir 'Output/']); end
if ~exist([user.maindir 'Output/tmp'], 'dir'); mkdir([user.maindir 'Output/tmp']); end

Nproteins = 200; % per replicate
Ninteractions = 10; % pairwise interactions
Ndifferentcomparisons = 3;
x = 1:user.Nfraction;



%% Make fake Chromatograms (chrom) to use as templates

chrom = zeros(Ninteractions,user.Nfraction);
Ngauss_fake = poissrnd(1,Ninteractions,1) + 1;
Ngauss_fake(Ngauss_fake>5) = 5;
for ii = 1:Ninteractions
  
  % chromatogram, shared between interactors and replicates
  C = rand(Ngauss_fake(ii),1)*(user.Nfraction - 10) + 5; % [5 50]
  H = rand(Ngauss_fake(ii),1)*15 + 3; % [3 18]
  W = rand(Ngauss_fake(ii),1)*6 + 2; % [2 8]
  for jj = 1:Ngauss_fake(ii)
    chrom(ii,:) = chrom(ii,:) + H(jj)*exp(-((x-C(jj))/W(jj)).^2);
  end
end

figure,
hist(Ngauss_fake,1:7)


%% Copy chrom into each replicate and channel

Chromatograms = cell(length(user.silacratios),1);
Proteins = cell(length(user.silacratios),1);
replicate = cell(length(user.silacratios),1);

for ii = 1:length(user.silacratios)
  Chromatograms{ii} = nan(Nproteins*user.Nreplicate,user.Nfraction);
  replicate{ii} = nan(Nproteins*user.Nreplicate, 1);
  Proteins{ii} = cell(Nproteins,1);
  
  cc = 0; % row counter in Chromatograms
  
  for mm = 1:Ninteractions
    
    compare_multiplier = 1;
    if mm <= Ndifferentcomparisons && ii == 1
      compare_multiplier = 0.5;
    end
    
    % "A" protein
    cc = cc+1;
    for jj = 1:user.Nreplicate
      kk = cc+(jj-1)*Nproteins;
      Proteins{ii}{kk} = ['A' num2str(mm)];
      replicate{ii}(kk) = jj;
      Chromatograms{ii}(kk,:) = chrom(mm,:);
    end
    
    % "B" protein
    cc = cc+1;
    for jj = 1:user.Nreplicate
      kk = cc+(jj-1)*Nproteins;
      Proteins{ii}{kk} = ['B' num2str(mm)];
      replicate{ii}(kk) = jj;
      Chromatograms{ii}(kk,:) = chrom(mm,:);
    end
  end
  
  % "Z" proteins, i.e. all the rest
  for nn = cc+1:Nproteins
    for jj = 1:user.Nreplicate
      kk = nn+(jj-1)*Nproteins;
      Proteins{ii}{kk} = ['Z' num2str(ii)];
      replicate{ii}(kk) = jj;
    end
  end
  
end



%% Dirty up the chromatograms

for ii = 1:length(user.silacratios)
  for cc = 1:size(Chromatograms{ii},1)
    
    % add noise
    Chromatograms{ii}(cc,:) = Chromatograms{ii}(cc,:) + rand(size(Chromatograms{ii}(cc,:)))*max(Chromatograms{ii}(cc,:))*.1;
    
    % add up to 15% nans
    Nnan = floor(rand * user.Nfraction * 0.15);
    I = randsample(user.Nfraction,Nnan);
    Chromatograms{ii}(I) = nan;
    
  end
end



%% Comparison: change some chromatograms in channel 1

ii = 1; % channel

cc = 0; % chromatogram counter
for mm = 1:Ninteractions
  
  compare_multiplier = 1;
  if mm <= Ndifferentcomparisons
    compare_multiplier = 0.5;
  end
  
  % "A" protein
  cc = cc+1;
  for jj = 1:user.Nreplicate
    kk = cc+(jj-1)*Nproteins;
    Chromatograms{ii}(kk,:) = Chromatograms{ii}(kk,:) * compare_multiplier;
  end
  
  % "B" protein
  cc = cc+1;
  for jj = 1:user.Nreplicate
    kk = cc+(jj-1)*Nproteins;
    Chromatograms{ii}(kk,:) = Chromatograms{ii}(kk,:) * compare_multiplier;
  end
  
end



%% Quick visualization

figure
for ii = 1:length(Chromatograms)
  subplot(length(Chromatograms),1,ii)
  imagesc(Chromatograms{ii})
  pause(.001)
end



%% Write fake data

% Chromatogram file(s)
for jj = 1:length(user.MQfiles)
  chromfile = user.MQfiles{jj};
  chromid = fopen(chromfile,'wt');
  % Write header
  fprintf(chromid,'%s,%s, ', 'Major protein group', 'Replicate');
  for ii = 1:user.Nfraction
    s = ['Ratio M/L ' num2str(ii)];
    fprintf(chromid,'%s, ', s);
  end
  fprintf(chromid,'\n');
  for ii = 1:size(Chromatograms{jj},1)
    fprintf(chromid,'%s, %6.4f,', Proteins{jj}{ii}, replicate{jj}(ii));
    fprintf(chromid,'%6.4g,',Chromatograms{jj}(ii,:)); %Chromatogram information
    fprintf(chromid,'\n');
  end
  fclose(chromid);
end

% Write the interactions in Corum pairwise file
%   - skip one TP
%   - add one FP
corumfile = user.corumpairwisefile;
corumfile = fopen(corumfile,'wt');
for ii = 1:round(Ninteractions) - 1
  sA = ['A' num2str(ii)];
  sB = ['B' num2str(ii)];
  fprintf(corumfile,'%s,%s, \n', sA, sB);
end
% Write one FP
sA = 'A1';
sB = 'A2';
fprintf(corumfile,'%s,%s, \n', sA, sB);
% Write 100 junk/filler interactions
 symbols = ['a':'z' 'A':'Z' '0':'9'];
for ii = 1:100
  s1 = randsample(symbols,8);
  s2 = randsample(symbols,8);
  fprintf(corumfile,'%s,%s, \n', s1, s2);
end
fclose(corumfile);

% Chromatogram file(s)
mpgfile = user.majorproteingroupsfile;
mpgid = fopen(mpgfile,'wt');
% Write header
fprintf(mpgid,'%s,\n', 'Majority protein IDs');
for ii = 1:Nproteins
  fprintf(mpgid,'%s,\n', Proteins{1}{ii});
end
fclose(mpgid);



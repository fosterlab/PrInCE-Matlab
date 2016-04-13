


%% 0. Initialize
tic
fprintf('\n    0. Initialize')

% Define folders, i.e. define where everything lives.
datadir = [user.maindir 'Data/']; % where data files live
% Make folders if necessary
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist([datadir '/Enrichment/'], 'dir'); mkdir([datadir '/Enrichment/']); end

% #### HACK. FIX THIS ####
clear ComplexesIn
% ComplexesIn{1} = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/Data/Complexes/Final_complexes_rep1.csv';
% ComplexesIn{2} = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/Data/Complexes/Final_complexes_rep2.csv';
% ComplexesIn{3} = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/Data/Complexes/Final_complexes_rep3.csv';
% ComplexesIn{4} = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/Data/Complexes/Final_complexes_rep4.csv';
% ComplexesIn{5} = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/Data/Complexes/Final_complexes_rep5.csv';
% ComplexesIn{6} = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/Data/Complexes/Final_complexes_rep6.csv';
ComplexesIn{1} = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/Data/Complexes/Final_complexes_precision47.csv';

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)




%% 1. Read input data
tic
fprintf('\n    1. Read input data')

predComplexes = cell(size(ComplexesIn));
uniqueProteins = [];
Ncomplex = zeros(size(ComplexesIn));
for ii = 1:length(ComplexesIn)
  data(ii) = importdata(ComplexesIn{ii});
  
  % sanity check
  if ~isequal(data(ii).textdata{1,2},'Predicted complex')
    error('Enrichment: Final_complexes file badly formatted, 1')
  end
  if size(data(ii).data,2)~=2
    error('Enrichment: Final_complexes file badly formatted, 2')
  end
  if (length(data(ii).textdata) - length(data(ii).data)) == 1
  elseif (length(data(ii).textdata) - length(data(ii).data)) > 1
    nn0 = length(data(ii).data);
    nn = (length(data(ii).textdata) - nn0) - 1;
    for jj = nn0+1 : nn0+nn
      data(ii).data(jj,:) = nan;
    end
  else
    error('Enrichment: Final_complexes file badly formatted, 3')
  end
  
  % make predComplexes, uniqueProteins
  tmp2 = data(ii).textdata(2:end,2);
  predComplexes{ii} = cell(size(tmp2));
  for jj = 1:length(tmp2)
    predComplexes{ii}{jj} = strsplit(tmp2{jj},' ');
    uniqueProteins = unique([predComplexes{ii}{jj} uniqueProteins]);
  end
  Ncomplex(ii) = length(tmp2);
  
end

Nprot = length(uniqueProteins);


% Make predComplexes2, indices of uniqueProteins
predComplexes2 = cell(size(predComplexes));
for ii = 1:length(ComplexesIn)
  predComplexes2{ii} = cell(size(predComplexes{ii}));
  for jj = 1:Ncomplex(ii)
    predComplexes2{ii}{jj} = find(ismember(uniqueProteins,predComplexes{ii}{jj}));
  end
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)




%% 2. Attach disease labels to proteins
% Find each protein in the fasta file, and obtain the gene name.
% Find that gene name in the omim file, and attach the disease label.

tic
fprintf('    2. Attach disease labels to proteins')

% load whole fastafile to variable
fastatext = fileread(user.fastafile);
Inewline = find(fastatext == 13 | fastatext == 10); % find all newline characters in the fastafile
Ign = strfind(fastatext, 'GN='); % find all Gene Name declarations in the fastafile

% load whole omimfile to variable
omimtext = fileread(user.omimfile);
inewline = find(omimtext == 13 | omimtext == 10); % find all newline characters in the fastafile

% DLsummvect will will summarize label attaching
% 1 = Disease Label attached
% 2 = protein not found in fasta file
% 3 = protein in fasta, but no gene name
% 4 = protein in fasta, gene name in fasta, but gene name not in omim
DLsummvect = nan(Nprot,1);

geneName = cell(size(uniqueProteins));
diseaseLabel = zeros(size(uniqueProteins,1),10);
for ii = 1:Nprot
  
  protName = uniqueProteins{ii};
  s = ['|' protName '|'];
  
  % look for protName in the fasta file
  Iprot = strfind(fastatext, s);
  
  % Search the fasta file for this protein's gene name
  % find the Ign between Iprot and Inewline
  if isempty(Iprot)
    DLsummvect(ii) = 2;
    continue;
  end
  Inextline = Inewline(find(Inewline>Iprot,1,'first'));
  Istart = Ign(Ign>Iprot & Ign<Inextline); % Gene Name occurring immediately after protName
  if isempty(Istart)
    DLsummvect(ii) = 3;
    continue;
  end
  Iend1 = strfind(fastatext(Istart:Istart+20),' '); % look for white space
  Iend2 = find(fastatext(Istart:Istart+20) == 13 | fastatext(Istart:Istart+20) == 10); % look for newline
  Iend = min([Iend1 Iend2]) + Istart -2;
  geneName{ii} = fastatext(Istart+3:Iend);
  
  % Search the omim file for this gene's disease label(s)
  s2 = [char(9) geneName{ii} char(9)];
  igenename = strfind(lower(omimtext),lower(s2));
  if isempty(igenename)
    DLsummvect(ii) = 4;
    continue;
  end
  for jj = 1:length(igenename)
    iprevline = inewline(find(inewline<igenename(jj),1,'last'));
    iend = find(omimtext(iprevline:igenename(jj))==9) + iprevline;
    iend = iend(1);
    diseaseLabel(ii,jj) = str2double(omimtext(iprevline+1:iend-2));
  end
  
  DLsummvect(ii) = 1;
end
uniqueDiseases = unique(diseaseLabel(:));
uniqueDiseases(uniqueDiseases==0) = [];
Ndis = length(uniqueDiseases);

% Make truncated disease labels
DL = cell(5,1);
for ii = 1:5
  DL{ii} = zeros(Nprot,1);
  for jj = 1:Nprot
    if diseaseLabel(jj,1)==0
      continue;
    end
    tmp = num2str(diseaseLabel(jj,1));
    DL{ii}(jj) = str2double(tmp(1:ii));
  end
  UD{ii} = unique(DL{ii});
  UD{ii}(UD{ii}==0) = [];
  ND(ii) = length(UD{ii});
end
DL{6} = diseaseLabel; % all digits!
UD{6} = uniqueDiseases;
ND(6) = Ndis;

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 3. Run enrichment analysis
%       i) Look at every Complex-Disease pair.
%      ii) For each pair, make a contingency table T.
%     iii) Calculate p(ii) = fisherexact(T).
%      iv) Run B-H correction on all p.

tic
fprintf('    3. Run enrichment analysis')

dd = 6; % digits in omim disease labels

p_enrich = cell(size(ComplexesIn));
pcutoff = nan(size(ComplexesIn));
Enriched = cell(size(ComplexesIn));
for cc = 1:length(ComplexesIn) % precision/replicate index
  cc
  
  % Fisher exact test
  p_enrich{cc} = nan(Ncomplex(cc),ND(dd));
  tmax = zeros(Ncomplex(cc),ND(dd));
  for ii = 1:Ncomplex(cc)
    ii
    I = predComplexes2{cc}{ii};  % index of proteins in this complex
    I0 = ~ismember(1:Nprot,predComplexes2{cc}{ii});  % index of proteins NOT in this cluster
    
    for jj = 1:ND(dd)
      dl = UD{dd}(jj); % current disease label
      dlcomp = DL{dd}(I,:); % all disease labels for this complex
      dlcomp0 = DL{dd}(I0,:); % all disease labels for proteins NOT in this complex
      
      % Make contingency table
      a = sum(dlcomp(:) == dl & dlcomp(:)~=0);  % in cluster       AND disease label
      b = sum(dlcomp(:) ~= dl & dlcomp(:)~=0);  % in cluster       AND NOT disease label
      c = sum(dlcomp0(:)== dl & dlcomp0(:)~=0); % NOT in cluster   AND disease label
      d = sum(dlcomp0(:)~= dl & dlcomp0(:)~=0); % NOT in cluster   AND NOT disease label
      T = [a b; c d];
      tmax(ii,jj) = max([a c]);
      [~,p_enrich{cc}(ii,jj)] = fishertest(T);
      
    end
  end
  
  % Benjamini-Hochberg
  p_asc = sort(p_enrich{cc}(:),'ascend');
  if size(p_asc,1)>size(p_asc,2)
    p_asc = p_asc';
  end
  m = length(p_asc);
  I = find( (1:m)/m*user.fdr - p_asc >= 0, 1, 'last');
  pcutoff(cc) = p_asc(I);
  
  % Format for output
  Ienriched = find(p_enrich{cc}(:)<pcutoff(cc));
  [I1,I2] = ind2sub(size(p_enrich{cc}),Ienriched);
  Enriched{cc} = cell(Ncomplex(cc),2);
  for ii = 1:length(I1)
    if isempty(Enriched{cc}{I1(ii),1})
      Enriched{cc}{I1(ii),1} = num2str(uniqueDiseases(I2(ii)));
    else
      Enriched{cc}{I1(ii),1} = strjoin({Enriched{cc}{I1(ii),1} num2str(uniqueDiseases(I2(ii)))}, ' ');
    end
    if isempty(Enriched{cc}{I1(ii),2})
      Enriched{cc}{I1(ii),2} = num2str(p_enrich{cc}(I1(ii),I2(ii)));
    else
      Enriched{cc}{I1(ii),2} = strjoin({Enriched{cc}{I1(ii),2} num2str(p_enrich{cc}(I1(ii),I2(ii)))}, '  ');
    end
  end
  
  figure
  imagesc(p_enrich{cc}<=pcutoff(cc))
end


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)





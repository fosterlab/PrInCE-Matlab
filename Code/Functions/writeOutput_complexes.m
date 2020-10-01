%WRITEOUTPUT_COMPLEXES Writes output tables for the PRINCE Complexes module.

Precision_values = round(desiredPrecision * 100);

Names = cell(size(csplit,1),1);
for ii = 1:size(csplit,1)
  if csplit(ii,1)==0 && csplit(ii,2)==0
    name = 'All_interactions_';
  else
    if csplit(ii,1)==0
      name1 = '';
    else
      rep = num2str(csplit(ii,1));
      name1 = ['Replicate_' rep '_'];
    end
    if csplit(ii,2)==0
      name2 = '';
    else
      rep = user.silacratios{csplit(ii,2)};
      name2 = ['Channel_' num2str(rep) '_'];
    end
    name = [name1  name2];
  end
  name = name(1:end-1);
  Names{ii} = name;
end


%% Complexes_*.csv

for ci = 1:size(csplit,1)
  fn = strcat([datadir 'Complexes_' Names{ci} '.csv']);
  
  fid3 = fopen(fn,'w');
  %Write Header
  fprintf (fid3,'%s,%s,%s,%s,%s,%s,%s,%s,%s,\n',...
    'ID','Predicted complex', 'Complex size', 'Complex density', 'Novel complex?', ...
    'Best CORUM match', 'Overlapping proteins', 'N overlap', 'CORUM coverage');
  
  for ii = 1:size(CL(ci).Members,1)    
    predComplex = strjoin(uniqueProteins(CL(ci).Members{ii}), ' ');
    sizePredComplex = length(CL(ci).Members{ii});

    Icorummatch = find(corumMatches{ci}(:,1)==ii & corumMatches{ci}(:,3)==1);
    corumComplex = '';
    overlap = '';
    novel = 1;
    if ~isempty(Icorummatch)
      novel = 0;
      Icorummatch = Icorummatch(1);
      corumI = corumMatches{ci}(Icorummatch,2);
      corumComplex = strjoin(uniqueProteins(corumComplex2{corumI}), ' ');
      overlap = strjoin(uniqueProteins(intersect(CL(ci).Members{ii},corumComplex2{corumI})));
    end
    
    fprintf (fid3,'%d, %s, %d, %6.3f, %d, %s, %s, %d, %6.3f,\n',...
      ii, predComplex, sizePredComplex, CL(ci).Density(ii), novel, corumComplex, overlap, ...
      corumMatches{ci}(Icorummatch,4), corumMatches{ci}(Icorummatch,5) );
  end
  fclose(fid3);
end



%% Corum_best_match_*prec.csv

for ci = 1:size(csplit,1)
  fn = strcat([datadir 'More_corum_matches_' Names{ci} '.csv']);
  
  fid3 = fopen(fn,'w');
  %Write Header
  fprintf (fid3,'%s,%s,%s,%s,%s,%s,%s\n',...
    'Predicted complex number','Predicted complex', 'CORUM match',...
    'Overlapping proteins', 'Best or 2nd best match?', 'N overlap', 'CORUM coverage');
  
  for ii = 1:size(corumMatches{ci},1)
    cmplxI = corumMatches{ci}(ii,1);
    corumI = corumMatches{ci}(ii,2);
    predComplex = strjoin(uniqueProteins(CL(ci).Members{cmplxI}), ' ');
    corumComplex = strjoin(uniqueProteins(corumComplex2{corumI}), ' ');
    overlap = strjoin(uniqueProteins(intersect(CL(ci).Members{cmplxI},corumComplex2{corumI})));
    if corumMatches{ci}(ii,3)==1
      s = 'Best';
    elseif corumMatches{ci}(ii,3)==2
      s = 'Second best';
    else
      s = 'error100';
    end
    
    % sanity check
    if length(intersect(CL(ci).Members{cmplxI},corumComplex2{corumI})) ~= corumMatches{ci}(ii,4)
      error('writeOutput_complexes: badly formatted corumMatches')
    end
    
    fprintf (fid3,'%d,%s,%s, %s, %s,%d,%6.3f,\n',...
      cmplxI, predComplex, corumComplex, overlap, s,...
      corumMatches{ci}(ii,4), corumMatches{ci}(ii,5) );
  end
  fclose(fid3);
end




%% Corum_best_match_*prec.csv

fn = strcat([datadir 'Summary_complexes.csv']);

fid3 = fopen(fn,'w');
%Write Header
fprintf (fid3,'%s,%s,%s,%s,%s,%s,\n',...
  'Name', 'Number of complexes', 'Median complex size', 'Avg complex density',...
  'Geometric accuracy', 'Matching ratio');

for ii = 1:size(csplit,1)
  
  name = Names{ii};
  Nmembers = length(CL(ii).Members);
  tmp = nan(Nmembers,1);
  for jj = 1:Nmembers
    tmp(jj) = length(CL(ii).Members{jj});
  end
  medsize = median(tmp);
  avgdens = mean(CL(ii).Density);
  
  fprintf (fid3,'%s,%d,%d,%6.3f,%6.3f,%6.3f,\n',...
    name, Nmembers, medsize, avgdens, ...
    CL(ii).GeomAcc,CL(ii).MatchRatio);
end
fclose(fid3);


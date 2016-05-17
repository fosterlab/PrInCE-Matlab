
Precision_values = round(user.desiredPrecision * 100);


%% Final_complexes_precision*.csv

for ci = 1:countPrec
  fn = strcat([datadir 'Complexes/Final_complexes_precision' num2str(Precision_values(ci)) '.csv']);
  
  fid3 = fopen(fn,'w');
  %Write Header
  fprintf (fid3,'%s,%s,%s,%s,%s,%s,%s,\n',...
    'Predicted complex number','Predicted complex', 'Novel complex?', 'Best CORUM match',...
    'Overlapping proteins', 'N overlap', 'CORUM coverage');
  
  for ii = 1:size(ComplexList(ci).Members,1)    
    predComplex = strjoin(uniqueProteins(ComplexList(ci).Members{ii}), ' ');

    Icorummatch = find(corumMatches{ci}(:,1)==ii & corumMatches{ci}(:,3)==1);
    corumComplex = '';
    overlap = '';
    novel = 1;
    if ~isempty(Icorummatch)
      novel = 0;
      Icorummatch = Icorummatch(1);
      corumI = corumMatches{ci}(Icorummatch,2);
      corumComplex = strjoin(uniqueProteins(corumComplex2{corumI}), ' ');
      overlap = strjoin(uniqueProteins(intersect(ComplexList(ci).Members{ii},corumComplex2{corumI})));
    end
    
    fprintf (fid3,'%d, %s, %d, %s, %s, %d, %6.3f,\n',...
      ii, predComplex, novel, corumComplex, overlap, ...
      corumMatches{ci}(Icorummatch,4), corumMatches{ci}(Icorummatch,5) );
  end
  fclose(fid3);
end


%% Corum_best_match_*prec.csv

for ci = 1:countPrec
  fn = strcat([datadir 'Complexes/Corum_best_match_precision' num2str(Precision_values(ci)) '.csv']);
  
  fid3 = fopen(fn,'w');
  %Write Header
  fprintf (fid3,'%s,%s,%s,%s,%s,%s,%s\n',...
    'Predicted complex number','Predicted complex', 'CORUM match',...
    'Overlapping proteins', 'Best or 2nd best match?', 'N overlap', 'CORUM coverage');
  
  for ii = 1:size(corumMatches{ci},1)
    cmplxI = corumMatches{ci}(ii,1);
    corumI = corumMatches{ci}(ii,2);
    predComplex = strjoin(uniqueProteins(ComplexList(ci).Members{cmplxI}), ' ');
    corumComplex = strjoin(uniqueProteins(corumComplex2{corumI}), ' ');
    overlap = strjoin(uniqueProteins(intersect(ComplexList(ci).Members{cmplxI},corumComplex2{corumI})));
    if corumMatches{ci}(ii,3)==1
      s = 'Best';
    elseif corumMatches{ci}(ii,3)==2
      s = 'Second best';
    else
      s = 'error100';
    end
    
    % sanity check
    if length(intersect(ComplexList(ci).Members{cmplxI},corumComplex2{corumI})) ~= corumMatches{ci}(ii,4)
      error('writeOutput_complexes: badly formatted corumMatches')
    end
    
    fprintf (fid3,'%d,%s,%s, %s, %s,%d,%6.3f,\n',...
      cmplxI, predComplex, corumComplex, overlap, s,...
      corumMatches{ci}(ii,4), corumMatches{ci}(ii,5) );
  end
  fclose(fid3);
end
%WRITEOUTPUT_FOLDCHANGES Writes output tables for the PRINCE FoldChanges
%   module

fn = [datadir '/Fold_change_atElutionPeaks_' user.comparisonpairs{1} '_vs_' user.comparisonpairs{2} '.csv'];
fid = fopen(fn,'w');
fprintf (fid,'%s,%s,%s,%s,%s\n','Protein ID','Replicate','Fraction (Gauss center)','log2 Fold change','log2 Fold change (normalized)');
for ii = 1:length(Finalised_Master_Gaussian_list.Protein_name)
  prot = Finalised_Master_Gaussian_list.Protein_name{ii,1};
  if isempty(prot); continue; end
  I = find(ismember(Combined_Gaussians.Protein_name,prot));
  reps = Combined_Gaussians.Replicate(I);
  centers = Combined_Gaussians.Center(I);
  foldChange = Combined_Gaussians.log2_of_gaussians(I);
  foldChange_norm = Combined_Gaussians.log2_normalised_gaussians(I);
  %reps = Finalised_Master_Gaussian_list.Replicate(ii,:);
  %centers = Finalised_Master_Gaussian_list.Center(ii,:);
  %foldChange = Finalised_Master_Gaussian_list.foldChange(ii,:);
  %foldChange_norm = Finalised_Master_Gaussian_list.foldChange_normalized(ii,:);
  for jj = 1:length(I)
    fprintf(fid,'%s,%d,%6.4f,%6.4f,%6.4f\n',...
      prot,reps(jj), centers(jj), foldChange(jj), foldChange_norm(jj));
  end
end


try
  % Average fold change over all fractions
  fn = [datadir '/Fold_change_Average_' user.comparisonpairs{1} '_vs_' user.comparisonpairs{2} '.csv'];
  fid = fopen(fn,'w');
  fprintf(fid,'Protein ID,Average fold change,');
  reps = unique(Combined_Gaussians.Replicate);
  for ii = 1:length(reps)
    fprintf(fid,'%s,',['Replicate ' num2str(reps(ii))]);
  end
  fprintf(fid,'P-value\n');
  for ii = 1:length(Finalised_Master_Gaussian_list.Protein_name)
    prot = Finalised_Master_Gaussian_list.Protein_name{ii,1};
    
    % Hack hack hack hack hack hack hack hack hack hack hack hack hack
    % This should be in FoldChanges.m
    I = find(ismember(txt_val{1}(:,2),prot));
    fc = nan(size(I));
    for jj = 1:length(fc)
      tmp = num_val{2}(I(jj),2:end) ./ num_val{1}(I(jj),2:end);
      fc(jj) = log2(nanmean(tmp));
    end
    [x1,pp] = ttest(fc);
    % Hack hack hack hack hack hack hack hack hack hack hack hack hack
    
    fprintf(fid,'%s,%6.4f,',prot,nanmean(fc));
    for jj = 1:length(fc)
      fprintf(fid,'%6.4f,',fc(jj));
    end
    fprintf(fid,'%6.4f\n',pp);
  end
  fclose all;
end


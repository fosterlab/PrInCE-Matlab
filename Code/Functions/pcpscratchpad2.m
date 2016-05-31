% May 31
%
% Think I just fixed THE BUG, and that it had to do with my Corum files + TP_Matrix error.

%% Once I reduce to just A-B interactions, are my corum files the same as Nick's?

% Apoptosis data
ff{1} = '/Users/Mercy/Academics/Foster/NickCodeData/Old runs/GregPCP_20160517/Data/Corum_correctly_formated_Uniprot_IDs.csv';
ff{2} = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Output/tmp/Corum_pairwise.csv';
corList_apt = cell(2,1);
for ii = 1:2
  fid=fopen(ff{ii}, 'rt');    %Input corum data base.
  Corum_Import= textscan(fid,'%s\t', 'Delimiter',',');
  fclose(fid);
  No=length(Corum_Import{1})/2;
  tmp = reshape(Corum_Import{1,1},2,No)';
  corList_apt{ii} = cell(size(tmp,1),1);
  for jj = 1:size(tmp,1)
    tmp(jj,:) = sort(tmp(jj,:));
    corList_apt{ii}{jj} = strjoin(tmp(jj,:),'_');
  end
  
  % reduce to unique
  %[~,idx]=unique(cell2mat(corList_apt{ii}),'rows');
  %corList_apt{ii} =  corList_apt{ii}(idx,:);
  corList_apt{ii} = unique(corList_apt{ii});
end
I = intersect(corList_apt{1},corList_apt{2});
figure,hold on
myVenn2([length(corList_apt{1}) length(corList_apt{2})],length(I))

%%
% Tissue data
ff{1} = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/Input/Mapped_mouse_Corum_list_20150109.csv';
ff{2} = '/Users/Mercy/Academics/Foster/Tissue_PCPSILAC/PCPSILAC_analysis/Output/tmp/Corum_pairwise.csv';
corList_tis = cell(2,1);
for ii = 1:2
  fid=fopen(ff{ii}, 'rt');    %Input corum data base.
  Corum_Import= textscan(fid,'%s\t', 'Delimiter',',');
  fclose(fid);
  No=length(Corum_Import{1})/2;
  tmp = reshape(Corum_Import{1,1},2,No)';
  corList_tis{ii} = cell(size(tmp,1),1);
  for jj = 1:size(tmp,1)
    tmp(jj,:) = sort(tmp(jj,:));
    corList_tis{ii}{jj} = strjoin(tmp(jj,:),'_');
  end
  
  % reduce to unique
  %[~,idx]=unique(cell2mat(corList_apt{ii}),'rows');
  %corList_apt{ii} =  corList_apt{ii}(idx,:);
  corList_tis{ii} = unique(corList_tis{ii});
end
I = intersect(corList_tis{1},corList_tis{2});
figure,hold on
myVenn2([length(corList_tis{1}) length(corList_tis{2})],length(I))



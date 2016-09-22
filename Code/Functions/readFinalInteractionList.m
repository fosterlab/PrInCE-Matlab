function [interaction_list, head] = readFinalInteractionList(fn)

% Quick function for reading Final_Interactions_list_... produced by ROC.
% NB, make sure this function exactly corresponds with writeOutput_ROC


column_names1 = {'Unique interactions','Protein A','Protein B','Center A',...
  'Center B','Replicate','Channel','Delta Center','R^2',...
  'Euc. Dist', 'Both proteins in Corum', 'Interaction in Corum','Interaction score',...
  'Score rank (low = better)','Score percentile (low = better)','Precision level',...
  'Global?'};

column_names2 = {'Unique interactions','Protein A', 'Protein B','Replicate',...
  'Channel','Both proteins in Corum','Interaction in Corum','Interaction score',...
  'Score rank (low = better)','Score percentile (low = better)','Precision level'};

column_names3 = {'Unique interactions','Protein A', 'Protein B',...
  'Channel','Both proteins in Corum','Interaction in Corum','Interaction score',...
  'Score rank (low = better)','Score percentile (low = better)','Precision level'};


fid = fopen(fn);

header = fgetl(fid);
header = strsplit(header,',');
header(cellfun('isempty',header)) = [];

if isequal(header,column_names1)
  Idata = [11 12 13 14 15 16 17];
  Itext = [1 2 3 4 5 6 7 8 9 10];
elseif isequal(header,column_names2)
  Idata = [6 7 8 9 10 11];
  Itext = [1 2 3 4 5];
elseif isequal(header,column_names3)
  Idata = [5 6 7 8 9 10];
  Itext = [1 2 3 4];
else
  error('Incorrectly formatted Final_Interactions_list')
end
head.data = header(Idata);
head.text = header(Itext);

interaction_list.data = nan(200000,length(Idata));
interaction_list.text = cell(200000,length(Itext));
cc = 0;
while ~feof(fid)
  
  current_line = fgetl(fid);
  current_line = strsplit(current_line,',');
  cc = cc+1;
  
  for ii = 1:length(Idata)
    interaction_list.data(cc,ii) = str2double(current_line{Idata(ii)});
  end
  for ii = 1:length(Itext)
    interaction_list.text{cc,ii} = current_line{Itext(ii)};
  end
end

interaction_list.data = interaction_list.data(1:cc,:);
interaction_list.text = interaction_list.text(1:cc,:);

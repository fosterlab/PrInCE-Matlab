function interaction_list = readFinalInteractionList(fn)

% Quick function for reading Final_Interactions_list_... produced by ROC.
% NB, make sure this function exactly corresponds with writeOutput_ROC


column_names = {'Unique interactions','Protein A','Protein B','Center A',...
  'Center B','Replicate','Delta Center','R^2',...
  'Euc. Dist', 'Both proteins in Corum', 'Interaction in Corum','Interaction score',...
  'Score rank (low = better)','Score percentile (low = better)','Precision level',...
  'Global?'};


fid = fopen(fn);

head = fgetl(fid);
head = strsplit(head,',');
head(cellfun('isempty',head)) = [];

% sanity check
if ~isequal(head,column_names)
  error('Incorrectly formatted Final_Interactions_list')
end

interaction_list = cell(100000,length(column_names));
cc = 0;
while ~feof(fid)
  
  current_line = fgetl(fid);
  current_line = strsplit(current_line,',');
  current_line(cellfun('isempty',current_line)) = [];
  
  cc = cc+1;
  
  interaction_list(cc,:) = current_line;
end

interaction_list = interaction_list(1:cc,:);
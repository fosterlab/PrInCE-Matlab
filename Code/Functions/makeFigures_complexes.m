%MAKEFIGURES_COMPLEXES Makes figures for the PRINCE Complexes module.

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



%% Colour scheme:
% black - reference only
% purple - recovered reference
% orange - novel
%
% orange: #f58d62, [245 141 98]
% purple: #8f6ca9, [143 108 169]
% light blue: #5ec5e2, [94 197 226]
% dark blue: #1482b2, [20 130 178]


%% Single complex, corum / predicted overlap

h = figure;
for ii = 1:size(csplit,1)
  for jj = 1:length(CL(ii).Members)
    
    predComplex_members = CL(ii).Members{jj};
    predComplex_connections = CL(ii).Connections{jj};
    if isempty(predComplex_members) | isempty(predComplex_connections)
      continue;
    end
    
    % find closest reference complex
    refComplex = [];
    sizeCorum = length(refComplex);
    if not(isempty(corumMatches{ii}))
      I = find(corumMatches{ii}(:,1)==jj & corumMatches{ii}(:,3)==1);
      if ~isempty(I)
        refComplex = corumComplex2{corumMatches{ii}(I,2)};
      end
      sizeCorum = length(refComplex);
    end
    
    set(0, 'CurrentFigure', h);
    clf reset;
    
    p = plotcomplex(predComplex_members,predComplex_connections,refComplex,uniqueProteins);
    p.MarkerSize = 7;
    
    sizePred = length(CL(ii).Members{jj});
    if sizeCorum+sizePred < 20
      set(gcf,'paperunits','inches','paperposition',[.25 .25 3 3],...
        'units','inches','position',[.25 .25 8 8])
    elseif sizeCorum+sizePred > 20 & sizeCorum+sizePred < 40
      set(gcf,'paperunits','inches','paperposition',[.25 .25 5 5],...
        'units','inches','position',[.25 .25 8 8])
    else
      set(gcf,'paperunits','inches','paperposition',[.25 .25 7 7],...
        'units','inches','position',[.25 .25 15 15])
    end
    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gca,'Visible','off')
    axis square
    sf=[figdir '/Predicted_vs_corum_' Names{ii} '_Complex_' num2str(jj)];
    saveas(gcf, sf, 'svg');
    saveas(gcf, sf, 'png');
    
  end
end



%% Hairball figure 2

ii = 1;

% make a square connection matrix
connectionMatrix = zeros(5000,5000);
conMatrixProteins = nan(5000,1);
Ipred_all = zeros(5000,1);
Iref_all = zeros(5000,1);

% fill in with predicted complexes
cc = 0;
for jj = 1:length(CL(ii).Members)
  predComplex_members = CL(ii).Members{jj};
  predComplex_connections = CL(ii).Connections{jj};
  
  % find closest reference complex
  refComplex = [];
  if not(isempty(corumMatches{ii}))
    I = find(corumMatches{ii}(:,1)==jj & corumMatches{ii}(:,3)==1);
    if ~isempty(I)
      refComplex = corumComplex2{corumMatches{ii}(I,2)};
    end
  end
  
  allProteins = unique([refComplex predComplex_members]);
  Nproteins1 = length(allProteins);
  I = cc+1 : cc+Nproteins1;
  
  Ipred = (ismember(allProteins,predComplex_members)); % indices of predicted complex members
  Iref = (ismember(allProteins,refComplex));% indices of reference complex members
  %Iref_only = Iref(~ismember(Iref,Ipred));
  %Ipred_only = Ipred(~ismember(Ipred,Iref));
  
  cm = zeros(Nproteins1,Nproteins1);
  cm(Ipred,Ipred) = predComplex_connections;
  cm(Iref,Iref) = 1;
  
  Ipred_all(I) = Ipred;
  Iref_all(I) = Iref;
  connectionMatrix(I,I) = cm;
  conMatrixProteins(I) = allProteins;
  cc = cc+Nproteins1;
end
connectionMatrix = connectionMatrix(1:cc,1:cc);
conMatrixProteins = conMatrixProteins(1:cc);

Ipred_all = Ipred_all(1:cc);
Iref_all = Iref_all(1:cc);
Iref_only = Iref_all==1 & Ipred_all==0;
Ipred_only = Iref_all==0 & Ipred_all==1;

% plot the network
G = graph(connectionMatrix,'OmitSelfLoops');
LWidths = G.Edges.Weight;
LWidths = (LWidths - min(LWidths) + 0.1).^2;
LWidths = sqrt(LWidths / max(LWidths)) * .5;

if not(isempty(G.Edges))
  figure
  p = plot(G,'Layout','force','LineWidth',LWidths);
  p.MarkerSize = 0.5;
  
  % Highlight Edges
  Edges = table2array(G.Edges);
  % turn prediction-to-anything edges purple
  Ipred_edge = ismember(Edges(:,1),find(Ipred_all)) | ismember(Edges(:,2),find(Ipred_all));
  highlight(p,Edges(Ipred_edge,1),Edges(Ipred_edge,2),'EdgeColor',[143 108 169]/255)
  % turn reference-only-to-anything edges black
  Iref_edge = ismember(Edges(:,1),find(Iref_only)) | ismember(Edges(:,2),find(Iref_only));
  highlight(p,Edges(Iref_edge,1),Edges(Iref_edge,2),'EdgeColor','k')
  % turn prediction-only edges orange
  Ipredonly = ismember(Edges(:,1),find(Ipred_only)) | ismember(Edges(:,2),find(Ipred_only));
  highlight(p,Edges(Ipredonly,1),Edges(Ipredonly,2),'EdgeColor',[245 141 98]/255)
  
  
  % Highlight Nodes
  % turn reference-only nodes black
  highlight(p,find(Iref_only),'NodeColor',[.3 .3 .3])
  % turn prediction-only nodes orange
  highlight(p,find(Ipred_only),'NodeColor',[245 141 98]/255)
  % turn overlapping (reference+prediction) nodes purple
  I = 1:cc;
  Iboth = I(~ismember(I,[find(Iref_only); find(Ipred_only)]));
  highlight(p,Iboth,'NodeColor',[143 108 169]/255)
  
  set(gca,'LooseInset',get(gca,'TightInset'));
  set(gcf,'paperunits','inches','paperposition',[.25 .25 7 7],...
    'units','inches','position',[.25 .25 15 15])
  axis square
  %axis([0.6688   54.5981    1.0965   51.4244])
  %axis([0 0 78 75])
  set(gca,'Visible','off','xtick',[],'ytick',[])
  sf=[figdir '/Hairball2_redicted_vs_corum'];
  saveas(gcf, sf, 'svg');
  print([sf '.png'], '-dpng', '-r1000');
end

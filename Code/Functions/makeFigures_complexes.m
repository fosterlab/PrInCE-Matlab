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

for ii = 1:size(csplit,1)
  for jj = 1:length(CL(ii).Members)
    
    predComplex_members = CL(ii).Members{jj};
    predComplex_connections = CL(ii).Connections{jj};
    
    % find closest reference complex
    I = find(corumMatches{ii}(:,1)==jj & corumMatches{ii}(:,3)==1);
    refComplex = [];
    if ~isempty(I)
      refComplex = corumComplex2{corumMatches{ii}(I,2)};
    end
    sizeCorum = length(refComplex);
    
    figure
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
    saveas(gcf, sf, 'epsc');
    close all
  end
end



%% Hairball figure 1
% 
% % make a square connection matrix
% connectionMatrix = zeros(Nproteins,Nproteins);
% ii = 1;
% 
% % fill in with predicted complexes
% for jj = 1:length(CL(ii).Members)
%   I = CL(ii).Members{jj};
%   connectionMatrix(I,I) = CL(ii).Connections{jj};
% end
% degreePred = sum(connectionMatrix>0);
% Ipred = find(degreePred>0);
% 
% % fill in with reference complexes
% for jj = 1:length(corumMatches{ii});
%   nc = corumMatches{ii}(jj,2);
%   I = corumComplex2{nc};
%   connectionMatrix(I,I) = 1;
% end
% 
% % reduce matrix to only proteins with degree>0
% degree = sum(connectionMatrix>0);
% Igood = find(degree>0);
% connectionMatrix = connectionMatrix(Igood,Igood);
% Ipred = find(ismember(Igood,Ipred));
% 
% % plot the network
% G = graph(connectionMatrix,'OmitSelfLoops');
% LWidths = G.Edges.Weight;
% LWidths = (LWidths - min(LWidths) + 0.1).^2;
% LWidths = LWidths / max(LWidths) * 2;
% 
% % Which edges are reference-only, recovered, or novel?
% Edges = table2array(G.Edges);
% Iref = [];
% Iref_edge = zeros(size(Edges(:,1)));
% Ipred_edge = zeros(size(Edges(:,1)));
% Irecovered_edge = zeros(size(Edges(:,1)));
% Iref_node = [];
% for jj = 1:size(corumMatches{ii},1)  
%   nc = corumMatches{ii}(jj,2);
%   Iref = find(ismember(Igood,corumComplex2{nc}));
%   Iref_node = unique([Iref_node Iref]);
%   
%   Iref_edge = Iref_edge | (ismember(Edges(:,1),Iref) & ismember(Edges(:,2),Iref));
%   Ipred_edge = Ipred_edge | (ismember(Edges(:,1),Ipred) | ismember(Edges(:,2),Ipred));
%   Irecovered_edge = Irecovered_edge | (Iref_edge & Ipred_edge);
% end
% 
% Irefonly = Iref_edge & ~Ipred_edge;
% Ipredonly = ~Iref_edge & Ipred_edge;
% Irecovered = Irecovered_edge;
% 
% 
% 
% figure
% p = plot(G,'Layout','force','LineWidth',LWidths);
% 
% % Highlight Edges
% % turn reference-to-reference edges black
% highlight(p,Edges(Iref_edge,1),Edges(Iref_edge,2),'EdgeColor','k')
% 
% % turn prediction-to-anything edges purple
% highlight(p,Edges(Ipred_edge,1),Edges(Ipred_edge,2),'EdgeColor',[143 108 169]/255)
% 
% % turn prediction-only edges orange
% highlight(p,Edges(Ipredonly,1),Edges(Ipredonly,2),'EdgeColor',[245 141 98]/255)
% 
% 
% Iref_only = Iref_node(~ismember(Iref_node,Ipred));
% Ipred_only = Ipred(~ismember(Ipred,Iref_node));
% 
% % Highlight Nodes
% % turn reference-only nodes black
% highlight(p,Iref_only,'NodeColor','k')
% 
% % turn prediction-only nodes orange
% highlight(p,Ipred_only,'NodeColor',[245 141 98]/255)
% 
% % turn overlapping (reference+prediction) nodes purple
% I = 1:length(Igood);
% Iboth = I(~ismember(I,[Iref_only Ipred_only]));
% highlight(p,Iboth,'NodeColor',[143 108 169]/255)
% 
% 
% set(gcf,'paperunits','inches','paperposition',[.25 .25 7 7],...
%   'units','inches','position',[.25 .25 15 15])
% axis square
% set(gca,'Visible','off','xtick',[],'ytick',[])
% sf=[figdir '/Hairball1_redicted_vs_corum'];
% saveas(gcf, sf, 'epsc');



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
  % find closest reference complex
  I = find(corumMatches{ii}(:,1)==jj & corumMatches{ii}(:,3)==1);
  refComplex = [];
  if ~isempty(I)
    refComplex = corumComplex2{corumMatches{ii}(I,2)};
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
axis([0 0 78 75])
set(gca,'Visible','off','xtick',[],'ytick',[])
sf=[figdir '/Hairball2_redicted_vs_corum'];
saveas(gcf, sf, 'epsc');
print([sf '.png'], '-dpng', '-r1000');


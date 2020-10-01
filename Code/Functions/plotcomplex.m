function [p,G] = plotcomplex(predComplex_members,predComplex_connections,refComplex,protNames)
%   [p,G]=PLOTCOMPLEX(Pc,conn,Rc,names) plots a force-directed diagram of
%   the predicted complex Pc and (partially) overlapping reference complex
%   Rc. Pc and Rc are vectors of indices that correspond to the protein
%   names in the cell array names. conn is the NxN weighted connection 
%   matrix of Pc. names is a cell array of protein IDs. The order of
%   elements in names corresponds to indices in Pc and Rc.
%
%   G is the output of GRAPH(X), where X is the square connection matrix
%   conn augmented with nodes and edges from the reference complex. All
%   reference edges are assumed to be of weight 1.
%
%   p is the output of plot(G).
%
%   Colour scheme:
%     Orange    Members of Pc not in Rc (novel)
%     Black     Members of Rc not in Pc (reference-only)
%     Purple    Members in Pc and Rc (recovered)
%
%   See also GRAPH.


% G = digraph(1,2:5);
% G = addedge(G,2,6:15);
% G = addedge(G,15,16:20);
% plot(G,'Layout','force')

% orange: #f58d62, [245 141 98]
% purple: #8f6ca9, [143 108 169]
% light blue: #5ec5e2, [94 197 226]
% dark blue: #1482b2, [20 130 178]
%
% Node colour scheme:
%   black - reference only
%   purpe - reference and predicted
%   orange - predicted only
%
% Edge colour scheme:
%   black - reference-only nodes to reference-only nodes nodes
%   orange - prediction-only nodes to anything
%   light blue - all other edges

if iscell(predComplex_members)
  if isstring(predComplex_members{1})
    
  end
end

if nargin<4
  protNames = [];
end

if ~isvector(predComplex_members) && ~isempty(refComplex)
  error('predComplex_members must be a numerical vector of predicted complex members.')
end
if iscolumn(predComplex_members)
  predComplex_members = predComplex_members';
end

if ~isvector(refComplex) && ~isempty(refComplex)
  error('refComplex must be a numerical vector of predicted complex members.')
end
if iscolumn(refComplex)
  refComplex = refComplex';
end

if size(predComplex_connections,1)~=size(predComplex_connections,2)
  error('predComplex_connections must be a square connection matrix.')
end

if size(predComplex_connections,1)~=length(predComplex_members)
  error('size(predComplex_connections,1) must equal length(predComplex_members).')
end


% Make a square connection matrix by combining predComplex + refComplex
allProteins = unique([refComplex predComplex_members]);
Nproteins = length(allProteins);

Ipred = find(ismember(allProteins,predComplex_members)); % indices of predicted complex members
Iref = find(ismember(allProteins,refComplex));% indices of reference complex members
Iref_only = Iref(~ismember(Iref,Ipred));
Ipred_only = Ipred(~ismember(Ipred,Iref));

connectionMatrix = zeros(Nproteins,Nproteins);
connectionMatrix(Ipred,Ipred) = predComplex_connections;
connectionMatrix(Iref,Iref) = 1;

G = graph(connectionMatrix,'OmitSelfLoops');
LWidths = G.Edges.Weight;
LWidths = (LWidths - min(LWidths) + 0.1).^2;
LWidths = LWidths / max(LWidths) * 4;
if ~isempty(protNames)
  p = plot(G,'Layout','force','LineWidth',LWidths,'NodeLabel',protNames(allProteins));
else
  p = plot(G,'Layout','force','LineWidth',LWidths);
  p.NodeLabel = {};
end
pause(.0001)

% Adjust Edges
Edges = table2array(G.Edges);
% turn prediction-to-anything edges purple
Ipred_edge = ismember(Edges(:,1),Ipred) | ismember(Edges(:,2),Ipred);
highlight(p,Edges(Ipred_edge,1),Edges(Ipred_edge,2),'EdgeColor',[143 108 169]/255)

% turn reference-only-to-anything edges black
Iref_edge = ismember(Edges(:,1),Iref_only) | ismember(Edges(:,2),Iref_only);
highlight(p,Edges(Iref_edge,1),Edges(Iref_edge,2),'EdgeColor','k')

% turn prediction-only-to-prediction-only edges orange
Ipredonly = ismember(Edges(:,1),Ipred_only) | ismember(Edges(:,2),Ipred_only);
highlight(p,Edges(Ipredonly,1),Edges(Ipredonly,2),'EdgeColor',[245 141 98]/255)


% Adjust Nodes
% turn reference-only nodes black
highlight(p,Iref_only,'NodeColor','k')

% turn prediction-only nodes orange
highlight(p,Ipred_only,'NodeColor',[245 141 98]/255)

% turn overlapping (reference+prediction) nodes purple
I = 1:Nproteins;
Iboth = I(~ismember(I,[Iref_only Ipred_only]));
highlight(p,Iboth,'NodeColor',[143 108 169]/255)


set(gca,'xtick',[],'ytick',[])

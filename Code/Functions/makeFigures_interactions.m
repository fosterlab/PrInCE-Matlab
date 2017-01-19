%MAKEFIGURES_INTERACTIONS Makes figures for the PRINCE Interactions module.

precPlot = num2str(round(desiredPrecision(di)*100));


%% How many interactions were found in at least N replicates?

sf = [figdir 'Number_interactions_per_channel.png'];

jj = 1;
PA = Precision_array;

if size(Precision_array,1)==1
  bar_x = 1:(number_of_channels);
  bar_y = Precision_array;
  bar_x(2) = 0;
  bar_y(2,:) = 0;
  x0 = [.5 1.5];
else
  x0 = [bar_x(1)-1 bar_x(end)+1];
  bar_x = 1:(number_of_channels);
  bar_y = Precision_array;
  prec_fig = bar_y(:,2) ./ (bar_y(:,1) + bar_y(:,2));
end

myC= [30/255 144/255 255/255
  255/255 215/255 0/255
  178/255 34/255 34/255
  193/255 205/255 193/255];

figure
subplot(2,1,1)
bar(bar_x, bar_y, 0.6, 'stack');
y = ylim;
for ii = 1:length(prec_fig)
  if isnan(prec_fig(ii))
    s = '--';
  else
    s = [num2str(round(prec_fig(ii)*1000)/10) '%'];
  end
  texty = sum(bar_y(ii,:));
  text(ii-0.25,texty+diff(y)*.04,s)
end
legend('Inter-complex (FP)', 'Intra-complex (TP)',  'Novel');
xlim(x0)
grid on
ylabel('Number of interactions','FontSize', 8);
xlabel('N conditions','FontSize', 8);
title('Number of interactions in AT LEAST N conditions','FontSize', 12);

% How many interactions were found in each channel?
kk = 2;
if kk == 1 % global precision
  I = interaction_final.global;
else % replicate precision
  I = ones(size(interaction_final.global));
end
[Ia,Ib] = find(interaction_final.replicate_numbers>0 & repmat(I,1,size(interaction_final.replicate_numbers,2))>0 );
N_intperrep = zeros(length(Ia), 1);
itype = zeros(length(Ia), 1);
for ii = 1:length(Ia)
  N_intperrep(ii) = interaction_final.replicate_numbers(Ia(ii),Ib(ii));
  if interaction_final.proteinInCorum(Ia(ii)) && interaction_final.interactionInCorum(Ia(ii))
    itype(ii) = 1;
  elseif interaction_final.proteinInCorum(Ia(ii)) && ~interaction_final.interactionInCorum(Ia(ii))
    itype(ii) = 2;
  else
    itype(ii) = 3;
  end
end

x = 1:number_of_replicates*number_of_channels;
h_tp = hist(N_intperrep(itype==1),x);
h_fp = hist(N_intperrep(itype==2),x);
h_nov = hist(N_intperrep(itype==3),x);
hh = [h_fp; h_tp; h_nov];
xlim_flag = 0;
if length(x)==1
  x(2) = x(1)+1;
  hh(:,2) = zeros(length(hh),1);
  xlim_flag = 1;
end

xlab = cell(number_of_replicates*number_of_channels,1);
for ii = 1:number_of_replicates*number_of_channels
  sr = nan(length(user.silacratios),1);
  for jj = 1:length(user.silacratios)
    tmp = findstr(user.silacratios{jj},GaussIn{ii});
    if ~isempty(tmp)
      sr(jj) = tmp;
    end
  end
  Irep = findstr('rep',GaussIn{ii});
  rep = GaussIn{ii}(Irep+3);
  xlab{ii} = ['Rep' rep '-' user.silacratios{sr>0}];
end

subplot(2,1,2),hold on
f2 = bar(x,hh', 0.6,'stacked');
xlim(x0)
ylabel('Number of interactions','FontSize', 10);
xlabel('Replicate - Isotopologue channel','FontSize', 10);
title('Number of interactions in EACH condition','FontSize', 12);
set(gca,'xtick',x,'xticklabel',xlab,'fontsize',9)
grid on

saveas(gcf, sf);


%% Final precision-recall, ROC curves

[~,I] = sort(interaction_final.precisionDropout(:),'descend');
tmp = interaction_final.precisionDropoutavg(:);
tmp(tmp==0) = nan;

cols = rand(number_of_channels,3);
cols(1,:) = [1 0 0];
if number_of_channels>1
  cols(2,:) = [0 1 0];
end

I2 = nan(length(I),1);
for ii = 1:length(I)
  tmp2 = interaction_final.channel{ii};
  if length(tmp2)>1
    I2(ii) = 0;
  else
    I2(ii) = tmp2;
  end
end
I2 = I2(I);

xi = 1:length(tmp);
figure,hold on
plot(xi,tmp(I),'k')
grid on
ylabel('Interaction score','fontsize',12)
xlabel('Number of interactions','fontsize',12)
% Save figure
set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
sf=[figdir 'Final_Precision_' precPlot];
saveas(gcf, sf, 'png');


%%


tpfp = classList;
tpfp(possList==0) = nan;
[~,I] = sort(score,'descend');
tpfp = tpfp(I);   % order  by score

I = find(tpfp==1);
tp2plot = nan(length(I)*3,2);
for ii = 1:length(I)
  tp2plot((ii-1)*3 + (1:2),1) = I(ii);
  tp2plot((ii-1)*3 + (1:2),2) = [0 1]';
end

I = find(tpfp==0);
fp2plot = nan(length(I)*3,2);
for ii = 1:length(I)
  fp2plot((ii-1)*3 + (1:2),1) = I(ii);
  fp2plot((ii-1)*3 + (1:2),2) = [0 1]';
end


figure

subplot(4,1,1:3)
loglog(cumsum(tpfp==1),'g','linewidth',2)
hold on
loglog(cumsum(tpfp==0),'r','linewidth',2)
xlim([.99 5000])
legend('Intra-complex (TP)','Inter-complex (FP)','location','northwest')
ylabel('Cumulative number of interactions, TP or FP')
grid on

subplot(4,1,4)
semilogx(fp2plot(:,1),fp2plot(:,2),'r');
hold on
semilogx(tp2plot(:,1),tp2plot(:,2),'g');
ylim([-.25 1.25])
set(gca,'ytick',[])
xlim([.99 5000])
xlabel('Interaction number, ranked by interaction score')

% Save figure
set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
sf=[figdir 'Ranked_interactions_cumulative_TP_FP' precPlot];
saveas(gcf, sf, 'png');

% %% Histogram of score for class=0, class=1
%
% Nx = max([15 length(classList)/3000]);
% x = linspace(min(score),max(score),Nx);
% h1 = hist(score(classList==1),x);
% h0 = hist(score(classList==0),x);
% 
% figure,hold on
% p0 = patch([0 x x(end)],[0 h0/sum(h0) 0],'r');
% p1 = patch([0 x x(end)],[0 h1/sum(h1) 0],'g');
% ax = axis;
% plot([1 1].*xcutoff(di),[ax(3) ax(4)],'--r')
% h = text(xcutoff(di)-diff(ax(1:2))*.02,ax(3)+diff(ax(3:4))*.65,...
%   ['Desired precision = ' num2str(round(desiredPrecision(di)*100)) '%']);
% set(h, 'rotation', 90)
% axis(ax)
% p0.FaceAlpha = 0.4;
% p1.FaceAlpha = 0.4;
% p0.EdgeColor = 'none';
% p1.EdgeColor = 'none';
% grid on
% hl = legend('Known non-interaction','Known interaction','location','northwest');
% set(hl,'fontsize',12)
% xlabel('Interaction score','fontsize',12)
% ylabel('Count, number of interactions, normalized','fontsize',12)
% title('Histogram of interaction scores for protein pairs in CORUM','fontsize',12)
% 
% % Save figure
% set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
% sf=[figdir 'ScoreHistogram_' precPlot '_GlobalPrecision'];
% saveas(gcf, sf, 'png');



%% Plot best 200 interactions

f1 = [user.maindir 'Output/Figures/Interactions/BestInteractions/'];
if not(exist(f1,'dir'))
  mkdir(f1)
end

% get chromatograms and proteins
chroms = cell(number_of_channels,length(rep2channel));
prots = cell(number_of_channels,length(rep2channel));
nfrac = 0;
for channel_counter = 1:number_of_channels
  replicatesThisChannel = find(rep2channel == channel_counter);
  for replicate_counter = replicatesThisChannel
    sf = [maindir '/Output/tmp/' 'data_rep' num2str(replicate_counter) '_chan' num2str(channel_counter) '.mat'];
    load(sf)
    chroms{channel_counter,replicate_counter} = Chromatograms_raw;
    prots{channel_counter,replicate_counter} = Protein;
    nfrac = max([nfrac size(Chromatograms_raw,2)]);
    clear TP_Matrix possibleInts Protein inverse_self Chromatograms Chromatograms_raw Dist
  end
end

[tmp,I] = sort(interaction_final.precisionDropout,'descend');
nplots = min([200  length(I)]);
for ii = 1:nplots
  intI = I(ii);
  chans = interaction_final.channel{intI};
  protA = interaction_final.proteinA{intI};
  protB = interaction_final.proteinB{intI};
  
  if isempty(chans)
    continue
  elseif length(chans)==1
    subplots = [1 1];
  elseif length(chans)==2
    subplots = [2 1];
  elseif length(chans)==3
    subplots = [2 2];
  elseif length(chans)==4
    subplots = [2 2];
  elseif length(chans)==5
    subplots = [3 2];
  elseif length(chans)==6
    subplots = [3 2];
  elseif length(chans)>6 && length(chans)<=0
    subplots = [3 3];
  elseif length(chans)>9
    subplots = [ceil(length(chans)/4) 4];
  end
  
  figure
  for jj = 1:length(chans)
    for kk = 1:length(rep2channel)
      if isempty(prots{jj,kk})
        continue;
      end
      Ia = find(ismember(prots{jj,kk}.Isoform,protA));
      Ib = find(ismember(prots{jj,kk}.Isoform,protB));
      
      subplot(subplots(1),subplots(2),jj), hold on
      plot(1:length(chroms{jj,kk}(Ia,:)),chroms{jj,kk}(Ia,:),'k','linewidth',2)
      plot(1:length(chroms{jj,kk}(Ib,:)),chroms{jj,kk}(Ib,:),'color',[.75 .75 .75],'linewidth',2)
      legend(protA,protB,'location','best')
    end
    ylabel('Protein amount')
    xlabel('Fraction number')
    title(user.silacratios{jj})
    xlim([0 nfrac+1])
  end
  ax = axis;
  s = ['Prec=' num2str(round(interaction_final.precisionDropout(intI)*100)) '%'];
  text(ax(1)+(ax(2)-ax(1))*.2,ax(3)+(ax(4)-ax(3))*.5,s)
  % Save figure
  set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
  sf=[f1 'Interactions_' num2str(ii) '_' protA '_' protB];
  saveas(gcf, sf, 'png');
  close all
end


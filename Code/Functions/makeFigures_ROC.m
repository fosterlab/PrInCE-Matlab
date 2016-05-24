

precPlot = num2str(round(desiredPrecision(precision_write_out_counter)*100));


%% How many interactions were found in at least N replicates?

sf = cell(2,1);
sf{1} = [figdir 'Interactions_in_multiple_replicates_' precPlot '_GlobalPrecision.png'];
sf{2} = [figdir 'Interactions_in_multiple_replicates_' precPlot '_ReplicatePrecision.png'];
titles = cell(2,1);
titles{1} =   ['Desired precision = ' precPlot '% (applied globally)'];
titles{2} =   ['Desired precision = ' precPlot '% (applied at replicate)'];


for jj = 1:2
  PA = Precision_array{jj};
  
  %Interaction_not_in_corum=(PA(:,3)-PA(:,1));
  %Interaction_in_corum_not_detected=(PA(:,1)-PA(:,2));
  bar_x = 1:(number_of_replicates*number_of_channels);
  %bar_y = [Interaction_in_corum_not_detected(:,1) (PA(:,2)) Interaction_not_in_corum(:,1)];
  bar_y = Precision_array{jj};
  prec_fig = bar_y(:,2) ./ (bar_y(:,1) + bar_y(:,2));
  
  myC= [30/255 144/255 255/255
    255/255 215/255 0/255
    178/255 34/255 34/255
    193/255 205/255 193/255];
  
  figure
  subplot(2,1,1)
  xlim_flag = 0;
  if length(bar_x)==1
    bar_x(2) = 0;
    bar_y(2,:) = zeros(size(bar_y));
    xlim_flag = 1;
  end
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
  legend('Proteins in corum (FP)', 'Interaction in corum (TP)',  'Proteins/Interactions not in corum');
  if xlim_flag
    x = [0.4 1.6];
    xlim(x)
  end
  grid on
  ylabel('Number of interactions','FontSize', 8);
  xlabel('Interaction observed in at least N channels/replicates','FontSize', 8);
  title(titles{jj})
  
  f1=subplot(2,1,2);
  bar(bar_x, bar_y, 0.6, 'stack');
  I2use = ceil(size(PA,1)/2);
  y2 = sum(bar_y(I2use,:))*1.2;
  if y2==0
    y2 = y(2)/2;
  end
  ylim([0,y2]);
  if xlim_flag
    x = [0.4 1.6];
    xlim(x)
  end
  y = ylim;
  for ii = 1:length(prec_fig)
    if isnan(prec_fig(ii))
      s = '--';
    else
      s = [num2str(round(prec_fig(ii)*1000)/10) '%'];
    end
    texty = sum(bar_y(ii,:));
    if texty<y(2)
      text(ii-0.25,texty+diff(y)*.04,s)
    end
  end
  grid on
  ylim(y)
  ylabel('Number of interactions','FontSize', 8);
  xlabel('Interaction observed in at least N channels/replicates','FontSize', 8);
  
  % Save figure
  saveas(gcf, sf{jj});
  
end



%% How many interactions were found in each replicate?

sf = cell(2,1);
sf{1} = [figdir 'Number_interactions_per_replicate' precPlot '_GlobalPrecision.png'];
sf{2} = [figdir 'Number_interactions_per_replicate' precPlot '_ReplicatePrecision.png'];
titles = cell(2,1);
titles{1} =   ['Desired precision = ' precPlot '% (applied globally)'];
titles{2} =   ['Desired precision = ' precPlot '% (applied at replicate)'];

for kk = 1:2
  
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
    x(2) = 0;
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
  
  figure,hold on
  f2 = bar(x,hh', 0.8,'stacked');
  legend('Proteins in corum (FP)', 'Interaction in corum (TP)',  'Proteins/Interactions not in corum',...
    'location','northwest');
  ylabel('Number of interactions','FontSize', 10);
  xlabel('Replicate - Isotopologue channel','FontSize', 10);
  title(titles{kk},'FontSize', 12);
  if xlim_flag
    x = xlim;
    x(2) = 1.6;
    xlim(x)
  end
  set(gca,'xtick',x,'xticklabel',xlab,'fontsize',9)
  grid on
  
  saveas(gcf, sf{kk});
end



%% Final precision-recall, ROC curves

Recall_from_Ninteract = Ninteract / all_positives;
I = Recall_from_Ninteract > max(Recall_plot)*1.02;
g = interaction_final.global==1;

x1 = min(Recall_from_Ninteract);
x2 = max(Recall_from_Ninteract);

figure,hold on

% dummy plotting for legend
scatter(-100,-100,10,'k','filled')
scatter(-100,-100,10,[.75 .75 .75],'filled')
plot([-100 -101], [-100 -101],'k')

scatter(Recall_plot(g), interaction_final.precisionDropout(g),10,'k','filled')
scatter(Recall_plot(~g), interaction_final.precisionDropout(~g),10,[.75 .75 .75],'filled')
plot(Recall_from_Ninteract,precRange,'k')
xlabel({'Recall' ['as fraction of ' num2str(all_positives) ' interactions in CORUM']},'fontsize',12)
ylabel('Precision','fontsize',12)
for ii = 1:length(desiredPrecision)
  plot([x1 x2],[1 1].*desiredPrecision(ii),'--k')
  plot([1 1].*final_Recall(1),[0 desiredPrecision(ii)],'--k')
  text(x2*.7,desiredPrecision(ii)+.02,['Global precision = ' num2str(round(desiredPrecision(ii)*100)) '%'])
end
axis([x1 x2 0 1])
grid on
legend('Interactions assessed globally','Interactions assessed per replicate','Estimate','location','northeast')

% Save figure
set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
sf=[figdir 'Final_PrecisionRecall_' precPlot '_GlobalPrecision'];
saveas(gcf, sf, 'png');


figure,hold on

% dummy plotting for legend
scatter(-100,-100,10,'k','filled')
scatter(-100,-100,10,[.75 .75 .75],'filled')
plot([-100 -101], [-100 -101],'k')

[~,I1] = min(abs(tprRange -final_Recall(1)));
[~,I2] = min(abs(tprRange -final_Recall(2)));
plot(fprRange,tprRange,'k')
scatter(fprRange(I1),tprRange(I1),50,'k')
scatter(fprRange(I2),tprRange(I2),50,[.6 .6 .6])
plot([0 1],[0 1],'--r')
xlabel('False positive rate, FP/(FP+TN)','fontsize',12)
ylabel('True positive rate, TP/(TP+FN)','fontsize',12)
axis([-0.01 1 0 1])
legend('Interactions assessed globally','Interactions assessed per replicate','location','southeast')
grid on
% Save figure
set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
sf=[figdir 'Final_ROC_' precPlot '_GlobalPrecision'];
saveas(gcf, sf, 'png');



%% Histogram of score for class=0, class=1

Nx = max([15 length(class)/3000]);
x = linspace(min(score),max(score),Nx);
h1 = hist(score(class==1),x);
h0 = hist(score(class==0),x);

figure,hold on
p0 = patch([0 x x(end)],[0 h0 0],'r');
p1 = patch([0 x x(end)],[0 h1 0],'g');
ax = axis;
plot([1 1].*xcutoff(di),[ax(3) ax(4)],'--r')
h = text(xcutoff(di)-diff(ax(1:2))*.02,ax(3)+diff(ax(3:4))*.65,...
  ['Global precision = ' num2str(round(desiredPrecision(di)*100)) '%']);
set(h, 'rotation', 90)
for ii = 1:size(xcutoff_rep,1)
  plot([1 1].*xcutoff_rep(ii,di),[ax(3) ax(4)],'linestyle','--','color',[.6 .6 .6])
end
axis(ax)
p0.FaceAlpha = 0.4;
p1.FaceAlpha = 0.4;
p0.EdgeColor = 'none';
p1.EdgeColor = 'none';
grid on
hl = legend('Known non-interaction','Known interaction','Global threshold','Replicate thresholds','location','northwest');
set(hl,'fontsize',12)
xlabel('Interaction score','fontsize',12)
ylabel('Count, number of interactions, normalized','fontsize',12)
title('Histogram of interaction scores for protein pairs in CORUM','fontsize',12)

% Save figure
set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
sf=[figdir 'ScoreHistogram_' precPlot '_GlobalPrecision'];
saveas(gcf, sf, 'png');



%% Normalized histogram of score for class=0, class=1
x = linspace(min(score),max(score),Nx);
h1 = hist(score(class==1),x);
h0 = hist(score(class==0),x);

figure,hold on
p0 = patch([0 x x(end)],[0 h0/sum(h0) 0],'r');
p1 = patch([0 x x(end)],[0 h1/sum(h1) 0],'g');
ax = axis;
plot([1 1].*xcutoff(di),[ax(3) ax(4)],'--r')
h = text(xcutoff(di)-diff(ax(1:2))*.02,ax(3)+diff(ax(3:4))*.65,...
  ['Global precision = ' num2str(round(desiredPrecision(di)*100)) '%']);
set(h, 'rotation', 90)
for ii = 1:size(xcutoff_rep,1)
  plot([1 1].*xcutoff_rep(ii,di),[ax(3) ax(4)],'linestyle','--','color',[.6 .6 .6])
end
axis(ax)
p0.FaceAlpha = 0.4;
p1.FaceAlpha = 0.4;
p0.EdgeColor = 'none';
p1.EdgeColor = 'none';
grid on
hl = legend('Known non-interaction','Known interaction','Global threshold','Replicate thresholds','location','northwest');
set(hl,'fontsize',12)
xlabel('Interaction score','fontsize',12)
ylabel('Count, number of interactions, normalized','fontsize',12)
title('Histogram of interaction scores for protein pairs in CORUM','fontsize',12)

% Save figure
set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
sf=[figdir 'ScoreHistogram_normalized_' precPlot '_GlobalPrecision'];
saveas(gcf, sf, 'png');


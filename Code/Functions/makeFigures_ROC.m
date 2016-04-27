
if ~exist('firstFlag','var')
  firstFlag = 1;
end

%% How many interactions were found in at least N replicates?

% #PARSE_HERE
% Make two figures: one global, one replicate-specific

Interaction_not_in_corum=(Precision_array(:,3)-Precision_array(:,1));
Interaction_in_corum_not_detected=(Precision_array(:,1)-Precision_array(:,2));

myC= [30/255 144/255 255/255
  255/255 215/255 0/255
  178/255 34/255 34/255
  193/255 205/255 193/255];

figure
subplot(2,1,1)
bar_x = 1:(number_of_replicates*number_of_channels);
bar_y = [Interaction_in_corum_not_detected(:,1) (Precision_array(:,2)) Interaction_not_in_corum(:,1)];
prec_fig = bar_y(:,2) ./ (bar_y(:,1) + bar_y(:,2));
xlim_flag = 0;
if length(bar_x)==1
  bar_x(2) = 0;
  bar_y(2,:) = zeros(size(bar_y));
  xlim_flag = 1;
end
b1 = bar(bar_x, bar_y, 0.6, 'stack');
y = ylim;
for ii = 1:length(prec_fig)
  if isnan(prec_fig(ii))
    s = '--';
  else
    s = [num2str(round(prec_fig(ii)*1000)/10) '%'];
  end
  texty = sum(bar_y(ii,:));
  text(ii-0.25,texty+diff(y)*.025,s)
end
legend('Proteins in corum (FP)', 'Interaction in corum (TP)',  'Proteins/Interactions not in corum');
ylim([0,(Interaction_in_corum_not_detected(1,1)+Precision_array(1,2)+Interaction_not_in_corum(1,1))*1.1]);
if xlim_flag
  x = [0.4 1.6];
  xlim(x)
end
title({'Interactions observed across isotoplogue channels' ['Desired precision = ' num2str(round(desiredPrecision(pri)*100)) '%']},'FontSize', 12);
ylabel('Number of interactions','FontSize', 8);
xlabel('Interaction observed in at least N channels/replicates','FontSize', 8);

f1=subplot(2,1,2);
b2 = bar(bar_x, bar_y, 0.6, 'stack');
legend('Proteins in corum (FP)', 'Interaction in corum (TP)',  'Proteins/Interactions not in corum');
I2use = ceil(size(Precision_array,2)/2);
y2 = max(sum(Precision_array(:,I2use)))*1.2;
if y2==0
  y2 = (Interaction_in_corum_not_detected(1,1)+Precision_array(1,2)+Interaction_not_in_corum(1,1))*1.1;
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
  text(ii-0.25,texty+diff(y)*.025,s)
end
ylim(y)
title('Interactions observed across isotoplogue channels (Zoom)','FontSize', 12);
ylabel('Number of interactions','FontSize', 8);
xlabel('Interaction observed in at least N channels/replicates','FontSize', 8);

% Save figure
Final_Interaction_figure=[figdir1 'Interactions_in_multiple_replicates_' mat2str(Precision_values(precision_write_out_counter)) '_precision.png'];
saveas(gcf, Final_Interaction_figure);



%% How many interactions were found in each replicate?

% #PARSE_HERE
% Make two figures: one global, one replicate-specific

[Ia,Ib] = find(~cellfun('isempty',interaction_final.replicate_numbers));

N_intperrep = zeros(length(Ia), 1);
itype = zeros(length(Ia), 1);
for ii = 1:length(Ia)
  N_intperrep(ii) = interaction_final.replicate_numbers{Ia(ii),Ib(ii)};
  if interaction_final.proteinInCorum{Ia(ii)} && interaction_final.interactionInCorum{Ia(ii)}
    itype(ii) = 1;
  elseif interaction_final.proteinInCorum{Ia(ii)} && ~interaction_final.interactionInCorum{Ia(ii)}
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
    sr(jj) = findstr(user.silacratios{1},GaussIn{1});
  end
  Irep = findstr('rep',GaussIn{ii});
  rep = GaussIn{ii}(Irep+3);
  xlab{ii} = ['Rep' rep '-' user.silacratios{sr>0}];
end

figure,hold on
f2 = bar(x,hh', 0.8,'stacked');
for k=1:3
  set(f2(k),'facecolor', myC(k,:), 'EdgeColor', 'k' );
end
legend('Proteins in corum (FP)', 'Interaction in corum (TP)',  'Proteins/Interactions not in corum',...
  'location','northwest');
ylabel('Number of interactions','FontSize', 10);
xlabel('Replicate - Isotopologue channel','FontSize', 10);
title(['Desired precision = ' num2str(round(desiredPrecision(pri)*100)) '%'],'FontSize', 12);
if xlim_flag
  x = xlim;
  x(2) = 1.6;
  xlim(x)
end
set(gca,'xtick',x,'xticklabel',xlab,'fontsize',9)

sf = [figdir1 'Number_interactions_per_replicate' mat2str(Precision_values(precision_write_out_counter)) '_precision.png'];
saveas(gcf, sf);



%% Final precision-recall, ROC curves

%#PARSE_HERE
% Make two figures, or can I put both global and replicate-specific here?
if firstFlag ==1
  % convert original Recall to fraction of possible CORUM interactions
%   TP = Ninteract;
%   allPositivesinCorum = final_FN + final_TP;
%   Recall = TP / allPositivesinCorum;
%   calcrec_final = final_TP / (final_FN + final_TP);
%   
%   recall_ratio = calcrec(di) / calcrec_final;
%   
%   x1 = min(Recall);
%   x2 = max(Recall);

  all_positives = final_TP + final_FN;
  Recall_from_Ninteract = Ninteract / all_positives;
  Recall_this_precision = calcrec * sum(class==1) / all_positives;
  
  figure,hold on
  plot(Recall_from_Ninteract,precRange)
  %plot(recRange,precRange)
  xlabel({'Recall' 'as fraction of all possible CORUM interactions'},'fontsize',12)
  ylabel('Precision','fontsize',12)
  for ii = 1:length(desiredPrecision)
    plot([x1 x2],[1 1].*desiredPrecision(ii),'--r')
    %plot([1 1].*calcrec(ii),[0 desiredPrecision(ii)],'--k')
    plot([1 1].*Recall_this_precision(ii),[0 desiredPrecision(ii)],'--k')
    text(.8,desiredPrecision(ii)+.02,['precision = ' num2str(round(desiredPrecision(ii)*100)) '%'])
  end
  axis([x1 x2 0 1])
  
  % Save figure
  set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
  sf=[figdir1 'Final_PrecisionRecall'];
  saveas(gcf, sf, 'png');
  
  
  figure,hold on
  plot(fprRange,tprRange)
  plot([0 1],[0 1],'--r')
  xlabel('False positive rate, FP/(FP+TN)','fontsize',12)
  ylabel('True positive rate, TP/(TP+FN)','fontsize',12)
  axis([0 1 0 1])
  
  % Save figure
  set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
  sf=[figdir1 'Final_ROC'];
  saveas(gcf, sf, 'png');
end


%% Histogram of score for class=0, class=1

if firstFlag ==1
  
  x = linspace(min(score),max(score),length(class)/50);
  h1 = hist(score(class==1),x);
  h0 = hist(score(class==0),x);
  
  figure,hold on
  p0 = patch([0 x x(end)],[0 h0 0],'r');
  p1 = patch([0 x x(end)],[0 h1 0],'g');
  ax = axis;
  for ii = 1:length(desiredPrecision)
    plot([1 1].*xcutoff(ii),[ax(3) ax(4)],'--r')
    h = text(xcutoff(ii)-diff(ax(1:2))*.02,ax(3)+diff(ax(3:4))*.75,...
      ['precision = ' num2str(round(desiredPrecision(ii)*100)) '%']);
    set(h, 'rotation', 90)
  end
  axis(ax)
  p0.FaceAlpha = 0.4;
  p1.FaceAlpha = 0.4;
  grid on
  hl = legend('Known non-interaction','Known interaction','location','northwest');
  set(hl,'fontsize',12)
  xlabel('Interaction score','fontsize',12)
  ylabel('Count, number of interactions','fontsize',12)
  title('Histogram of interaction scores for protein pairs in CORUM','fontsize',12)
  
  % Save figure
  set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
  sf=[figdir1 'ScoreHistogram'];
  saveas(gcf, sf, 'png');
end


%% Normalized histogram of score for class=0, class=1
if firstFlag ==1
  x = linspace(min(score),max(score),length(class)/50);
  h1 = hist(score(class==1),x);
  h0 = hist(score(class==0),x);
  
  figure,hold on
  p0 = patch([0 x x(end)],[0 h0/sum(h0) 0],'r');
  p1 = patch([0 x x(end)],[0 h1/sum(h1) 0],'g');
  ax = axis;
  for ii = 1:length(desiredPrecision)
    plot([1 1].*xcutoff(ii),[ax(3) ax(4)],'--r')
    h = text(xcutoff(ii)-diff(ax(1:2))*.02,ax(3)+diff(ax(3:4))*.75,...
      ['precision = ' num2str(round(desiredPrecision(ii)*100)) '%']);
    set(h, 'rotation', 90)
  end
  axis(ax)
  p0.FaceAlpha = 0.4;
  p1.FaceAlpha = 0.4;
  grid on
  hl = legend('Known non-interaction','Known interaction','location','northwest');
  set(hl,'fontsize',12)
  xlabel('Interaction score','fontsize',12)
  ylabel('Count, number of interactions, normalized','fontsize',12)
  title('Histogram of interaction scores for protein pairs in CORUM','fontsize',12)
  
  % Save figure
  set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
  sf=[figdir1 'ScoreHistogram_normalized'];
  saveas(gcf, sf, 'png');
end
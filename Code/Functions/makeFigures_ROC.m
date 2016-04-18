%Plot Precision calculations
Interaction_not_in_corum=(Precision_array(:,3)-Precision_array(:,1));
Interaction_in_corum_not_detected=(Precision_array(:,1)-Precision_array(:,2));

myC= [30/255 144/255 255/255
  255/255 215/255 0/255
  178/255 34/255 34/255
  193/255 205/255 193/255];

figure
subplot(2,1,1);
f1_figure=bar(1:(number_of_replicates*number_of_channels), [Interaction_in_corum_not_detected(:,1) (Precision_array(:,2)) Interaction_not_in_corum(:,1)], 0.6, 'stack');

for k=1:3
  set(f1_figure(k),'facecolor', myC(k,:), 'EdgeColor', 'k' );
end
legend('Proteins in corum (FP)', 'Interaction in corum (TP)',  'Proteins/Interactions not in corum');
ylim([0,(Interaction_in_corum_not_detected(1,1)+Precision_array(1,2)+Interaction_not_in_corum(1,1))*1.1]);
title({'Interactions observed across isotoplogue channels' ['Precision = ' num2str(round(desiredPrecision(pri)*100)) '%']},'FontSize', 12);
ylabel('Number of interactions','FontSize', 8);
xlabel('isotoplogue channels','FontSize', 8);

f1=subplot(2,1,2);
f2_figure=bar(1:(number_of_replicates*number_of_channels), [Interaction_in_corum_not_detected(:,1) (Precision_array(:,2)) Interaction_not_in_corum(:,1)], 0.6, 'stack');

for k=1:3
  set(f2_figure(k),'facecolor', myC(k,:), 'EdgeColor', 'k' );
end
legend('Proteins in corum (FP)', 'Interaction in corum (TP)',  'Proteins/Interactions not in corum');
I2use = ceil(size(Precision_array,2)/2);
y2 = max(sum(Precision_array(:,I2use)))*1.2;
if y2==0
  y2 = (Interaction_in_corum_not_detected(1,1)+Precision_array(1,2)+Interaction_not_in_corum(1,1))*1.1;
end
ylim([0,y2]);
title('Interactions observed across isotoplogue channels (Zoom)','FontSize', 12);
ylabel('Number of interactions','FontSize', 8);
xlabel('isotoplogue channels','FontSize', 8);

% Save figure
Final_Interaction_figure=[figdir1 'Observed_interactions_across_replicates_' mat2str(Precision_values(precision_write_out_counter)) '_precision.png'];
saveas(gcf, Final_Interaction_figure);



%% Final precision-recall, ROC curves

if firstFlag ==1
  figure,hold on
  plot(recRange,precRange)
  xlabel('Recall','fontsize',12)
  ylabel('Precision','fontsize',12)
  for ii = 1:length(desiredPrecision)
    plot([0 1],[1 1].*desiredPrecision(ii),'--r')
    plot([1 1].*calcrec(ii),[0 desiredPrecision(ii)],'--k')
    text(.8,desiredPrecision(ii)+.02,['precision = ' num2str(round(desiredPrecision(ii)*100)) '%'])
  end
  axis([0 1 0 1])
  
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
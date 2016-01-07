%Plot Precision calculations
Interaction_not_in_corum=(Precision_array(:,3)-Precision_array(:,1));
Interaction_in_corum_not_detected=(Precision_array(:,1)-Precision_array(:,2));

myC= [30/255 144/255 255/255
  255/255 215/255 0/255
  178/255 34/255 34/255
  193/255 205/255 193/255];

figure
f1=subplot(2,1,1);
f1_figure=bar(1:(number_of_replicates*number_of_channels), [Interaction_in_corum_not_detected(:,1) (Precision_array(:,2)) Interaction_not_in_corum(:,1)], 0.6, 'stack');

for k=1:3
  set(f1_figure(k),'facecolor', myC(k,:), 'EdgeColor', 'k' );
end
legend('Proteins in corum (FP)', 'Interaction in corum (TP)',  'Proteins/Interactions not in corum');
ylim([0,(Interaction_in_corum_not_detected(1,1)+Precision_array(1,2)+Interaction_not_in_corum(1,1))*1.1]);
title('Interactions observed across isotoplogue channels','FontSize', 12);
ylabel('Number of interactions','FontSize', 8);
xlabel('isotoplogue channels','FontSize', 8);

f1=subplot(2,1,2);
f2_figure=bar(1:(number_of_replicates*number_of_channels), [Interaction_in_corum_not_detected(:,1) (Precision_array(:,2)) Interaction_not_in_corum(:,1)], 0.6, 'stack');

for k=1:3
  set(f2_figure(k),'facecolor', myC(k,:), 'EdgeColor', 'k' );
end
legend('Proteins in corum (FP)', 'Interaction in corum (TP)',  'Proteins/Interactions not in corum');
ylim([0,(Interaction_in_corum_not_detected(3,1)+Precision_array(3,2)+Interaction_not_in_corum(3,1))*1.2]);
title('Interactions observed across isotoplogue channels (Zoom)','FontSize', 12);
ylabel('Number of interactions','FontSize', 8);
xlabel('isotoplogue channels','FontSize', 8);

% Save figure
Final_Interaction_figure=strcat('Observed_interactions_across_replicates_',mat2str(Precision_values(precision_write_out_counter)),'_precision.png');
saveas(f1, Final_Interaction_figure);
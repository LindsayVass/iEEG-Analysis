function confusionplot(data,labels,figuretitle,rotateXlabels)
% Make a visual plot of a confusion matrix

% Plot the data
imagesc(data);

% Set axis labels
set(gca,'XTick',[1:1:size(data,1)]);
set(gca,'YTick',[1:1:size(data,1)]);
set(gca,'Ticklength',[0 0]);
set(gca,'XTickLabel',labels);
set(gca,'YTickLabel',labels);

% Make gridlines
for x = 1.5:1:size(data,1);
    line([x x],ylim,'linestyle','-','color','black');
    line(ylim,[x x],'linestyle','-','color','black');
end

% Make it big
%set(gcf,'PaperPositionMode','auto');
%set(gcf,'Position',[1 1 1600 1600]);

% Set the title
title(figuretitle,'FontSize',12);
if(rotateXlabels == 1)
    xticklabel_rotate;
end

end


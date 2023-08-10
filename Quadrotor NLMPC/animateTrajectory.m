
fig = figure('Name','Trajectory Animation');
hold on
for i=1:size(xHistory,1)
    clf;
    animate(time(i), xHistory(i,:));
    exportgraphics(fig, "F.gif", Append=true);   % to save gif of the fig
    pause(Ts);    
end
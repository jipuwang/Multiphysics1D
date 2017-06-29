
%% For scalar flux
fontSizeLarge=14;
% quantity='temperature';
quantity='scalar flux';
for i=1:3
  figure(i)
  if i==1
    title('noFB\_LSMoC\_const\_cubic\_phi0\_RMS','FontSize',fontSizeLarge,'FontWeight','bold')
  end
  if i==2
    title('linear\_LSMoC\_const\_cubic\_phi0\_RMS','FontSize',fontSizeLarge,'FontWeight','bold')
  end
  if i==3
    title('sqrt\_LSMoC\_const\_cubic\_phi0\_RMS','FontSize',fontSizeLarge,'FontWeight','bold')
  end
  xlabel('mesh size [cm]','FontSize',fontSizeLarge,'FontWeight','bold')
  set(gca,'FontSize',fontSizeLarge)
  ylabel(quantity,'FontSize',fontSizeLarge,'FontWeight','bold');
  set(gca,'FontSize',fontSizeLarge,'FontWeight','bold')
end

quantity='temperature';
for i=4:6
  figure(i)
  if i==4
    title('noFB\_LSMoC\_const\_cubic\_T\_RMS','FontSize',fontSizeLarge,'FontWeight','bold')
  end
  if i==5
    title('linear\_LSMoC\_const\_cubic\_T\_RMS','FontSize',fontSizeLarge,'FontWeight','bold')
  end
  if i==6
    title('sqrt\_LSMoC\_const\_cubic\_T\_RMS','FontSize',fontSizeLarge,'FontWeight','bold')
  end
  xlabel('mesh size [cm]','FontSize',fontSizeLarge,'FontWeight','bold')
  set(gca,'FontSize',fontSizeLarge)
  ylabel(quantity,'FontSize',fontSizeLarge,'FontWeight','bold');
  set(gca,'FontSize',fontSizeLarge,'FontWeight','bold')
end
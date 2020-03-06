function divergence_settings
  view(2);
  set(gca,'fontsize',50);
  colormap('redblue');
  ylim([-90,90]);
  xlim([-180,180]);
  xticks([-180:60:180])
  yticks([-90:30:90]);
  hcb = colorbar('southoutside');
  set(hcb,"fontsize",50)
  caxis([-1.2,1.2]);
  set(hcb,'XTick',[-1.2:0.3:1.2])
endfunction

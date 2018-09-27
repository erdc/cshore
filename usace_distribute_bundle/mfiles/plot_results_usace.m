if isfield(results,'run_success')
  if ~results.run_success;
    disp('Due to failed execution, no plotting is attempted for run')
    disp(results.header)
    return
  end
end

if ~exist('iplotheader')|iplotheader==0
  dumtitle = '';
else
  dumtitle = results.header(2:end);
end

if ~exist('iplotbc');iplotbc=1;end
if ~exist('iplothydro');iplothydro=1;end
if ~exist('iplotrunup');iplotrunup=1;end


%%%%%%%%%%%%%%%%%%%%find the indices for plotting%%%%%%%%%%%%%%%%%%%%

num_of_snaps = 10; %change this to suit the plot
num_of_snaps = max(2,num_of_snaps);
if length(results.morpho)>1
  di = (length(results.morpho)-1)/(num_of_snaps-1);
  inds_morph = unique(round(di*round([0:length(results.morpho)-1]./di))) + 1;
  inds_hydro = unique([1 inds_morph(2:end)-1]);
else
  inds_morph=1;
  if length(results.hydro)>1
    di = (length(results.hydro)-1)/(num_of_snaps-1);
    inds_hydro = unique(round(di*round([0:length(results.hydro)-1]./di))) + 1;
  else
    inds_hydro=1;
  end
end

%%%%%%%%%%%%%%%%%%%%plot boundary conditions%%%%%%%%%%%%%%%%%%%%
if iplotbc
  figure;
  subplot(3,1,1)
  plot(results.bc.time_offshore,results.bc.Hrms_offshore);hold all
  plot(results.bc.time_offshore,results.bc.Hrms_offshore,'.');hold all
  ylabel('H_{rms}[m]','fontname','times','fontsize',14,'fontangle','italic')
  title('Wave Conditions at Seaward Boundary','fontname','times','fontsize',14,'fontangle','italic')
  subplot(3,1,2)
  plot(results.bc.time_offshore,results.bc.Tp_offshore);hold all
  plot(results.bc.time_offshore,results.bc.Tp_offshore,'.');hold all
  ylabel('T_{p}[m]','fontname','times','fontsize',14,'fontangle','italic')  
  subplot(3,1,3)
  plot(results.bc.time_offshore,results.bc.angle_offshore);hold all
  plot(results.bc.time_offshore,results.bc.angle_offshore,'.');hold all
  ylabel('angle[deg]','fontname','times','fontsize',14,'fontangle','italic')  
  xlabel('T[s]','fontname','times','fontsize',14,'fontangle','italic')  
  figure;
  subplot(3,1,1)
  plot(results.bc.time_offshore,results.bc.wave_setup_offshore);hold all
  plot(results.bc.time_offshore,results.bc.wave_setup_offshore,'.');hold all
  ylabel('wave setup[m]','fontname','times','fontsize',14,'fontangle','italic')
  title('Water Level Conditions at Seaward Boundary','fontname','times','fontsize',14,'fontangle','italic')
  subplot(3,1,2)
  plot(results.bc.time_offshore,results.bc.strm_tide_offshore);hold all
  plot(results.bc.time_offshore,results.bc.strm_tide_offshore,'.');hold all
  ylabel('water level[m]','fontname','times','fontsize',14,'fontangle','italic')
  subplot(3,1,3)
  d = results.bc.strm_tide_offshore+results.bc.strm_tide_offshore-results.morpho(1).zb(1);
  plot(results.bc.time_offshore,d);hold all
  plot(results.bc.time_offshore,d,'.');hold all
  ylabel('water depth[m]','fontname','times','fontsize',14,'fontangle','italic')
end




%%%%%%%%%%%%%%%%%%%%plot morphology%%%%%%%%%%%%%%%%%%%%
if exist('iplotmorpho')
  if iplotmorpho
    clear hh hlabs
    figure;
    if results.params.iveg
      ind= find(results.veg.n>1);
      hh123=fill([results.morpho(1).x(ind); flipud(results.morpho(1).x(ind))], ...
                 [results.morpho(1).zb(ind); flipud(results.morpho(1).zb(ind))+3],[.8 1 .8]);hold on
    end
    
    cnt = 0;
    for j= 1:length(inds_morph)
      i = inds_morph(j);
      cnt = cnt+1;
      hh(cnt)=plot(results.morpho(i).x,results.morpho(i).zb,'linewidth',2);hold all
      hlabs{cnt} = [num2str(results.morpho(i).time), ' s'];
      if results.params.iveg
        if isfield(results.morpho,'ivegitated')
          if ~isempty(results.morpho(i).ivegitated)
            ind = find(results.morpho(i).ivegitated);
            plot(results.morpho(i).x(ind),results.morpho(i).zb(ind),'+','linewidth',4,'color',get(hh(cnt),'Color'));
          end
        end
      end
      if results.params.isedav==1
        if isfield(results.morpho,'zb_p')
          plot(results.morpho(i).x,results.morpho(i).zb_p,'linewidth',3,'color',[ .8  .5  .0])
        end
      end
    end
    
    legend(hh,hlabs,'Location','Northwest')
    title('Bottom Position','fontname','times','fontsize',14,'fontangle','italic')
    ylabel('z[m]','fontname','times','fontsize',14,'fontangle','italic')
    xlabel('x[m]','fontname','times','fontsize',14,'fontangle','italic')
  end
end
%%%%%%%%%%%%%%%%%%%%plot hydrodynamics%%%%%%%%%%%%%%%%%%%%

if iplothydro
  clear hh hlabs
  figure;
  cnt = 0;
  t = [0 [results.hydro.time_end]];
  for j = 1:length(inds_hydro)
    i = inds_hydro(j);
    cnt = cnt+1;
    hh(cnt)=plot(results.hydro(i).x,results.hydro(i).Hrms,'linewidth',2);hold all
    hlabs{cnt} = [num2str(t(i)),' to ' num2str(t(i+1)), ' s'];
  end
  legend(hh,hlabs,1)
  title(vertcat('Root-mean-square wave height',dumtitle),...
        'fontname','times','fontsize',14,'fontangle','italic','interpreter','none')
  ylabel('H_{rms}[m]','fontname','times','fontsize',14,'fontangle','italic')
  xlabel('x[m]','fontname','times','fontsize',14,'fontangle','italic')
end

%%%%%%%%%%%%%%%%%%%%plot runup%%%%%%%%%%%%%%%%%%%%

if iplotrunup
  if isfield(results.hydro,'runup_2_percent')
    clear hh hlabs
    figure;
    cnt = 0;
    t = [0 [results.hydro.time_end]];
    for j = 1:length(inds_hydro)
      i = inds_hydro(j);
      if results.params.iprofl
	x = results.morpho(i).x;
	zb = results.morpho(i).zb;
      else
	x = results.morpho(1).x;
	zb = results.morpho(1).zb;
      end 
      z0 = zb-results.hydro(i).runup_2_percent;      
      crosspt = find(z0(2:end).*z0(1:end-1)<0);
      if ~isempty(crosspt)
	cnt = cnt+1;
	crosspt = crosspt(1);
	xrunup = interp1(z0(crosspt-1:crosspt+1),x(crosspt-1:crosspt+1),0);
	hh(cnt)=plot(x,zb);hold all
	plot(xrunup,results.hydro(i).runup_2_percent,'ko',...
             'markerfacecolor',get(hh(cnt),'color'))
	hlabs{cnt} = [num2str(t(i)),' to ' num2str(t(i+1)),...
                      ' s, R_{2%} = ',sprintf('%3.2f',results.hydro(i).runup_2_percent),' m'];
      end
    end
    if exist('hh')
      legend(hh,hlabs,'Location','Northwest')
    end
    title(vertcat('Runup Position',dumtitle),'fontname','times','fontsize',14,'fontangle','italic','interpreter','none')
    ylabel('z_{b}[m]','fontname','times','fontsize',14,'fontangle','italic')
    xlabel('x[m]','fontname','times','fontsize',14,'fontangle','italic')
  end
end

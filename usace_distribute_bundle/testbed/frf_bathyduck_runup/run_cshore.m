exec_cmd = '../../bin/CSHORE_USACE_LINUX.out';

cnt  = 0;
for j = 2:3 
  cnt = cnt+1;
  system([exec_cmd,' ',num2str(j)]);
  allresults(cnt).results= load_results_usace(num2str(j));
  load(['g',num2str(j)]);
  allresults(cnt).g = g;
end
save allresults allresults 

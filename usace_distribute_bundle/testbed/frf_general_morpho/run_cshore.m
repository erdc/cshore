exec_cmd = '../../bin/CSHORE_USACE_LINUX.out';

cnt  = 0;
fnames = dir('*.infile');


for j = 1:length(fnames)
  fname = fnames(j).name(1:3);
  system([exec_cmd,' ',fname]);
  allresults(j).results= load_results_usace(fname);
  load(['g_in',fname]);
  allresults(j).g = g;
  allresults(j).in = in;
end
save allresults allresults 

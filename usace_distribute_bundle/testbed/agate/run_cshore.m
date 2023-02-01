exec_cmd = '../../bin/CSHORE_USACE_LINUX.out';
addpath ../../mfiles/
fnames = dir('*.infile');

for j = 1:length(fnames);
  fn = fnames(j).name;
  fn = fn(1:length(fn)-7);
  system([exec_cmd,' ',fn]);
  allresults(j).results= load_results_usace(fn);
  allresults(j).exp_case = fn;
  % load(['g',num2str(j)]);
  % allresults(cnt).g = g;
end
save allresults allresults 

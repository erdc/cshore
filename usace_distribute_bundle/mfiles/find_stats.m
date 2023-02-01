function stats = find_stats(t,eta,numsegs)
if nargin<3;numsegs=5;end


stats.dt = (t(end)-t(1))/(length(t)-1);
stats.fs = 1/stats.dt;
nfft = round(length(eta)/numsegs);
wind = round(nfft/4);
[stats.Petaeta,stats.f] = pwelch(eta,nfft,stats.fs,wind);
stats.df = stats.f(2)-stats.f(1);
stats.mo=stats.df*sum(stats.Petaeta) ;
stats.Hmo = 4*sqrt(stats.mo);
[j1 j2] = max(stats.Petaeta);
stats.Tp = 1/stats.f(j2);
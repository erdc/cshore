function [kh] = dispersion(depth, period)
%[kh] = dispersion(depth, period)

  grav = 9.806 ; 
  omeg = 2*pi/period ; 
  options=optimset('Algorithm','Levenberg-Marquardt','Display','off');
  
  fun   = @(k) omeg^2 - grav*k.*tanh(k*depth) ; 
  wvnum = fsolve (fun, 1.0, options) ; 
  kh    = wvnum*depth ; 
  
end
function [SV1,SV2,SV3 ] = IRF( t, N ,T, par,B )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% State Vector
SV1   = zeros(N,2);
SV2   = zeros(N,2);
SV3   = zeros(N,2);

for i = 1: N
  if t(i)<T
   % Standard Equilibrium   ----------------------------------
     ir = par.in(1,1)-(-par.r(1,1));
     A = [par.rho;1]*ir;
     X = A - B*[exp(par.delta*(t(i)-T)); exp(par.lambda*(t(i)-T))]*ir; 
     SV1(i,:) = X';
   % Backward Stable------------------------------------------------
%      pi_T = lambda/(lambda-delta)*ir;
     X = A +(par.delta/(par.lambda-par.delta))*[par.delta;1]*exp(par.lambda*(t(i)-T))*ir;
     SV2(i,:) = X';
     % No-Jump Equilibrium
     pi_T = exp(par.delta*T)*((exp(-par.delta*T)-1)*par.lambda - (exp(-par.lambda*T)-1)*par.delta)...
         /(par.lambda-par.delta)*ir;
     X = A - B*[exp(par.delta*(t(i)-T)); exp(par.lambda*(t(i)-T))]*ir...
         + [par.lambda;1]*pi_T*exp(par.delta*(t(i)-T)); 
     SV3(i,:) = X';     
  else % for t>T
%    Standard
     pi_T = 0;
%      ir = in(1,2) - r(1,2);
     SV1(i,:) = pi_T ;
%      Backward Stable
     pi_T = par.lambda/(par.lambda-par.delta);
     X = [par.lambda;1]*pi_T*par.r(1,1)*exp(par.delta*(t(i)-T));
     SV2(i,:) = X';
% No-Jump Equilibium
     pi_T = exp(par.delta*T)*((exp(-par.delta*T)-1)*par.lambda - (exp(-par.lambda*T)-1)*par.delta)...
         /(par.lambda-par.delta)*par.r(1,1);
     X = [par.lambda;1]*pi_T*exp(par.delta*(t(i)-T));
     SV3(i,:) = X';     

  end
end
end


function [t,y,ssTR,ssRT] = delayQN_ODE(X0,MU,NC,dt,TF)
import Gillespie.*

% Make sure vector components are doubles
X0 = double(X0);
MU = double(MU);

% Make sure all vectors are row vectors
if(iscolumn(X0))
    X0 = X0';
end
if(iscolumn(MU))
    MU = MU';
end
if(iscolumn(NC))
    NC = NC';
end

p.MU = MU;
p.NC = NC;
p.delta = 10^5; % context switch rate (super fast)

%states name

%task ordering


% Jump matrix
jump=[-1,  +1,  +0,  +0,  +0;
    +1,  -1,  +0,  +0,  +0;
    +0,  +0,  +1,  -1,  +0;
    +0,  +0,  -1,  +0,  -1;
    ];
T = @(X)propensities_2state(X,p);
opts = odeset('Events',@(t,y)eventfun(t,y,jump,T));
[t,y]=ode15s(@(t,y) jump'*T(y),[0,TF], X0,opts);

ssTR=T(y(end,:)');
ssRT=[];
end

% Propensity rate vector (CTMC)
function Rate = propensities_2state(X, p)
Rate = [p.MU(1)*X(1);
    p.MU(2)*min(X(2),X(3));
    p.MU(4)*X(4);
    p.MU(5)*X(5);];
Rate(isnan(Rate))=0;
end

function [x,isterm,dir] = eventfun(t,y,jump,T)
dy = jump'*T(y);
x = norm(dy) - 1e-6;
isterm = 1;
dir = 0;
end
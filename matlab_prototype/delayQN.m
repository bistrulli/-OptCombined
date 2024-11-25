function X = delayQN(X0,MU,NC,TF,rep,dt)
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
stoich_matrix=[-1,  +1,  +0,  +0,  +0;
    +1,  -1,  +0,  +0,  +0;
    +0,  +0,  +1,  -1,  +0;
    +0,  +0,  -1,  +0,  -1;
    ];

tspan = [0, TF];
pfun = @propensities_2state;

X = zeros(length(X0), ceil(TF/dt) + 1, rep);
for i = 1:rep
    [t, x] = firstReactionMethod(stoich_matrix, pfun, tspan, X0, p,[],10000000);
    tsin = timeseries(x,t);
    tsout = resample(tsin, linspace(0, TF, ceil(TF/dt)+1), 'zoh');
    X(:, :, i) = tsout.Data';
end
end

% Propensity rate vector (CTMC)
function Rate = propensities_2state(X, p)
Rate = [p.MU(1)*X(1);
    p.MU(2)*min(X(2),X(3));
    p.MU(4)*X(4);
    p.MU(5)*X(5);];
Rate(isnan(Rate))=0;
end
%% Authorship
% David E. Acuña-Ureta, Ph.D.
% Assistant Professor
% Department of Mechanical and Metallurgical Engineering
% Pontificia Universidad Católica de Chile
% E-mail: david.acuna@ing.puc.cl

%% Definitions
% Table 1: Simulation parameters
Ts = 1.0;           % sample time
Ecrit = 1389900;    % battery storage capacity
sigma_w = 1e-6;     % standard deviation process noise
sigma_eta = 0.9;    % (set 0.9 or 0.1) standard deviation measurement noise
v0 = 41.405;        % open-circuit voltage model parameter
vL = 33.481;        % open-circuit voltage model parameter
alpha = 5.319e-3;   % open-circuit voltage model parameter
beta = 11.505;      % open-circuit voltage model parameter
gamma = 1.5538;     % open-circuit voltage model parameter
R = 0.26;           % internal resistance

% state-space functions
Voc = @(x) vL + (v0-vL)*exp((x-1)*gamma) + alpha*vL*(x-1) + ...
    (1-alpha)*vL*(exp(-beta)-exp(-beta*sqrt(x)));               % Eq. (3): open-circuit voltage model
V = @(x,u) Voc(x) - R*u + randn*sigma_eta;                      % Eq. (2): measurement equation
state_trans = @(x,u) x - u*V(x,u)*(Ts/Ecrit) + randn*sigma_w;   % Eq. (1): state transition equation

% Exogenous inputs (current discharge) modeled as a Gaussian mixture
w = [0.3, 0.5, 0.2];	% Eqs. (43)-(45): weights of the mixture
u = [3, 6, 9];          % Eqs. (40)-(42): expected values
sig = 1;              	% standard deviation (same for all Gaussians)

% Threshold
v_cutoff = 33;  % cut-off voltage

% Uncertain event likelihood function (a.k.a. "uncertain hazard zone")
PEv = @(v,u,sigma) 1 - (1/2)*(1+erf(( v - (v_cutoff+R*u) )/(sqrt(2)*sigma))); % uncertain hazard zone

% Other simulation parameters
kp = 800;   % time at which prognostics begin
kh = 5000;  % time at which prognostics finish
x0 = 0.8;   % State-of-Charge at time 'kp'

%% Method 1 (Standard): Certain event approach
N = 1e3;                                                    % number of Monte Carlo simulations
histogram_bins = zeros(1,kh-kp+1);                          % count of frequency at each histogram bin
start_time = cputime;                                       % starting time
for traj=1:N                                                % simulation of trajectories
    X = x0;                                                 % set initial state
    U = u(randsrc(1,1,[(1:length(u));w])) + randn*sig;      % set initial exogenous input
    for k=2:kh-kp+1
        X = state_trans(X,U);                               % sample x_{k+1} ~ p(x_{k+1}|x_k,u_k)
        U = u(randsrc(1,1,[(1:length(u));w])) + randn*sig;	% sample u_{k+1} ~ p(u_{k+1})
        if (V(X,U)<=v_cutoff)                               % if cut-off voltage threshold is hit
            histogram_bins(k) = histogram_bins(k) + 1;      % sum contribution to histogram
            break;                                          % stop simulating this trajectory
        end
    end
end
computation_time = cputime - start_time;                % Eq. (34): computation time
if sigma_eta == 0.9                                     % load ground truth probability distribution
    load('ground_truth_eta09.mat');	
elseif sigma_eta == 0.1
    load('ground_truth_eta01.mat');
end
EoD_time_PMF = histogram_bins./sum(histogram_bins);    	% Eq. (22): EoD time probability distribution
l1_distance = norm(EoD_ground_truth-EoD_time_PMF,1);    % Eq. (33): l1-distance to ground truth

% plot results
figure
hold on
h1 = area(kp:kh,EoD_ground_truth);
h1.FaceColor = [0.6350, 0.0780, 0.1840];
h1.FaceAlpha = 0.4;
h1.EdgeColor = [0.6350, 0.0780, 0.1840];
h2 = area(kp:kh,EoD_time_PMF);
h2.FaceColor = [0, 0.4470, 0.7410];
h2.FaceAlpha = 0.4;
h2.EdgeColor = [0, 0.4470, 0.7410];
xlabel('Time (s)')
ylabel('Probability')
legend('Ground truth',sprintf('Method 1: N=%d',N),'location','northwest')
if sigma_eta == 0.9
    axis([800 3300 0 2.5e-3])
elseif sigma_eta == 0.1
    axis([3400 4900 0 4e-3])
end
fprintf('l1-distance (efficacy): %f\n',l1_distance)                       % display efficacy
fprintf('Computation time (efficiency): %f seconds\n',computation_time)   % display efficiency

%% Method 2 (Proposed): Uncertain event approach
N = 1;                                              % number of Monte Carlo simulations
histogram_bins = zeros(1,kh-kp+1);                	% count of frequency at each histogram bin
start_time = cputime;                               % starting time
rho_min = 0e-7;                                     % parameter criteria for stopping simulations
for traj=1:N                                        % simulation of trajectories
    X = x0;                                                             % set initial state
    W = 1;                                                              % product in Eq. (26)
    U = u(randsrc(1,1,[(1:length(u));w])) + randn*sig;                  % set initial exogenous input
    for k=2:kh-kp+1
        X = state_trans(X,U);                                           % sample x_{k+1} ~ p(x_{k+1}|x_k,u_k)
        U = u(randsrc(1,1,[(1:length(u));w])) + randn*sig;              % sample u_{k+1} ~ p(u_{k+1})
        Expected_Voc = w(1)*PEv(Voc(X),u(1),sqrt(sigma_eta^2+R^2*sig^2)) + ...
            w(2)*PEv(Voc(X),u(2),sqrt(sigma_eta^2+R^2*sig^2)) + ...
            w(3)*PEv(Voc(X),u(3),sqrt(sigma_eta^2+R^2*sig^2));          % Eq. (51)
        histogram_bins(k) = histogram_bins(k) + Expected_Voc*W;         % sum contribution to histogram
        W = W*(1-PEv(Voc(X),U,sigma_eta));                              % update product in Eq. (26)
        if W < rho_min                                                  % if the product is negligible
            break;                                                      % stop simulating this trajectory
        end 
    end
end
computation_time = cputime - start_time;                % Eq. (34): computation time
if sigma_eta == 0.9                                     % load ground truth probability distribution
    load('ground_truth_eta09.mat');	
elseif sigma_eta == 0.1
    load('ground_truth_eta01.mat');
end
EoD_time_PMF = histogram_bins./sum(histogram_bins);     % Eq. (22): EoD time probability distribution
l1_distance = norm(EoD_ground_truth-EoD_time_PMF,1);    % Eq. (33): l1-distance to ground truth

% plot results
figure
hold on
h1 = area(kp:kh,EoD_ground_truth);
h1.FaceColor = [0.6350, 0.0780, 0.1840];
h1.FaceAlpha = 0.4;
h1.EdgeColor = [0.6350, 0.0780, 0.1840];
h2 = area(kp:kh,EoD_time_PMF);
h2.FaceColor = [0, 0.5, 0];
h2.FaceAlpha = 0.4;
h2.EdgeColor = [0, 0.5, 0];
xlabel('Time (s)')
ylabel('Probability')
legend('Ground truth',sprintf('Method 2: N=%d',N),'location','northwest')
if sigma_eta == 0.9
    axis([900 3300 0 1.4e-3])
elseif sigma_eta == 0.1
    axis([3500 4800 0 4.0e-3])
end
fprintf('l1-distance (efficacy): %.4f\n',l1_distance)                       % display efficacy
fprintf('Computation time (efficiency): %.4f seconds\n',computation_time)   % display efficiency
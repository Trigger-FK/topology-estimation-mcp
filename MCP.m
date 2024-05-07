close all;
clear;

%% Load data
sample_duration = 1000;
load(sprintf('Data/%ddata.mat', sample_duration));
if sample_duration == 1000
    lambda = 0.8;
    gamma = 0.1;
elseif sample_duration == 200
    lambda = 1.8;
    gamma = 0.1;
elseif sample_duration == 100
    lambda = 0.80;
    gamma = 0.05;
end

%% Set initial state and stopping criterion in DCA
q0 = ones([m, 1])/n;
max_iter = 100;
tol = 1e-6; 

%% make the correct weight vector
omega = [];
for i = 1:n-1
    omega = [omega; abs(A(i, i+1:n).')];
end

%% Define the convex function
g = @(u) -log(det(E*diag(u)*E.' + e*e.'/n)) + 2/(sigma^2)*trace(S*(E*diag(u)*E.')) + lambda*norm(u,1);
h = @(u) lambda*ME(u, gamma);

%% DCA
q_pre2 = zeros([m,1]);
q_pre = q0;
tic
for k = 1:max_iter
    cvx_begin
        variable q(m)

        g_q = -log_det(E*diag(q)*E' + e*e'/n) + 2/(sigma^2)*trace(S*E*diag(q)*E') + lambda*norm(q,1);
        dh = lambda*dME(q_pre, gamma);

        minimize (g_q - dh'*q)
        subject to
            lambda_min(E*diag(q)*E.' + e*e.'/n) > 0; % positive definite matrix
            q >= 0;
    cvx_end

    if norm((g(q)-h(q))-(g(q_pre)-h(q_pre)),2) < tol || norm((g(q)-h(q))-(g(q_pre2)-h(q_pre2)),2) < tol
        break;
    end
    q_pre2 = q_pre;
    q_pre = q;
end
toc

%% Display the result
A_vec = [];
edge = [];
count = 0;
for i = 1:n
    for j=i+1:n
        count = count+1;
        edge = [edge, abs(A(i,j))];
        if abs(A(i,j)) > 0
            A_vec = [A_vec, count];
        end
    end
end

disp(q)
q_active = (abs(q)>=0.1);
find(q_active)'
A_vec

L2_error = norm((g(q)-h(q))-(g(omega)-h(omega)),2)
relative_error = norm(omega-q,2)/norm(omega,2)

figure;
stem(q, 'LineWidth',4);
hold on
stem(edge, 'LineWidth',4);
grid on
xlabel('edge');
xticks(1:m)
xlim([0 m+1])
ylim([-0.1 2.2])
legend('estimate', 'true', 'Location', 'northeastoutside');
title('solution q');

%% Save the result
if not(exist(sprintf("EstimateResult/MCP_%d",sample_duration),'dir'))
    mkdir(sprintf("EstimateResult/MCP_%d",sample_duration))
end

date=datetime('now','Format','MM-dd HH-mm-ss');
save(sprintf("EstimateResult/MCP_%d/MCP_%d %s",sample_duration,sample_duration,date),'omega','gamma','lambda','edge','q',"L2_error","relative_error")

%% ----------The Moreau Envelope----------
function y = ME(q, gamma)
    y = 0;
    for i = 1:length(q)
        if q(i) >= gamma
            y = y + (q(i) - gamma/2);
        elseif 0 <= q(i) && q(i) < gamma
            y = y + (power(q(i),2) / (2*gamma));
        end
    end
end

%% ----------The Gradient of Moreau Envelope----------
function y = dME(q, gamma)
    y = zeros(size(q));
    for i = 1:length(q)
        if q(i) >= gamma
            y(i) = 1;
        elseif 0 <= q(i) && q(i) < gamma
            y(i) = q(i)/gamma;
        end
    end
end

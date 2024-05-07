close all;
clear;

%% Load data
sample_duration = 1000;
load(sprintf('Data/%ddata.mat',sample_duration));

%% Set parameters
if sample_duration == 1000
    a = 0;
elseif sample_duration == 200
    a = 0.01;
elseif sample_duration == 100
    a = 0.01 ;
end

%% Set the correct weight vector
omega = [];
for i = 1:n-1
    omega = [omega; abs(A(i,i+1:n).')];
end

%% L1
tic
cvx_begin
    variable q(m)
    J = -log_det(E*diag(q)*E' + e*e'/n) + 2*trace(S*E*diag(q)*E')/sigma^2;
    g = a*norm(q,1);

    minimize (J + g)
    subject to
        lambda_min(E*diag(q)*E' + e*e.'/n) > 0; 
        q >= 0
cvx_end
toc

%% Define the convex function
J = @(u) -log(det(E*diag(u)*E.' + e*e.'/n)) + 2*trace(S*(E*diag(u)*E.'))/sigma^2;
g = @(u) a*norm(u,1);

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

L2_error = norm((J(q) + g(q))-(J(omega) + g(omega)),2)
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
if not(exist(sprintf("EstimateResult/L1_%d",sample_duration),'dir'))
    mkdir(sprintf("EstimateResult/L1_%d",sample_duration))
end

date=datetime('now','Format','MM-dd HH-mm-ss');
save(sprintf("EstimateResult/L1_%d/L1_%d %s",sample_duration,sample_duration,date),'omega','a','edge','q',"L2_error","relative_error")
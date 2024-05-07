close all;
clear;

%% Set sample duration
sample_duration = 1000;

%% Set initial parameter
rng(1);
n = 100;
m = n*(n-1)/2;
t_max = 50;
h = 0.001;
t = 0:h:t_max;
sigma = 0.10;

%% Define matrices
I = eye(n);
e = ones(n,1);

A = GridAdjMatrix(n);
D = diag(sum(A, 2));
L = D - A;

E = zeros([n m]);
num = 1;
for i = 1:n
    for j = i+1:n
        temp = zeros([n 1]);
        temp(i) = 1;
        temp(j) = -1;
        E(:,num) = temp.';
        num = num + 1;
    end
end

%% Set buffer values
S_all = zeros(n,n);
mu_row = n;
mu_col = length(t);

realization_num = 100;
for realization = 1:realization_num
    mu = randn(n, 1);
    mu_graph = zeros(mu_row, mu_col);
    mu_graph(:, 1) = mu;
    mu_tild = zeros(mu_row, mu_col);
    mu_tild(:, 1) = (I - e*e'/n) * mu;
    
    %% Control sequence
    for k = 1:mu_col-1
        mu = mu + h*(-L*mu) + sqrt(h)*sigma*randn(n, 1);
        mu_graph(:, k+1) = mu;
        mu_tild(:, k+1) = (I - ones(n,1)*ones(n,1)'/n) * mu;
    end
    
    %% Calculate the sample covariance matrix
    sample_data = mu_tild(:, end-sample_duration+1:end);
    S_all = S_all + sample_data*sample_data'/sample_duration;
end
S = S_all / realization_num;

%% Save the variables
if not(exist("Data",'dir'))
    mkdir("Data")
end
save(sprintf('Data/%ddata.mat', sample_duration))

%% Define the GridAdjMatrix Function
function A = GridAdjMatrix(N)
    n = sqrt(N);
    if mod(n, 1) ~= 0
        error('the square root is not natural number: sqrt(%d) = %f', N, n);
    end
    a = zeros(N,N);
    for i=1:N
        if mod(i, n) ~= 0 && i+1 <= N
            a(i,i+1) = 1;
        end
        if i+n <= N
            a(i,i+n) = 1;
        end
    end
    weight = 0.01 + (2 - 0.01) * rand(N,N);
    a = weight .* a;
    A = a + a';
end

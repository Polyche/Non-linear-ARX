load('iddata-07.mat');
u = id.u;
y = id.y;
Ts = id.Ts;
u_val = val.u;
y_val = val.y;
N_val = length(y_val);
na_max = 4;
nb_max = 4;
m_max  = 3;
nk = 1;
N = length(y);
results = [];
for na = 1:na_max
    for nb = 1:nb_max
        for m = 1:m_max
            Phi = [];
            for k = max(na, nb + nk):N
                delayed_y = zeros(1, na);
                for j = 1:na
                    if (k - j) > 0
                        delayed_y(j) = y(k - j);
                    end
                end
                delayed_u = zeros(1, nb);
                for j = 1:nb
                    if (k - nk - j + 1) > 0
                        delayed_u(j) = u(k - nk - j + 1);
                    end
                end
                delayed_vars = [delayed_y, delayed_u];
                poly_terms = [1, monomial_terms(delayed_vars, m)];
                Phi = [Phi; poly_terms];
            end
            y_regress = y(max(na, nb + nk):end);
            theta = Phi \ y_regress;
            y_pred = zeros(N,1);
            for k = max(na, nb + nk):N
                delayed_y = zeros(1, na);
                for j = 1:na
                    if (k - j) > 0
                        delayed_y(j) = y(k - j);
                    end
                end
                delayed_u = zeros(1, nb);
                for j = 1:nb
                    if (k - nk - j + 1) > 0
                        delayed_u(j) = u(k - nk - j + 1);
                    end
                end
                delayed_vars = [delayed_y, delayed_u];
                poly_terms = [1, monomial_terms(delayed_vars, m)];
                y_pred(k) = poly_terms * theta;
            end
            y_sim = zeros(N,1);
            for k = max(na, nb + nk):N
                delayed_y = zeros(1, na);
                for j = 1:na
                    if (k - j) > 0
                        delayed_y(j) = y_sim(k - j);
                    end
                end
                delayed_u = zeros(1, nb);
                for j = 1:nb
                    if (k - nk - j + 1) > 0
                        delayed_u(j) = u(k - nk - j + 1);
                    end
                end
                delayed_vars = [delayed_y, delayed_u];
                poly_terms = [1, monomial_terms(delayed_vars, m)];
                y_sim(k) = poly_terms * theta;
            end
            mse_pred = mean((y - y_pred).^2);
            mse_sim = mean((y - y_sim).^2);
            y_pred_val = zeros(N_val,1);
            for k = max(na, nb + nk):N_val
                delayed_y = zeros(1, na);
                for j = 1:na
                    if (k - j) > 0
                        delayed_y(j) = y_val(k - j);
                    end
                end
                delayed_u = zeros(1, nb);
                for j = 1:nb
                    if (k - nk - j + 1) > 0
                        delayed_u(j) = u_val(k - nk - j + 1);
                    end
                end
                delayed_vars = [delayed_y, delayed_u];
                poly_terms = [1, monomial_terms(delayed_vars, m)];
                y_pred_val(k) = poly_terms * theta;
            end
            y_sim_val = zeros(N_val,1);
            for k = max(na, nb + nk):N_val
                delayed_y = zeros(1, na);
                for j = 1:na
                    if (k - j) > 0
                        delayed_y(j) = y_sim_val(k - j);
                    end
                end
                delayed_u = zeros(1, nb);
                for j = 1:nb
                    if (k - nk - j + 1) > 0
                        delayed_u(j) = u_val(k - nk - j + 1);
                    end
                end
                delayed_vars = [delayed_y, delayed_u];
                poly_terms = [1, monomial_terms(delayed_vars, m)];
                y_sim_val(k) = poly_terms * theta;
            end
            mse_pred_val = mean((y_val - y_pred_val).^2);
            mse_sim_val = mean((y_val - y_sim_val).^2);
            results = [results; na, nb, m, mse_pred, mse_sim, mse_pred_val, mse_sim_val];
        end
    end
end
fig = uifigure('Name', 'MSE Results', 'Position', [100, 100, 1000, 400]);
uitable(fig, 'Data', results, ...
    'ColumnName', {'na','nb','m','MSE y_{pred} (ID)','MSE y_{sim}  (ID)','MSE y_{pred} (VAL)','MSE y_{sim}  (VAL)'}, ...
    'Position', [25, 50, 950, 300], 'FontSize', 12);
na = 2; nb = 1; m = 3;
Phi = [];
for k = max(na, nb + nk):N
    delayed_y = zeros(1, na);
    for j = 1:na
        if (k - j) > 0
            delayed_y(j) = y(k - j);
        end
    end
    delayed_u = zeros(1, nb);
    for j = 1:nb
        if (k - nk - j + 1) > 0
            delayed_u(j) = u(k - nk - j + 1);
        end
    end
    delayed_vars = [delayed_y, delayed_u];
    poly_terms = [1, monomial_terms(delayed_vars, m)];
    Phi = [Phi; poly_terms];
end
y_regress = y(max(na, nb + nk):end);
theta = Phi \ y_regress;
y_pred = zeros(N,1);
for k = max(na, nb + nk):N
    delayed_y = zeros(1, na);
    for j = 1:na
        if (k - j) > 0
            delayed_y(j) = y(k - j);
        end
    end
    delayed_u = zeros(1, nb);
    for j = 1:nb
        if (k - nk - j + 1) > 0
            delayed_u(j) = u(k - nk - j + 1);
        end
    end
    delayed_vars = [delayed_y, delayed_u];
    poly_terms = [1, monomial_terms(delayed_vars, m)];
    y_pred(k) = poly_terms * theta;
end
y_sim = zeros(N,1);
for k = max(na, nb + nk):N
    delayed_y = zeros(1, na);
    for j = 1:na
        if (k - j) > 0
            delayed_y(j) = y_sim(k - j);
        end
    end
    delayed_u = zeros(1, nb);
    for j = 1:nb
        if (k - nk - j + 1) > 0
            delayed_u(j) = u(k - nk - j + 1);
        end
    end
    delayed_vars = [delayed_y, delayed_u];
    poly_terms = [1, monomial_terms(delayed_vars, m)];
    y_sim(k) = poly_terms * theta;
end
u_val = val.u;
y_val = val.y;
N_val = length(y_val);
y_pred_val = zeros(N_val,1);
for k = max(na, nb + nk):N_val
    delayed_y = zeros(1, na);
    for j = 1:na
        if (k - j) > 0
            delayed_y(j) = y_val(k - j);
        end
    end
    delayed_u = zeros(1, nb);
    for j = 1:nb
        if (k - nk - j + 1) > 0
            delayed_u(j) = u_val(k - nk - j + 1);
        end
    end
    delayed_vars = [delayed_y, delayed_u];
    poly_terms = [1, monomial_terms(delayed_vars, m)];
    y_pred_val(k) = poly_terms * theta;
end
y_sim_val = zeros(N_val,1);
for k = max(na, nb + nk):N_val
    delayed_y = zeros(1, na);
    for j = 1:na
        if (k - j) > 0
            delayed_y(j) = y_sim_val(k - j);
        end
    end
    delayed_u = zeros(1, nb);
    for j = 1:nb
        if (k - nk - j + 1) > 0
            delayed_u(j) = u_val(k - nk - j + 1);
        end
    end
    delayed_vars = [delayed_y, delayed_u];
    poly_terms = [1, monomial_terms(delayed_vars, m)];
    y_sim_val(k) = poly_terms * theta;
end
figure('Name','Identification Results','NumberTitle','off');
plot(y,'b','LineWidth',1.2); hold on; plot(y_pred,'r','LineWidth',1.2);
plot(y_sim,'g','LineWidth',1.2);
legend('Actual Output','One-Step-Ahead Prediction','Simulation','Location','Best');
title(sprintf('Identification Dataset: na=%d, nb=%d, nk=%d, m=%d',na,nb,nk,m));
xlabel('Sample'); ylabel('Output');
figure('Name','Validation Results','NumberTitle','off');
plot(y_val,'b','LineWidth',1.2); hold on; plot(y_pred_val,'r','LineWidth',1.2);
plot(y_sim_val,'g','LineWidth',1.2);
legend('Actual Output','One-Step-Ahead Prediction','Simulation','Location','Best');
title(sprintf('Validation Dataset: na=%d, nb=%d, nk=%d, m=%d',na,nb,nk,m));
xlabel('Sample'); ylabel('Output');

function terms = monomial_terms(vars, degree)
num_vars = length(vars);
terms = [];
for d = 1:degree
    for i = 1:num_vars
        terms = [terms, vars(i)^d];
    end
end
end
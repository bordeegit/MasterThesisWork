close all

options = optimoptions('fsolve');
options.MaxFunctionEvaluations = 2000;
options.Display = 'off';

N_opt = 1000;

W_vec = zeros(N_opt,3);

for ind = 1:N_opt
        
    W_meas               = W_log.signals.values(ind,:)';
    
    x_meas               = states.signals.values(ind,:)';
    x_meas(1)            = pi/2 - x_meas(1);
    x_meas(3)            = -x_meas(3);
    Psi_meas             = Psi.signals.values(ind);
    Fline_meas           = [0; 0; Forces.signals.values(ind,end)];
    pm.A                 = area;
    pm.m                 = mass;
    pm.rho               = rho;
    pm.g                 = g;
    pm.n_l               = n_line;
    pm.d_l               = Line_diameter;
    pm.CD_l              = CD_Line;
    pm.rho_l             = Line_density;
    pm.alpha_0           = alpha_0;
    pm.alpha_max         = alpha_max;
    pm.alpha_min         = alpha_min;
    pm.alpha_var         = alpha_var;
    pm.Cl_var            = CL_var;
    pm.Cd_var            = CD_var;

    xdot_meas            = statesdot.signals.values(ind,:)';
    
    fun = @(W)accDiff(x_meas,Psi_meas,Fline_meas,W,pm, xdot_meas);
    W0 = W_meas; 
    %W0 = [0;0;0];
    W = fsolve(fun, W0, options);
    W_vec(ind,:) = W'; 
    
    fprintf("Iter %d, Estimated Wind is [%.2f,%.2f,%.2f]\n\t\t  Actual Wind is [%.2f,%.2f,%.2f]\n", ...
        ind, W(1), W(2), W(3), W_meas(1),W_meas(2),W_meas(3));

end

figure(1)
plot(1:N_opt, W_log.signals.values(1:N_opt,1), 1:N_opt, W_vec(:,1));
figure(2)
plot(1:N_opt, W_log.signals.values(1:N_opt,2), 1:N_opt, W_vec(:,2));
figure(3)
plot(1:N_opt, W_log.signals.values(1:N_opt,3), 1:N_opt, W_vec(:,3));


function F = accDiff(x_meas,Psi_meas,Fline_meas,W,pm, xdot_meas)
    
    xdot = kiteEq(x_meas,Psi_meas,Fline_meas,W,pm);
    F = xdot_meas(1:3) - xdot(1:3);

end

close all

options = optimoptions('fsolve');
options.MaxFunctionEvaluations = 6000;
options.MaxIterations = 2000;
options.Display = 'off';
options.Algorithm = 'levenberg-marquardt';

N_opt = 500;

W0 = W_log.signals.values(1,:)';
W_vec = [zeros(N_opt,3)];

x_noiseLvl = 0;
Psi_noiseLvl = 0;
Fline_noiseLvl = 0;

% Generating Parameters

W_meas               = W_log.signals.values(1,:)';
    
x_meas               = states.signals.values(1,:)';
x_meas(1)            = pi/2 - x_meas(1);
x_meas(3)            = -x_meas(3);
x_meas               = x_meas + x_noiseLvl*x_meas.*randn(size(x_meas)); %adding noise to measurement
Psi_meas             = Psi.signals.values(1);
Psi_meas             = Psi_meas + Psi_noiseLvl*Psi_meas.*randn(size(Psi_meas));
Fline_meas           = [0; 0; Forces.signals.values(1,end)];
Fline_meas           = Fline_meas + Fline_noiseLvl*Fline_meas.*randn(size(Fline_meas));
xdot_meas            = statesdot.signals.values(1,:)';

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

%codegen kiteEq -args {x_meas,Psi_meas,Fline_meas,W_meas,pm} -lang:c++

tic
for ind = 2:N_opt
    
    fun = @(W)accDiff(x_meas,Psi_meas,Fline_meas,W,pm, xdot_meas);
    [W,~,exit,out] = fsolve(fun, W0, options);
    W_vec(ind-1,:) = W'; 

    fprintf("Step %d, Estimated Wind is [%-2.2f,%-2.2f,%-2.2f], Actual Wind is [%-2.2f,%-2.2f,%-2.2f]", ...
        ind-1, W(1), W(2), W(3), W_meas(1),W_meas(2),W_meas(3));

    if exit == 0
        fprintf(" Exceeded, iter: %d, fncCount: %d", out.iterations, out.funcCount);
    end
    if exit == -2 || exit == -3
        fprintf(" Eq not solved, exit: %d, msg: \n %s \n", exit, out.message);
    end
    fprintf("\n");
    
    W0 = W;   % soft start

    W_meas               = W_log.signals.values(ind,:)';
    
    x_meas               = states.signals.values(ind,:)';
    x_meas(1)            = pi/2 - x_meas(1);
    x_meas(3)            = -x_meas(3);
    x_meas               = x_meas + x_noiseLvl*x_meas.*randn(size(x_meas)); %adding noise to measurement
    Psi_meas             = Psi.signals.values(ind);
    Psi_meas             = Psi_meas + Psi_noiseLvl*Psi_meas.*randn(size(Psi_meas));
    Fline_meas           = [0; 0; Forces.signals.values(ind,end)];
    Fline_meas           = Fline_meas + Fline_noiseLvl*Fline_meas.*randn(size(Fline_meas));
    xdot_meas            = statesdot.signals.values(ind,:)';
    
end
toc

figure(1)
plot(1:N_opt, W_log.signals.values(1:N_opt,1), 1:N_opt, W_vec(:,1));
figure(2)
plot(1:N_opt, W_log.signals.values(1:N_opt,2), 1:N_opt, W_vec(:,2));
figure(3)
plot(1:N_opt, W_log.signals.values(1:N_opt,3), 1:N_opt, W_vec(:,3));


function F = accDiff(x_meas,Psi_meas,Fline_meas,W,pm, xdot_meas)
    
    xdot = kiteEq_mex(x_meas,Psi_meas,Fline_meas,W,pm);
    F = xdot_meas - xdot;

end

%% ASEN 3128 - Lab 5 - Main
% 
%
% Author: Cole MacPherson
% Collaborators: R. Briones, E. Cervantes, A. Filsoof 
% Date: 26th March 2021

%% Housekeeping
clc;
clear;
close all;
tic

%% Part 1
% 1.1 (Change Case II 747 tables to SI units)

    g = 9.81; % gravity [m/s^2]
    
    % 747-100 data
    S = convlength(5500,'ft','m'); % wing surface area [m^2]
    b = 195.68/3.281; % span [m]
    c_bar = 27.31/3.281; % average chord [m]
    CG_pos = 0.25; % CG location [%c]
    altitude = 20000/3.281; % altitude [m]
    M = 0.5; % mach number
    V = 518/3.281; % velocity [m/s]
    W = (6.366e5)*4.448; % weight [N]
    m = W/g; % mass [kg]
    I_x = 1.82e7*1.3558179619; % moment of interia about the x-axis [kgm^2]
    I_y = 3.31e7*1.3558179619; % moment of interia about the y-axis [kgm^2]
    I_z = 4.97e7*1.3558179619; % moment of interia about the z-axis [kgm^2]
    I_zx = 9.7e5*1.3558179619; % moment of inertia about the z-x [kgm^2]
    gamma = 1.4; % specific heat of air
    R = 287.05; % gas constant [J/kgK]
    theta_0 = 0; % [rad]
    [T_inf, a, P, rho] = atmoscoesa(altitude);
    zeta = deg2rad(-6.8); % [rad]
    C_L_0 = 0.654; % 
    C_D_0 = 0; %
    C_w_0 = (2*W)/(rho*V^2*S); % coefficient of weight
    C_D = 0.04; % coefficient of drag

    % dimensional derivatives
    Xu = (-4.883e1)*4.448*3.281; % [Ns/m]
    Xw = (1.546e3)*4.448*3.281; % [Ns/m]
    Xq = 0; %[Ns/m]
    Xw_dot = 0; % [Ns/ft]
    Xdelta_e = (3.994e4)*4.448; % [N/rad]
    Zu = (-1.342e3)*4.448*3.281; % [Ns/m]
    Zw = (-8.561e3)*4.448*3.281; % [Ns/m]
    Zq = (-1.263e5)*4.448; % [Ns/rad]
    Zw_dot = (3.104e2)*4.448*3.281; % [N(s^2)/m]
    Zdelta_e = (-3.341e5)*4.448; % [N/rad]
    Mu = (8.176e3)*4.448; % [Ns]
    Mw = (-5.627e4)*4.448; % [Ns}
    Mq = (-1.394e7)*(4.448/3.281); % [Nsm/rad]
    Mw_dot = (-4.138e3)*4.448; % [Ns^2]
    Mdelta_e = (-3.608e7)*(4.448/3.281); % [Nm/rad]
    
% 1.2
    % 1.2a (convert derivative table into stability frame)
    SM = stabFrameMatrix(Xu, Zu, Mu, Xw, Zw, Mw, Xq, Zq, Mq, Zw_dot, Mw_dot, zeta);
    
    % 1.2b (non-dimensionalized derivatives)
    ND = [(SM(1,1)-(rho*V*S*C_w_0*sin(theta_0)))/(0.5*rho*V*S) (SM(1,2)+(rho*V*S*C_w_0*cos(theta_0)))/(0.5*rho*V*S) SM(1,3)/(0.5*rho*V*c_bar*S);...
        SM(2,1)/(0.5*rho*V*S) SM(2,2)/(0.5*rho*V*S) SM(2,3)/(0.5*rho*V*c_bar*S);...
        SM(3,1)/(0.25*rho*V*c_bar*S) SM(3,2)/(0.25*rho*V*c_bar*S) SM(3,3)/(0.25*rho*V*(c_bar^2)*S);...
        SM(4,1)/(0.25*rho*c_bar*S) SM(4,2)/(0.25*rho*c_bar*S) SM(4,3)/(0.25*rho*(c_bar^2)*S)];
    
% 1.3 (find the A matrix)
    A = nonDimensionalize(SM,m,g,V,I_y);
    
% 1.4 (find eigenvectors/values and the damping ratio and natural frequencies)
    fprintf('Part 1 Problem 4:\n\t');
    [eigvecs,eigvals] = eig(A);
    w_n1 = sqrt(real(eigvals(3,3))^2+imag(eigvals(3,3))^2);
    zeta1 = -real(eigvals(3,3))/w_n1;
    
    fprintf(['Mode 1 (Phugoid mode):\n\t\tEigenvalue: %f ' char(177) ' %fi\n'],real(eigvals(3,3)),imag(eigvals(3,3)));
    fprintf(['\t\tNatural Frequency: ' char(177) ' %f [Hz]\n\t\tDamping Ratio: ' char(177) ' %f\n'],abs(w_n1),abs(zeta1));
    
    w_n2 = sqrt(real(eigvals(1,1))^2+imag(eigvals(1,1))^2);
    zeta2 = -real(eigvals(1,1))/w_n2;
    
    fprintf(['\tMode 2 (Short-period mode):\n\t\tEigenvalue: %f ' char(177) ' %fi\n'],real(eigvals(1,1)),imag(eigvals(1,1)));
    fprintf(['\t\tNatural Frequency: ' char(177) ' %f [Hz]\n\t\tDamping Ratio: ' char(177) ' %f\n\n'],abs(w_n2),abs(zeta2));
    
% 1.5 (Compare above results with Lanchester approximation)
    fprintf('Part 1 Problem 5:\t');
    
    B = [Mq/I_y (V*((Mw*cos(zeta))))/I_y;...
        1 0];
    eigval3 = eig(B);
    
    fprintf(['\n\tShort period mode approximation:\n\t\tEigenvalue: %f ' char(177) ' %fi\n'],real(eigval3(1)),imag(eigval3(1)));
    
    T_phugoid = (2*pi)/w_n1;
    T_lanchester = pi*sqrt(2)*(V/g);
    
    fprintf('\n\tThe oscillation periods: \n\t\tPhugiod Mode: %f [sec]\n\t\tLanchester approximation: %f [sec]\n\n',T_phugoid,T_lanchester);
    
% 1.6 (Simulate the linearized longitudinal dynamics)
    % 1.6.a (verify trim state)
    state_vec = [0;0;0;0];
    t_span = [0 1000];
    [t_6a, state_vec_6a] = ode45(@(t_6a, state_vec_6a) LLDodefunc(t_6a,state_vec_6a,A), t_span, state_vec);
    
    plotLLD(t_6a,state_vec_6a,'trimmed','6.a'); % plot results
    
    % 1.6.b.i (Delta u = 10 [m/s])
    state_vec = [10;0;0;0];
    [t_6bi, state_vec_6bi] = ode45(@(t_6bi, state_vec_6bi) LLDodefunc(t_6bi,state_vec_6bi,A), t_span, state_vec);
    
    plotLLD(t_6bi,state_vec_6bi,'\Deltau^{E} = 10 [m/s]','6.b.i'); % plot results
    
    % 1.6.b.ii (Delta w = 10 [m/s])
    state_vec = [0;10;0;0];
    [t_6bii, state_vec_6bii] = ode45(@(t_6bii, state_vec_6bii) LLDodefunc(t_6bii,state_vec_6bii,A), t_span, state_vec);
    
    plotLLD(t_6bii,state_vec_6bii,'\Delta\omega^{E} = 10 [m/s]','6.b.ii'); % plot results
    
    % 1.6.b.iii (Delta q = 0.1 [rad/s])
    state_vec = [0;0;0.1;0];
    [t_6biii, state_vec_6biii] = ode45(@(t_6biii, state_vec_6biii) LLDodefunc(t_6biii,state_vec_6biii,A), t_span, state_vec);
    
    plotLLD(t_6biii,state_vec_6biii,'\Deltaq = 0.1 [rad/s]','6.b.iii'); % plot results
    
    % 1.6.b.iv (Delta theta = 0.1 [rad])
    state_vec = [0;0;0;0.1];
    [t_6biv, state_vec_6biv] = ode45(@(t_6biv, state_vec_6biv) LLDodefunc(t_6biv,state_vec_6biv,A), t_span, state_vec);
    
    plotLLD(t_6biv,state_vec_6biv,'\Delta\theta = 0.1 [rad]','6.b.iv'); % plot results
    
    % 1.6.c (Which initial deviations are good at exciting the Short Period mode? The Phugoid?)
    %{
        The pitch rate and y velocity deviations are good at exciting the
        Short Period mode and the picth and x velocity are better at
        exciting the Phugoid mode
    %}
     
%% Part 2
% 2.7 (Add elevator control derivatives)
    fprintf('Part 2 Problem 7:\n')
    
    X_deltae = (-3.818e-6)*(1/2)*rho*V^2*S;
    Z_deltae = (-0.3648)*(1/2)*rho*V^2*S;
    M_deltae = (-1.444)*(1/2)*rho*V^2*S*c_bar;
    X_deltae = convforce(X_deltae,'lbf','N');
    Z_deltae = convforce(Z_deltae,'lbf','N');
    M_deltae = M_deltae*1.3558179619;
    
    fprintf(['\tX_' char(948) '_e = ' sprintf('%f',X_deltae) ' [N]\n']);
    fprintf(['\tZ_' char(948) '_e = ' sprintf('%f',Z_deltae) ' [N]\n']);
    fprintf(['\tM_' char(948) '_e = ' sprintf('%f',M_deltae) ' [N]\n\n']);
    
% 2.8 (Elevator feedback control law)
    % 2.8.a (increased pitch stiffness in the short period mode)        
        ks = 1:0.01:3;
        k1s = (Mq*(1-sqrt(ks)))/M_deltae;
        k2s = (V*Mw*(1-ks))/M_deltae;
        A_sp = [Mq/I_y V*Mw/I_y; 1 0];
        B_sp = [M_deltae/I_y; 0];
        A8a = @(k1,k2) A_sp + B_sp*[-k1 -k2];
        eigvals8a = arrayfun(@(k1,k2) eig(A8a(k1,k2))',k1s,k2s,'UniformOutput',false);
        eigvals8a = [eigvals8a{:}]';

        % plot
        figure
        hold on
        plot(eigvals8a,'x');
        yline(0);
        title('8.a Locus of Eigenvalues');
        xlabel('Real');
        ylabel('Im');
        grid on
        grid minor
        
        figure
        hold on
        plot(k1s,real(eigvals8a(1:2:end)),'x');
        yline(0);
        title('8.a Locus of Eigenvalues');
        xlabel('Real');
        ylabel('Im');
        grid on
        grid minor
        
    % 2.8.b (Implement full linearized longitudinal dynamics)
        B8b = [X_deltae/m 0; Z_deltae/(m-Zw_dot) 0; M_deltae/I_y + (Mw_dot/I_y)*(Z_deltae/(m-Zw_dot)) 0; 0 0];
        A_tot = @(k1,k2) (A + B8b*[0 0 -k1 -k2; 0 0 0 0]);
        eigvals8b = arrayfun(@(k1,k2) eig(A_tot(k1,k2))',k1s,k2s,'UniformOutput',false);
        eigvals8b = [eigvals8b{:}]';
        
        % plot
        figure
        hold on
        plot(eigvals8b,'x');
        xline(0);
        yline(0);
        title('8.b Locus of Eigenvalues');
        xlabel('Real');
        ylabel('Im');
        grid on
        grid minor
        
    % 2.8.c (Simulate the closed loop behavior)
        
        k = 1;
        k1 = (Mq*(1-sqrt(k)))/M_deltae;
        k2 = (V*Mw*(1-k))/M_deltae;
        A_tot = A + B8b*[0 0 -k1 -k2; 0 0 0 0];
        state_vec8c = [0; 0; 0; 0.1];
        [t_k1, state_vec_k1] = ode45(@(t_k1, state_vec_k1) LLDodefunc(t_k1,state_vec_k1,A_tot), t_span, state_vec8c);
        
        plotLLD(t_k1,state_vec_k1,'\Delta\theta = 0.1 [rad] (k=-1)','8.c'); % plot results
            
        k = 2;
        t_span = [0 100];
        k1 = (Mq*(1-sqrt(k)))/M_deltae;
        k2 = (V*Mw*(1-k))/M_deltae;
        A_tot = A + B8b*[0 0 -k1 -k2; 0 0 0 0];
        [t_k2, state_vec_k2] = ode45(@(t_k2, state_vec_k2) LLDodefunc(t_k2,state_vec_k2,A_tot), t_span, state_vec8c);
        
        plotLLD(t_k2,state_vec_k2,'\Delta\theta = 0.1 [rad] (k=-2)','8.c'); % plot results
        
    % 2.8.d
        %{
            Determine the reason for the effects on phugoid induced by the 
            short period control. Is this primarily due to the natural 
            coupling in the A matrix between the rotational and 
            translational states, or due to the additional control effects 
            in the B matrix (not modeled in the short period control 
            design) of elevator actuation on the translational states?
    
            The B matrix dominates

        %}
        

%% End Housekeeping
toc
    
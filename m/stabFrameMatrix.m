function stabilityMatrix = stabFrameMatrix(Xu, Zu, Mu, Xw, Zw, Mw, Xq, Zq, Mq, Zw_dot, Mw_dot, zeta)

    Xu_prime = Xu*cos(zeta)^2 - (Xw+Zu)*sin(zeta)*cos(zeta) + Zw*sin(zeta)^2;
    Xw_prime = Xw*cos(zeta)^2 + (Xu-Zw)*sin(zeta)*cos(zeta) - Zu*sin(zeta)^2;
    Xq_prime = Xq*cos(zeta) - Zq*sin(zeta);
    Xudot_prime = Zw_dot*sin(zeta)^2;
    Xwdot_prime = -Zw_dot*sin(zeta)*cos(zeta);
    
    Zu_prime = Zu*cos(zeta)^2 - (Zw-Xu)*sin(zeta)*cos(zeta) - Xw*sin(zeta)^2;
    Zw_prime = Zw*cos(zeta)^2 + (Zu+Xw)*sin(zeta)*cos(zeta) + Xu*sin(zeta)^2;
    Zq_prime = Zq*cos(zeta) + Zq*sin(zeta);
    Zudot_prime = -Zw_dot*sin(zeta)*cos(zeta);
    Zwdot_prime = Zw_dot*cos(zeta)^2;
    
    Mu_prime = Mu*cos(zeta) - Mw*sin(zeta);
    Mw_prime = Mw*cos(zeta) + Mu*sin(zeta);
    Mq_prime = Mq;
    Mudot_prime = -Mw_dot*sin(zeta);
    Mwdot_prime = Mw_dot*cos(zeta);

    X = [Xu_prime; Xw_prime; Xq_prime; Xudot_prime; Xwdot_prime];
    Z = [Zu_prime; Zw_prime; Zq_prime; Zudot_prime; Zwdot_prime];
    M = [Mu_prime; Mw_prime; Mq_prime; Mudot_prime; Mwdot_prime];
    
    stabilityMatrix = [X, Z, M]; % stability derivatives matrix in the stability frame
    
end
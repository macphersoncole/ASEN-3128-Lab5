function plotLLD(t,state_vec,state,problem)

    figure
    sgtitle([problem ' Linearized Longitudinal Dynamics ' state]);
    subplot(2,2,1);
    plot(t, state_vec(:,1),'linewidth',2);
    grid on
    grid minor
    xlabel('time [s]')
    ylabel('\Delta u^{E} [m/s]')
    title('Change in u versus time')
    xlim([0 max(t)])

    subplot(2,2,2);
    plot(t, state_vec(:,2),'linewidth',2);
    grid on
    grid minor
    xlabel('time [s]')
    ylabel('\Delta w^{E} [m/s]')
    title('Change in w versus time')
    xlim([0 max(t)])

    subplot(2,2,3);
    plot(t, state_vec(:,3),'linewidth',2);
    grid on
    grid minor
    xlabel('time [s]')
    ylabel('\Delta q [rad/s]')
    title('Change in pitch rate versus time')
    xlim([0 max(t)])

    subplot(2,2,4);
    plot(t, state_vec(:,4),'linewidth',2);
    grid on
    grid minor
    xlabel('time [s]')
    ylabel('\Delta \theta [rad]')
    title('Change in pitch versus time')
    xlim([0 max(t)])

end
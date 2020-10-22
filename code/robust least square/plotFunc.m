function plotFunc(param, scope, results, ifPlot)
    SGDA_results = results.SGDA;
    AGDA_results = results.AGDA;
    SSGDA_results = results.SSGDA;
    SAGDA_results = results.SAGDA;
    
    iter = param.iter;
    u0 = param.u0;
    
    info_str = strcat(', Iter: ', num2str(iter), ', Initial: (',...
                num2str(u0(1)), ', ', num2str(u0(2)), ')');
    
    MS = 16;
    LW = 2;
    LegFont = 12;
    
    %% 1. Objective
    figure
    axes('XScale', 'linear', 'YScale', 'linear');
    hold on
    %h_SGDA_obj = loglog(0:iter, SGDA_results.hist_obj, 'r-', 'LineWidth', LW);
    %h_AGDA_obj = loglog(0:iter, AGDA_results.hist_obj, 'b--', 'LineWidth', LW);
    h_SSGDA_obj = loglog(0:1000, SSGDA_results.hist_obj, 'r--', 'LineWidth', LW);
    h_SAGDA_obj = loglog(0:1000, SAGDA_results.hist_obj, 'b-', 'LineWidth', LW);
    hold off
    set(gca, 'fontsize', 18)
    l = legend([h_SSGDA_obj, h_SAGDA_obj], 'SGDA', 'AGDA', 'Location','northeast');
    l.FontSize = LegFont;
    %title(strcat('Objective Values', info_str), 'fontsize', 18)
    xlabel('Iterations', 'fontsize', 18)
    ylabel('Objective Value', 'fontsize', 18)
    
% %     filename = strcat(dir, '/Compare_Obj.fig');
%     if ifPlot
%         savefig(filename)
%     end
    %% 2.1 Stationary Point
    figure
    axes('XScale', 'linear', 'YScale', 'log');
    hold on
    h_SGDA_meas = loglog(0:iter, SGDA_results.meas, 'r-', 'LineWidth', LW);
    h_AGDA_meas = loglog(0:iter, AGDA_results.meas, 'b--', 'LineWidth', LW);
    hold off
    set(gca, 'fontsize', 18)
    l = legend([h_SGDA_meas, h_AGDA_meas], 'SGDA', 'AGDA', 'Location','northeast');
    l.FontSize = LegFont;
    %title('(b) Convergence of deterministic GDA', 'fontsize', 18)
    xlabel('Iterations', 'fontsize', 18)
    ylabel('$\Vert x_t-x^*\Vert^2+\Vert y_t-y^*\Vert^2$','Interpreter','latex', 'fontsize', 18)
    
%     filename = strcat(dir, '/Compare_Meas.fig');
%     if ifPlot
%         savefig(filename)
%     end




    %% 2.2. Stationary Point
    figure
    axes('XScale', 'linear', 'YScale', 'log');
    hold on
    h_SSGDA_meas = loglog(0:1000, SSGDA_results.meas, 'g-', 'LineWidth', LW);
    h_SAGDA_meas = loglog(0:1000, SAGDA_results.meas, 'm--', 'LineWidth', LW);
    hold off
    set(gca, 'fontsize', 18)
    l = legend([h_SSGDA_meas, h_SAGDA_meas], 'Stoc-SGDA', 'Stoc-AGDA', 'Location','northeast');
    l.FontSize = LegFont;
    %title('(c) Convergence of stochastic GDA', 'fontsize', 18)
    xlabel('Iterations', 'fontsize', 18)
    ylabel('$\Vert x_t-x^*\Vert^2+\Vert y_t-y^*\Vert^2$','Interpreter','latex', 'fontsize', 18)
    
    
    
    %% 3. Trajectory
    figure
    contour(scope.x, scope.y, scope.z, 100)
    hold on
    h_init = plot(param.u0(1), param.u0(2), 'm.', 'MarkerSize', MS*1.5);
    h_SGDA = plot(SGDA_results.hist(1,:), SGDA_results.hist(2,:),...
                'r-', 'LineWidth', LW);
    h_final_SGDA = plot(SGDA_results.final(1), SGDA_results.final(2), 'ro', ...
                        'LineWidth', LW, 'MarkerSize', MS);
    
    h_AGDA = plot(AGDA_results.hist(1,:), AGDA_results.hist(2,:),...
                'b--', 'LineWidth', LW);
    h_final_AGDA = plot(AGDA_results.final(1), AGDA_results.final(2), 'bo', ...
                    'LineWidth', LW, 'MarkerSize', MS);
                
    index1_SSGDA = abs(SSGDA_results.hist(1,:))<100;
    index2_SSGDA = abs(SSGDA_results.hist(2,:))<100;
    index_SSGDA = min(index1_SSGDA, index2_SSGDA);
    h_SSGDA = plot(SSGDA_results.hist(1,index_SSGDA), SSGDA_results.hist(2,index_SSGDA),...
                'g-', 'LineWidth', LW);
    h_final_SSGDA = plot(SSGDA_results.final(1), SSGDA_results.final(2), 'yo', ...
                        'LineWidth', LW, 'MarkerSize', MS);
    
    index1_SAGDA = abs(SAGDA_results.hist(1,:))<100;
    index2_SAGDA = abs(SAGDA_results.hist(2,:))<100;
    index_SAGDA = min(index1_SAGDA, index2_SAGDA);
    h_SAGDA = plot(SAGDA_results.hist(1,index_SAGDA), SAGDA_results.hist(2,index_SAGDA),...
                'm--', 'LineWidth', LW);
    h_final_SAGDA = plot(SAGDA_results.final(1), SAGDA_results.final(2), 'mo', ...
                    'LineWidth', LW, 'MarkerSize', MS);            
    hold off

    l = legend([h_SGDA, h_AGDA, h_SSGDA, h_SAGDA, h_init],...
            'SGDA', 'AGDA','Stoc-SGDA', 'Stoc-AGDA', 'Initial',...
            'Location','southeast');
    l.FontSize = LegFont;

    xlabel('x', 'fontsize', 18);
    ylabel('y', 'fontsize', 18);
    %title('(d) Trajactories of deterministic GDA', 'fontsize', 18)
    
%     filename = strcat(dir, '/Compare_Traj.fig');
%     if ifPlot
%         savefig(filename)
%     end
end
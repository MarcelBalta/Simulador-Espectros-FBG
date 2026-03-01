function Simulacion_FBG()
    % OBTENER TAMAÑO DE PANTALLA PARA CENTRAR VENTANAS
    screen = get(0, 'ScreenSize');
    sw = screen(3); sh = screen(4);

    % VENTANA PRINCIPAL CENTRADA
    fig_w = 1150; fig_h = 750; 
    fig = uifigure('Name', 'Simulación FBG - Control Dinámico', ...
        'Position', [(sw-fig_w)/2, (sh-fig_h)/2, fig_w, fig_h]);
    
    % VARIABLES DE DATOS
    x_exp = linspace(1545, 1555, 1000)'; 
    y_exp = zeros(size(x_exp));         

    % PANEL DE CONTROL
    panel_width = 330; 
    panel = uipanel(fig, 'Title', 'PARÁMETROS FÍSICOS', ...
        'TitlePosition', 'centertop', 'FontWeight', 'bold', 'Position', [20 50 panel_width 660]);
    
    ancho_botones = 260; ancho_edit = 180;    
    offset_botones = (panel_width - ancho_botones - 20) / 2; 
    offset_edit = (panel_width - ancho_edit - 20) / 2;

    function h = create_latex_input(y_pos, latex_text, start_val)
        uilabel(panel, 'HorizontalAlignment', 'center', 'Position', [0 y_pos panel_width-20 25], ...
            'Text', latex_text, 'Interpreter', 'latex', 'FontSize', 15);
        h = uieditfield(panel, 'numeric', 'HorizontalAlignment', 'center', ...
             'Position', [offset_edit y_pos-42 ancho_edit 35], ...
             'Value', start_val, 'ValueDisplayFormat', '%.6f', 'FontSize', 14);
    end

    % BLOQUE DE INPUTS
    edit_L = create_latex_input(600, 'Longitud $L$ (mm)', 3.32);
    edit_L.ValueDisplayFormat = '%.3f'; 
    edit_Lambda = create_latex_input(530, 'Periodo $\Lambda$ (nm)', 530.342);
    edit_dn = create_latex_input(460, '$\Delta n_{eff}$', 0.00045);
    edit_dn.ValueDisplayFormat = '%.8f'; 
    edit_neff = create_latex_input(390, 'Índice efectivo $n_{eff}$', 1.44729);

    % BLOQUE DE BOTONES
    y_base = 280; espaciado = 55; alto_btn = 45;

    uibutton(panel, 'Text', '1. CARGAR DATOS TXT', 'Position', [offset_botones y_base ancho_botones alto_btn], ...
        'BackgroundColor', [0.1 0.1 0.3], 'FontColor', 'w', 'FontWeight', 'bold', 'FontSize', 13, 'ButtonPushedFcn', @(btn,e) cargar_datos_callback());
    uibutton(panel, 'Text', '2. AUTO-AJUSTE', 'Position', [offset_botones y_base-espaciado ancho_botones alto_btn], ...
        'BackgroundColor', [0.4 0.2 0.6], 'FontColor', 'w', 'FontWeight', 'bold', 'FontSize', 13, 'ButtonPushedFcn', @(btn,e) run_robust_fit());
    uibutton(panel, 'Text', '3. ACTUALIZAR VISTA', 'Position', [offset_botones y_base-espaciado*2 ancho_botones alto_btn], ...
        'BackgroundColor', [0.2 0.6 1.0], 'FontColor', 'w', 'FontWeight', 'bold', 'FontSize', 13, 'ButtonPushedFcn', @(btn,e) update_main_plot());
    uibutton(panel, 'Text', '4. COMPARAR SIMULACIONES', 'Position', [offset_botones y_base-espaciado*3 ancho_botones alto_btn], ...
        'BackgroundColor', [0.8 0.4 0.2], 'FontColor', 'w', 'FontWeight', 'bold', 'FontSize', 13, 'ButtonPushedFcn', @(btn,e) open_comparison_tool());
    uibutton(panel, 'Text', '5. BARRIDO / EVOLUCIÓN', 'Position', [offset_botones y_base-espaciado*4 ancho_botones alto_btn], ...
        'BackgroundColor', [0.1 0.5 0.1], 'FontColor', 'w', 'FontWeight', 'bold', 'FontSize', 13, 'ButtonPushedFcn', @(btn,e) open_sweep_tool());

    % Ejes principales
    ax1 = uiaxes(fig, 'Position', [380 400 730 310], 'FontSize', 14);
    ax2 = uiaxes(fig, 'Position', [380 60 730 310], 'FontSize', 14);

    % --- FUNCIONES CORE ---
    function [x, y] = calculate_fbg(L_mm, Lam_nm, dn, n_eff)
        x = x_exp; lambda_v = x * 1e-9; 
        lambda_B = 2 * n_eff * (Lam_nm * 1e-9);
        lambda_centro = lambda_B * (1 + (dn / n_eff));
        kappa = (pi * dn) / lambda_B;
        sigma = 2 * pi * n_eff * ((1./lambda_v) - (1/lambda_centro));
        ratio_sq = (sigma ./ kappa).^2;
        y = zeros(size(lambda_v));
        idx_s = ratio_sq > 1;
        if any(idx_s)
            s_g = sqrt(sigma(idx_s).^2 - kappa^2);
            y(idx_s) = (kappa^2 * sin(s_g * (L_mm*1e-3)).^2) ./ (sigma(idx_s).^2 - kappa^2 * cos(s_g * (L_mm*1e-3)).^2);
        end
        idx_h = ratio_sq <= 1;
        if any(idx_h)
            h_g = sqrt(kappa^2 - sigma(idx_h).^2);
            y(idx_h) = (kappa^2 * sinh(h_g * (L_mm*1e-3)).^2) ./ (kappa^2 * cosh(h_g * (L_mm*1e-3)).^2 - sigma(idx_h).^2);
        end
        y = max(0, min(1, y)); 
    end

    function cargar_datos_callback()
        [file, path] = uigetfile('*.txt', 'Selecciona el espectro experimental');
        if isequal(file,0), return; end
        data = importdata(fullfile(path, file));
        if isstruct(data), x_raw = data.data(:,1); y_raw = data.data(:,2);
        else, x_raw = data(:,1); y_raw = data(:,2); end
        x_exp = x_raw; y_exp = y_raw; update_main_plot();
    end

    function open_comparison_tool()
        cw = 950; ch = 680;
        comp_fig = uifigure('Name', 'Comparador', 'Position', [(sw-cw)/2, (sh-ch)/2, cw, ch]);
        c_ax = uiaxes(comp_fig, 'Position', [70 230 830 400], 'FontSize', 14); grid(c_ax, 'on');
        p_comp = uipanel(comp_fig, 'Position', [50 20 850 180], 'Title', 'Ajuste de Comparación');
        
        drop = uidropdown(p_comp, 'Position', [30 50 160 35], 'Items', {'Delta n', 'Longitud L', 'Periodo'}, 'Value', 'Delta n', 'FontSize', 14, 'ValueChangedFcn', @(d,e) update_vals(d.Value));
        uilabel(p_comp, 'Text', 'Sim A', 'Position', [210 90 170 25], 'FontSize', 14, 'FontWeight', 'bold', 'FontColor', [0 0 0.8], 'HorizontalAlignment', 'center');
        uilabel(p_comp, 'Text', 'Sim B', 'Position', [400 90 170 25], 'FontSize', 14, 'FontWeight', 'bold', 'FontColor', [0.8 0 0], 'HorizontalAlignment', 'center');
        val_A = uieditfield(p_comp, 'numeric', 'Position', [210 50 170 35], 'Value', edit_dn.Value, 'FontSize', 14);
        val_B = uieditfield(p_comp, 'numeric', 'Position', [400 50 170 35], 'Value', edit_dn.Value * 1.1, 'FontSize', 14);
        uibutton(p_comp, 'Text', 'SUPERPONER', 'Position', [620 45 180 55], 'BackgroundColor', [0.8 0.4 0.2], 'FontColor', 'w', 'FontSize', 14, 'FontWeight', 'bold', 'ButtonPushedFcn', @(btn,e) run_comp());

        function update_vals(choice)
            switch choice
                case 'Delta n', val_A.Value=edit_dn.Value; val_B.Value=edit_dn.Value*1.1;
                case 'Longitud L', val_A.Value=edit_L.Value; val_B.Value=edit_L.Value*1.2;
                case 'Periodo', val_A.Value=edit_Lambda.Value; val_B.Value=edit_Lambda.Value+0.1;
            end
        end
        function run_comp()
            LA=edit_L.Value; LamA=edit_Lambda.Value; dnA=edit_dn.Value; LB=edit_L.Value; LamB=edit_Lambda.Value; dnB=edit_dn.Value;
            lblA = 'A'; lblB = 'B';
            switch drop.Value
                case 'Delta n', dnA=val_A.Value; dnB=val_B.Value; lblA = sprintf('\\Deltan = %.6f', dnA); lblB = sprintf('\\Deltan = %.6f', dnB);
                case 'Longitud L', LA=val_A.Value; LB=val_B.Value; lblA = sprintf('L = %.3f mm', LA); lblB = sprintf('L = %.3f mm', LB);
                case 'Periodo', LamA=val_A.Value; LamB=val_B.Value; lblA = sprintf('\\Lambda = %.3f nm', LamA); lblB = sprintf('\\Lambda = %.3f nm', LamB);
            end
            [xA, yA] = calculate_fbg(LA, LamA, dnA, edit_neff.Value);
            [xB, yB] = calculate_fbg(LB, LamB, dnB, edit_neff.Value);
            plot(c_ax, xA, yA, 'Color', [0 0 0.8], 'LineWidth', 2, 'DisplayName', lblA); hold(c_ax, 'on');
            plot(c_ax, xB, yB, 'Color', [0.8 0 0], 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', lblB); hold(c_ax, 'off');
            legend(c_ax, 'show', 'Location', 'northeast', 'FontSize', 11); grid(c_ax, 'on'); 
            xlabel(c_ax, 'Longitud de onda (nm)'); ylabel(c_ax, 'Reflectividad');
            ylim(c_ax, [0 1]);
        end
    end
    
    function open_sweep_tool()
        bw = 950; bh = 680;
        sweep_fig = uifigure('Name', 'Barrido Espectral', 'Position', [(sw-bw)/2, (sh-bh)/2, bw, bh]);
        s_ax = uiaxes(sweep_fig, 'Position', [80 230 820 400], 'FontSize', 14); grid(s_ax, 'on');
        p_sweep = uipanel(sweep_fig, 'Position', [50 20 850 180], 'Title', 'Configuración de Pasadas');
        
        s_drop = uidropdown(p_sweep, 'Position', [40 55 180 35], 'Items', {'Delta n', 'Longitud L'}, 'Value', 'Delta n', 'FontSize', 14, 'ValueChangedFcn', @(d,e) update_s_val());
        uilabel(p_sweep, 'Text', 'Valor Final', 'Position', [250 95 160 25], 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        uilabel(p_sweep, 'Text', 'Nº de Saltos', 'Position', [440 95 150 25], 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        
        % Inicialización automática con el valor exacto de la ventana principal
        s_end = uieditfield(p_sweep, 'numeric', 'Position', [250 55 160 35], 'Value', edit_dn.Value, 'FontSize', 14);
        s_steps = uieditfield(p_sweep, 'numeric', 'Position', [440 55 150 35], 'Value', 8, 'Limits', [1 100], 'FontSize', 14);
        
        uibutton(p_sweep, 'Text', 'EJECUTAR', 'Position', [630 45 180 60], 'BackgroundColor', [0.1 0.5 0.1], 'FontColor', 'w', 'FontSize', 14, 'FontWeight', 'bold', 'ButtonPushedFcn', @(btn,e) run_sweep());
        
        function update_s_val()
            if strcmp(s_drop.Value, 'Delta n'), s_end.Value = edit_dn.Value; else, s_end.Value = edit_L.Value; end
        end
        
        function run_sweep()
            cla(s_ax); hold(s_ax, 'on');
            N = s_steps.Value; 
            v_final = s_end.Value;
            % Lógica original: empezar desde el primer incremento hasta llegar al final
            vals = linspace(v_final/N, v_final, N); 
            colors = jet(N);
            for i = 1:N
                curr_L = edit_L.Value; curr_dn = edit_dn.Value;
                if strcmp(s_drop.Value, 'Delta n')
                    curr_dn = vals(i); lbl = sprintf('\\Deltan = %.6f', curr_dn);
                else
                    curr_L = vals(i); lbl = sprintf('L = %.3f mm', curr_L);
                end
                [~, y_loop] = calculate_fbg(curr_L, edit_Lambda.Value, curr_dn, edit_neff.Value);
                plot(s_ax, x_exp, y_loop, 'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', lbl);
            end
            if N <= 15, legend(s_ax, 'show', 'Location', 'northeastoutside', 'FontSize', 10); end
            ylim(s_ax, [0 1]); grid(s_ax, 'on');
            title(s_ax, ['Evolución de ', s_drop.Value], 'FontSize', 16); 
            xlabel(s_ax, 'Longitud de onda (nm)'); ylabel(s_ax, 'Reflectividad'); 
            hold(s_ax, 'off');
        end
    end

    function run_robust_fit()
        if all(y_exp == 0), uialert(fig, 'Cargue datos primero.', 'Aviso'); return; end
        [~, idx_pico] = max(y_exp); lambda_pico_nm = x_exp(idx_pico);
        dn_semilla = 0.0004; lambda_est_nm = lambda_pico_nm / (2 * (edit_neff.Value + dn_semilla));
        f_model = @(p, x) get_y_only_optimized(p(1), p(2), p(3), edit_neff.Value);
        options = optimoptions('lsqcurvefit', 'Display', 'none', 'TolFun', 1e-8);
        try
            p_opt = lsqcurvefit(f_model, [3, lambda_est_nm, 0.0004], x_exp, y_exp, [0.1, lambda_est_nm*0.99, 1e-6], [80.0, lambda_est_nm*1.01, 1e-2], options);
            edit_L.Value = p_opt(1); edit_Lambda.Value = p_opt(2); edit_dn.Value = p_opt(3);
        catch, end
        update_main_plot();
    end

    function y = get_y_only_optimized(L_mm, Lam_nm, dn, neff), [~, y] = calculate_fbg(L_mm, Lam_nm, dn, neff); end
    
    function update_main_plot()
        [~, y_sim] = calculate_fbg(edit_L.Value, edit_Lambda.Value, edit_dn.Value, edit_neff.Value);
        r_sq = 1 - sum((y_exp - y_sim).^2)/sum((y_exp - mean(y_exp)).^2);
        plot(ax1, x_exp, y_exp, 'Color', [0 0.2 0.8], 'LineWidth', 1.2, 'DisplayName', 'Exp.'); hold(ax1, 'on');
        plot(ax1, x_exp, y_sim, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Sim.'); hold(ax1, 'off');
        title(ax1, 'Ajuste FBG'); grid(ax1, 'on'); legend(ax1, 'show');
        xlabel(ax1, 'Longitud de onda (nm)'); ylabel(ax1, 'Reflectividad');
        ylim(ax1, [0 1]);
        delete(findall(ax1, 'Tag', 'r2text'));
        text(ax1, 'Units', 'normalized', 'Position', [0.98, 0.45], 'String', sprintf('$R^2$=%.4f', r_sq), 'Interpreter', 'latex', 'Tag', 'r2text', 'FontSize', 20, 'HorizontalAlignment', 'right', 'BackgroundColor', [1 1 1 0.7]);
        semilogy(ax2, x_exp, y_exp + 1e-5, 'Color', [0 0.2 0.8]); hold(ax2, 'on');
        semilogy(ax2, x_exp, y_sim + 1e-5, 'r--'); hold(ax2, 'off');
        title(ax2, 'Ajuste Logarítmico'); grid(ax2, 'on');
        xlabel(ax2, 'Longitud de onda (nm)'); ylabel(ax2, 'Reflectividad (log)');
        ylim(ax2, [1e-4 1.2]); 
    end

    update_main_plot();
end
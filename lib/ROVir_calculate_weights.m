function [eigenvec, SIR] = ROVir_calculate_weights(img_orig, S, I, SIRplot)
    % ROVir coils implementation for 3D data
    % Based on https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.28706
    % 
    % [eigenvec, pct_signal, pct_interf, SIR] = ROVir_calculate_weights(k_orig, S, I, SIRplot)
    %   k_orig = original k-space data for ROVir weight calculation
    %               dimensions: [nx, ny, nz, nc]
    %   S, I = binary arrays (same size as k_orig), indicating which areas
    %               are signal/interference, respectively ([nx, ny, nz])
    %   SIRplot = boolean: plot Signal, Interference, SIR values (fig 5)
    % 
    % matthijs.debuck@ndcn.ox.ac.uk
    % updates 20230821, e.m.schrauben@amsterdamumc.nl

    if ~exist('SIRplot', 'var')
        SIRplot = 0;
    end

    [nx, ny, nz, nc] = size(img_orig);

    %% Select Signal (S) and Interference (I) regions
    imgs_origS = img_orig.*S;
    imgs_origI = img_orig.*I;

    %% 1) Form matrices A and B (eqs 7 and 10)
    g_S = reshape(imgs_origS, [nx*ny*nz, nc]);  %eq. 8 for signal
    clear imgs_origS
    g_S = g_S(S(:), :);
    A = g_S' * g_S; clear g_S

    g_I = reshape(imgs_origI, [nx*ny*nz, nc]);  %eq. 8 for interference
    clear imgs_origI
    g_I = g_I(I(:), :);
    B = g_I' * g_I; clear g_I
    
    
    %% 2) Generalized eigenvalue decomposition
    [eigenvec, eigenval] = eig(A, B, 'vector'); 

    % sort: highest eigenvalues first
    [~,ind] = sort(eigenval, 'descend');
    eigenvec = eigenvec(:,ind);
    
    % normalize eigenvectors
    eigenvec = bsxfun(@rdivide, eigenvec, vecnorm(eigenvec,2));
    
    %% 3) Apply Gram-Schmidt orthonormalization
    eigenvec = GramSchmidt(eigenvec);
    
    %% Compute SIR, before sorting SIR values again 
    %SIR changes due to orthonormalization; calculate and reorder SIR
    %accordingly
    
    signal_ch_preSort = zeros(nc,1,'single'); interf_ch_preSort = zeros(nc,1,'single'); 
    for ch_j = 1:nc
        signal_ch_preSort(ch_j) = abs(eigenvec(:,ch_j)' * A * eigenvec(:,ch_j)); %eq. 6
        interf_ch_preSort(ch_j) = abs(eigenvec(:,ch_j)' * B * eigenvec(:,ch_j)); %eq. 9
    end

    SIR_preSort = signal_ch_preSort./interf_ch_preSort;  %signal-interference ratio (eq 11)

    %% Sort eigenvectors again (note: not used in ROVir paper)
    [~,ind] = sort(SIR_preSort, 'descend');
    eigenvec = eigenvec(:,ind);
    
    %% Calculate SIR again (recalculated i.o. reordered as an extra check)
    signal_ch = zeros(nc,1,'single'); interf_ch = zeros(nc,1,'single'); 
    for ch_j = 1:nc
        signal_ch(ch_j) = abs(eigenvec(:,ch_j)' * A * eigenvec(:,ch_j)); %eq. 6
        interf_ch(ch_j) = abs(eigenvec(:,ch_j)' * B * eigenvec(:,ch_j)); %eq. 9
    end

    SIR = signal_ch./interf_ch;  %signal-interference ratio (eq 11)
    
    %% Plot Signal, Interference, SIR
    if SIRplot
        figure('Position', [300 200 1300 350])

        peak_SIplot = max(max(signal_ch), max(interf_ch));
        subplot(1,3,1)
        plot(1:nc, signal_ch/peak_SIplot, 'LineWidth', 2)
        title('Signal')
        ylabel('a.u.'), xlabel('Virtual coil (j)')
        ylim([0 1])
        set(gca, 'FontSize', 15)

        subplot(1,3,2)
        plot(1:nc, interf_ch/peak_SIplot, 'LineWidth', 2)
        title('Interference')
        xlabel('Virtual coil (j)')
        ylim([0 1])
        set(gca, 'FontSize', 15)

        subplot(1,3,3)
        plot(1:nc, SIR, 'LineWidth', 2)
        title('SIR')
        xlabel('Virtual coil (j)')
        set(gca, 'FontSize', 15)
        
    end
end
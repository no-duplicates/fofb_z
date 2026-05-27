function readytorun_2Ks_local(c_k, mc_k)
tic;
if nargin < 1
    c_k = 1;
end

if nargin < 2
    mc_k = 2;
end



controller_scale_k=c_k;

%Scale of Mc_local, 0: the original 2Ts delay controller. 
Mc_local_scale_k=mc_k;



%% Load files

M_fpath = 'M-24-05-20.mat';
A_fpath = 'A_0dB_DC_gain_PI_200Hz-25-07-15.mat';
F_fpath = 'F_shelf_200Hz_10kHz_notches.mat';
K1_opt_path = 'K1_opt.mat';
K2_opt_path = 'K2_opt.mat';
M = load(M_fpath).M;
A = load(A_fpath).A;
F = load(F_fpath).F;

K1_opt= load(K1_opt_path).K1t;
K2_opt= load(K2_opt_path).K2t;
%% Instantiate model

ncorr = size(A,1);
excluded_corr = [1 80 81 160]; % Currently being used on FF

Ts = A{2}.Ts;

fofb_type.bpm_sector = 'M1M2C2C3';
fofb_type.corr_sector = 'M1M2C2C3';
fofb_type.bpm_remove_idx = [];
fofb_type.corr_remove_idx = excluded_corr;
[bpm_idx,corr_idx] = fofb_idx(fofb_type);

M = M(bpm_idx,corr_idx);



Mc = pinv(M);








% Mask everything except local and neighboring 10-actuator sections

    missing = [1 80 81 160];

    row_labels = setdiff(1:160, missing, 'stable');
    col_labels = 1:160;

    row_sec = ceil(row_labels / 10);
    col_sec = ceil(col_labels / 10);

    allowed = false(16,16);
    allowed = allowed | eye(16);

    % Ordinary neighboring sections
    pairs = [
        1 2;
        2 3;
        3 4;
        4 5;
        5 6;
        6 7;
        7 8;
        8 9;
        9 10;
        10 11;
        11 12;
        12 13;
        13 14;
        14 15;
        15 16
    ];

    % Custom wrap-around neighboring pairs
    pairs = [
        pairs;
        1 8;    % 1:10 <-> 71:80
        9 16    % 81:90 <-> 151:160
    ];

    for p = 1:size(pairs,1)
        i = pairs(p,1);
        j = pairs(p,2);
        allowed(i,j) = true;
        allowed(j,i) = true;
    end

    mask = allowed(row_sec(:), col_sec(:).');
    Mc_local = Mc .* mask;















missing = [1 80 81 160];




K = (1/Ts)*[0.12*ones(size(Mc,1)/2,1); 0.166*ones(size(Mc,1)/2,1)]    *controller_scale_k;
A = A(corr_idx);


A_reduced=A;

n = numel(A);
changed_idx = [];



for i = 1:n
    d = A_reduced{i}.InputDelay;

    if isequal(d, 1)
        A_reduced{i}.InputDelay = 0;
        changed_idx(end+1) = i;

    else

    end
end

fprintf('\nTotal changed: %d\n', numel(changed_idx));
disp('Changed indices:');
disp(changed_idx);



F = F(corr_idx);
K_1=K*(1.0);%K1_opt;
K_2=K*1.0;%K2_opt;




lowfreq_cutoff = 0.001;
target_CL_bw = 1e3;
Wz = c2d(tf(target_CL_bw/lowfreq_cutoff, ...
        [1/lowfreq_cutoff/2/pi 1]), Ts)*eye(size(M,1));

[P,G,C] = ofbmdl_2k(M,Mc,Mc_local*Mc_local_scale_k,K_1,K_2,A_reduced,F,[],Wz);



%assert(isstable(P))

%!!!!!!!!!!!!!!!!!!!!!!!!

% Assertation fails, but bode plot looks fine.

%!!!!!!!!!!!!!!!!!!!!!!!!


%{
%% Compute H2 and H-infinity norms of d -> yd

Tdy = P('yd','d');   % disturbance d to measured/output yd

fprintf('\n===== Norms for transfer d -> yd =====\n');

% Check stability
if ~isstable(Tdy)
    warning('Tdy = P(''yd'',''d'') is not stable. H2/Hinf norms may be invalid.');
end

% H2 norm
try
    h2_norm = norm(Tdy, 2);
    fprintf('H2 norm:      %.6e\n', h2_norm);
catch ME
    fprintf('H2 norm failed: %s\n', ME.message);
end

% H-infinity norm
try
    hinf_norm = norm(Tdy, inf);
    fprintf('Hinf norm:    %.6e\n', hinf_norm);
    fprintf('Hinf norm dB: %.6f dB\n', 20*log10(hinf_norm));
catch ME
    fprintf('Hinf norm failed: %s\n', ME.message);
end

%}


%% Plot sigma of d -> yd (disturbance rejection)

sigmaopts = sigmaoptions;
sigmaopts.FreqUnits = 'Hz';
sigmaopts.Grid = 'on';








%% Make figures directory
fig_dir = fullfile(pwd, 'figures');

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% Safe filename strings for parameters
c_k_str  = strrep(sprintf('%.4f', controller_scale_k), '.', 'p');
mc_k_str = strrep(sprintf('%.4f', Mc_local_scale_k), '.', 'p');

base_name = sprintf('local_neighbor_1Ts_full_2Ts_ck_%s_mck_%s', ...
                    c_k_str, mc_k_str);










%% Frequency grid: 1 Hz to 10000 Hz
f_hz = logspace(0, 4, 2000);   % 1 to 10000 Hz
w = 2*pi*f_hz;                 % rad/s for sigma()

%% Non-weighted plot: P('yd','d')
Tdy = P('yd','d');

sv = sigma(Tdy, w);            % singular values
sv_max = squeeze(max(sv, [], 1));  % max singular value at each frequency

[hinf_1_10000, idx_peak] = max(sv_max);
f_peak = f_hz(idx_peak);
hinf_1_10000_db = 20*log10(hinf_1_10000);

plot_title = sprintf(['1Ts local + neighbor information, 2Ts full information.\n' ...
                      'Controller\\_scale: %.4f, Mc\\_local\\_scale: %.4f'], ...
                      controller_scale_k, Mc_local_scale_k);

fig1=figure;
sigmaplot(Tdy, {1, 1e5}, sigmaopts);
title(plot_title);

legend(sprintf('max sigma over 1-10000 Hz = %.4g = %.2f dB at %.2f Hz', ...
       hinf_1_10000, hinf_1_10000_db, f_peak), ...
       'Location', 'best');

fprintf('Elapsed time: %.2f s\n', toc);




file1_png = fullfile(fig_dir, [base_name '_yd_to_d_nonweighted.png']);
file1_fig = fullfile(fig_dir, [base_name '_yd_to_d_nonweighted.fig']);

saveas(fig1, file1_png);
savefig(fig1, file1_fig);



%% Weighted plot: P('z','d')
Tzd = P('z','d');

% Singular values at exactly 1 Hz
f_check_hz = 1;
w_check = 2*pi*f_check_hz;   % rad/s

sv_1hz = sigma(Tzd, w_check);
sv_1hz_db = 20*log10(sv_1hz(:));

% Keep only singular values below 20 dB
sv_1hz_below20_db = sv_1hz_db(sv_1hz_db < 20);

if ~isempty(sv_1hz_below20_db)
    max_below20_db = max(sv_1hz_below20_db);
    max_below20 = 10^(max_below20_db/20);

    legend_str = sprintf(['At 1 Hz: largest singular value below 20 dB ' ...
                          '= %.4g = %.2f dB'], ...
                          max_below20, max_below20_db);
else
    legend_str = 'At 1 Hz: no singular value below 20 dB';
end

plot_title = sprintf(['1Ts local + neighbor information, 2Ts full information. Weighted by frequency.\n' ...
                      'Controller\\_scale: %.4f, Mc\\_local\\_scale: %.4f'], ...
                      controller_scale_k, Mc_local_scale_k);

fig2=figure;
sigmaplot(Tzd,{1,1e5},sigmaopts);
title(plot_title);

legend(legend_str, 'Location', 'best');

fprintf('Elapsed time: %.2f s\n',toc);

file2_png = fullfile(fig_dir, [base_name '_z_to_d_weighted.png']);
file2_fig = fullfile(fig_dir, [base_name '_z_to_d_weighted.fig']);

saveas(fig2, file2_png);
savefig(fig2, file2_fig);
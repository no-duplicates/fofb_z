

clc; clear; tic;

remove_delay=true;
mask_all=false;
mask_far=true;

%% Load files

M_fpath = 'M-24-05-20.mat';
A_fpath = 'A_0dB_DC_gain_PI_200Hz-25-07-15.mat';
F_fpath = 'F_shelf_200Hz_10kHz_notches.mat';

M = load(M_fpath).M;
A = load(A_fpath).A;
F = load(F_fpath).F;

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






if mask_all % mask all information from other sections
    % M is your 156x160 matrix

missing = [1 80 81 160];

% Row actuator labels after removing missing actuators
row_labels = setdiff(1:160, missing, 'stable');   % length = 156

% Column actuator labels are still full 1:160
col_labels = 1:160;

% Section index for each actuator: 1,...,8
row_sec = ceil(row_labels / 20);
col_sec = ceil(col_labels / 20);

% Keep only entries whose row and column are in the same 20-actuator section
mask = (row_sec(:) == col_sec(:).');   % 156x160 logical matrix

Mc = Mc .* mask;
end




if mask_far%mask everything except local and neighbours

missing = [1 80 81 160];

row_labels = setdiff(1:160, missing, 'stable');
col_labels = 1:160;

% now 10 actuators per section -> 16 sections
row_sec = ceil(row_labels / 10);
col_sec = ceil(col_labels / 10);

allowed = false(16,16);
allowed = allowed | eye(16);

% trivial neighboring sections
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

% custom neighboring pairs
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
Mc = Mc.* mask;

end


K = (1/Ts)*[0.12*ones(size(Mc,1)/2,1); 0.166*ones(size(Mc,1)/2,1)]    *1;
A = A(corr_idx);



if remove_delay%remove the delay
n = numel(A);
changed_idx = [];

for i = 1:n
    d = A{i}.InputDelay;

    if isequal(d, 1)
        A{i}.InputDelay = 0;
        changed_idx(end+1) = i;

    else

    end
end

fprintf('\nTotal changed: %d\n', numel(changed_idx));
disp('Changed indices:');
disp(changed_idx);

end

F = F(corr_idx);

[P,G,C] = ofbmdl(M,Mc,K,A,F);
assert(isstable(P))



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

%% Plot sigma of d -> yd (disturbance rejection)

sigmaopts = sigmaoptions;
sigmaopts.FreqUnits = 'Hz';
sigmaopts.Grid = 'on';




% Build title string
if remove_delay
    ts_str = '1Ts';
else
    ts_str = '2Ts';
end

if mask_all
    info_str = 'local information';
elseif mask_far
    info_str = 'local + neighbor information';
else
    info_str = 'full information';
end

plot_title = sprintf('%s, %s', ts_str, info_str);


figure;
sigmaplot(P('yd','d'),{1,1e5},sigmaopts);

title(plot_title);
fprintf('Elapsed time: %.2f s\n',toc);
%% cases.m
% Physical parameters for running the paper's NS
% cases (Figs 10-15), using the validated Keller-box boundary-layer
% solution (lambda, delta0, Lscr, Mbar for n=0.76, Pr=0.72, gamma=1.4).
%
% No freestream condition is changed anywhere in this script (Minf,
% T0inf, Tinf, Lstar, ReL all match the paper exactly) -- every fix
% below works by choosing R_scale (and, for Fig.14, which R values are
% geometrically viable at all) rather than perturbing the freestream.
%
% GEOMETRIC CONSTRAINT (applies to any expansion-corner case): the body
% surface reaches r=0 (the axis) at a downstream distance
%   reach_mm = r*_mm / tan(|alpha*|)
% which, expressed in scaled triple-deck x-units, reduces to the
% freestream-INDEPENDENT ratio
%   room_TD_units = R_scale / tan(|alpha*|)
% This ratio can only be changed by changing R_scale (if R_scale is not
% itself the tested variable) or alpha_scale (if it is). No freestream
% adjustment can alter it.
%
%
% FIG.14 (R-sweep, alpha=-12, adiabatic): R_scale IS the tested
% variable, so it cannot be arbitrarily changed. Only R=4,6,8,10 are
% included. R=4 is therefore also used as the default expansion radius
% for Figs 13 and 15.

clear; clc

%% ---------------- Assumptions (UNCHANGED from the paper, throughout) ----------------
gamma = 1.4;
n     = 0.76;
Pr    = 0.72;

Minf   = 6.0;
T0inf  = 600.0;     % K
Tinf   = 73.2;      % K
Lstar  = 0.1;       % m
ReL    = 1.08e5;

lambda = 0.72164;
delta0 = 0.25078;
Lscr   = -0.16308;
Mbar   = 1.69104;

hinf = 1/((gamma-1)*Minf^2);
Re0  = ReL*hinf^n;
chi  = Minf^2*Re0^(-0.5);
K    = (gamma-1)^(-0.5)*lambda^1.25*Minf^0.5*Re0^(-0.125);

r_factor = sqrt(Pr);
Taw      = Tinf*(1 + r_factor*(gamma-1)/2*Minf^2);
gw_aw    = (Taw/Tinf)*hinf;

fprintf('Re0=%.6g, chi=%.4f, Taw=%.2fK, g_w_aw=%.4f  (all UNCHANGED from the paper)\n\n', ...
    Re0, chi, Taw, gw_aw);

x_domain_mm        = [-30.0, 30.0];
x_downstream_mm    = x_domain_mm(2);
R_expansion        = 4.0;
domain_safety_frac = 0.65;   % suggested domain half-width, as a fraction of reach-axis, for Fig.14
R_cylinder_R4_mm = R_expansion * ...
    denomR(gw_aw,gamma,lambda,Minf,n,Re0)*Lstar*1e3;

fprintf('Expansion x-domain = [%.1f, %.1f] mm\n', x_domain_mm(1), x_domain_mm(2));
fprintf('Adiabatic R=4 physical radius = %.6f mm\n\n', R_cylinder_R4_mm);

%% ---------------- Compression cases (Figs 10-12): UNCHANGED ----------------
compression_cases = {
    'Fig.10 (compression)', 4.0,   0.0
    'Fig.10 (compression)', 4.0,  -2.5
    'Fig.10 (compression)', 4.0,  -5.0
    'Fig.10 (compression)', 4.0,  -7.5
    'Fig.10 (compression)', 4.0, -10.0
    'Fig.11 (compression)', 6.0,   0.0
    'Fig.11 (compression)', 6.0,  -2.5
    'Fig.11 (compression)', 6.0,  -5.0
    'Fig.11 (compression)', 6.0,  -7.5
    'Fig.11 (compression)', 6.0, -10.0
    'Fig.12 (compression)', 9.0,   0.0
    'Fig.12 (compression)', 9.0,  -2.5
    'Fig.12 (compression)', 9.0,  -5.0
    'Fig.12 (compression)', 9.0,  -7.5
    'Fig.12 (compression)', 9.0, -10.0
};
print_SL_case_table(compression_cases, ones(size(compression_cases,1),1), ...
    lambda, Re0, K, Lscr, n, gamma, Minf, Lstar, Tinf, hinf, gw_aw, Taw, ...
    'Compression (R_scale=1, unchanged)');

%% ---------------- Expansion alpha-sweep (Fig.13): common R_scale=4 ----------------
fig13_cases = {
    'Fig.13 (expansion, alpha-sweep)', -11.5, 0.0
    'Fig.13 (expansion, alpha-sweep)', -12.0, 0.0
    'Fig.13 (expansion, alpha-sweep)', -12.3, 0.0
    'Fig.13 (expansion, alpha-sweep)', -12.4, 0.0
};
R_scale_fig13 = R_expansion*ones(size(fig13_cases,1),1);
print_SL_case_table(fig13_cases, R_scale_fig13, ...
    lambda, Re0, K, Lscr, n, gamma, Minf, Lstar, Tinf, hinf, gw_aw, Taw, ...
    'Fig.13 (common R_scale=4)');

%% ---------------- Expansion R-sweep (Fig.14): LARGE R ONLY ----------------
fig14_alpha = -12.0;
fig14_R_all = [4,6,8,10];

ang14 = alpha_to_deg(fig14_alpha, lambda, Re0);
fprintf('--- Fig.14 (R-sweep, alpha=%.1f, adiabatic): angle=%.2f deg, Tw*=%.2fK ---\n', ...
    fig14_alpha, ang14, Taw);
fprintf('%8s %10s %14s %14s %14s %12s\n', ...
    'R_scale','r*(mm)','reach-axis(mm)','TD-units room','domain +-(mm)','status');

fig14_R_kept = [];
for i = 1:numel(fig14_R_all)
    R = fig14_R_all(i);
    r_mm    = R*denomR(gw_aw,gamma,lambda,Minf,n,Re0)*Lstar*1e3;
    reach   = r_mm/tand(abs(ang14));
    td_room = R/tand(abs(ang14));
    if td_room < 1.5
        status = 'EXCLUDED';
        domain_str = '--';
    elseif td_room < 3.5
        status = 'borderline';
        domain_str = sprintf('%.2f', domain_safety_frac*reach);
        fig14_R_kept(end+1) = R; %#ok<SAGROW>
    else
        status = 'recommended';
        domain_str = sprintf('%.2f', domain_safety_frac*reach);
        fig14_R_kept(end+1) = R; %#ok<SAGROW>
    end
    fprintf('%8.0f %10.3f %14.3f %14.2f %14s %12s\n', R, r_mm, reach, td_room, domain_str, status);
end
fprintf('Note: R=1,2 EXCLUDED (< 1.5 TD-units of room); this is an intrinsic\n');
fprintf('geometric/theoretical constraint at this alpha, independent of any\n');
fprintf('freestream or meshing choice -- see discussion in the text.\n\n');

%% ---------------- Expansion SL-sweep (Fig.15): common physical radius at R=4 ----------------
fig15_alpha = -12.3;
fig15_SL    = [0.0, -2.5, -5.0, -7.5, -10.0];

R_scale_fig15 = zeros(size(fig15_SL));
for i = 1:numel(fig15_SL)
    SL = fig15_SL(i);
    if SL == 0
        gw = gw_aw;
    else
        gw = gw_from_SL(SL, K, Lscr, n);
    end
    r_at_R1 = denomR(gw,gamma,lambda,Minf,n,Re0)*Lstar*1e3;
    R_scale_fig15(i) = R_cylinder_R4_mm/r_at_R1;
end

fig15_cases = cell(numel(fig15_SL),3);
for i = 1:numel(fig15_SL)
    fig15_cases{i,1} = 'Fig.15 (expansion, SL-sweep)';
    fig15_cases{i,2} = fig15_alpha;
    fig15_cases{i,3} = fig15_SL(i);
end
fprintf('Fig.15 common physical cylinder radius = %.6f mm\n', R_cylinder_R4_mm);
fprintf('Required downstream domain = %.1f mm\n\n', x_downstream_mm);
print_SL_case_table(fig15_cases, R_scale_fig15, ...
    lambda, Re0, K, Lscr, n, gamma, Minf, Lstar, Tinf, hinf, gw_aw, Taw, ...
    'Fig.15 (R=4 reference; common physical DNS radius)');

%% ================= local functions =================
function ang = alpha_to_deg(alpha_scale, lambda, Re0)
    ang = alpha_scale * lambda^0.5 * Re0^(-0.25) * 180/pi;
end

function gw = gw_from_SL(SL, K, Lscr, n)
    gw = (SL/(K*Lscr))^(-1/(n+0.5));
end

function d = denomR(gw, gamma, lambda, Minf, n, Re0)
    d = (gamma-1)^0.5 * lambda^(-1.25) * Minf^1.5 * gw^(n+0.5) * Re0^(-0.375);
end

function print_SL_case_table(cases, R_scale_vec, lambda, Re0, K, Lscr, n, gamma, Minf, Lstar, Tinf, hinf, gw_aw, Taw, title_str)
    fprintf('--- %s ---\n', title_str);
    fprintf('%-32s %8s %8s %8s %10s %10s %12s\n', 'Figure','alpha','angle','R_scale','Tw*(K)','r*(mm)','reach(mm)');
    nRows = size(cases,1);
    for i = 1:nRows
        alpha_s = cases{i,2};
        SL      = cases{i,3};
        R_s     = R_scale_vec(i);
        ang = alpha_to_deg(alpha_s, lambda, Re0);
        if SL == 0
            gw = gw_aw; Tw = Taw;
        else
            gw = gw_from_SL(SL, K, Lscr, n);
            Tw = gw*Tinf/hinf;
        end
        r_mm = R_s*denomR(gw,gamma,lambda,Minf,n,Re0)*Lstar*1e3;
        reach = r_mm/tand(abs(ang));
        fprintf('%-32s %8.2f %8.2f %8.2f %10.2f %10.4f %12.3f\n', ...
            cases{i,1}, alpha_s, ang, R_s, Tw, r_mm, reach);
    end
    fprintf('\n');
end

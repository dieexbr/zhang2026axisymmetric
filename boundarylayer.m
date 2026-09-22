%% boundarylayer.m
% Complete determination of the upstream hypersonic boundary layer for
% arbitrary n, Pr, gamma, and freestream conditions:
%   Keller-box solve -> lambda, delta0, script-L, Mbar, M0(eta)
%   -> g_w, Re0, chi, S, S*L, alpha, R
%
% Requires: pack, unpack, residual_vec, numeric_jacobian, solve_kellerbox
% (as defined in lambda_kellerbox_general.m)

clear; clc; clear all

gamma = 1.4;
n  = 0.76;
Pr = 0.72;
TH_FLOOR = 1e-10;

%% 1) Converge the boundary-layer profile via continuation
N = 400;
[phi, Y] = solve_kellerbox(N, 1.0, 1.0, gamma, TH_FLOOR, [], []);   % baseline

n_path = linspace(1.0, n, 12);
for k = 2:numel(n_path)
    [phi, Y] = solve_kellerbox(N, n_path(k), 1.0, gamma, TH_FLOOR, Y, phi);
end

Pr_path = linspace(1.0, Pr, 8);
for k = 2:numel(Pr_path)
    [phi, Y] = solve_kellerbox(N, n, Pr_path(k), gamma, TH_FLOOR, Y, phi);
end

[F,S,Th,w] = unpack(Y);
lambda = S(1);
fprintf('lambda(n=%.2f,Pr=%.2f) = %.6f\n', n, Pr, lambda);

%% 2) Recover eta(phi), delta0, M0(eta), script-L, Mbar
% deta/dphi = Theta^n / S  -- evaluate at box MIDPOINTS to avoid the
% exact 0/0 singularity at phi=1 (both Theta and S vanish there).
h    = diff(phi);
Thm  = 0.5*(Th(1:end-1)+Th(2:end));
Sm   = 0.5*(S(1:end-1)+S(2:end));
phim = 0.5*(phi(1:end-1)+phi(2:end));

integrand_eta = Thm.^n ./ Sm;
eta = [0; cumsum(integrand_eta.*h)];
delta0 = eta(end);

M0m = phim ./ sqrt(gamma*(gamma-1)*Thm);   % M0 at eta(midpoints)

Lscr_integrand = (1.0./M0m.^2 - 1.0) .* integrand_eta;
Lscr = sum(Lscr_integrand.*h);
Mbar = 1.0/sqrt(1 + Lscr/delta0);

fprintf('delta0 = %.6f\n', delta0);
if Lscr < 0
    regime_str = 'supercritical, Mbar>1';
else
    regime_str = 'subcritical, Mbar<1';
end
fprintf('script-L = %.6f  (%s)\n', Lscr, regime_str);

fprintf('Mbar = %.6f\n', Mbar);

figure
subplot(1,3,1)
plot(phim, Thm, phim, phim, 'LineWidth', 1.3); grid on
xlabel('\phi'); legend('\Theta','\phi','Location','southeast')
title('Profile vs \phi')
subplot(1,3,2)
plot(eta(2:end), M0m, 'LineWidth', 1.3); grid on
xlabel('\eta = Y'); ylabel('M_0'); title('Mach-number distribution')
subplot(1,3,3)
plot(eta(2:end), Lscr_integrand, 'LineWidth', 1.3); grid on
xlabel('\eta = Y'); ylabel('1/M_0^2 - 1'); title('script-L integrand')

%% 3) Freestream nondimensionalization (Chuvakhov et al. 2021 example)
Minf  = 6.0;
T0inf = 600;    % K
Tinf  = 73.2;   % K
Lstar = 0.1;    % m
ReL   = 1.08e5;

hinf = 1.0/((gamma-1)*Minf^2);
Re0  = ReL*hinf^n;
chi  = Minf^2*Re0^(-0.5);
fprintf('\nh_inf'' = %.5f, Re0 = %.5g, chi = %.4f\n', hinf, Re0, chi);

%% 4) Wall-cooling and geometric parameters, for an example wall temperature
Tw_ex = 150;    % K
gw_ex = (Tw_ex/Tinf)*hinf;
S_ex  = (gamma-1)^(-0.5)*lambda^1.25*Minf^0.5*gw_ex^(-(n+0.5))*Re0^(-0.125);
SL_ex = S_ex*Lscr;
fprintf('Tw*=%.0fK -> g_w=%.4f, S=%.4f, S*L=%.4f\n', Tw_ex, gw_ex, S_ex, SL_ex);

alpha_phys_deg = 15;
alpha_scale = deg2rad(alpha_phys_deg)*lambda^(-0.5)*Re0^0.25;
fprintf('alpha(15 deg) = %.4f\n', alpha_scale);

Rstar_phys = 1.53e-3;   % m
Rprime = Rstar_phys/Lstar;
Mexp_R = 1.5;
denomR = (gamma-1)^0.5*lambda^(-1.25)*Minf^Mexp_R*gw_ex^(n+0.5)*Re0^(-0.375);
R_scale = Rprime/denomR;
fprintf('R(r*=1.53mm) = %.4f\n', R_scale);

%% ===================== local functions =====================
function Y = pack(F,S,Th,w)
    Y = zeros(4*numel(F),1);
    Y(1:4:end)=F; Y(2:4:end)=S; Y(3:4:end)=Th; Y(4:4:end)=w;
end

function [F,S,Th,w] = unpack(Y)
    F=Y(1:4:end); S=Y(2:4:end); Th=Y(3:4:end); w=Y(4:4:end);
end

function R = residual_vec(Y, phi, h, phim, n, Pr, gamma, TH_FLOOR)
    N = numel(phi)-1; Ntot = 4*(N+1);
    [F,S,Th,w] = unpack(Y);
    Fm  = 0.5*(F(1:end-1)+F(2:end));
    Sm  = 0.5*(S(1:end-1)+S(2:end));
    Thm = 0.5*(Th(1:end-1)+Th(2:end));
    wm  = 0.5*(w(1:end-1)+w(2:end));

    Thm_safe = max(Thm, TH_FLOOR);                    % positivity floor
    Sm_safe  = Sm;
    tiny = abs(Sm) < 1e-12;
    Sm_safe(tiny) = sign(Sm(tiny) + 1e-300)*1e-12;     % avoid /0

    dF  = diff(F)./h;  dS = diff(S)./h;
    dTh = diff(Th)./h; dw = diff(w)./h;

    R1 = dF - phim.*Thm_safe.^(n-1)./((gamma-1)*Sm_safe);
    R2 = dS + Fm/2;
    R3 = dTh - wm;
    R4 = dw - ( (Fm/2).*(1-Pr).*wm - Pr*Sm_safe )./Sm_safe;

    R = zeros(Ntot,1);
    R(1:4:4*N)=R1; R(2:4:4*N)=R2; R(3:4:4*N)=R3; R(4:4:4*N)=R4;
    R(Ntot-3)=F(1); R(Ntot-2)=Th(1); R(Ntot-1)=Th(end); R(Ntot)=S(end);
end

function [J,R0] = numeric_jacobian(Y,phi,h,phim,n,Pr,gamma,TH_FLOOR)
    Ntot = numel(Y);
    R0 = residual_vec(Y,phi,h,phim,n,Pr,gamma,TH_FLOOR);
    rows=[]; cols=[]; vals=[];
    eps_=1e-7;
    for k = 1:Ntot
        dYk = eps_*max(1.0,abs(Y(k)));
        Yp = Y; Yp(k)=Yp(k)+dYk;
        Rp = residual_vec(Yp,phi,h,phim,n,Pr,gamma,TH_FLOOR);
        col = (Rp-R0)/dYk;
        nz = find(col ~= 0);
        rows=[rows; nz]; cols=[cols; k*ones(numel(nz),1)]; vals=[vals; col(nz)]; %#ok<AGROW>
    end
    J = sparse(rows,cols,vals,Ntot,Ntot);
end

function [phi,Y] = solve_kellerbox(N, n, Pr, gamma, TH_FLOOR, Y_init, phi_init)
    j = (0:N)';
    phi = 0.5*(1 - cos(pi*j/N));       % clustered at both ends
    h = diff(phi);
    phim = 0.5*(phi(1:end-1)+phi(2:end));

    if isempty(Y_init)
        lam_guess = 0.525;
        Sg  = lam_guess*(1-phi);
        Thg = 0.5*phi.*(1-phi);
        wg  = 0.5 - phi;
        integrand = phi.*max(Thg,TH_FLOOR).^(n-1)./((gamma-1)*max(Sg,1e-8));
        Fg = [0; cumsum(0.5*(integrand(1:end-1)+integrand(2:end)).*h)];
        Y = pack(Fg,Sg,Thg,wg);
    else
        [Fp,Sp,Thp,wp] = unpack(Y_init);
        F  = interp1(phi_init,Fp,phi,'linear','extrap');
        S  = interp1(phi_init,Sp,phi,'linear','extrap');
        Th = interp1(phi_init,Thp,phi,'linear','extrap');
        w  = interp1(phi_init,wp,phi,'linear','extrap');
        Y = pack(F,S,Th,w);
    end

    tol = 1e-12; maxit = 50;
    for it = 1:maxit
        [J,R] = numeric_jacobian(Y,phi,h,phim,n,Pr,gamma,TH_FLOOR);
        if norm(R) < tol, break; end
        dY = J\(-R);
        maxstep = max(abs(dY)./(abs(Y)+1e-3));
        if maxstep > 0.4, alpha = 0.4/maxstep; else, alpha = 1.0; end
        Y = Y + alpha*dY;
    end
end
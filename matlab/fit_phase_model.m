function out = fit_phase_model(fdMHz, dphiDeg, fcMHz, varargin)
%FIT_PHASE_MODEL  Подгонка моделей разности фаз канала: M1 / M2 / M3.
%
%   Модели (fd — МГц со знаком, Δτ — нс, фазы — градусы):
%     M1:  Δφ = C                            — только офсет;
%     M2:  Δφ = −0.36·fd·Δτ                  — только задержка (модель
%                                               прежней версии анализатора);
%     M3:  Δφ = C(fc) − 0.36·fd·Δτ           — офсет + задержка.
%   M3 строится СОВМЕСТНО по всем несущим: одна задержка Δτ на всё
%   множество записей и свой офсет C(fc) на каждой fc. Тем самым сразу
%   проверяется постоянство задержки при перенастройке и получается
%   калибровочная таблица (fc, C, Δτ) для компенсации:
%       x_k  <-  x_k · exp(−j·(C_k(fc) − 0.36·f·Δτ_k)),   f — МГц от fc.
%
%   СИНТАКСИС:
%     out = fit_phase_model(fd, dphi, fc)
%     out = fit_phase_model(fd, dphi, fc, 'TauStep', 0.02, 'Verbose', true)
%
%   ВХОДЫ (векторы равной длины, только записи с детектированным тоном):
%     fdMHz   — смещение тона со знаком, МГц (0 допустим: точка несёт
%               только информацию об офсете);
%     dphiDeg — разность фаз Δφ_k − φ_ref, град (в [−180, 180));
%     fcMHz   — несущая, МГц.
%
%   ВЫХОД out (структура):
%     out.M1.C, out.M1.rms               — офсет (круговое среднее), град;
%                                          СКО невязки, град;
%     out.M2.tau, out.M2.rms, .resid     — задержка M2, нс; СКО невязки;
%                                          невязки точек, град;
%     out.M3.tau, out.M3.rms, .resid     — совместная задержка M3, нс;
%     out.M3.fc, out.M3.C                — векторы: fc (МГц) и C(fc), град;
%     out.M3.Cmap                        — офсеты по fc (имена полей — валидные
%                                          идентификаторы): Cmap.fc700, Cmap.fc1000;
%     out.ambNs                          — период неоднозначности M3 (см. ниже);
%     out.ambM2Ns                        — период неоднозначности M2, нс;
%     out.nPoints, out.nFdDistinct       — число записей / различных |fd|≠0;
%     out.hasPairs                       — у каждого |fd| есть обе ветви ±;
%     out.separable                      — true, если данные определяют и C,
%                                          и Δτ (есть ≥2 различных fd хотя
%                                          бы на одной несущей);
%     out.dof                            — nPoints − (nFc + 1): избыточность
%                                          фита (≤0 — фит точный по
%                                          построению, модель не проверить);
%     out.warnM2, out.warnM3             — тексты предупреждений ('' = ок).
%
%   НЕОДНОЗНАЧНОСТЬ. У M3 (офсеты C(fc) свободны) пара (C, Δτ+T) даёт те
%   же предсказания, если 0.36·fd_i·T ≡ 0.36·fd_j·T (mod 360) для всех
%   точек одной fc, т.е. T·(fd_i−fd_j) ≡ 0 (mod 1000). Отсюда
%       ambNs = 1000 / D,   D = НОД разностей |fd_i − fd_j| внутри каждой fc.
%   Для полной пары ±fd: D = 2|fd| -> период 500/|fd| (та же величина, что
%   у формулы пары ΔΣ = −(Δφ+−Δφ−)/(0.72·fd)). У M2 (офсета нет)
%   неоднозначность шире: T·fd ≡ 0 (mod 1000) для всех точек ->
%       ambM2Ns = 1000 / НОД(|fd|).
%   Для {±5, ±10, ±20} МГц: оба периода 200 нс; взаимно простые |fd|
%   (например, 7 и 13) снимают неоднозначность почти полностью.
%   Если внутри каждой fc только одно значение fd, D = 0: задержка M3
%   не определяется вовсе (офсет всё забирает) -> tau = NaN.
%
%   МЕТОД: сетка по Δτ (шаг TauStep) в пределах ±amb/2 (свой диапазон для
%   M2 и M3 — один период неоднозначности); офсет на каждом шаге —
%   круговое среднее после вычитания задержечного члена с заворотом фазы
%   (корректно при любых заворотах); сетка считается векторно сразу для
%   всех Δτ; уточнение золотым сечением.
%
%   Совместимо с MATLAB (R2016b и новее — используется неявное расширение
%   размерностей) и GNU Octave (без тулбоксов).

% ---- разбор аргументов ------------------------------------------------------
if nargin < 3, error('fit_phase_model: требуется fd, dphi, fc'); end
opts = struct('TauStep', 0.02, 'Verbose', true);
for ii = 1:2:numel(varargin)
    if ~isfield(opts, varargin{ii})
        error('fit_phase_model: неизвестный параметр %s', varargin{ii});
    end
    opts.(varargin{ii}) = varargin{ii+1};
end

x  = double(fdMHz(:));
y  = double(dphiDeg(:));
fv = double(fcMHz(:));
if numel(x) ~= numel(y) || numel(x) ~= numel(fv)
    error('fit_phase_model: fd, dphi, fc разной длины');
end
if isempty(x), error('fit_phase_model: нет данных для фита'); end
if any(~isfinite(x)) || any(~isfinite(y)) || any(~isfinite(fv))
    error('fit_phase_model: нечисловые значения (NaN/Inf) во входах');
end

% ---- диагностика набора fd --------------------------------------------------
nPts = numel(x);
fcList = unique(fv);                     % несущие по возрастанию
nFc = numel(fcList);
xAbs = unique(abs(x(x ~= 0)));           % различные |fd| без нулей, МГц
nFdDistinct = numel(xAbs);
hasPairs = true;
for ii = 1:nFdDistinct
    if ~any(x == xAbs(ii)) || ~any(x == -xAbs(ii))
        hasPairs = false;                % у этого |fd| нет пары ±fd
    end
end

% D — НОД разностей fd внутри каждой fc (неоднозначность и разделимость M3)
D = 0;
for jf = 1:nFc
    fdc = x(fv == fcList(jf));
    for ii = 1:numel(fdc)
        for jj = ii+1:numel(fdc)
            D = gcd_local(D, abs(fdc(ii) - fdc(jj)));
        end
    end
end
separable = (D > 0);
dof = nPts - nFc - 1;                    % точки − (офсеты + задержка)

if nFdDistinct >= 1
    g = xAbs(1);
    for ii = 2:nFdDistinct
        g = gcd_local(g, xAbs(ii));
    end
    ambM2 = 1e3 / g;
else
    ambM2 = NaN;
end
if nFdDistinct == 0
    ambM3 = NaN;
elseif D > 0
    ambM3 = 1e3 / D;
else
    ambM3 = Inf;                         % задержка не определяется
end

% ---- предупреждения ---------------------------------------------------------
warnM2 = '';
warnM3 = '';
if nPts < 2
    warnM3 = 'менее 2 записей — параметры фита недостоверны';
    warnM2 = warnM3;
end
if nFdDistinct == 0
    if nPts >= 2
        warnM3 = 'все fd = 0 — задержка не определяется (доступен только офсет C)';
        warnM2 = warnM3;
    end
else
    if D == 0
        w = ['внутри каждой fc одно значение fd — задержка M3 не определяется ' ...
             '(офсет забирает весь сигнал); нужны ≥ 2 различных fd хотя бы ' ...
             'на одной несущей'];
    else
        w = '';
        if dof <= 0
            w = sprintf(['избыточность фита %d (точек не больше, чем параметров) — ' ...
                         'невязка обнулится при любой верной модели, ' ...
                         'адекватность не проверить'], dof);
        end
        if ~hasPairs
            if isempty(w), w = ''; else, w = [w, '; ']; end
            w = [w, 'не у всех |fd| есть пары ±fd — офсет и задержка ' ...
                    'разделяются только по наклону (менее устойчиво, чем по парам)'];
        end
    end
    if ~isempty(w)
        if isempty(warnM3), warnM3 = w; else, warnM3 = [warnM3, '; ', w]; end
    end
end

% ---- M1: только офсет -------------------------------------------------------
C1  = circmean_local(y);
r1  = rms_local(wrapdeg_local(y - C1));
out.M1 = struct('C', C1, 'rms', r1);

% ---- M2: только задержка (модель прежнего анализатора) ---------------------
if nFdDistinct == 0
    out.M2 = struct('tau', NaN, 'rms', NaN, 'resid', NaN(numel(y), 1));
else
    T2 = make_grid_local(ambM2, opts.TauStep);
    % невязка y − (−0.36·fd·τ) с заворотом, сразу для всех τ сетки (G×N):
    Zg = repmat(y(:)', numel(T2), 1) + 0.36 * (x(:)' .* T2);
    Zg = mod(Zg + 180, 360) - 180;       % заворот к [−180, 180)
    r2s = sqrt(mean(Zg.^2, 2));
    [~, j2] = min(r2s);
    lo = T2(j2) - opts.TauStep;  hi = T2(j2) + opts.TauStep;
    for it = 1:40                         % золотое сечение вокруг оптимума
        m1 = lo + (hi - lo) * 0.382;  m2 = hi - (hi - lo) * 0.382;
        if rms_local(wrapdeg_local(y + 0.36 * x * m1)) < ...
           rms_local(wrapdeg_local(y + 0.36 * x * m2))
            hi = m2;
        else
            lo = m1;
        end
    end
    t2 = (lo + hi) / 2;
    resid2 = wrapdeg_local(y + 0.36 * x * t2);
    out.M2 = struct('tau', t2, 'rms', rms_local(resid2), 'resid', resid2);
end

% ---- M3: офсет C(fc) на каждой несущей + общая задержка --------------------
if ~isfinite(ambM3)
    out.M3 = struct('tau', NaN, 'rms', NaN, 'resid', NaN(numel(y), 1), ...
                    'fc', fcList(:), 'C', NaN(nFc, 1), 'Cmap', struct());
else
    T3 = make_grid_local(ambM3, opts.TauStep);
    masks = false(numel(y), nFc);
    fcIdx = zeros(numel(y), 1);           % номер fc каждой точки
    for jf = 1:nFc
        masks(:, jf) = (fv == fcList(jf));
        fcIdx(masks(:, jf)) = jf;
    end
    % офсет C(fc) — круговое среднее по точкам каждой fc, для всех τ сразу:
    Zg = repmat(y(:)', numel(T3), 1) + 0.36 * (x(:)' .* T3);
    Zg = mod(Zg + 180, 360) - 180;
    Ez = exp(1j * deg2rad(Zg));           % G×N
    Cg = zeros(numel(T3), nFc);
    for jf = 1:nFc
        Cg(:, jf) = rad2deg(angle(mean(Ez(:, masks(:, jf)), 2)));
    end
    Rg = mod(Zg - Cg(:, fcIdx) + 180, 360) - 180;   % невязки M3, G×N
    r3s = sqrt(mean(Rg.^2, 2));
    [~, j3] = min(r3s);
    lo = T3(j3) - opts.TauStep;  hi = T3(j3) + opts.TauStep;
    for it = 1:40                         % золотое сечение вокруг оптимума
        m1 = lo + (hi - lo) * 0.382;  m2 = hi - (hi - lo) * 0.382;
        if m3_rms_local(x, y, masks, m1) < m3_rms_local(x, y, masks, m2)
            hi = m2;
        else
            lo = m1;
        end
    end
    t3 = (lo + hi) / 2;
    z = wrapdeg_local(y + 0.36 * x * t3);
    C3 = zeros(nFc, 1);
    resid3 = z;
    for jf = 1:nFc
        C3(jf) = circmean_local(z(masks(:, jf)));
        resid3(masks(:, jf)) = wrapdeg_local(z(masks(:, jf)) - C3(jf));
    end
    Cmap = struct();
    for jf = 1:nFc
        % имя поля struct обязано быть валидным идентификатором MATLAB
        % (не может начинаться с цифры и содержать точку), поэтому
        % префикс fc + замена недопустимых символов: Cmap.fc700, Cmap.fc1000
        fld = regexprep(sprintf('fc%g', fcList(jf)), '[^0-9a-zA-Z_]', '_');
        Cmap.(fld) = C3(jf);
    end
    out.M3 = struct('tau', t3, 'rms', rms_local(resid3), 'resid', resid3, ...
                    'fc', fcList(:), 'C', C3, 'Cmap', Cmap);
end

out.ambNs = ambM3;
out.ambM2Ns = ambM2;
out.nPoints = nPts;
out.nFdDistinct = nFdDistinct;
out.hasPairs = hasPairs;
out.separable = separable;
out.dof = dof;
out.warnM2 = warnM2;
out.warnM3 = warnM3;

% ---- печать -----------------------------------------------------------------
if opts.Verbose
    fprintf('  M1 (офсет):           C = %+8.2f град, СКО невязки %7.2f град\n', ...
            out.M1.C, out.M1.rms);
    if isnan(out.M2.tau)
        fprintf('  M2 (задержка): не оценивается (%s)\n', warnM2);
    else
        fprintf('  M2 (только задержка): Δτ = %+8.2f нс, СКО невязки %7.2f град\n', ...
                out.M2.tau, out.M2.rms);
    end
    if isnan(out.M3.tau)
        fprintf('  M3 (офсет+задержка): не оценивается (%s)\n', warnM3);
    else
        fprintf('  M3 (офсет+задержка):  Δτ = %+8.2f нс (общая на %d fc), СКО %7.2f град\n', ...
                out.M3.tau, numel(out.M3.fc), out.M3.rms);
        for jf = 1:numel(out.M3.fc)
            fprintf('      C(%g МГц) = %+8.2f град\n', out.M3.fc(jf), out.M3.C(jf));
        end
        fprintf('      период неоднозначности M3 = %.0f нс', ambM3);
        if isfinite(ambM2) && abs(ambM2 - ambM3) > 1e-9
            fprintf(' (M2: %.0f нс)', ambM2);
        end
        fprintf('\n');
    end
    if ~isempty(warnM3)
        fprintf('      ВНИМАНИЕ: %s\n', warnM3);
    end
end
end % function fit_phase_model

%% ======================= локальные функции =================================
function T = make_grid_local(ambNs, step)
%make_grid_local  Сетка Δτ в пределах ±amb/2 (один период неоднозначности).
tauMax = ceil(ambNs / 2);
T = (-tauMax : step : tauMax);
T = T(:);
end

function r = m3_rms_local(x, y, masks, t)
%m3_rms_local  СКО невязки M3 при задержке t (офсеты — по своим fc).
z = wrapdeg_local(y + 0.36 * x * t);
res = z;
for jf = 1:size(masks, 2)
    res(masks(:, jf)) = wrapdeg_local(z(masks(:, jf)) - ...
                                      circmean_local(z(masks(:, jf))));
end
r = rms_local(res);
end

function y = wrapdeg_local(x)
%wrapdeg_local  Приведение углов (град) к [-180, 180).
y = mod(x + 180, 360) - 180;
end

function c = circmean_local(a)
%circmean_local  Круговое среднее углов (град) — через средний вектор.
c = rad2deg(angle(mean(exp(1j * deg2rad(a)))));
end

function r = rms_local(x)
r = sqrt(mean(x.^2));
end

function g = gcd_local(a, b)
%gcd_local  НОД двух чисел (fd может быть дробным после str2double).
a = abs(a); b = abs(b);
while b > 1e-9
    t = mod(a, b);
    a = b; b = t;
end
g = a;
end

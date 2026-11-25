%% ANALYZE_PHASE_DIFFS  Анализ разности фаз каналов по файлам data_<...>.bin.
%
% Скрипт (пункт 4 задания). Для каждого файла data_<DDMMHH>_<f>_<fd>.bin
% (fd может быть ОТРИЦАТЕЛЬНЫМ и записываться как 5, так и -5):
%   - читает IQ-кадр (формат PELENGATOR IQ CAPTURE V1);
%   - оценивает фазы каналов на частоте тона fd (из имени файла, со знаком);
%   - строит матрицу разностей фаз 4x4 (теплокарта с числами) и кривые
%     Δφ_k - Δφ_1 относительно канала 1 (ch0);
%   - подгоняет модели Δφ_k(fd) по всем записям с тоном (fit_phase_model):
%       M2: Δφ = −0.36·fd·Δτ            — только задержка (прежняя версия);
%       M3: Δφ = C_k(fc) − 0.36·fd·Δτ   — офсет + задержка, СОВМЕСТНО по
%                                         всем fc: одна Δτ_k на всё
%                                         множество записей + офсеты C_k(fc);
%     и выводит калибровочную таблицу (fc, C_k, Δτ_k) для компенсации
%     каналов: x_k ← x_k·exp(−j·(C_k(fc) − 0.36·f·Δτ_k)), f — МГц от fc;
%   - диагностирует адекватность: СКО невязки фита сравнивается с СКО
%     измерения (dphi_*_std_deg); если невязка на порядок выше шума,
%     печатается предупреждение «данные не описываются чисто задержечной
%     моделью», а Δτ(M2) помечается неадекватным (флаг в CSV/MAT);
%   - сводный график Δφ(fd) с линиями фита M3 и оценками Δτ по файлам.
%
% ВЫВОД РЕЗУЛЬТАТОВ (в папку данных):
%   phase_diff_results.csv       — по строке на файл: fc, fd, SNR, фазы,
%                                  Δφ_k-ch0 и их СКО (формат прежний);
%   phase_diff_calibration.csv   — калибровочная таблица M3: по строке на
%                                  (канал, fc) — C, Δτ, период
%                                  неоднозначности, СКО невязок M2/M3,
%                                  СКО измерения, флаги разделимости и
%                                  адекватности;
%   phase_diff_results.mat       — то же + полные матрицы Δφ каждого файла
%                                  и все результаты фитов;
%   консольные таблицы + предупреждения.
%
% НАСТРОЙКА: путь к папке с данными — в блоке НАСТРОЙКИ.
% Запуск: analyze_phase_diffs (файл-функция — параметры внутри, локальные
% функции ниже; скрипт-локальные функции Octave не поддерживает).
% Опционально: analyze_phase_diffs(struct('DataDir', '...')) — переопределить
% настройки (DataDir, FileMask, SNRMin, MakeFigures, SaveCsv, SaveMat, SavePng).
% Совместимо с MATLAB (R2016b и новее) и GNU Octave (не требует тулбоксов).

function analyze_phase_diffs(overrides)

if nargin < 1, overrides = struct(); end

%% ------------------------- НАСТРОЙКИ -----------------------------------
DataDir    = '251126';    % папка с файлами data_<time>_<f>_<fd>.bin
                           % (можно задать полный путь)
FileMask   = 'data_*_*.bin';% маска файлов анализа
SNRMin     = 20;            % мин. SNR тона, дБ (ниже — файл помечается)
MakeFigures= true;          % false -> только расчёт и CSV/MAT (headless-запуск)
SaveCsv    = true;          % выгрузка CSV
SaveMat    = true;          % выгрузка .mat
SavePng    = true;          % сохранение окон в PNG (требует MakeFigures)

% применение переопределений (analyze_phase_diffs(struct('DataDir',...)))
ovFields = fieldnames(overrides);
for iv = 1:numel(ovFields)
    if isfield(overrides, ovFields{iv})
        eval([ovFields{iv}, ' = overrides.', ovFields{iv}, ';']);
    end
end
%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(thisDir);

if ~exist(DataDir, 'dir')
    cand = {fullfile(thisDir, DataDir), ...
            fullfile(thisDir, '..', DataDir), ...
            fullfile(thisDir, '..', 'testdata')};
    ok = '';
    for i = 1:numel(cand)
        if exist(cand{i}, 'dir'), ok = cand{i}; break; end
    end
    if isempty(ok)
        error('Папка с данными не найдена: %s. Укажите путь в DataDir.', DataDir);
    end
    DataDir = ok;
end

files = dir(fullfile(DataDir, FileMask));
if ~isempty(files)
    files = sort({files.name});
end
if isempty(files)
    error('В папке %s нет файлов по маске %s.', DataDir, FileMask);
end
fprintf('Найдено файлов: %d в %s\n', numel(files), DataDir);

%% --- разбор имён и расчёт ---------------------------------------------------
nChMax = 0;
R = struct('file', {}, 'fcMHz', {}, 'fdMHz', {}, 'snrDb', {}, ...
           'tone', {}, 'phaseDeg', {}, 'dphiDeg', {}, 'dphiStdDeg', {}, ...
           'dphiMatrix', {}, 'freqUsed', {}, 'nCh', {});

for i = 1:numel(files)
    fname = files{i};
    [p, fext] = fileparts(fname);
    assert(isempty(p)); %#ok<ASRTUC>  % файлы плоские
    parts = regexp(fext, '^data_(\d+)_([+-]?\d+)_([+-]?\d+)$', 'tokens', 'once');
    [iq, hdr] = read_iq_capture(fullfile(DataDir, fname));
    if isempty(parts)
        fprintf('[%2d/%2d] %s: имя не разобрано, fc=%g МГц, fd ищу по спектру\n', ...
                i, numel(files), fname);
        fdHint = [];
        fcMHz = hdr.fc / 1e6;
        fdMHz = NaN;
    else
        fcMHz = str2double(parts{2});
        fdMHz = str2double(parts{3});      % знак сохраняется (-5 == тон на -5 МГц)
        fdHint = fdMHz * 1e6;
        % сверка с заголовком
        if abs(hdr.fc - fcMHz*1e6) > 1
            fprintf(['[%2d/%2d] %s: ВНИМАНИЕ fc в имени (%g) != заголовку (%g)\n'], ...
                    i, numel(files), fname, fcMHz, hdr.fc/1e6);
        end
    end
    est = estimate_tone_phase(iq, hdr.fs, 'FreqHint', fdHint, 'SNRMin', SNRMin);
    if isempty(fdHint)
        fdMHz = est.freqUsed / 1e6;
    end

    R(i).file = fname;
    R(i).fcMHz = fcMHz;
    R(i).fdMHz = fdMHz;
    R(i).snrDb = est.snrDb;
    R(i).tone = est.toneDetected;
    R(i).phaseDeg = est.phaseDeg;
    R(i).dphiDeg = est.dphiDeg;          % относительно канала 1 (ch0)
    R(i).dphiStdDeg = est.dphiStdDeg;
    R(i).freqUsed = est.freqUsed;
    R(i).nCh = size(iq, 1);
    % полная матрица 4x4
    nCh = size(iq, 1);
    D = est.phaseDeg * ones(1, nCh) - ones(nCh, 1) * est.phaseDeg';
    D = mod(D + 180, 360) - 180;
    R(i).dphiMatrix = D;
    nChMax = max(nChMax, nCh);

    if est.toneDetected
        fprintf(['[%2d/%2d] %s: fc=%g МГц, fd=%g МГц, SNR=%.1f дБ, ' ...
                 'Δφ(k-1) = %s град\n'], i, numel(files), fname, ...
                fcMHz, fdMHz, est.snrDb, mat2str(round(est.dphiDeg'*10)/10));
    else
        fprintf(['[%2d/%2d] %s: fc=%g МГц, fd=%g МГц, ТОН НЕ ОБНАРУЖЕН ' ...
                 '(SNR=%.1f дБ < %d) — Δφ не достоверны!\n'], ...
                i, numel(files), fname, fcMHz, fdMHz, est.snrDb, SNRMin);
    end
end

%% --- окно 1: матрицы Δφ 4x4 для каждого файла -------------------------------
nFiles = numel(files);
fig1 = [];
if MakeFigures
fig1 = figure('Name', 'Матрицы разности фаз', 'Color', 'w', ...
              'Position', [40 40 1200 720]);
ncols = ceil(sqrt(nFiles));
nrows = ceil(nFiles / ncols);
for i = 1:nFiles
    subplot(nrows, ncols, i);
    D = R(i).dphiMatrix;
    imagesc(D);
    set(gca, 'YDir', 'reverse');
    colormap(gca, hsv_colormap());
    hold on;
    nCh = R(i).nCh;
    for r = 1:nCh
        for c = 1:nCh
            txt = sprintf('%.1f', D(r, c));
            plot(c, r, 'r.', 'MarkerSize', 4); %#ok<*NOFP>  % точки в центрах
            text(c, r - 0.28, txt, 'HorizontalAlignment', 'center', ...
                 'FontSize', 7, 'Color', 'w', 'FontWeight', 'bold');
        end
    end
    hold off;
    xticks_local(1:nCh, arrayfun(@(k) sprintf('ch%d', k-1), 1:nCh, ...
                 'UniformOutput', false));
    yticks_local(1:nCh, arrayfun(@(k) sprintf('ch%d', k-1), 1:nCh, ...
                 'UniformOutput', false));
    axis equal tight;
    caxis_local([-180 180]);
    cb = colorbar();
    ylabel(cb, 'Δφ, град');
    ttl = sprintf('%s  fd=%g МГц', R(i).file, R(i).fdMHz);
    if ~R(i).tone
        ttl = [ttl, '  (нет тона!)'];
    end
    title(ttl, 'FontSize', 8);
end
if exist('sgtitle', 'file') || exist('sgtitle', 'builtin')
    sgtitle('Матрицы Δφ = φ_i − φ_j (град), строки минус столбцы');
end
end % MakeFigures

%% --- фиты моделей M2/M3 по всем записям (доработка по п. 6.1 отчёта) --------
% M2: Δφ = −0.36·fd·Δτ           — только задержка (модель прежней версии);
% M3: Δφ = C_k(fc) − 0.36·fd·Δτ  — офсет + задержка, СОВМЕСТНЫЙ фит по всем
%                                  несущим: одна Δτ_k на всё множество
%                                  записей + офсет C_k(fc) на каждой fc.
% Адекватность: СКО невязки фита сравнивается с СКО измерения
% (dphi_*_std_deg из CSV); невязка на порядок выше шума означает, что
% данные не описываются чисто задержечной моделью — печатается
% предупреждение, и Δτ(M2) помечается неадекватным (флаг в CSV/MAT).
fprintf('\n=== Фиты моделей разности фаз (M2 — только задержка, M3 — офсет+задержка) ===\n');
fcList = [];
for i = 1:nFiles
    if R(i).tone && isfinite(R(i).fcMHz)
        fcList = unique([fcList, R(i).fcMHz]);
    end
end
nFc = numel(fcList);
tausM2     = nan(nChMax, 1);    % M2: задержка, нс
tausStdM2  = nan(nChMax, 1);    % СКО невязки M2, град
tausM3     = nan(nChMax, 1);    % M3: совместная задержка (все fc), нс
tausStdM3  = nan(nChMax, 1);    % СКО невязки M3, град
Ck         = nan(nChMax, nFc);  % офсеты C_k(fc), град
tauAmb     = nan(nChMax, 1);    % период неоднозначности M3, нс
                                      % (1000/НОД разностей fd внутри fc)
tauAmbM2   = nan(nChMax, 1);    % период неоднозначности M2, нс (1000/НОД|fd|)
noiseStd   = nan(nChMax, 1);    % СКО измерения Δφ (медиана), град
m2Adequate = true(nChMax, 1);   % false = чисто задержечная модель неприменима
separableK = true(nChMax, 1);   % данные определяют и C, и Δτ
hasPairsK  = true(nChMax, 1);   % у всех |fd| есть пары ±fd
warnTxt    = cell(nChMax, 1);   % предупреждения фита M3

for k = 2:nChMax
    x = []; y = []; fv = []; sv = [];
    for i = 1:nFiles
        if R(i).tone && numel(R(i).dphiDeg) >= k ...
                && isfinite(R(i).fcMHz) && isfinite(R(i).fdMHz)
            x(end+1)  = R(i).fdMHz;         %#ok<*AGROW>
            y(end+1)  = R(i).dphiDeg(k);    %#ok<*AGROW>
            fv(end+1) = R(i).fcMHz;         %#ok<*AGROW>
            sv(end+1) = R(i).dphiStdDeg(k); %#ok<*AGROW>
        end
    end
    if isempty(x)
        fprintf('ch%d−ch0: нет записей с детектированным тоном — фит пропущен\n', k-1);
        m2Adequate(k) = false; separableK(k) = false; hasPairsK(k) = false;
        continue;
    end
    fprintf('ch%d−ch0 (%d записей, %d несущих):\n', k-1, numel(x), numel(unique(fv)));
    f = fit_phase_model(x, y, fv, 'Verbose', true);
    if isfinite(f.M2.tau)
        tausM2(k) = f.M2.tau;  tausStdM2(k) = f.M2.rms;
    end
    if isfinite(f.M3.tau)
        tausM3(k) = f.M3.tau;  tausStdM3(k) = f.M3.rms;
        for jf = 1:nFc
            jj = find(abs(f.M3.fc - fcList(jf)) < 1e-9, 1);
            if ~isempty(jj), Ck(k, jf) = f.M3.C(jj); end
        end
    end
    tauAmb(k)   = f.ambNs;
    tauAmbM2(k) = f.ambM2Ns;
    if ~isempty(sv), noiseStd(k) = median(sv); end
    separableK(k) = f.separable;
    hasPairsK(k)  = f.hasPairs;
    warnTxt{k}    = f.warnM3;

    % --- диагностика адекватности: невязка M2 против шума измерения ---
    if isfinite(tausStdM2(k)) && isfinite(noiseStd(k)) && noiseStd(k) > 0 ...
            && tausStdM2(k) / noiseStd(k) > 10
        m2Adequate(k) = false;
        fprintf(['  ВНИМАНИЕ (ch%d−ch0): СКО невязки M2 = %.1f град при шуме ' ...
                 'измерения %.3f град — данные НЕ описываются чисто ' ...
                 'задержечной моделью (доминирует офсет): Δτ(M2) = %+.2f нс ' ...
                 'не является задержкой канала, используйте M3.\n'], ...
                k-1, tausStdM2(k), noiseStd(k), tausM2(k));
    end
    fprintf('\n');
end

%% --- калибровочная таблица M3 -------------------------------------------------
if any(isfinite(tausM3))
    fprintf(['\n=== Калибровочная таблица M3: офсеты C_k(fc), град; ' ...
             'Δτ_k общая на все fc, нс ===\n']);
    fprintf('%10s', 'fc, МГц');
    for k = 2:nChMax
        fprintf(' %12s', sprintf('C(ch%d−ch0)', k-1));
    end
    fprintf('\n');
    for jf = 1:nFc
        fprintf('%10g', fcList(jf));
        for k = 2:nChMax
            if isfinite(Ck(k, jf)), fprintf(' %12.2f', Ck(k, jf));
            else, fprintf(' %12s', '—'); end
        end
        fprintf('\n');
    end
    fprintf('%10s', 'Δτ_k, нс');
    for k = 2:nChMax
        if isfinite(tausM3(k)), fprintf(' %12.2f', tausM3(k));
        else, fprintf(' %12s', '—'); end
    end
    fprintf('\n');
    fprintf('%10s', 'амб, нс');
    for k = 2:nChMax
        if isfinite(tauAmb(k)), fprintf(' %12.0f', tauAmb(k));
        elseif isnan(tauAmb(k)), fprintf(' %12s', '—');
        else, fprintf(' %12s', 'Inf'); end
    end
    fprintf('   (период неоднозначности M3)\n');
    fprintf(['Компенсация канала: x_k ← x_k·exp(−j·(C_k(fc) − 0.36·f·Δτ_k)), ' ...
             'f — МГц от fc (модель валидна в полосе снятых |fd|).\n']);
end

%% --- окно 2: Δφ_k - φ_1 по файлам и Δφ(fd) ----------------------------------
fig2 = [];
if MakeFigures
fig2 = figure('Name', 'Разности фаз и закономерность', 'Color', 'w', ...
              'Position', [80 80 1150 760]);

% 2а: Δφ(k-1) по файлам (группированные столбцы)
subplot(2, 2, [1 2]);
cols = lines_local(nChMax);
hold on;
for i = 1:nFiles
    dp = R(i).dphiDeg;
    nc = numel(dp);
    for k = 2:nc
        patch([i-0.32+(k-2)*0.16, i-0.24+(k-2)*0.16, ...
               i-0.24+(k-2)*0.16, i-0.32+(k-2)*0.16], ...
              [0 0 dp(k) dp(k)], cols(k, :), 'EdgeColor', 'none');
    end
    if ~R(i).tone
        plot(i, 190, 'rx', 'LineWidth', 2, 'MarkerSize', 9);
        text(i, 196, 'нет тона', 'HorizontalAlignment', 'center', ...
             'FontSize', 7, 'Color', 'r');
    end
end
plot([0 nFiles+1], [0 0], 'k-');
grid on; hold off;
set(gca, 'XTick', 1:nFiles, 'XTickLabel', ...
    arrayfun(@(i) sprintf('%s\\newlinefd=%g', R(i).file, R(i).fdMHz), ...
             1:nFiles, 'UniformOutput', false), 'XTickLabelRotation', 0);
xlim([0 nFiles+1]); ylim([-200 220]);
ylabel('Δφ_k − φ_{ch0}, град');
title('Разности фаз каналов относительно ch0 по файлам');
legend(arrayfun(@(k) sprintf('ch%d−ch0', k-1), 2:nChMax, ...
       'UniformOutput', false), 'Location', 'southeast');

% 2б: Δφ против fd — закономерность
subplot(2, 2, 3);
hold on;
for k = 2:nChMax
    x = []; y = [];
    for i = 1:nFiles
        if R(i).tone && numel(R(i).dphiDeg) >= k
            x(end+1) = R(i).fdMHz; %#ok<*AGROW>
            y(end+1) = R(i).dphiDeg(k); %#ok<*AGROW>
        end
    end
    if ~isempty(x)
        [xs, si] = sort(x);  ys = y(si);
        plot(xs, ys, 'o-', 'Color', cols(k, :), 'LineWidth', 1.2, ...
             'MarkerFaceColor', cols(k, :));
    end
end
% линии совместного фита M3: Δφ = C_k(fc) − 0.36·fd·Δτ_k (по каждой fc)
for k = 2:nChMax
    if ~isfinite(tausM3(k)), continue; end
    lineCol = 0.55 * cols(k, :) + 0.45;    % светлее цвета канала
    for jf = 1:nFc
        fdl = [];
        for i = 1:nFiles
            if R(i).tone && isfinite(R(i).fcMHz) ...
                    && abs(R(i).fcMHz - fcList(jf)) < 1e-9 ...
                    && isfinite(R(i).fdMHz)
                fdl(end+1) = R(i).fdMHz; %#ok<*AGROW>
            end
        end
        if isempty(fdl) || ~isfinite(Ck(k, jf)), continue; end
        xd = linspace(min(fdl), max(fdl), 200);
        ym = wrapdeg_local(Ck(k, jf) - 0.36 * xd * tausM3(k));
        ym([false, abs(diff(ym)) > 180]) = NaN;   % разрывы на заворотах ±360°
        plot(xd, ym, ':', 'Color', lineCol, 'LineWidth', 1.3);
    end
end
grid on; hold off;
xlabel('fd, МГц'); ylabel('Δφ_k − φ_{ch0}, град');
title('Δφ от fd: точки и совместный фит M3 (пунктир, по fc)');
ylim([-190 190]);

% 2в: оценка задержек Δτ_k = Δφ / (−2π fd) — накопление по fd
subplot(2, 2, 4);
hold on;
for k = 2:nChMax
    x = []; y = [];
    for i = 1:nFiles
        if R(i).tone && numel(R(i).dphiDeg) >= k && R(i).fdMHz ~= 0
            tau = -R(i).dphiDeg(k) / (0.36 * R(i).fdMHz); % нс (знак fd учтён)
            x(end+1) = R(i).fdMHz; %#ok<*AGROW>
            y(end+1) = tau; %#ok<*AGROW>
        end
    end
    if ~isempty(x)
        plot(x, y, 's', 'Color', cols(k, :), 'MarkerFaceColor', cols(k, :), ...
             'MarkerSize', 7);
    end
end
% совместные оценки Δτ_k (M3) — горизонтальные пунктирные линии
fdAll = [];
for i = 1:nFiles
    if R(i).tone && isfinite(R(i).fdMHz)
        fdAll = [fdAll, R(i).fdMHz]; %#ok<*AGROW>
    end
end
if ~isempty(fdAll)
    xr = [min(fdAll) - 1, max(fdAll) + 1];
    for k = 2:nChMax
        if isfinite(tausM3(k))
            plot(xr, [tausM3(k) tausM3(k)], ':', 'Color', cols(k, :), ...
                 'LineWidth', 1.5);
        end
    end
end
grid on; hold off;
xlabel('fd, МГц'); ylabel('Δτ_k = −Δφ/(0.36·fd), нс');
title('Оценка задержки по файлам; пунктир — совместная Δτ (M3)');
end % MakeFigures

%% --- консольная таблица ------------------------------------------------------
fprintf('\n=== Сводная таблица (Δφ относительно ch0, град) ===\n');
fprintf('%-28s %8s %8s %8s  ', 'файл', 'fc,МГц', 'fd,МГц', 'SNR,дБ');
fprintf('%s\n', strjoin(arrayfun(@(k) sprintf('%9s', ...
        sprintf('ch%d-ch0', k-1)), 2:nChMax, 'UniformOutput', false), ''));
for i = 1:nFiles
    fprintf('%-28s %8g %8g %8.1f  ', R(i).file, R(i).fcMHz, R(i).fdMHz, R(i).snrDb);
    for k = 2:nChMax
        if R(i).tone && numel(R(i).dphiDeg) >= k
            fprintf('%9.1f', R(i).dphiDeg(k));
        else
            fprintf('%9s', '—');
        end
    end
    fprintf('\n');
end

%% --- выгрузка результатов -----------------------------------------------------
outDir = DataDir;
if SaveCsv
    csvFile = fullfile(outDir, 'phase_diff_results.csv');
    fid = fopen(csvFile, 'w');
    fprintf(fid, 'file,fc_MHz,fd_MHz,snr_dB,tone_detected');
    for k = 2:nChMax
        fprintf(fid, ',dphi_ch%d_ch0_deg', k-1);
    end
    for k = 2:nChMax
        fprintf(fid, ',dphi_ch%d_ch0_std_deg', k-1);
    end
    fprintf(fid, ',freq_used_MHz\n');
    for i = 1:nFiles
        fprintf(fid, '%s,%g,%g,%.2f,%d', R(i).file, R(i).fcMHz, ...
                R(i).fdMHz, R(i).snrDb, R(i).tone);
        for k = 2:nChMax
            if R(i).tone && numel(R(i).dphiDeg) >= k
                fprintf(fid, ',%.3f', R(i).dphiDeg(k));
            else
                fprintf(fid, ',');
            end
        end
        for k = 2:nChMax
            if R(i).tone && numel(R(i).dphiDeg) >= k
                fprintf(fid, ',%.3f', R(i).dphiStdDeg(k));
            else
                fprintf(fid, ',');
            end
        end
        fprintf(fid, ',%.6f\n', R(i).freqUsed / 1e6);
    end
    fclose(fid);
    fprintf('CSV сохранён: %s\n', csvFile);
    % калибровочная таблица M3: по строке на (канал, fc)
    calFile = fullfile(outDir, 'phase_diff_calibration.csv');
    fid = fopen(calFile, 'w');
    fprintf(fid, ['ch,fc_MHz,C_deg,dtau_ns,dtau_amb_ns,dtau_amb_m2_ns,' ...
                  'rms_m2_deg,rms_m3_deg,noise_std_deg,separable,' ...
                  'has_pm_pairs,m2_adequate\n']);
    for k = 2:nChMax
        for jf = 1:nFc
            if isfinite(Ck(k, jf))
                fprintf(fid, 'ch%d,%g,%.3f,%.3f,%.1f,%.1f,%.3f,%.3f,%.4f,%d,%d,%d\n', ...
                        k-1, fcList(jf), Ck(k, jf), tausM3(k), tauAmb(k), ...
                        tauAmbM2(k), tausStdM2(k), tausStdM3(k), noiseStd(k), ...
                        separableK(k), hasPairsK(k), m2Adequate(k));
            end
        end
    end
    fclose(fid);
    fprintf('CSV калибровки сохранён: %s\n', calFile);
end
if SaveMat
    matFile = fullfile(outDir, 'phase_diff_results.mat');
    save(matFile, 'R', 'files', 'nChMax', 'DataDir', 'SNRMin', ...
         'fcList', 'tausM2', 'tausStdM2', 'tausM3', 'tausStdM3', 'Ck', ...
         'tauAmb', 'tauAmbM2', 'noiseStd', 'm2Adequate', 'separableK', ...
         'hasPairsK', 'warnTxt');
    fprintf('MAT сохранён: %s\n', matFile);
end
if SavePng && MakeFigures && ~isempty(fig1)
    saveas(fig1, fullfile(outDir, 'phase_diff_matrices.png'));
    saveas(fig2, fullfile(outDir, 'phase_diff_summary.png'));
    fprintf('PNG сохранены в %s\n', outDir);
end
end % function analyze_phase_diffs

%% ======================= локальные функции =================================
function y = wrapdeg_local(x)
    %WRAPDEG_LOCAL  Приведение углов (град) к [-180, 180).
    y = mod(x + 180, 360) - 180;
end

function cmap = hsv_colormap()
    % дiverging-карта красный-белый-синий для фаз [-180,180]
    n = 256;
    half = floor(n/2);
    top = [linspace(0.10, 0.35, half)', linspace(0.30, 0.55, half)', linspace(0.85, 0.95, half)'];
    bot = flipud(top);
    cmap = [top; [1 1 1]; bot];
end

function xticks_local(vals, labels)
    set(gca, 'XTick', vals, 'XTickLabel', labels);
end

function yticks_local(vals, labels)
    set(gca, 'YTick', vals, 'YTickLabel', labels);
end

function caxis_local(v)
    if exist('caxis', 'file') || exist('caxis', 'builtin')
        caxis(v);
    else
        clim(v);
    end
end

function cols = lines_local(n)
    base = [0 0.447 0.741; 0.850 0.325 0.098; 0.929 0.694 0.125; ...
            0.494 0.184 0.556; 0.466 0.674 0.188; 0.301 0.745 0.933; ...
            0.635 0.078 0.184];
    cols = repmat(base, ceil(n / size(base, 1)), 1);
    cols = cols(1:n, :);
end

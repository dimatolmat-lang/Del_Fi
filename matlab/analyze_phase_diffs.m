%% ANALYZE_PHASE_DIFFS  Анализ разности фаз каналов по файлам data_<...>.bin.
%
% Скрипт (пункт 4 задания). Для каждого файла data_<DDMMHH>_<f>_<fd>.bin:
%   - читает IQ-кадр (формат PELENGATOR IQ CAPTURE V1);
%   - оценивает фазы каналов на частоте тона fd (из имени файла);
%   - строит матрицу разностей фаз 4x4 (теплокарта с числами) и кривые
%     Δφ_k - Δφ_1 относительно канала 1 (ch0);
%   - оценивает межканальные задержки Δτ = Δφ / (2*pi*fd);
%   - сводный график Δφ(fd): если задержки каналов постоянны, точки для
%     каждого канала лягут на прямую Δφ_k(fd) = -360*fd*Δτ_k (mod 360).
%
% ВЫВОД РЕЗУЛЬТАТОВ:
%   phase_diff_results.csv  — по строке на файл: fc, fd, SNR, фазы, Δφ_k-1,
%                             Δτ_k-1, СКО Δφ;
%   phase_diff_results.mat  — то же + полные матрицы Δφ каждого файла;
%   консольная таблица + заключение о линейности Δφ(fd).
%
% НАСТРОЙКА: путь к папке с данными — в блоке НАСТРОЙКИ.
% Совместимо с MATLAB и GNU Octave (не требует тулбоксов).

%% ------------------------- НАСТРОЙКИ -----------------------------------
DataDir    = 'testdata';    % папка с файлами data_<time>_<f>_<fd>.bin
                           % (можно задать полный путь)
FileMask   = 'data_*_*.bin';% маска файлов анализа
SNRMin     = 20;            % мин. SNR тона, дБ (ниже — файл помечается)
SaveCsv    = true;          % выгрузка CSV
SaveMat    = true;          % выгрузка .mat
SavePng    = true;          % сохранение окон в PNG
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
    parts = regexp(fext, '^data_(\d+)_(\d+)_(\d+)$', 'tokens', 'once');
    [iq, hdr] = read_iq_capture(fullfile(DataDir, fname));
    if isempty(parts)
        fprintf('[%2d/%2d] %s: имя не разобрано, fc=%g МГц, fd ищу по спектру\n', ...
                i, numel(files), fname);
        fdHint = [];
        fcMHz = hdr.fc / 1e6;
        fdMHz = NaN;
    else
        fcMHz = str2double(parts{2});
        fdMHz = str2double(parts{3});
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
fig1 = figure('Name', 'Матрицы разности фаз', 'Color', 'w', ...
              'Position', [40 40 1200 720]);
nFiles = numel(files);
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

%% --- окно 2: Δφ_k - φ_1 по файлам и Δφ(fd) ----------------------------------
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
grid on; hold off;
xlabel('fd, МГц'); ylabel('Δφ_k − φ_{ch0}, град');
title('Δφ от частоты смещения (разворот на 360° — ожидание)');
ylim([-190 190]);

% 2в: оценка задержек Δτ_k = Δφ / (−2π fd) — накопление по fd
subplot(2, 2, 4);
hold on;
for k = 2:nChMax
    x = []; y = [];
    for i = 1:nFiles
        if R(i).tone && numel(R(i).dphiDeg) >= k
            tau = -R(i).dphiDeg(k) / (0.36 * R(i).fdMHz); % нс
            x(end+1) = R(i).fdMHz; %#ok<*AGROW>
            y(end+1) = tau; %#ok<*AGROW>
        end
    end
    if ~isempty(x)
        plot(x, y, 's', 'Color', cols(k, :), 'MarkerFaceColor', cols(k, :), ...
             'MarkerSize', 7);
    end
end
grid on; hold off;
xlabel('fd, МГц'); ylabel('Δτ_k = −Δφ/(0.36·fd), нс');
title('Оценка межканальной задержки по каждому файлу');

%% --- оценка суммарной задержки (фит Δφ = −0.36·fd·Δτ, fd в МГц, τ в нс) ----
fprintf('\n=== Оценка межканальных задержек (фит Δφ(fd) = −0.36·fd·Δτ) ===\n');
taus = nan(nChMax, 1);
tausStd = nan(nChMax, 1);
for k = 2:nChMax
    x = []; y = [];
    for i = 1:nFiles
        if R(i).tone && numel(R(i).dphiDeg) >= k
            x(end+1) = R(i).fdMHz;  %#ok<*AGROW>
            y(end+1) = R(i).dphiDeg(k); %#ok<*AGROW>
        end
    end
    if numel(x) >= 2
        % Grid-search по Δτ: невязка каждой точки с заворотом на 360°,
        % поэтому метод корректен при любом разбросе fd и любых заворотах фазы.
        % ВАЖНО: Δτ и Δτ + 1/gcd(fd) неразличимы. Если все fd кратны 5 МГц,
        % задержка определена с точностью до 200 нс — результат приводится
        % к диапазону [-1/(2*шаг), 1/(2*шаг)] нс. Для однозначности
        % записывайте файлы со взаимно простыми fd (например, 10 и 13 МГц).
        fdStep = min(diff(sort(unique(x))));
        if numel(unique(x)) > 1
            tauAmb = 1e3 / fdStep;           % нс — период неоднозначности (fdStep в МГц)
        else
            tauAmb = 1e3 / x(1);
        end
        [xg, si] = sort(x);  yg = y(si);
        tauGrid = (-ceil(tauAmb/2):0.02:ceil(tauAmb/2));   % нс
        bestRes = inf; bestTau = NaN;
        for tau = tauGrid
            % модель: Δφ_model = -0.36 * fd[МГц] * τ[нс] (т.е. -2π·f·τ)
            res = wrapdeg_local(yg - (-0.36 * xg * tau));
            r2 = rms_local(res);
            if r2 < bestRes, bestRes = r2; bestTau = tau; end
        end
        % уточнение: золотое сечение вокруг bestTau с шагом сетки
        step = 0.02;
        lo = bestTau - step; hi = bestTau + step;
        for it = 1:40
            m1 = lo + (hi - lo) * 0.382;  m2 = hi - (hi - lo) * 0.382;
            r1 = rms_local(wrapdeg_local(yg - (-0.36 * xg * m1)));
            r2 = rms_local(wrapdeg_local(yg - (-0.36 * xg * m2)));
            if r1 < r2, hi = m2; else, lo = m1; end
        end
        bestTau = (lo + hi) / 2;
        bestRes = rms_local(wrapdeg_local(yg - (-0.36 * xg * bestTau)));
        taus(k) = bestTau;     % нс
        tausStd(k) = bestRes;
        fprintf('ch%d-ch0: Δτ = %8.2f нс  (СКО невязки %.2f град, %d файлов, неоднозначность ±%.0f нс)\n', ...
                k-1, bestTau, bestRes, numel(x), tauAmb/2);
    elseif numel(x) == 1
        taus(k) = -y(1) / (0.36 * x(1));
        fprintf('ch%d−ch0: Δτ = %8.2f нс  (единственный файл)\n', k-1, taus(k));
    end
end

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
end
if SaveMat
    matFile = fullfile(outDir, 'phase_diff_results.mat');
    save(matFile, 'R', 'taus', 'tausStd', 'files', 'nChMax', ...
         'DataDir', 'SNRMin');
    fprintf('MAT сохранён: %s\n', matFile);
end
if SavePng
    saveas(fig1, fullfile(outDir, 'phase_diff_matrices.png'));
    saveas(fig2, fullfile(outDir, 'phase_diff_summary.png'));
    fprintf('PNG сохранены в %s\n', outDir);
end

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

function r = rms_local(x)
    r = sqrt(mean(x.^2));
end

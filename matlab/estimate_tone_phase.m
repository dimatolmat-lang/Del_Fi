function res = estimate_tone_phase(iq, fs, varargin)
%ESTIMATE_TONE_PHASE  Оценка частоты и фаз тонового сигнала в каналах.
%
%   RES = ESTIMATE_TONE_PHASE(IQ, FS)
%   RES = ESTIMATE_TONE_PHASE(IQ, FS, 'FreqHint', FD, 'SNRMin', 20, 'NSeg', 8)
%
%   Вход:
%     IQ — матрица nChannels x nSamples (комплексные отсчёты каналов);
%     FS — частота дискретизации, Гц.
%   Пары имя-значение:
%     'FreqHint' — известная частота смещения fd, Гц (положительная).
%                  Если задана, фазы считаются точно на ней.
%     'SNRMin'   — минимальный SNR тона (дБ), ниже которого тон считается
%                  необнаруженным. По умолчанию 20.
%     'NSeg'     — число сегментов для оценки стабильности фазы (8).
%
%   Выход RES — структура:
%     .toneDetected  — true/false: есть ли уверенный тон;
%     .snrDb         — SNR тона (дБ), средний по каналам;
%     .freqUsed      — частота, на которой считались фазы, Гц;
%     .phaseDeg      — фазы каналов, град (nChannels x 1);
%     .dphiDeg       — разности фаз phi_k - phi_1, град, в [-180,180);
%     .dphiStdDeg    — круговое СКО разностей фаз по сегментам, град;
%     .amp           — амплитуда тона в каналах (единицы отсчётов АЦП);
%     .method        — 'freqHint' | 'fft'.
%
%   Метод: окно Ханна по всей записи; частота — по пику суммарного спектра
%   каналов с параболическим уточнением; фаза канала — угол взвешенного
%   DFT-отклика на частоте тона. Разности фаз каналов почти не зависят от
%   малой ошибки частоты, поэтому Δφ устойчива.
%
%   Совместимо с MATLAB и GNU Octave (не требует тулбоксов).

    %% --- разбор аргументов ------------------------------------------------
    FreqHint = [];
    SNRMin = 20;
    NSeg = 8;
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~ischar(name)
            error('estimate_tone_phase:arg', ...
                  'Имя параметра должно быть строкой (пара имя-значение).');
        end
        if k == numel(varargin)
            error('estimate_tone_phase:arg', ...
                  'Параметру ''%s'' не задано значение.', name);
        end
        val = varargin{k+1};
        switch lower(name)
            case 'freqhint',  FreqHint = val;
            case 'snrmin',    SNRMin = val;
            case 'nseg',      NSeg = val;
            otherwise
                error('estimate_tone_phase:arg', ...
                      'Неизвестный параметр: %s', name);
        end
    end

    [nCh, nS] = size(iq);
    if nCh < 2
        error('estimate_tone_phase:chans', 'Нужно не менее 2 каналов.');
    end
    if ~isempty(FreqHint) && (FreqHint <= 0 || FreqHint > fs/2)
        error('estimate_tone_phase:hint', ...
              'FreqHint должна быть в (0, fs/2] = (0, %g] Гц.', fs/2);
    end

    win  = hann_local(nS);          % столбец nS x 1
    fAxis = fftfreqs_local(nS, fs); % столбец, порядок fftshift

    %% --- оценка частоты тона ----------------------------------------------
    if ~isempty(FreqHint)
        freqUsed = FreqHint;
        method = 'freqHint';
    else
        % суммарный спектр каналов: тон общий, пик надёжнее
        Ssum = zeros(nS, 1);
        for c = 1:nCh
            Ssum = Ssum + abs(fftshift(fft(iq(c,:).' .* win))).^2;
        end
        [~, kpk] = max(Ssum);
        binW = fs / nS;
        df = 0;
        if kpk > 1 && kpk < nS
            ym = 10*log10(Ssum(kpk-1));
            y0 = 10*log10(Ssum(kpk));
            yp = 10*log10(Ssum(kpk+1));
            den = ym - 2*y0 + yp;
            if den < 0    % парабола вершиной вверх в лог-домене
                df = 0.5 * (ym - yp) / den * binW;
            end
        end
        % если пик в отрицательной области спектра — берём |f|: для общего
        % сигнала это тот же тон, фазы каналов инвертируются согласованно
        freqUsed = abs(fAxis(kpk) + df);
        method = 'fft';
    end

    %% --- фазы каналов на частоте тона --------------------------------------
    n = (0:nS-1).';
    e = exp(-1i * 2 * pi * freqUsed * n / fs);   % столбец
    phaseDeg = zeros(nCh, 1);
    amp = zeros(nCh, 1);
    for c = 1:nCh
        x = iq(c,:).';                 % столбец
        resp = sum(x .* win .* e);
        phaseDeg(c) = angle(resp) * 180 / pi;
        amp(c) = 2 * abs(resp) / sum(win);
    end

    %% --- SNR тона по каждому каналу (окно Ханна, вся запись) --------------
    % Сигнал: пиковый бин и по 2 соседних (протечка окна).
    % Шум: медиана остальных бинов, умноженная на ENBW Ханна (1.5 бина).
    snrCh = zeros(nCh, 1);
    for c = 1:nCh
        S = abs(fftshift(fft(iq(c,:).' .* win))).^2 / (sum(win)^2 / 4);
        [~, kpk] = max(S);
        sigBins = max(1, kpk-2) : min(nS, kpk+2);
        sigPow = sum(S(sigBins));
        noiseS = S;
        noiseS(sigBins) = [];
        noiseFloor = median(noiseS);
        snrCh(c) = 10 * log10(sigPow / (noiseFloor * 1.5 * numel(sigBins)));
    end
    snrDb = mean(snrCh);

    %% --- стабильность разностей фаз по сегментам ----------------------------
    nSeg = max(2, min(NSeg, floor(nS / 128)));
    L = floor(nS / nSeg);
    wseg = hann_local(L);
    segPhase = zeros(nCh, nSeg);
    for s = 1:nSeg
        idx = (s-1)*L+1 : s*L;
        eseg = exp(-1i * 2 * pi * freqUsed * (idx-1).' / fs);
        for c = 1:nCh
            segPhase(c, s) = angle(sum(iq(c, idx).' .* wseg .* eseg)) * 180/pi;
        end
    end
    segDphi = zeros(nCh, nSeg);
    for s = 1:nSeg
        segDphi(:, s) = wrapdeg(segPhase(:, s) - segPhase(1, s));
    end
    dphiStdDeg = zeros(nCh, 1);
    for c = 1:nCh
        dphiStdDeg(c) = circStdDeg(segDphi(c, :));
    end

    %% --- результат ----------------------------------------------------------
    dphiDeg = wrapdeg(phaseDeg - phaseDeg(1));

    res = struct('toneDetected', snrDb >= SNRMin, ...
                 'snrDb', snrDb, 'freqUsed', freqUsed, ...
                 'phaseDeg', phaseDeg, 'dphiDeg', dphiDeg, ...
                 'dphiStdDeg', dphiStdDeg, 'amp', amp, ...
                 'method', method);
end

%% --- локальные функции ------------------------------------------------------
function w = hann_local(L)
    %HANN_LOCAL  Окно Ханна, столбец L x 1 (без Signal Processing Toolbox).
    if L == 1
        w = 1;
    else
        w = 0.5 * (1 - cos(2*pi*(0:L-1)' / (L-1)));
    end
end

function f = fftfreqs_local(n, fs)
    %FFTFREQS_LOCAL  Частоты FFT в порядке fftshift: [-fs/2, fs/2).
    f = (-floor(n/2) : ceil(n/2)-1)' * (fs / n);
end

function y = wrapdeg(x)
    %WRAPDEG  Приведение углов (град) к [-180, 180).
    y = mod(x + 180, 360) - 180;
end

function s = circStdDeg(aDeg)
    %CIRCSTDDEG  Круговое СКО углов (град) через длину результирующего вектора.
    v = exp(1i * aDeg * pi / 180);
    R = abs(sum(v)) / numel(v);
    if R >= 1
        s = 0;
    else
        s = sqrt(-2 * log(R)) * 180 / pi;
    end
end

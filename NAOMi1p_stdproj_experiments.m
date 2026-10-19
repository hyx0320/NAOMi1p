function out = NAOMi1p_stdproj_experiments(varargin)
%NAOMI1P_STDPROJ_EXPERIMENTS  光学 / 背景波动 / 帧率 的调参实验 + std 投影对比
%
% 三块实验（不要合成一个大网格）
%   【optics】NA、像素采样                —— 看对 PSF 与成像的影响
%   【noise 】血管搏动 / 散射 / 血红蛋白吸收 / 读出噪声 sigma_p —— 压背景波动
%   【fps   】10 Hz vs 30 Hz              —— 等时长设计，指标用 max ΔF/F 分布
%
% 核心做法
%   optics / noise：复用参考运行的 vol_out.mat 与 S.mat（同一块脑组织、同一套活动），
%                   只重跑 simulate_1p_optical_propagation + scan_volume_1p（约 10 min/配置）
%   fps           ：两侧都用 rng(1) 重新生成 S.mat，且**等时长**
%                   （10Hz×300帧 == 30Hz×900帧，100Hz 粒度步数都是 3000 → 脉冲序列完全相同）
%   基线          ：直接读参考运行的 mov_w_bg.tiff 前 N 帧，不重跑
%
% 用法
%   NAOMi1p_stdproj_experiments                        % 默认只跑 optics+noise（约 2 h）
%   NAOMi1p_stdproj_experiments('dryrun', true)        % 只看计划
%   NAOMi1p_stdproj_experiments('group','noise')       % 只跑减背景那一块
%   NAOMi1p_stdproj_experiments('group','fps')         % 帧率专项（+40 min）
%   NAOMi1p_stdproj_experiments('group','all')         % 全部
%   NAOMi1p_stdproj_experiments('configs',{'scat0.5','sig50'})
%
% 输出（不覆盖参考运行目录）
%   D:\MesoTwins\NAOMi1p_stdproj\cfg_<名字>\       每配置的电影与 stdproj.mat
%   ...\stdproj_compare.png                        多面板 std 投影（统一色标）
%   ...\stdproj_diff.png                           基线 vs 最优 vs 差异
%   ...\experiments_summary.csv / .md              指标表 + 结论

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REF_DEFAULT = ['D:\MesoTwins\NAOMi1p_output\w_dend\res_0.80\' ...
               'vol_600_150_NA_0.35_Hz_10_exp_5_d_30k_pw_0.50\1'];

ip = inputParser;
ip.addParameter('refdir',  REF_DEFAULT, @(x) ischar(x) || isstring(x));
ip.addParameter('outroot', 'D:\MesoTwins\NAOMi1p_stdproj', @(x) ischar(x) || isstring(x));
ip.addParameter('frames',  200, @(x) isnumeric(x) && x >= 20);
ip.addParameter('group',   'main', @(x) any(strcmpi(x, {'main','optics','noise','fps','all'})));
ip.addParameter('configs', {}, @iscell);
ip.addParameter('dryrun',  false, @(x) islogical(x) || isnumeric(x));
ip.addParameter('resume',  true,  @(x) islogical(x) || isnumeric(x));
ip.addParameter('allowExpensive', false, @(x) islogical(x) || isnumeric(x));
ip.parse(varargin{:});
opt = ip.Results;
opt.refdir  = char(opt.refdir);
opt.outroot = char(opt.outroot);
nframes = opt.frames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0) 路径与输入检查
root = fileparts(mfilename('fullpath'));
addpath(genpath(root)); rehash path;
fprintf('\n========== NAOMi1p 调参实验（std 投影 / 背景波动 / 帧率）==========\n');
fprintf('脚本     : %s\n', [mfilename('fullpath') '.m']);
fprintf('参考运行 : %s\n', opt.refdir);
fprintf('输出目录 : %s\n', opt.outroot);
fprintf('基准帧数 : %d\n', nframes);

% ---- 静音一条已知无害的原版警告 ----
% genCorticalLightPathLite_1p:296-299 把 N3/2（N3 为奇数时是半整数）直接当索引用，
% 于是每个像素弹 4 条"冒号运算符需要整数操作数"。该段在一个 i×j 双重循环里，
% 600×600 网格约 9 万次迭代 → 最多 36 万条警告：刷屏、拖慢、淹没真正报错。
% MATLAB 对半整数索引本来就是四舍五入（裁剪窗口长度不变，仅整体偏移 1 个细网格点
% = 0.25 µm），所以静音这条 ID 不改变任何数值结果。跑完会恢复。
warn_state = warning('off', 'MATLAB:colon:nonIntegerIndex'); %#ok<WNOFF>

need = {'simulate_1p_optical_propagation','scan_volume_1p','fun_time_trace_generation', ...
        'masked_3DGP_v2','setCellFluoresence','RUSH_ai148d_config'};
for k = 1:numel(need)
    if isempty(which(need{k}))
        error('找不到函数 %s：请确认本脚本在项目根目录、各子目录完整。', need{k});
    end
end

ref_vol  = fullfile(opt.refdir, 'vol_out.mat');
ref_smat = fullfile(opt.refdir, 'firing_rate_0.001_smod_flag_other', 'S.mat');
ref_sopt = fullfile(opt.refdir, 'firing_rate_0.001_smod_flag_other', 'spikes_opts.mat');
ref_prm  = fullfile(opt.refdir, 'prams.mat');
ref_mov  = fullfile(opt.refdir, 'mov_w_bg.tiff');
ref_movc = fullfile(opt.refdir, 'mov_wo_bg.tiff');
ref_psf  = fullfile(opt.refdir, 'psf.tiff');
needed = {ref_vol, ref_smat, ref_sopt, ref_prm, ref_mov, ref_movc};
for k = 1:numel(needed)
    if exist(needed{k}, 'file') ~= 2, error('缺少参考文件：%s', needed{k}); end
end
if exist(opt.outroot, 'dir') ~= 7, mkdir(opt.outroot); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) 配置表
cfg = build_configs(nframes);
switch lower(opt.group)
    case 'main',   cfg = cfg(ismember({cfg.group}, {'optics','noise'}));
    case 'all',    % keep all
    otherwise,     cfg = cfg(strcmp({cfg.group}, lower(opt.group)));
end
if isempty(cfg), error('group=%s 下没有配置', opt.group); end
if ~isempty(opt.configs)
    cfg = cfg(ismember({cfg.name}, opt.configs));
    if isempty(cfg), error('指定配置名都不存在'); end
end

fprintf('\n--- 计划：%d 个配置（基线 A_base 直接读参考运行）---\n', numel(cfg));
fprintf('%-11s %-5s %-6s %-8s %-6s %-7s %-7s %-5s %-12s %s\n', ...
        'name','NA','px','sig_p','搏动','散射x','吸收x','fps','光学','说明');
for k = 1:numel(cfg)
    fprintf('%-11s %-5.2f %-6.2f %-8g %-6g %-7.2f %-7.2f %-5g %-12s %s\n', cfg(k).name, ...
        cfg(k).NA, cfg(k).pixel_size, cfg(k).sigma_p, cfg(k).dilate, ...
        cfg(k).scatter_scale, cfg(k).hemo_scale, cfg(k).frate, ...
        optics_label(cfg(k).optics_class), cfg(k).note);
end
% 输出/耗时提示
nrec = sum(strcmp({cfg.optics_class}, 'recompute'));
if nrec > 0
    if opt.allowExpensive
        fprintf('\n⚠ %d 个配置需要重算光学（每个约 4.6 h），已允许 → 总耗时可能 > %.0f 小时\n', ...
                nrec, nrec*4.6);
    else
        fprintf('\n提示：%d 个配置需要重算光学（每个约 4.6 h），默认跳过。\n', nrec);
        fprintf('      要跑请加：''allowExpensive'', true\n');
    end
end
if opt.dryrun
    fprintf('\n[dryrun] 未实际运行。\n'); out = struct('cfg', cfg); return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) 载入参考运行（复用）
fprintf('\n--- 载入参考运行（只读）---\n');
t0 = tic;
fprintf('  prams.mat ... ');       PR = load(ref_prm);  fprintf('ok\n');
fprintf('  S.mat ... ');           SS = load(ref_smat); fprintf('ok  (%d×%d)\n', size(SS.S.soma,1), size(SS.S.soma,2));
fprintf('  spikes_opts.mat ... '); SO = load(ref_sopt); fprintf('ok\n');
fprintf('  vol_out.mat ... ');     VV = load(ref_vol);  fprintf('ok  (%.0f s)\n', toc(t0));

if size(SS.S.soma,2) < nframes + 100
    error('S.mat 只有 %d 列，不足 %d；请把 frames 调小。', size(SS.S.soma,2), nframes+100);
end
kc = nframes + 100;
base_act = SS.S;
base_act.soma = base_act.soma(:,1:kc);
if isfield(base_act,'dend') && ~isempty(base_act.dend), base_act.dend = base_act.dend(:,1:kc); end
if isfield(base_act,'bg')   && ~isempty(base_act.bg),   base_act.bg   = base_act.bg(:,  1:kc); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) 基线
fprintf('\n--- 基线 A_base：读参考运行前 %d 帧 ---\n', nframes);
base = rep_template();
base.name='A_base'; base.NA=0.35; base.pixel_size=0.8; base.sigma_p=200;
base.dilate=10; base.exp_level=5; base.scatter_scale=1; base.hemo_scale=1;
base.frate=10; base.note='参考运行（原版代码）';

[mu_b, sd_b, ~, nb] = tiff_stats(ref_mov, nframes);
[~, ~, mx_b] = tiff_stats(ref_movc, nframes);
base.frames = nb;
base.image_size = [size(sd_b,1), size(sd_b,2)];
base.metrics = compute_metrics(sd_b, mx_b, mu_b, SS.S.soma(:,1:kc));
base.sd_map = sd_b;
if exist(ref_psf,'file') == 2
    base.psf_lat_um = psf_lateral_fwhm(read_tiff_stack(ref_psf));
end
fprintf('  bg/sig=%.3f | std对比度=%.2f | PSF横向FWHM=%.1f um | maxΔF/F p90=%.3f\n', ...
        base.metrics.bg_over_sig, base.metrics.std_contrast, base.psf_lat_um, base.metrics.maxdff_p90);

vmax = pct(sd_b(:), 99.5); if ~(vmax > 0), vmax = 1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4) 逐个跑
R = base;
for k = 1:numel(cfg)
    c = cfg(k);
    outdir = fullfile(opt.outroot, ['cfg_' c.name]);
    cache  = fullfile(outdir, 'stdproj.mat');
    if opt.resume && exist(cache,'file')
        fprintf('\n[%d/%d] %s 已完成，跳过\n', k, numel(cfg), c.name);
        D = load(cache); R(end+1) = D.rep; %#ok<AGROW>
        continue;
    end
    if strcmp(c.optics_class,'recompute') && ~opt.allowExpensive
        fprintf('\n[%d/%d] %s 跳过 —— 需要重算光学（约 4.6 h）\n', k, numel(cfg), c.name);
        fprintf('        要跑请加参数：''allowExpensive'', true\n');
        continue;
    end
    fprintf('\n[%d/%d] %s —— %s\n', k, numel(cfg), c.name, c.note);
    tk = tic;
    try
        rep = run_one(c, PR, VV, base_act, SO.spike_opts, outdir, opt);
        rep.elapsed_s = toc(tk);
        fprintf('  ✔ %.0f s | bg/sig=%.3f | std对比度=%.2f | PSF=%.1f um | maxΔF/F p90=%.3f\n', ...
            rep.elapsed_s, rep.metrics.bg_over_sig, rep.metrics.std_contrast, ...
            rep.psf_lat_um, rep.metrics.maxdff_p90);
        save(cache, 'rep');
        R(end+1) = rep; %#ok<AGROW>
    catch ME
        fprintf('  ✘ 失败：%s\n', ME.message);
        fid = fopen(fullfile(opt.outroot,'errors.log'), 'a');
        if fid > 0
            fprintf(fid, '=== %s | %s ===\n%s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'), c.name, ME.message);
            for ii = 1:min(6,numel(ME.stack))
                fprintf(fid, '    at %s (line %d)\n', ME.stack(ii).name, ME.stack(ii).line);
            end
            fclose(fid);
        end
        fprintf('  详见 %s\n', fullfile(opt.outroot,'errors.log'));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5) 出图与报告
warning(warn_state);                       % 恢复警告设置
fprintf('\n--- 生成图与报告 ---\n');
make_plots(R, vmax, opt.outroot);
write_summary(R, opt.outroot, opt.refdir);
fprintf('\n=============== 完成 ===============\n');
fprintf('  %s\n  %s\n  %s\n  %s\n\n', ...
    fullfile(opt.outroot,'stdproj_compare.png'), fullfile(opt.outroot,'stdproj_diff.png'), ...
    fullfile(opt.outroot,'experiments_summary.csv'), fullfile(opt.outroot,'experiments_summary.md'));
out = struct('reports', R, 'vmax', vmax, 'outroot', opt.outroot);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ======================= 记录模板（字段顺序必须一致）=======================
function r = rep_template()
m0 = struct('sig_std',NaN,'bg_std',NaN,'bg_over_sig',NaN,'std_contrast',NaN, ...
            'mean_pixel',NaN,'max_pixel',NaN,'range_used',NaN,'sat_frac',NaN, ...
            'maxdff_p50',NaN,'maxdff_p90',NaN);
r = struct('name','', 'NA',NaN, 'pixel_size',NaN, 'sigma_p',NaN, 'dilate',NaN, ...
           'exp_level',NaN, 'scatter_scale',NaN, 'hemo_scale',NaN, 'frate',NaN, ...
           'note','', 'frames',NaN, 'image_size',[NaN NaN], 'psf_lat_um',NaN, ...
           'metrics',m0, 'sd_map',[], 'elapsed_s',NaN);
end

%% ======================= 配置表 =======================
function C = build_configs(nframes)
C = [];
% ---- optics：光学参数 ----
C = add(C,'optics','NA0.5','NA',0.50, 'note','光学：NA 0.35→0.50');
C = add(C,'optics','NA0.7','NA',0.70, 'note','光学：NA 0.35→0.70');
C = add(C,'optics','px1.6','pixel_size',1.6, 'note','采样：0.8→1.6 um/px');
% ---- noise：压背景波动 ----
C = add(C,'noise','vd2',    'dilate',2,         'note','血管搏动 10→2');
C = add(C,'noise','vd0',    'vessel_off',true,  'note','关闭血管调制（vessel_img=1，静态遮挡也一起去掉）');
C = add(C,'noise','sig100', 'sigma_p',100,      'note','读出噪声 200→100');
C = add(C,'noise','sig50',  'sigma_p',50,       'note','读出噪声 200→50');
C = add(C,'noise','scat0.5','scatter_scale',0.5,'note','散射强度 x0.5（压弥漫背景）');
C = add(C,'noise','scat2',  'scatter_scale',2.0,'note','散射强度 x2（对照，应更差）');
C = add(C,'noise','hemo0',  'hemo_scale',0,     'note','血红蛋白吸收归零（去黑阴影）');
C = add(C,'noise','combo',  'combo',true,       'note','搏动2+sigma100+exp15（含散射会多花4.6h，故不含）');
% ---- fps：帧率（等时长；两侧都重生成 S）----
C = add(C,'fps','fr10r','frate',10,'frames',nframes,'regen_S',true, ...
        'note', sprintf('%d 帧 @10Hz（重生成 S，%.0f s）', nframes, (nframes+100)/10));
nf30 = round((nframes+100)/10*30) - 100;      % 等时长换算
C = add(C,'fps','fr30r','frate',30,'frames',nf30,'regen_S',true, ...
        'note', sprintf('%d 帧 @30Hz（等时长 %.0f s，约 3x 成本）', nf30, (nf30+100)/30));
end

function C = add(C, grp, name, varargin)
c = struct('name',name,'group',grp,'NA',0.35,'pixel_size',0.8,'sigma_p',200, ...
           'dilate',10,'exp_level',5,'vessel_off',false,'scatter_scale',1, ...
           'hemo_scale',1,'frate',10,'regen_S',false,'frames',200,'note','', ...
           'optics_class','reuse');
isCombo = false; k = 1;
while k <= numel(varargin)
    f = varargin{k}; v = varargin{k+1};
    if strcmp(f,'combo')
        isCombo = v;
    else
        c.(f) = v;
    end
    k = k + 2;
end
if isCombo
    % 注意：故意不含散射 —— 改 scatter_wt 会强制重算光学（约 4.6 h），
    % 而散射已由 scat0.5 单独测过。这里只组合"免费"的旋钮。
    c.dilate = 2; c.sigma_p = 100; c.exp_level = 15;
end
% ---- 判断本配置是否需要重算光学 ----
% psf / psfT / psfB / mask 同时受 NA 与 scatter_wt 影响 → 必须重算
% colmask 只是 10.^(-A*hemoabs) → 换 hemoabs 等价于对 colmask 取幂（精确、免重算）
if (abs(c.NA - 0.35) > 1e-9) || (abs(c.scatter_scale - 1) > 1e-9)
    c.optics_class = 'recompute';
elseif abs(c.hemo_scale - 1) > 1e-9
    c.optics_class = 'hemo';
else
    c.optics_class = 'reuse';
end
if isempty(C), C = c; else, C(end+1) = c; end %#ok<AGROW>
end

function s = optics_label(cls)
switch cls
    case 'reuse', s = '复用缓存';
    case 'hemo',  s = '变换colmask';
    otherwise,    s = '重算≈4.6h';
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ======================= 跑一个配置 =======================
function rep = run_one(c, PR, VV, base_act, spike_opts_in, outdir, opt)
if exist(outdir,'dir') ~= 7, mkdir(outdir); end

% 每个配置从同一随机状态出发 → 散射介质、相机噪声的"随机实现"在各配置间一致，
% 差异只来自被改的那个参数（散射强度不同时介质形状相同、仅幅度不同）。
rng(1, 'twister');

vol_params = struct; neur_params = struct; vasc_params = struct;
dend_params = struct; bg_params = struct; axon_params = struct;
psf_params = struct; scan_params = struct; noise_params = struct;
wdm_params = struct; spike_opts = struct; vessel_mod = struct;
RUSH_ai148d_config;

vol_params = PR.vol_params;                        % 与 vol_out 完全一致

% ---- 光学覆盖 ----
psf_params.NA    = c.NA;
psf_params.objNA = c.NA;
if isfield(psf_params,'scatter_wt') && ~isempty(psf_params.scatter_wt)
    psf_params.scatter_wt = psf_params.scatter_wt * c.scatter_scale;
end
if isfield(psf_params,'hemoabs') && ~isempty(psf_params.hemoabs)
    psf_params.hemoabs = psf_params.hemoabs * c.hemo_scale;
end

% ---- 传感器 / 成像覆盖 ----
noise_params.sigma_p = c.sigma_p;
pixel_size = c.pixel_size;
exp_level  = c.exp_level;

vessel_mod.flag = true;
vessel_mod.frate = c.frate;
vessel_mod.vessel_dynamics = 2;
vessel_mod.FOV_N = 2;
vessel_mod.max_dilate_amp = max(c.dilate, 1);
vessel_mod.seed = 10;
if c.vessel_off, vessel_mod.flag = false; end

% ---- 活动时序：帧率实验需重生成（同种子 + 等时长 → 脉冲序列一致）----
if c.regen_S
    fprintf('  重生成时序（%g Hz, %d 帧）... ', c.frate, c.frames);
    t0 = tic;
    rng(1, 'twister');
    fun_time_trace_generation(VV.vol_out, c.frames, c.frate, outdir);
    sp = dir(fullfile(outdir,'firing_rate_*'));
    if isempty(sp), error('未生成 firing_rate_* 目录'); end
    Sx = load(fullfile(outdir, sp(1).name, 'S.mat'));
    neur_act = Sx.S;
    so = load(fullfile(outdir, sp(1).name, 'spikes_opts.mat'));
    spike_opts = so.spike_opts;
    fprintf('%.0f s（%d 列）\n', toc(t0), size(neur_act.soma,2));
else
    neur_act   = base_act;
    spike_opts = spike_opts_in;
end

% ---- ① 光学：复用 / 代数变换 / 重算（见 get_psf_struct）----
fprintf('  光学 ... '); t1 = tic;
PSF_struct = get_psf_struct(c, psf_params, vol_params, VV.vol_out, opt);
fprintf('%.0f s\n', toc(t1));

% ---- ② 扫描 ----
fprintf('  扫描成像 ... '); t2 = tic;
scan_volume_1p(VV.vol_out, PSF_struct, neur_act, vol_params, scan_params, ...
               noise_params, spike_opts, wdm_params, vessel_mod, ...
               pixel_size, exp_level, outdir);
fprintf('%.0f s\n', toc(t2));

% ---- ③ std 投影与指标 ----
fprintf('  统计 std 投影 ... ');
[mu_w, sd_w, ~, nw] = tiff_stats(fullfile(outdir,'mov_w_bg.tiff'),  []);
[~,    ~,    mx_c] = tiff_stats(fullfile(outdir,'mov_wo_bg.tiff'), []);
fprintf('ok（%d 帧）\n', nw);

rep = rep_template();
rep.name = c.name;  rep.NA = c.NA;  rep.pixel_size = c.pixel_size;
rep.sigma_p = c.sigma_p;  rep.dilate = c.dilate;  rep.exp_level = c.exp_level;
rep.scatter_scale = c.scatter_scale;  rep.hemo_scale = c.hemo_scale;
rep.frate = c.frate;  rep.note = c.note;  rep.frames = nw;
rep.image_size = [size(sd_w,1), size(sd_w,2)];
rep.psf_lat_um = psf_lateral_fwhm(PSF_struct.psf);
rep.metrics = compute_metrics(sd_w, mx_c, mu_w, neur_act.soma);
rep.sd_map = sd_w;

for sub = {'debug','debug_gt'}
    dd = fullfile(outdir, sub{1});
    if exist(dd,'dir') == 7 && ~isempty(strfind(dd, outdir))
        try, rmdir(dd,'s'); catch, end %#ok<CTCH>
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ======================= 光学结果：复用 / 变换 / 重算 =======================
function PSF_struct = get_psf_struct(c, psf_params, vol_params, vol_out, opt)
% 为什么能复用：
%   PSF_struct 里 scan_volume_1p 只用 6 个字段 ——
%     psf / psfT / psfB / mask / colmask / blur
%   其中 psf*、mask 只由 (NA, scatter_wt, 几何) 决定；
%   colmask = 10.^(-A*hemoabs) 只由 hemoabs 决定，且是纯代数关系。
%   pixel_size / sigma_p / vessel_mod / exp_level 完全不进入光学计算。
% 所以：
%   NA、scatter_wt 改变   → 必须重算（~4.6 h/次）
%   只改 hemoabs          → colmask.^k 精确等价（秒级）
%   其它                  → 直接读缓存（秒级）
keyname = sprintf('PSF_NA%g_scat%g.mat', c.NA, c.scatter_scale);
cachedir = fullfile(opt.outroot, '_psf_cache');
if exist(cachedir,'dir') ~= 7, mkdir(cachedir); end
cachefile = fullfile(cachedir, keyname);

switch c.optics_class
    case {'reuse','hemo'}
        if exist(cachefile,'file') == 2
            L = load(cachefile); PSF_struct = L.PSF_struct;
            fprintf('(缓存 %s) ', keyname);
        elseif abs(c.NA-0.35) < 1e-9 && abs(c.scatter_scale-1) < 1e-9 ...
                && exist(fullfile(opt.refdir,'PSF_struct.mat'),'file') == 2
            % 首次：从参考运行借用（同一套 NA=0.35 / 散射=1 的光学结果）
            L = load(fullfile(opt.refdir,'PSF_struct.mat')); PSF_struct = L.PSF_struct;
            PSF_struct_cache_save(cachefile, PSF_struct);
            fprintf('(借用参考运行 PSF_struct.mat) ');
        else
            error(['缺少光学缓存 %s。\n' ...
                   '  该配置的分类是 %s，本应可复用，但没有对应缓存。\n' ...
                   '  请先跑一次 NA=0.35 / 散射=1 的配置生成缓存，或加 ''allowExpensive'',true 重算。'], ...
                   cachefile, c.optics_class);
        end
        if strcmp(c.optics_class,'hemo')
            % colmask_new = (10^(-A*hemo_old))^(hemo_new/hemo_old)
            PSF_struct.colmask = PSF_struct.colmask .^ c.hemo_scale;
            fprintf('(colmask^%g) ', c.hemo_scale);
        end
    otherwise   % 'recompute'
        fprintf('(重算，约 4.6 h) ');
        PSF_struct = simulate_1p_optical_propagation(vol_params, psf_params, vol_out);
        PSF_struct_cache_save(cachefile, PSF_struct);
end
end

function PSF_struct_cache_save(f, PSF_struct) %#ok<INUSD>
try
    save(f, 'PSF_struct', '-v7.3');
catch ME
    warning('光学缓存写入失败（不影响本次运行）：%s', ME.message);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ======================= 指标 =======================
function m = compute_metrics(sd, mx_clean, mu_noisy, soma_traces)
v = sort(mx_clean(:));
thr = v(max(1, round(0.90*numel(v))));
sig_mask = mx_clean > thr;

m = struct();
m.sig_std      = median(sd(sig_mask));
m.bg_std       = median(sd(~sig_mask));
m.bg_over_sig  = m.bg_std / max(m.sig_std, 1e-9);          % ← 越小越好
m.std_contrast = pct(sd(:), 99) / max(median(sd(:)), 1e-9);
m.mean_pixel   = mean(mu_noisy(:));
m.max_pixel    = max(mu_noisy(:));
m.range_used   = double(max(mu_noisy(:))) / 65535;
m.sat_frac     = mean(mu_noisy(:) >= 65000);

% max ΔF/F 分布（对应论文图里的 Max ΔF/F histogram；帧率实验的关键指标）
m.maxdff_p50 = NaN; m.maxdff_p90 = NaN;
if nargin >= 4 && ~isempty(soma_traces)
    F  = double(soma_traces);
    F0 = median(F, 2);
    ok = F0 > 1e-9;
    if any(ok)
        dff = (F(ok,:) - F0(ok)) ./ F0(ok);
        mx  = max(dff, [], 2);
        m.maxdff_p50 = pct(mx, 50);
        m.maxdff_p90 = pct(mx, 90);
    end
end
end

function w = psf_lateral_fwhm(psf)
psf = double(psf);
prof = squeeze(sum(sum(psf,2),1)).';
if isempty(prof) || max(prof) <= 0, w = NaN; return; end
prof = prof / max(prof);
[~, zpk] = max(prof);
pl = double(psf(:,:,zpk));
[~, imax] = max(pl(:));
[ry, cx] = ind2sub(size(pl), imax);
wx = fwhm1(pl(ry,:));  wy = fwhm1(pl(:,cx).');
w = mean([wx wy], 'omitnan');
if ~isfinite(w), w = NaN; end
end

function w = fwhm1(v)
v = double(v(:).');
if isempty(v) || max(v) <= 0, w = NaN; return; end
v = v / max(v);
g = v > 0.5;  f = find(diff(g));
if numel(f) < 2, w = NaN; return; end
s1 = (v(f(1)+1)-0.5)/(v(f(1)+1)-v(f(1)));
s2 = (v(f(2))-0.5)/(v(f(2))-v(f(2)+1));
w = s1 + s2 + sum(g);
end

function q = pct(x, p)
x = sort(double(x(:)));  x = x(~isnan(x));
if isempty(x), q = NaN; return; end
q = x(min(numel(x), max(1, round(p/100*numel(x)))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ======================= TIFF 工具 =======================
function [mu, sd, mx, n] = tiff_stats(fname, nframes)
% 单次读盘同时算 均值/标准差/最大值，避免把整部电影读进内存
t = Tiff(fname,'r');  c = onCleanup(@() t.close());
s1=[]; s2=[]; mx=[]; k=0;
while true
    if k > 0
        if t.lastDirectory(), break; end
        t.nextDirectory();
    end
    if ~isempty(nframes) && k >= nframes, break; end
    img = double(t.read());  k = k + 1;
    if k == 1, s1 = zeros(size(img)); s2 = zeros(size(img)); mx = img; end
    s1 = s1 + img;  s2 = s2 + img.*img;  mx = max(mx, img);
end
n = k;
if k == 0, mu=[]; sd=[]; mx=[]; return; end
mu = s1/k;  sd = sqrt(max(s2/k - mu.*mu, 0));
end

function psf = read_tiff_stack(fname)
t = Tiff(fname,'r');  c = onCleanup(@() t.close());
psf = []; k = 0;
while true
    if k > 0
        if t.lastDirectory(), break; end
        t.nextDirectory();
    end
    img = t.read();  k = k + 1;
    if k == 1, psf = zeros([size(img,1),size(img,2),1],'like',img); end
    psf(:,:,k) = img;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ======================= 出图 =======================
function make_plots(R, vmax, outroot)
n = numel(R);  cols = ceil(sqrt(n));  rows = ceil(n/cols);  tgt = size(R(1).sd_map);

f1 = figure('Position',[60 60 430*cols 400*rows],'Color','w');
for k = 1:n
    subplot(rows,cols,k);
    img = resize_to(clamp_map(R(k).sd_map, vmax), tgt);
    imagesc(img); axis image off; colormap(parula); caxis([0 vmax]);
    title(sprintf('%s\nbg/sig=%.2f', R(k).name, R(k).metrics.bg_over_sig), ...
          'Interpreter','none','FontSize',9);
end
cb = colorbar('Position',[0.94 0.1 0.015 0.8]);  ylabel(cb,'时间 std (a.u.)');
sgtitle(sprintf('时间标准差（std）投影对比 — 统一色标 0~%.0f，各 %d 帧', vmax, R(1).frames), ...
        'FontSize',12,'FontWeight','bold');
p1 = fullfile(outroot,'stdproj_compare.png');
print(f1,'-dpng','-r150',p1);  close(f1);  fprintf('  已写 %s\n', p1);

ib = find(strcmp({R.name},'A_base'),1);  if isempty(ib), ib = 1; end
vals = arrayfun(@(x) x.metrics.bg_over_sig, R);  [~, ibest] = min(vals);
A = resize_to(clamp_map(R(ib).sd_map, vmax), tgt);
B = resize_to(clamp_map(R(ibest).sd_map, vmax), tgt);
D = A - B;  dm = max(abs(D(:)));  if ~(dm > 0), dm = 1; end

f2 = figure('Position',[80 80 1500 470],'Color','w');
subplot(1,3,1); imagesc(A); axis image off; colormap(parula); caxis([0 vmax]);
title(sprintf('基线 %s  bg/sig=%.2f', R(ib).name, R(ib).metrics.bg_over_sig), ...
      'Interpreter','none','FontSize',10); colorbar;
subplot(1,3,2); imagesc(B); axis image off; caxis([0 vmax]);
title(sprintf('最优 %s  bg/sig=%.2f', R(ibest).name, R(ibest).metrics.bg_over_sig), ...
      'Interpreter','none','FontSize',10); colorbar;
subplot(1,3,3); imagesc(D); axis image off; colormap(jet); caxis([-dm dm]);
title('差异（基线 − 最优）：蓝=波动被压低','FontSize',10); colorbar;
sgtitle('std 投影：基线 vs 最优 vs 差异图','FontSize',12,'FontWeight','bold');
p2 = fullfile(outroot,'stdproj_diff.png');
print(f2,'-dpng','-r150',p2);  close(f2);  fprintf('  已写 %s\n', p2);
end

function img = clamp_map(x, vmax)
if ndims(x) == 3, x = mean(x,3); end
img = min(double(x), vmax);
end

function Y = resize_to(X, sz)
if isequal(size(X), sz), Y = X; else, Y = imresize(X, sz); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ======================= 报告 =======================
function write_summary(R, outroot, refdir)
fname = fullfile(outroot,'experiments_summary.csv');
fid = fopen(fname,'w');
if fid < 0, warning('无法写 %s', fname); return; end
fprintf(fid, ['name,NA,pixel_size_um,sigma_p,max_dilate_amp,exp_level,' ...
    'scatter_scale,hemo_scale,frate,' ...
    'bg_over_sig,sig_std,bg_std,std_contrast,psf_lat_fwhm_um,' ...
    'mean_pixel,range_used,sat_frac,maxdff_p50,maxdff_p90,frames,elapsed_s,note\n']);
for k = 1:numel(R)
    r = R(k);
    fprintf(fid, ['%s,%.3g,%.3g,%g,%g,%g,%.3g,%.3g,%g,' ...
        '%.4f,%.4f,%.4f,%.4f,%.3f,%.1f,%.4f,%.5f,%.4f,%.4f,%d,%.0f,"%s"\n'], ...
        r.name, r.NA, r.pixel_size, r.sigma_p, r.dilate, r.exp_level, ...
        r.scatter_scale, r.hemo_scale, r.frate, ...
        r.metrics.bg_over_sig, r.metrics.sig_std, r.metrics.bg_std, ...
        r.metrics.std_contrast, r.psf_lat_um, r.metrics.mean_pixel, ...
        r.metrics.range_used, r.metrics.sat_frac, ...
        r.metrics.maxdff_p50, r.metrics.maxdff_p90, ...
        r.frames, r.elapsed_s, strrep(r.note, ',', ';'));
end
fclose(fid);  fprintf('  已写 %s\n', fname);

mname = fullfile(outroot,'experiments_summary.md');
fid = fopen(mname,'w');
fprintf(fid, '# NAOMi1p 调参实验小结\n\n');
fprintf(fid, '- 参考运行：`%s`\n', refdir);
fprintf(fid, '- optics / noise 块：同一份 `vol_out.mat` 与 `S.mat`（严格控制变量）\n');
fprintf(fid, '- fps 块：两侧均 `rng(1)` 重生成 S，且等时长（10Hz×300帧 == 30Hz×900帧，100Hz 粒度步数同为 3000）\n');
fprintf(fid, '- `bg/sig` = 背景区时间std中位数 ÷ 信号区时间std中位数，**越小越好**\n\n');
fprintf(fid, '| 配置 | NA | px | sigma_p | 搏动 | 散射x | 吸收x | fps | bg/sig | std对比度 | PSF(um) | maxΔF/F p90 |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for k = 1:numel(R)
    r = R(k);
    fprintf(fid, '| %s | %.2f | %.2f | %g | %g | %.2f | %.2f | %g | **%.3f** | %.2f | %.1f | %.3f |\n', ...
        r.name, r.NA, r.pixel_size, r.sigma_p, r.dilate, r.scatter_scale, ...
        r.hemo_scale, r.frate, r.metrics.bg_over_sig, r.metrics.std_contrast, ...
        r.psf_lat_um, r.metrics.maxdff_p90);
end
[~, ib] = min(arrayfun(@(x) x.metrics.bg_over_sig, R));
fprintf(fid, '\n**bg/sig 最低**：`%s`（%.3f）\n', R(ib).name, R(ib).metrics.bg_over_sig);

% 帧率专项
is_fps = strncmp({R.name}, 'fr', 2);
if any(is_fps)
    fprintf(fid, '\n### 帧率专项（等时长）\n\n');
    fprintf(fid, '| 配置 | fps | 帧数 | maxΔF/F p50 | maxΔF/F p90 | std对比度 |\n');
    fprintf(fid, '|---|---|---|---|---|---|\n');
    for k = 1:numel(R)
        if is_fps(k) || strcmp(R(k).name,'A_base')
            r = R(k);
            fprintf(fid, '| %s | %g | %d | %.4f | %.4f | %.2f |\n', ...
                r.name, r.frate, r.frames, r.metrics.maxdff_p50, ...
                r.metrics.maxdff_p90, r.metrics.std_contrast);
        end
    end
    fprintf(fid, '\n> 若 10 Hz 的 maxΔF/F 明显低于 30 Hz，说明**欠采样削低了峰值**（抓不到高频钙瞬变）。\n');
end
fclose(fid);  fprintf('  已写 %s\n', mname);
end

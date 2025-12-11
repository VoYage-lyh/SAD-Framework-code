function [analysis_params, sim_params] = ConfigAdapter(preConfig, identifiedParams)
% ConfigAdapter_v2 - 参数适配器（修正版）
%
% 工作流程：
%   1. preConfig: 来自GUI的预配置（拓扑、几何、质量）
%   2. identifiedParams: 来自参数识别的结果（刚度、阻尼、非线性）- 可选
%
% 输出：
%   analysis_params: 用于参数识别代码的参数
%   sim_params: 用于仿真代码的参数（包含识别结果）

    if nargin < 2
        identifiedParams = [];
    end
    
    %% 1. 严格性检查：必须提供参数识别结果
    if isempty(identifiedParams)
        error('ConfigAdapter:NoIdentifiedParams', ...
              '错误：未提供 identifiedParams（参数识别结果）。\n' + ...
              '为了避免硬编码数值，仿真必须基于 analyse_chibi_data.m 的识别结果运行。');
    end
    
    %% 2. 生成参数识别所需参数 (保持不变)
    analysis_params = struct();
    analysis_params.fs_target = preConfig.signal.fs_target;
    analysis_params.cutoff_freq = preConfig.signal.cutoff_freq;
    analysis_params.filter_order = preConfig.signal.filter_order;
    analysis_params.nfft = preConfig.signal.nfft;
    analysis_params.freq_range = [preConfig.signal.freq_range_min, preConfig.signal.freq_range_max];
    analysis_params.snr_threshold = preConfig.signal.snr_threshold;
    analysis_params.n_nodes = 3;
    analysis_params.n_dof = 6;
    analysis_params.node_labels = {'Root', 'Mid', 'Tip'};
    analysis_params.direction_labels = {'Y', 'Z'};
    analysis_params.topology = preConfig.topology;
    
    %% 3. 生成仿真所需参数
    sim_params = struct();
    
    % 基础设置
    sim_params.parallel_execution_max_workers = preConfig.basic.parallel_max_workers;
    sim_params.workFolder = preConfig.basic.workFolder;
    sim_params.model_name = preConfig.basic.modelName;
    sim_params.gravity_g = preConfig.basic.gravity_g;
    sim_params.use_parallel = preConfig.basic.useParallel;
    sim_params.config = preConfig.topology;
    
    % --- 构建主干参数 (严格依赖识别结果) ---
    sim_params.trunk = buildTrunkParams(preConfig, identifiedParams);
    
    % --- 构建分枝参数 (严格依赖识别结果 + 激活完整性验证) ---
    sim_params.predefined_params = generatePredefinedParams(preConfig, identifiedParams);
    
    % --- 构建果实配置 ---
    sim_params.fruit_config = struct();
    sim_params.fruit_config.attach_secondary_mid = preConfig.fruit.attach_secondary_mid;
    sim_params.fruit_config.attach_secondary_tip = preConfig.fruit.attach_secondary_tip;
    sim_params.fruit_config.attach_tertiary_mid = preConfig.fruit.attach_tertiary_mid;
    sim_params.fruit_config.attach_tertiary_tip = preConfig.fruit.attach_tertiary_tip;
    
    % [修复] 使用 GUI 配置的每节点果实数量，移除硬编码
    if isfield(preConfig.fruit, 'fruits_per_node')
        sim_params.fruit_config.fruits_per_node = preConfig.fruit.fruits_per_node;
    else
        sim_params.fruit_config.fruits_per_node = 1;
    end
    
    % 生成默认果实参数 (用于未特定指明位置的果实)
    sim_params.default_fruit_params = buildFruitParamsStrict(preConfig, identifiedParams);
    
    % 激励参数
    sim_params.excitation = preConfig.excitation;
    
    % [修复] 更新激励频率为第一阶固有频率 (优先从Trunk获取，适配聚合结构体)
    found_freq = false;
    if isfield(identifiedParams, 'branches') && isfield(identifiedParams.branches, 'Trunk') && ...
       isfield(identifiedParams.branches.Trunk, 'linear') && ...
       isfield(identifiedParams.branches.Trunk.linear, 'natural_freqs_x')
        freqs = identifiedParams.branches.Trunk.linear.natural_freqs_x;
        if ~isempty(freqs)
            sim_params.excitation.frequency_hz = freqs(1);
            found_freq = true;
        end
    end
    
    % 兼容旧版单一结构体
    if ~found_freq && isfield(identifiedParams, 'linear') && ...
       isfield(identifiedParams.linear, 'natural_freqs_x')
        freqs = identifiedParams.linear.natural_freqs_x;
        if ~isempty(freqs)
            sim_params.excitation.frequency_hz = freqs(1);
        end
    end

    % 仿真控制
    sim_params.sim_stop_time = preConfig.simulation.stop_time;
    sim_params.sim_fixed_step = preConfig.simulation.fixed_step;
    sim_params.has_identified_params = true;
    
    fprintf('ConfigAdapter: 参数适配完成。模型参数已基于实验数据生成。\n');
end

%% ==================== 构建主干参数 ====================
function trunk = buildTrunkParams(preConfig, identifiedParams)
    trunk = struct();
    
    % --- 几何与质量 (保持不变) ---
    trunk.length = preConfig.trunk.length;
    trunk.diameter_base = preConfig.trunk.diameter_base;
    trunk.diameter_tip = preConfig.trunk.diameter_tip;
    % [移除] trunk.z_factor = preConfig.trunk.z_factor; <--- z_factor 不再来源于配置
    
    m_total = preConfig.trunk.total_mass;
    m_dist = preConfig.trunk.mass_distribution;
    trunk.root.m = m_total * m_dist(1);
    trunk.mid.m = m_total * m_dist(2);
    trunk.tip.m = m_total * m_dist(3);
    
    % --- 刚度和阻尼 (严格数据驱动) ---
    
    % 1. 定位线性参数源
    if isfield(identifiedParams, 'branches') && isfield(identifiedParams.branches, 'Trunk')
        target_lin = identifiedParams.branches.Trunk.linear;
        fprintf('    主干：加载 Trunk 分枝识别参数。\n');
    elseif isfield(identifiedParams, 'linear')
        target_lin = identifiedParams.linear;
        fprintf('    主干：加载全局线性参数。\n');
    else
        error('ConfigAdapter:NoTrunkParams', ...
              '错误：未在 identifiedParams 中找到主干(Trunk)或全局(linear)参数。\n请先运行参数识别。');
    end
    
    % 2. 严格检查 Y 方向物理参数向量 (identified_params_x)
    % 格式应为: [k_g, c_g, k_rm, c_rm, k_mt, c_mt]
    if ~isfield(target_lin, 'identified_params_x') || isempty(target_lin.identified_params_x) || length(target_lin.identified_params_x) < 6
        error('ConfigAdapter:PhysicalParamsMissingY', ...
              '严重错误：参数识别结果中缺少有效的 Y 方向物理参数向量 "identified_params_x"。');
    end
    
    % 3. 赋值 Y 方向参数 (Y方向 = 主振动方向)
    p_vec_y = target_lin.identified_params_x;
    
    trunk.root.k_y_conn_to_base = p_vec_y(1); % k_g
    trunk.root.c_y_conn_to_base = p_vec_y(2); % c_g
    
    trunk.root.k_y_conn = p_vec_y(3); % k_rm
    trunk.mid.k_y_conn  = p_vec_y(5); % k_mt
    trunk.tip.k_y_conn  = 0;        
    
    trunk.root.c_y_conn = p_vec_y(4);
    trunk.mid.c_y_conn  = p_vec_y(6);
    trunk.tip.c_y_conn  = 0;

    fprintf('    [严格模式] Y 方向参数应用: k_g=%.1f, k_rm=%.1f。\n', p_vec_y(1), p_vec_y(3));

    
    % 4. Z 方向参数处理 (强制要求识别结果)
    
    if isfield(target_lin, 'identified_params_z') && length(target_lin.identified_params_z) >= 6
        % Z-A: 使用 Z 方向独立识别参数 (p_vec_z)
        p_vec_z = target_lin.identified_params_z;
        fprintf('    [严格模式] 应用 Z 方向独立识别参数。\n');
        
        % 4.1 自动计算 z_factor 并记录 (作为识别结果输出)
        if p_vec_y(1) ~= 0
            % 基于接地刚度 k_g 计算 z_factor
            z_fac = p_vec_z(1) / p_vec_y(1); 
            trunk.z_factor = z_fac; 
            fprintf('    [自动计算] Z-factor (k_g,z / k_g,y) = %.4f\n', z_fac);
        else
            % 特殊情况处理
            trunk.z_factor = 1.0; 
            warning('ConfigAdapter:ZeroKG', 'Y方向 k_g 接近零，Z-factor 设为 1.0。');
        end

        % 4.2 赋值 Z 方向参数
        trunk.root.k_z_conn_to_base = p_vec_z(1);
        trunk.root.c_z_conn_to_base = p_vec_z(2);
        
        trunk.root.k_z_conn = p_vec_z(3);
        trunk.mid.k_z_conn  = p_vec_z(5);
        trunk.tip.k_conn_z  = 0;
        
        trunk.root.c_z_conn = p_vec_z(4);
        trunk.mid.c_z_conn  = p_vec_z(6);
        trunk.tip.c_z_conn  = 0;
        
    else
        % Z-B: Z方向独立识别参数缺失 -> 报错中断
        error('ConfigAdapter:NoZParams', ...
              ['严重错误：缺少 Z 方向识别参数 (identified_params_z)。\n' ...
               '在严格模式下，Z方向参数必须由数据识别得出，禁止使用估算值或几何因子。\n' ...
               '请检查分析步骤是否包含 Z 方向识别。']);
    end
end

%% ==================== 验证预定义参数完整性 ====================
function validatePredefinedParams(predefined, preConfig)
    % 验证所有必需的分枝参数是否都已生成
    
    missingBranches = {};
    
    % 检查一级分枝
    num_p = preConfig.topology.num_primary_branches;
    for p = 1:num_p
        branch_id = sprintf('P%d', p);
        if ~isfield(predefined, branch_id)
            missingBranches{end+1} = branch_id;
        end
    end
    
    % 检查二级分枝
    for p = 1:num_p
        num_s = preConfig.topology.secondary_branches_count(p);
        for s = 1:num_s
            branch_id = sprintf('P%d_S%d', p, s);
            if ~isfield(predefined, branch_id)
                missingBranches{end+1} = branch_id;
            end
        end
    end
    
    % 检查三级分枝
    for p = 1:num_p
        if p <= length(preConfig.topology.tertiary_branches_count)
            tertiary_for_p = preConfig.topology.tertiary_branches_count{p};
            num_s = preConfig.topology.secondary_branches_count(p);
            for s = 1:num_s
                if s <= length(tertiary_for_p)
                    num_t = tertiary_for_p(s);
                    for t = 1:num_t
                        branch_id = sprintf('P%d_S%d_T%d', p, s, t);
                        if ~isfield(predefined, branch_id)
                            missingBranches{end+1} = branch_id;
                        end
                    end
                end
            end
        end
    end
    
    % 报告缺失的分枝
    if ~isempty(missingBranches)
        error('ConfigAdapter:MissingData', ...
              '以下分枝在预配置中缺少几何参数，请在GUI中完整配置:\n  %s', ...
              strjoin(missingBranches, ', '));
    end
end

%% ==================== 为分枝生成默认段参数 ====================
function predefined = generatePredefinedParams(preConfig, identifiedParams)
    % 生成预定义分枝参数
    % 包含：严格的数据依赖 + 激活完整性验证 + 几何参数校验
    
    % 基础检查
    if isempty(preConfig) || isempty(identifiedParams)
        error('ConfigAdapter:MissingData', '缺少配置或识别参数');
    end
    
    % [修复] 移除此处对 taper_factors 的全局提取
    % 原因：taper 是在 estimateStiffnessDamping 中针对每个分枝单独计算的
    % identified_taper = identifiedParams.linear.taper_factors; <--- DELETE
    
    predefined = struct();
    fruitConfig = preConfig.fruit;
    
    % --- 定义内部辅助函数：处理单个分枝的通用逻辑 ---
    function processBranch(name, type_struct)
        % 1. 确定等级
        lvl = determineBranchLevel(name);
        
        % 2. 获取几何
        if ~isfield(type_struct, name)
            error('缺少分枝 %s 的几何参数', name);
        end
        geom = type_struct.(name);
        validateBranchGeometry(geom, name);
        
        % 3. 获取基准刚度 AND 专属递减因子
        % 此时 estimateStiffnessDamping 返回三个值
        [k_b, c_b, specific_taper] = estimateStiffnessDamping(geom, identifiedParams, name);
        
        % 4. 使用专属递减因子生成分段参数
        % 这样 Root, Mid, Tip 都会精确等于实验识别值
        predefined.(name) = generateBranchSegmentParams(geom, k_b, c_b, specific_taper);
        
        predefined.(name).branch_level = lvl;
        
        % 5. 挂果逻辑
        if shouldAttachFruit(name, lvl, fruitConfig)
            if shouldAttachAtPosition(lvl, 'mid', fruitConfig)
                predefined.(name).fruit_at_mid = buildFruitParamsStrict(preConfig, identifiedParams);
            end
            if shouldAttachAtPosition(lvl, 'tip', fruitConfig)
                predefined.(name).fruit_at_tip = buildFruitParamsStrict(preConfig, identifiedParams);
            end
        end
    end

    % --- 遍历生成 ---
    % 一级分枝
    fprintf('    正在处理一级分枝...\n');
    for p = 1:preConfig.topology.num_primary_branches
        processBranch(sprintf('P%d', p), preConfig.primary);
    end
    
    % 二级分枝
    fprintf('    正在处理二级分枝...\n');
    for p = 1:preConfig.topology.num_primary_branches
        for s = 1:preConfig.topology.secondary_branches_count(p)
            processBranch(sprintf('P%d_S%d', p, s), preConfig.secondary);
        end
    end
    
    % 三级分枝
    fprintf('    正在处理三级分枝...\n');
    for p = 1:preConfig.topology.num_primary_branches
        if p <= length(preConfig.topology.tertiary_branches_count)
            tertiary_counts = preConfig.topology.tertiary_branches_count{p};
            for s = 1:length(tertiary_counts)
                for t = 1:tertiary_counts(s)
                    processBranch(sprintf('P%d_S%d_T%d', p, s, t), preConfig.tertiary);
                end
            end
        end
    end

    % === 完整性验证 ===
    validatePredefinedParams(predefined, preConfig);
    
    fprintf('  分枝参数生成完成 (已通过几何校验与完整性验证)。\n');
end

%% ==================== 生成单个分枝的段参数 ====================
function params = generateBranchSegmentParams(branchGeom, k_base, c_base, identified_taper)
    % 生成单个分枝的段参数 - 严格模式
    % 递减因子必须从实验数据获取，无数据则报错
    % 输入:
    %   branchGeom - 分枝几何结构体（来自GUI配置）
    %   k_base - 基础刚度（从实验识别）
    %   c_base - 基础阻尼（从实验识别）
    %   identified_taper - 从实验识别的递减因子结构体（必须提供）
    
   % 验证递减因子必须存在且有效
    if nargin < 4 || isempty(identified_taper)
        error('ConfigAdapter:MissingTaperFactors', ...
              ['严重错误：生成分枝参数时缺少 "identified_taper" (刚度/阻尼递减因子)。\n' ...
               'ConfigAdapter 无法推断分枝沿长度方向的刚度变化。\n' ...
               '请检查 analyse_chibi_data.m 是否成功识别并输出了 taper_factors。']);
    end
    
    if ~isfield(identified_taper, 'k') || ~isfield(identified_taper, 'c') || ...
       length(identified_taper.k) ~= 3 || length(identified_taper.c) ~= 3
        error('ConfigAdapter:InvalidTaperFactors', ...
              ['严重错误：提供的 "identified_taper" 格式不正确。\n' ...
               '必须包含 "k" 和 "c" 字段，且每个字段必须是长度为 3 的向量 (对应 Root/Mid/Tip)。']);
    end
    
    if ~isstruct(identified_taper)
        error('ConfigAdapter:InvalidData', ...
              'identified_taper必须是结构体');
    end
    
    if ~isfield(identified_taper, 'k') || ~isfield(identified_taper, 'c')
        error('ConfigAdapter:MissingData', ...
              'identified_taper必须包含k和c字段（从实验FRF识别的刚度阻尼递减因子）');
    end
    
    if length(identified_taper.k) ~= 3 || length(identified_taper.c) ~= 3
        error('ConfigAdapter:InvalidData', ...
              '递减因子必须是长度为3的向量 [root, mid, tip]');
    end
    
    k_taper = identified_taper.k;
    c_taper = identified_taper.c;
    
    params = struct();
    
    % Z方向因子
    if ~isfield(branchGeom, 'z_factor')
        error('ConfigAdapter:MissingData', ...
              '分枝几何缺少z_factor，请在GUI中配置');
    end
    z_factor = branchGeom.z_factor;
    
    % 质量分配
    if ~isfield(branchGeom, 'mass_dist')
        error('ConfigAdapter:MissingData', ...
              '分枝几何缺少mass_dist（质量分配），请在GUI中配置');
    end
    m_dist = branchGeom.mass_dist;
    
    if ~isfield(branchGeom, 'total_mass')
        error('ConfigAdapter:MissingData', ...
              '分枝几何缺少total_mass，请在GUI中配置');
    end
    m_total = branchGeom.total_mass;
    
    % Root段
    params.root.m = m_total * m_dist(1);
    params.root.k_y_conn = k_base * k_taper(1);
    params.root.c_y_conn = c_base * c_taper(1);
    params.root.k_z_conn = k_base * k_taper(1) * z_factor;
    params.root.c_z_conn = c_base * c_taper(1) * z_factor;
    
    % Mid段
    params.mid.m = m_total * m_dist(2);
    params.mid.k_y_conn = k_base * k_taper(2);
    params.mid.c_y_conn = c_base * c_taper(2);
    params.mid.k_z_conn = k_base * k_taper(2) * z_factor;
    params.mid.c_z_conn = c_base * c_taper(2) * z_factor;
    
    % Tip段
    params.tip.m = m_total * m_dist(3);
    params.tip.k_y_conn = k_base * k_taper(3);
    params.tip.c_y_conn = c_base * c_taper(3);
    params.tip.k_z_conn = k_base * k_taper(3) * z_factor;
    params.tip.c_z_conn = c_base * c_taper(3) * z_factor;
    
    % 保存几何信息（必须有）
    if ~isfield(branchGeom, 'length')
        error('ConfigAdapter:MissingData', '分枝几何缺少length');
    end
    if ~isfield(branchGeom, 'diameter_base')
        error('ConfigAdapter:MissingData', '分枝几何缺少diameter_base');
    end
    if ~isfield(branchGeom, 'diameter_tip')
        error('ConfigAdapter:MissingData', '分枝几何缺少diameter_tip');
    end
    
    params.geometry.length = branchGeom.length;
    params.geometry.diameter_base = branchGeom.diameter_base;
    params.geometry.diameter_tip = branchGeom.diameter_tip;
    params.taper_factors.k = k_taper;
    params.taper_factors.c = c_taper;
end


%% ==================== 辅助函数 ====================
function [k_base, c_base, branch_taper] = estimateStiffnessDamping(branchGeom, identifiedParams, branchName)
    % 获取分枝刚度和阻尼 - 严格独立数据驱动版
    
    target_linear_params = [];
    
    % 1. 【严格匹配】必须在 branches 结构体中找到对应名称的数据
    if isfield(identifiedParams, 'branches') && isfield(identifiedParams.branches, branchName)
        branch_data = identifiedParams.branches.(branchName);
        if isfield(branch_data, 'linear')
            target_linear_params = branch_data.linear;
        end
    end
    
    % 2. 【严格校验】无数据直接报错，严禁回退到全局平均值
    if isempty(target_linear_params)
        error('ConfigAdapter:MissingExperimentData', ...
              ['严重错误：分枝 "%s" 没有任何实验识别数据！\n' ...
               '违反了"严格数据驱动"原则。请在参数识别阶段加载该分枝的实验文件。\n' ...
               '不允许使用全局平均值或父级参数代替。'], ...
               branchName);
    end
    
    % 3. 数据完整性检查
    if ~isfield(target_linear_params, 'K') || ~isfield(target_linear_params, 'C')
        error('ConfigAdapter:InvalidData', '分枝 "%s" 的识别结果损坏，缺少 K 或 C 矩阵。', branchName);
    end
    
    % 4. 提取矩阵
    K_vals = diag(target_linear_params.K); 
    C_vals = diag(target_linear_params.C);
    
    % 5. 计算基础值 (Base) 和 递减因子 (Taper)
    k_base = max(K_vals);
    if k_base <= 0
        error('ConfigAdapter:BadData', '分枝 "%s" 的刚度识别值无效(全部<=0)。', branchName);
    end
    k_taper = K_vals / k_base; 
    
    c_base = max(C_vals);
    if c_base <= 0
        error('ConfigAdapter:WarnData', '分枝 "%s" 的阻尼识别值异常，将导致仿真不稳定。', branchName);
        % 即使阻尼异常也不建议给默认值，最好报错，但此处保留警告让用户决定
    end
    c_taper = C_vals / c_base; % 防止除零
    
    branch_taper = struct();
    branch_taper.k = k_taper;
    branch_taper.c = c_taper;
    
    fprintf('    [√] 分枝 %s 参数加载成功 (k_base=%.1f, c_base=%.3f)\n', branchName, k_base, c_base);
end

function level = determineBranchLevel(branchName)
    % 根据分枝名称确定层级
    % P1, P2, P3 -> 1 (一级)
    % P1_S1, P2_S1 -> 2 (二级)
    % P1_S1_T1 -> 3 (三级)
    
    underscores = strfind(branchName, '_');
    if isempty(underscores)
        level = 1;  % 一级分枝 (P1, P2, ...)
    elseif length(underscores) == 1
        level = 2;  % 二级分枝 (P1_S1, ...)
    else
        level = 3;  % 三级分枝 (P1_S1_T1, ...)
    end
end

function attach = shouldAttachFruit(branchName, level, fruitConfig)
    % 根据配置确定是否需要挂果
    % 一级分枝不直接挂果
    % 二级和三级根据配置决定
    
    if level == 1
        attach = false;
    elseif level == 2
        attach = fruitConfig.attach_secondary_mid || fruitConfig.attach_secondary_tip;
    else  % level == 3
        attach = fruitConfig.attach_tertiary_mid || fruitConfig.attach_tertiary_tip;
    end
end

function attach = shouldAttachAtPosition(level, position, fruitConfig)
    % 确定特定位置是否挂果
    
    if level == 2
        if strcmp(position, 'mid')
            attach = fruitConfig.attach_secondary_mid;
        else
            attach = fruitConfig.attach_secondary_tip;
        end
    elseif level == 3
        if strcmp(position, 'mid')
            attach = fruitConfig.attach_tertiary_mid;
        else
            attach = fruitConfig.attach_tertiary_tip;
        end
    else
        attach = false;
    end
end



function validateBranchGeometry(branchGeom, branchName)
    % 验证分枝几何参数完整性 - 严格模式
    
    required_fields = {'total_mass', 'length', 'diameter_base', 'diameter_tip', 'mass_dist'};
    
    for i = 1:length(required_fields)
        field = required_fields{i};
        if ~isfield(branchGeom, field)
            error('ConfigAdapter:MissingData', ...
                  '分枝 %s 几何参数缺少字段: %s，请在GUI中完整配置', branchName, field);
        end
        
        value = branchGeom.(field);
        
        if strcmp(field, 'mass_dist')
            if length(value) ~= 3
                error('ConfigAdapter:InvalidData', ...
                      '分枝 %s 的mass_dist必须是长度为3的向量 [root, mid, tip]', branchName);
            end
            if abs(sum(value) - 1) > 0.01
                error('ConfigAdapter:InvalidData', ...
                      '分枝 %s 的mass_dist之和必须为1，当前为%.4f', branchName, sum(value));
            end
        else
            if isempty(value) || (isnumeric(value) && value <= 0)
                error('ConfigAdapter:InvalidData', ...
                      '分枝 %s 的 %s 值无效(空或非正)', branchName, field);
            end
        end
    end
end

function fruit_params = buildFruitParamsStrict(preConfig, identifiedParams)
    % 功能：基于识别出的统计模型生成果实参数
    % 核心原则：严禁使用预设的固定值 (如 F=5N)，必须由模型预测得出
    
    fruit_params = struct();
    
    % 1. 物理属性 (来自 GUI 预配置，这些是几何输入)
    fruit_params.m = preConfig.fruit.mass;
    fruit_params.diameter = preConfig.fruit.diameter;
    % 检查 preConfig 是否包含新字段（兼容旧版配置）
    if isfield(preConfig.fruit, 'k_pedicel')
        k_val = preConfig.fruit.k_pedicel;
        c_val = preConfig.fruit.c_pedicel;
    else
        warning('ConfigAdapter:LegacyConfig', '预配置中缺少果柄动力学参数，使用默认值。请更新 GUI 配置。');
        k_val = 2000;
        c_val = 0.5;
    end
    % 应用到 Y/Z 两个方向 (假设各向同性，或者也可以在GUI分开配置)
    fruit_params.k_pedicel_y = k_val;
    fruit_params.c_pedicel_y = c_val;
    fruit_params.k_pedicel_z = k_val;
    fruit_params.c_pedicel_z = c_val;

    % 2. 核心：计算断裂力 F_break
    % 需要从 identifiedParams 中提取脱落模型
    % 假设所有分枝共享同一个脱落机制模型，我们取第一个有效分枝的模型即可
    det_model = [];
    % [修复] 适配聚合结构体
    if isfield(identifiedParams, 'branches')
        % 遍历寻找含有 detachment_model 的分枝
        fn = fieldnames(identifiedParams.branches);
        for i = 1:length(fn)
            if isfield(identifiedParams.branches.(fn{i}), 'detachment_model')
                det_model = identifiedParams.branches.(fn{i}).detachment_model;
                break;
            end
        end
    elseif isfield(identifiedParams, 'detachment_model')
        det_model = identifiedParams.detachment_model;
    end
    
    if ~isempty(det_model)
        % 提取 GUI 中配置的平均特征作为输入
        % 注意：这里将直径从 m 转换为 cm 以匹配标定数据的单位
        D_input_cm = (preConfig.fruit.diameter * 100); 
        
        % 假设场景：中等高度，末端挂果，无开裂
        H_input_m = 1.5; 
        P_input_idx = 1.0; % 末端
        S_input_crack = 0; % 无开裂
        
        % 调用模型的预测接口
        predicted_F = det_model.predict(H_input_m, P_input_idx, D_input_cm, S_input_crack);
        
        % 加上随机波动 (模拟果实间的差异，基于识别出的误差 sigma)
        % 如果希望确定性仿真，去掉 randn 部分
        final_F_break = predicted_F + det_model.sigma_epsilon * randn();
        
        fruit_params.F_break = max(1.0, final_F_break); % 确保力为正值
        
        fprintf('    [ConfigAdapter] 果实断裂力 F_break 已由数据模型计算: %.2f N (模型预测值)\n', fruit_params.F_break);
    else
        % 如果没有识别结果，为了不让程序崩溃，报错提示
        error('ConfigAdapter:NoDetachmentModel', ...
              '未找到果实脱落力模型。请确保在参数识别阶段成功运行了 Stage 4 (基于内置的20组数据)。');
    end
end
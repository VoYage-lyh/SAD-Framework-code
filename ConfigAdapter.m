function [analysis_params, sim_params] = ConfigAdapter(preConfig, identifiedParams)
% ConfigAdapter_v4 - 参数适配器（终极增强版）
%
% 修复记录：
%   1. [关键] 自动补全 topology 中的统计字段 (num_primary_branches 等)，防止仿真构建报错。
%   2. [关键] 自动将 trunk 几何信息注入 config，防止绘图函数报错。
%   3. 保留了之前修复的 Z 方向数据提取和并行参数名问题。

    if nargin < 2
        identifiedParams = [];
    end
    
    %% 1. 生成参数识别所需参数
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

    %% 2. 检查是否需要生成仿真参数
    if isempty(identifiedParams)
        sim_params = [];
        fprintf('ConfigAdapter: 仅返回分析参数，未生成仿真参数 (identifiedParams 为空)。\n');
        return; 
    end
    
    %% 3. 生成仿真所需参数
    sim_params = struct();
    
    % --- 基础设置 ---
    sim_params.parallel_max_workers = preConfig.basic.parallel_max_workers; 
    sim_params.workFolder = preConfig.basic.workFolder;
    
    if isfield(preConfig.basic, 'modelName')
        sim_params.model_name = preConfig.basic.modelName;
    else
        sim_params.model_name = 'MDOF_Hierarchical_Vibration_Sim'; 
    end
    
    sim_params.gravity_g = preConfig.basic.gravity_g;
    sim_params.use_parallel = preConfig.basic.useParallel;
    
    % =========================================================================
    % [核心修复] 拓扑结构自动补全与增强
    % =========================================================================
    sim_params.config = preConfig.topology;
    
    % 1. 补全 num_primary_branches 及计数数组
    if isfield(sim_params.config, 'structure')
        structure = sim_params.config.structure;
        num_p = length(structure);
        
        % 强制覆盖/添加统计字段
        sim_params.config.num_primary_branches = num_p;
        sim_params.config.secondary_branches_count = zeros(1, num_p);
        sim_params.config.tertiary_branches_count = cell(1, num_p);
        
        for i = 1:num_p
            vec = structure{i};
            % 处理空分枝标记 -1
            if isequal(vec, -1) || (length(vec)==1 && vec(1) == -1)
                sim_params.config.secondary_branches_count(i) = 0;
                sim_params.config.tertiary_branches_count{i} = [];
            else
                sim_params.config.secondary_branches_count(i) = length(vec);
                sim_params.config.tertiary_branches_count{i} = vec;
            end
        end
        fprintf('    [ConfigAdapter] 已自动补全拓扑统计数据 (一级分枝数: %d)。\n', num_p);
    end
    
    % 2. 补全 trunk 信息 (防止绘图函数 plot_tree_topology 报错)
    if isfield(preConfig, 'trunk')
        sim_params.config.trunk = preConfig.trunk;
    end
    % =========================================================================
    
    % --- 构建主干参数 ---
    sim_params.trunk = buildTrunkParams(preConfig, identifiedParams);
    
    % --- 构建分枝参数 ---
    sim_params.predefined_params = generatePredefinedParams(preConfig, identifiedParams);
    
    % --- 构建果实配置 ---
    sim_params.fruit_config = struct();
    sim_params.fruit_config.attach_secondary_mid = preConfig.fruit.attach_secondary_mid;
    sim_params.fruit_config.attach_secondary_tip = preConfig.fruit.attach_secondary_tip;
    sim_params.fruit_config.attach_tertiary_mid = preConfig.fruit.attach_tertiary_mid;
    sim_params.fruit_config.attach_tertiary_tip = preConfig.fruit.attach_tertiary_tip;
    
    if isfield(preConfig.fruit, 'fruits_per_node')
        sim_params.fruit_config.fruits_per_node = preConfig.fruit.fruits_per_node;
    else
        sim_params.fruit_config.fruits_per_node = 1;
    end
    
    % 生成默认果实参数
    sim_params.default_fruit_params = buildFruitParamsStrict(preConfig, identifiedParams);
    
    % 激励参数
    sim_params.excitation = preConfig.excitation;
    
    % 更新激励频率
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
    
    fprintf('ConfigAdapter: 参数适配完成。\n');
end

%% ==================== 构建主干参数 ====================
function trunk = buildTrunkParams(preConfig, identifiedParams)
    trunk = struct();
    
    % --- 1. 几何与质量 ---
    trunk.length = preConfig.trunk.length;
    trunk.diameter_base = preConfig.trunk.diameter_base;
    trunk.diameter_tip = preConfig.trunk.diameter_tip;
    
    m_total = preConfig.trunk.total_mass;
    m_dist = preConfig.trunk.mass_distribution;
    trunk.root.m = m_total * m_dist(1);
    trunk.mid.m = m_total * m_dist(2);
    trunk.tip.m = m_total * m_dist(3);
    
    % --- 2. 刚度和阻尼 (线性 + 非线性) ---
    target_lin = [];
    target_nonlin = [];
    
    % 尝试从 branches.Trunk 结构获取
    if isfield(identifiedParams, 'branches') && isfield(identifiedParams.branches, 'Trunk')
        branch_data = identifiedParams.branches.Trunk;
        if isfield(branch_data, 'linear'), target_lin = branch_data.linear; end
        if isfield(branch_data, 'nonlinear'), target_nonlin = branch_data.nonlinear; end
        fprintf('    [ConfigAdapter] 主干：加载 Trunk 分枝识别参数。\n');
    elseif isfield(identifiedParams, 'linear')
        target_lin = identifiedParams.linear;
        if isfield(identifiedParams, 'nonlinear'), target_nonlin = identifiedParams.nonlinear; end
        fprintf('    [ConfigAdapter] 主干：加载全局参数。\n');
    else
        error('ConfigAdapter:NoTrunkParams', ...
              '错误：未在 identifiedParams 中找到主干(Trunk)或全局(linear)参数。');
    end
    
    % 线性赋值 (Y)
    if ~isfield(target_lin, 'identified_params_x') || isempty(target_lin.identified_params_x)
        error('ConfigAdapter:MissingYData', '主干参数缺少 Y 方向数据 (identified_params_x)。');
    end
    p_vec_y = target_lin.identified_params_x;
    trunk.root.k_y_conn_to_base = p_vec_y(1); trunk.root.c_y_conn_to_base = p_vec_y(2);
    trunk.root.k_y_conn = p_vec_y(3);         trunk.root.c_y_conn = p_vec_y(4);
    trunk.mid.k_y_conn  = p_vec_y(5);         trunk.mid.c_y_conn  = p_vec_y(6);
    trunk.tip.k_y_conn  = 0;                  trunk.tip.c_y_conn  = 0;

    % 线性赋值 (Z)
    p_vec_z = [];
    if isfield(target_lin, 'identified_params_z')
        p_vec_z = target_lin.identified_params_z;
    elseif isfield(target_lin, 'K_z') 
        error('ConfigAdapter:TrunkFormat', '主干 Z 数据格式不兼容，期望 identified_params_z 向量。');
    end
    
    if isempty(p_vec_z)
        error('ConfigAdapter:NoZParams', '严重错误：主干缺少 Z 方向识别参数 (identified_params_z)。');
    end
    
    trunk.root.k_z_conn_to_base = p_vec_z(1); trunk.root.c_z_conn_to_base = p_vec_z(2);
    trunk.root.k_z_conn = p_vec_z(3);         trunk.root.c_z_conn = p_vec_z(4);
    trunk.mid.k_z_conn  = p_vec_z(5);         trunk.mid.c_z_conn  = p_vec_z(6);
    trunk.tip.k_conn_z  = 0;                  trunk.tip.c_z_conn  = 0;

    % 非线性赋值
    trunk.root.k3_y_conn_to_base = 0; trunk.root.c2_y_conn_to_base = 0;
    trunk.root.k3_y_conn = 0;         trunk.root.c2_y_conn = 0;
    trunk.mid.k3_y_conn  = 0;         trunk.mid.c2_y_conn  = 0;
    trunk.root.k3_z_conn_to_base = 0; trunk.root.c2_z_conn_to_base = 0; 
    trunk.root.k3_z_conn = 0;         trunk.root.c2_z_conn = 0;
    trunk.mid.k3_z_conn = 0;          trunk.mid.c2_z_conn = 0;

    if ~isempty(target_nonlin) && isfield(target_nonlin, 'k3_coeffs')
        k3 = target_nonlin.k3_coeffs; c2 = target_nonlin.c2_coeffs;
        if length(k3)>=1, trunk.root.k3_y_conn_to_base = k3(1); trunk.root.c2_y_conn_to_base = c2(1); end
        if length(k3)>=2, trunk.root.k3_y_conn = k3(2);         trunk.root.c2_y_conn = c2(2); end
        if length(k3)>=3, trunk.mid.k3_y_conn  = k3(3);         trunk.mid.c2_y_conn  = c2(3); end
        fprintf('    [ConfigAdapter] 主干 Y 非线性参数已应用。\n');
    end
end

%% ==================== 验证预定义参数完整性 ====================
function validatePredefinedParams(predefined, preConfig)
    missingBranches = {};
    incompleteParams = {};
    structure = preConfig.topology.structure;
    
    numP = length(structure);
    for p = 1:numP
        bid_p = sprintf('P%d', p);
        if ~isfield(predefined, bid_p), missingBranches{end+1} = bid_p;
        elseif ~checkBranchParamIntegrity(predefined.(bid_p)), incompleteParams{end+1} = bid_p; end
        
        vec = structure{p};
        if isequal(vec, -1) || (length(vec)==1 && vec(1) == -1), continue; end
        
        numS = length(vec);
        for s = 1:numS
            bid_s = sprintf('P%d_S%d', p, s);
            if ~isfield(predefined, bid_s), missingBranches{end+1} = bid_s;
            elseif ~checkBranchParamIntegrity(predefined.(bid_s)), incompleteParams{end+1} = bid_s; end
            
            numT = vec(s);
            for t = 1:numT
                bid_t = sprintf('P%d_S%d_T%d', p, s, t);
                if ~isfield(predefined, bid_t), missingBranches{end+1} = bid_t;
                elseif ~checkBranchParamIntegrity(predefined.(bid_t)), incompleteParams{end+1} = bid_t; end
            end
        end
    end
    
    if ~isempty(missingBranches)
        error('ConfigAdapter:MissingData', '以下分枝在预配置中缺少几何参数:\n  %s', strjoin(missingBranches, ', '));
    end
    if ~isempty(incompleteParams)
        warning('ConfigAdapter:IncompleteParams', '以下分枝参数不完整:\n  %s', strjoin(incompleteParams, ', '));
    end
end

function is_ok = checkBranchParamIntegrity(branch_data)
    is_ok = true;
    segments = {'root', 'mid', 'tip'};
    required_fields = {'m', 'k_y_conn', 'c_y_conn', 'k_z_conn', 'c_z_conn'};
    for i = 1:length(segments)
        seg = segments{i};
        if isfield(branch_data, seg)
            seg_data = branch_data.(seg);
            for k = 1:length(required_fields)
                if ~isfield(seg_data, required_fields{k})
                    is_ok = false; return;
                end
            end
        else
            is_ok = false; return;
        end
    end
end

%% ==================== 为分枝生成默认段参数 ====================
function predefined = generatePredefinedParams(preConfig, identifiedParams)
    
    if isempty(preConfig) || isempty(identifiedParams)
        error('ConfigAdapter:MissingData', '缺少配置或识别参数');
    end
    
    structure = preConfig.topology.structure;
    predefined = struct();
    fruitConfig = preConfig.fruit;
    
    function processBranch(name, type_struct)
        lvl = determineBranchLevel(name);
        if ~isfield(type_struct, name), error('缺少分枝 %s 的几何参数', name); end
        geom = type_struct.(name);
        validateBranchGeometry(geom, name);
        
        [k_b, c_b, k3_b, c2_b, specific_taper, z_params] = estimateStiffnessDamping(geom, identifiedParams, name);
        
        predefined.(name) = generateBranchSegmentParams(geom, k_b, c_b, k3_b, c2_b, specific_taper, z_params);
        
        predefined.(name).branch_level = lvl;
        if shouldAttachFruit(name, lvl, fruitConfig)
            if shouldAttachAtPosition(lvl, 'mid', fruitConfig)
                predefined.(name).fruit_at_mid = buildFruitParamsStrict(preConfig, identifiedParams);
            end
            if shouldAttachAtPosition(lvl, 'tip', fruitConfig)
                predefined.(name).fruit_at_tip = buildFruitParamsStrict(preConfig, identifiedParams);
            end
        end
    end

    numP = length(structure);
    fprintf('    [ConfigAdapter] 正在基于 Cell 拓扑生成参数 (智能双向模式)...\n');
    
    for p = 1:numP
        branch_name_p = sprintf('P%d', p);
        processBranch(branch_name_p, preConfig.primary);
        
        vec = structure{p};
        if isequal(vec, -1) || (length(vec)==1 && vec(1) == -1), continue; end
        
        numS = length(vec);
        for s = 1:numS
            branch_name_s = sprintf('P%d_S%d', p, s);
            processBranch(branch_name_s, preConfig.secondary);
            
            numT = vec(s);
            for t = 1:numT
                branch_name_t = sprintf('P%d_S%d_T%d', p, s, t);
                processBranch(branch_name_t, preConfig.tertiary);
            end
        end
    end

    validatePredefinedParams(predefined, preConfig);
    fprintf('  分枝参数库生成完成。\n');
end

%% ==================== 生成单个分枝的段参数 ====================
function params = generateBranchSegmentParams(branchGeom, k_base, c_base, k3_base, c2_base, identified_taper, z_params)
    
    if isempty(identified_taper), error('ConfigAdapter:MissingTaper', '缺少 Y 方向分布因子 identified_taper。'); end
    if isempty(z_params)
        error('ConfigAdapter:MissingZData', '严重错误：Z 方向实测数据为空。参数提取步骤可能失败。');
    end
    
    k_taper = identified_taper.k; c_taper = identified_taper.c;
    k3_taper = identified_taper.k3; c2_taper = identified_taper.c2;
    
    params = struct();
    m_total = branchGeom.total_mass;
    m_dist = branchGeom.mass_dist;
    
    % Root段
    params.root.m = m_total * m_dist(1);
    params.root.k_y_conn = k_base * k_taper(1);
    params.root.c_y_conn = c_base * c_taper(1);
    params.root.k3_y_conn = k3_base * k3_taper(1);
    params.root.c2_y_conn = c2_base * c2_taper(1);
    params.root.k_z_conn = z_params.k_base * z_params.k_taper(1);
    params.root.c_z_conn = z_params.c_base * z_params.c_taper(1);
    params.root.k3_z_conn = z_params.k3_base * z_params.k3_taper(1);
    params.root.c2_z_conn = z_params.c2_base * z_params.c2_taper(1);
    
    % Mid段
    params.mid.m = m_total * m_dist(2);
    params.mid.k_y_conn = k_base * k_taper(2);
    params.mid.c_y_conn = c_base * c_taper(2);
    params.mid.k3_y_conn = k3_base * k3_taper(2);
    params.mid.c2_y_conn = c2_base * c2_taper(2);
    params.mid.k_z_conn = z_params.k_base * z_params.k_taper(2);
    params.mid.c_z_conn = z_params.c_base * z_params.c_taper(2);
    params.mid.k3_z_conn = z_params.k3_base * z_params.k3_taper(2);
    params.mid.c2_z_conn = z_params.c2_base * z_params.c2_taper(2);
    
    % Tip段
    params.tip.m = m_total * m_dist(3);
    params.tip.k_y_conn = k_base * k_taper(3);
    params.tip.c_y_conn = c_base * c_taper(3);
    params.tip.k3_y_conn = k3_base * k3_taper(3);
    params.tip.c2_y_conn = c2_base * c2_taper(3);
    params.tip.k_z_conn = z_params.k_base * z_params.k_taper(3);
    params.tip.c_z_conn = z_params.c_base * z_params.c_taper(3);
    params.tip.k3_z_conn = z_params.k3_base * z_params.k3_taper(3);
    params.tip.c2_z_conn = z_params.c2_base * z_params.c2_taper(3);
    
    params.geometry.length = branchGeom.length;
    params.geometry.diameter_base = branchGeom.diameter_base;
    params.geometry.diameter_tip = branchGeom.diameter_tip;
    params.taper_factors = identified_taper;
end

%% ==================== 辅助函数: 提取基准值和分布 ====================
function [k_base, c_base, k3_base, c2_base, branch_taper, z_params] = estimateStiffnessDamping(branchGeom, identifiedParams, branchName)
    
    t_lin_y = []; t_nonlin_y = [];
    t_lin_z = []; t_nonlin_z = [];
    
    % 1. 尝试从分枝专用数据读取
    if isfield(identifiedParams, 'branches') && isfield(identifiedParams.branches, branchName)
        branch_data = identifiedParams.branches.(branchName);
        
        if isfield(branch_data, 'linear'), t_lin_y = branch_data.linear; end
        
        if isfield(branch_data, 'linear_z')
            t_lin_z = branch_data.linear_z;
        elseif isfield(branch_data, 'linear')
             if isfield(branch_data.linear, 'K_z')
                 t_lin_z = struct('K', branch_data.linear.K_z, 'C', branch_data.linear.C_z);
             elseif isfield(branch_data.linear, 'identified_params_z')
                 vec = branch_data.linear.identified_params_z;
                 k_d = [vec(1)+vec(3), vec(3)+vec(5), vec(5)];
                 c_d = [vec(2)+vec(4), vec(4)+vec(6), vec(6)];
                 t_lin_z = struct('K', diag(k_d), 'C', diag(c_d));
             end
        end
        if isempty(t_lin_z) && isfield(branch_data, 'identified_params_z')
             vec = branch_data.identified_params_z;
             k_d = [vec(1)+vec(3), vec(3)+vec(5), vec(5)];
             c_d = [vec(2)+vec(4), vec(4)+vec(6), vec(6)];
             t_lin_z = struct('K', diag(k_d), 'C', diag(c_d));
        end
        
        if isfield(branch_data, 'nonlinear'), t_nonlin_y = branch_data.nonlinear; end
        if isfield(branch_data, 'nonlinear_z'), t_nonlin_z = branch_data.nonlinear_z; end
    end
    
    % 2. 尝试从全局数据读取 (兜底)
    g_lin_y = []; g_lin_z = [];
    g_nonlin_y = []; g_nonlin_z = [];
    
    if isfield(identifiedParams, 'linear')
        lin = identifiedParams.linear;
        if isfield(lin, 'K_x')
            g_lin_y = struct('K', lin.K_x, 'C', lin.C_x);
        elseif isfield(lin, 'identified_params_x')
             vec = lin.identified_params_x;
             k_d = [vec(1)+vec(3), vec(3)+vec(5), vec(5)];
             c_d = [vec(2)+vec(4), vec(4)+vec(6), vec(6)];
             g_lin_y = struct('K', diag(k_d), 'C', diag(c_d));
        else
            g_lin_y = lin;
        end
        if isfield(lin, 'K_z')
            g_lin_z = struct('K', lin.K_z, 'C', lin.C_z);
        elseif isfield(lin, 'identified_params_z')
             vec = lin.identified_params_z;
             k_d = [vec(1)+vec(3), vec(3)+vec(5), vec(5)];
             c_d = [vec(2)+vec(4), vec(4)+vec(6), vec(6)];
             g_lin_z = struct('K', diag(k_d), 'C', diag(c_d));
        end
        
        if isfield(identifiedParams, 'nonlinear'), g_nonlin_y = identifiedParams.nonlinear; end
        if isfield(identifiedParams, 'nonlinear_z'), g_nonlin_z = identifiedParams.nonlinear_z; end
    end
    
    % 3. 参数融合
    target_lin = t_lin_y; if isempty(target_lin), target_lin = g_lin_y; end
    target_lin_z = t_lin_z; if isempty(target_lin_z), target_lin_z = g_lin_z; end
    target_nonlin = t_nonlin_y; if isempty(target_nonlin), target_nonlin = g_nonlin_y; end
    target_nonlin_z = t_nonlin_z; if isempty(target_nonlin_z), target_nonlin_z = g_nonlin_z; end
    
    % 4. 提取 Y 参数
    if isempty(target_lin) || ~isfield(target_lin, 'K') || ~isfield(target_lin, 'C')
        error('ConfigAdapter:MissingYData', '分枝 "%s" 缺少 Y 方向线性数据 (identified_params_x)。', branchName);
    end
    K_vals = diag(target_lin.K); C_vals = diag(target_lin.C);
    k_base = max(K_vals); k_taper = K_vals / k_base; 
    c_base = max(C_vals); c_taper = C_vals / c_base;
    [k3_base, c2_base, k3_taper, c2_taper] = extractNonlinearBaseTaper(target_nonlin);
    
    % 5. 提取 Z 参数
    if isempty(target_lin_z) || ~isfield(target_lin_z, 'K')
        error('ConfigAdapter:MissingZData', '分枝 "%s" 缺少 Z 方向线性数据！', branchName);
    end
    z_params = struct();
    K_vals_z = diag(target_lin_z.K); C_vals_z = diag(target_lin_z.C);
    z_params.k_base = max(K_vals_z); z_params.c_base = max(C_vals_z);
    z_params.k_taper = K_vals_z / z_params.k_base; z_params.c_taper = C_vals_z / z_params.c_base;
    [z_params.k3_base, z_params.c2_base, z_params.k3_taper, z_params.c2_taper] = extractNonlinearBaseTaper(target_nonlin_z);
    
    branch_taper = struct('k', k_taper, 'c', c_taper, 'k3', k3_taper, 'c2', c2_taper);
end

function [k3_base, c2_base, k3_taper, c2_taper] = extractNonlinearBaseTaper(nl_struct)
    if ~isempty(nl_struct) && isfield(nl_struct, 'k3_coeffs')
        k3_vals = nl_struct.k3_coeffs(:); c2_vals = nl_struct.c2_coeffs(:);
        if length(k3_vals) < 3, k3_vals(end+1:3) = 0; end
        if length(c2_vals) < 3, c2_vals(end+1:3) = 0; end
        [max_k3, idx_k3] = max(abs(k3_vals(1:3)));
        if max_k3 > 1e-12, k3_base = k3_vals(idx_k3); k3_taper = k3_vals(1:3)/k3_base;
        else, k3_base=0; k3_taper=[0;0;0]; end
        [max_c2, idx_c2] = max(abs(c2_vals(1:3)));
        if max_c2 > 1e-12, c2_base = c2_vals(idx_c2); c2_taper = c2_vals(1:3)/c2_base;
        else, c2_base=0; c2_taper=[0;0;0]; end
    else
        k3_base=0; c2_base=0; k3_taper=[0;0;0]; c2_taper=[0;0;0];
    end
end

function level = determineBranchLevel(branchName)
    underscores = strfind(branchName, '_');
    if isempty(underscores), level = 1;
    elseif length(underscores) == 1, level = 2;
    else, level = 3; end
end

function attach = shouldAttachFruit(branchName, level, fruitConfig)
    if level == 1, attach = false;
    elseif level == 2, attach = fruitConfig.attach_secondary_mid || fruitConfig.attach_secondary_tip;
    else, attach = fruitConfig.attach_tertiary_mid || fruitConfig.attach_tertiary_tip; end
end

function attach = shouldAttachAtPosition(level, position, fruitConfig)
    if level == 2
        if strcmp(position, 'mid'), attach = fruitConfig.attach_secondary_mid; else, attach = fruitConfig.attach_secondary_tip; end
    elseif level == 3
        if strcmp(position, 'mid'), attach = fruitConfig.attach_tertiary_mid; else, attach = fruitConfig.attach_tertiary_tip; end
    else, attach = false; end
end

function validateBranchGeometry(branchGeom, branchName)
    if ~isfield(branchGeom, 'mass_dist') || length(branchGeom.mass_dist)~=3
        error('ConfigAdapter:InvalidData', '分枝 %s 的mass_dist无效', branchName);
    end
end

function fruit_params = buildFruitParamsStrict(preConfig, identifiedParams)
    fruit_params = struct();
    fruit_params.m = preConfig.fruit.mass;
    fruit_params.diameter = preConfig.fruit.diameter;
    k_val = 2000; c_val = 0.5;
    if isfield(preConfig.fruit, 'k_pedicel'), k_val = preConfig.fruit.k_pedicel; c_val = preConfig.fruit.c_pedicel; end
    fruit_params.k_pedicel_y = k_val; fruit_params.c_pedicel_y = c_val;
    fruit_params.k_pedicel_z = k_val; fruit_params.c_pedicel_z = c_val;
    fruit_params.F_break = 5.0; 
    
    det_model = [];
    if isfield(identifiedParams, 'branches')
        fn = fieldnames(identifiedParams.branches);
        for i = 1:length(fn)
            if isfield(identifiedParams.branches.(fn{i}), 'detachment_model')
                det_model = identifiedParams.branches.(fn{i}).detachment_model; break;
            end
        end
    elseif isfield(identifiedParams, 'detachment_model')
        det_model = identifiedParams.detachment_model;
    end
    
    if ~isempty(det_model)
        try
            pred = det_model.predict(1.5, 1.0, preConfig.fruit.diameter*100, 0);
            fruit_params.F_break = max(1.0, pred);
        catch
            fprintf('    [Warn] 脱落模型预测失败，使用默认值。\n');
        end
    end
end
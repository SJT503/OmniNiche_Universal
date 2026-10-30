# 跨器官动态特异性漂移与空间基质迟滞效应主导临床试验成熟转化与联合治疗协同

## 摘要

慢性器官疾病与实体瘤的晚期临床试验失败率超过85%，即便获人类遗传学证据支持或在健康参考图谱中呈高细胞类型表达特异性的靶点亦难以幸免。既有靶点筛选体系依赖静态健康基线表达，且将行政或商业性停研混同于生物学无效。本研究构建通用跨病种平台 OmniNiche，在12个非中枢神经系统（non-CNS）疾病轴中将多智能体临床停研解卷积、原代人类单细胞动态特异性漂移与空间多细胞随机朗之万动力学相耦合。对41,396项II-IV期干预性临床试验（其中29,791项解析至1,155个人类靶基因）的审计表明，41.0%的提前终止源于行政或入组因素而非靶点失效。在去混杂后的11,843项科学结局试验与23,967项确定性非行政试验中，健康基线特异性仅呈近随机区分度（受试者工作特征曲线下面积 AUROC 分别为 0.5283 与 0.5405），原因在于六大人源单细胞队列中的组织重塑诱发21.5%至37.1%的靶基因发生主导谱系室切换。将动态特异性漂移与药物模态依赖的基质屏障相容性整合后，折外预测 AUROC 提升至 0.6728（每标准差比值比 OR = 1.86，P = 5.70e-55），最高与最低四分位 III/IV 期成熟比扩大至 2.56 倍。仅基于原代人肝纤维化活检标定并零样本迁移至肺与肾纤维化的14项 IIb/III 期试验，OmniNiche 取得 AUROC = 0.9583，阐明了人源 CCR2 塌缩悖论、17.43 h 单药时间窗悬崖及 F4 期肝硬化双室 Bliss 协同效应（+0.9819）。结合12张人类癌组织 10x Visium 肌成纤维细胞胶原排斥图谱，模型在32例严格治疗前基线肿瘤活检中盲法预测抗 PD-1 应答达 AUROC = 0.7551，并筛选出破除基质屏障的逆转耐药联合方案。

**关键词**：临床试验损耗；多智能体临床解卷积；单细胞转录组；动态特异性漂移；空间转录组；Waddington 准势能景观；泛器官纤维化；肿瘤免疫治疗

## 引言

临床新药研发长期受困于高昂的晚期损耗，近90%的候选药物在进入I期临床至监管获批之间折戟 [1, 2]。过去十年的回顾性研究表明，拥有人类胚系遗传学证据支持的药物靶点获批概率约为无遗传支持靶点的两倍 [3, 4]，并推动了 Open Targets Platform 等大规模靶点-适应证知识图谱的建立 [5]。与此同时，伴随 Human Protein Atlas（HPA）单细胞转录组等健康人体参考图谱的完善 [6]，以 Yanai 特异性指数（$\tau_0$）[7] 度量的健康基线细胞类型特异性被广泛用作最大化在靶效力并规避旁观器官毒性的核心筛选指标。

尽管如此，两大结构性盲区制约了静态靶点指标在非神经系统实体器官疾病中的临床转化预测力。首先，公开临床试验注册库将研发管线调整、入组缓慢或资金撤回等非生物学停研，与真实的药效学无效及剂量限制性器官毒性混为一谈 [8]，直接在未分层的注册状态上开展基准评测会引入严重的结局标签噪声。其次，代谢功能障碍相关脂肪性肝炎（MASH）、特发性肺纤维化（IPF）、慢性肾脏病（CKD）、实体癌及黏膜自身免疫病中的慢性组织损伤会剧烈重塑细胞组成与细胞外基质结构。在纤维化与促纤维增生重塑过程中，浸润髓系细胞与驻留间质前体细胞发生保守的状态跃迁，转化为瘢痕相关 `SPP1+ TREM2+` 巨噬细胞及致密肌成纤维细胞屏障环 [9]，在瓦解健康 $\tau_0$ 谱系归属的同时，对单克隆抗体与细胞毒性淋巴细胞构筑起大分子空间位阻屏障。

近期人工智能虚拟细胞与虚拟组织建模路线图指出，准确预测治疗应答亟需超越静态转录组快照，转向具备跨器官与跨药物模态泛化能力的多细胞动力学系统 [10]。与此同时，多智能体科学推理架构已展现出在大规模生物医学证据提取与结构化审计中的高精度优势 [11]。然而，目前尚无任何框架能够将宏观尺度的临床停研原因解卷积、原代人类单细胞动态特异性漂移以及空间参数化的随机吸引子动力学融为一体。

本研究开发并验证了覆盖12个非中枢神经系统（non-CNS）疾病轴的端到端转化医学平台 OmniNiche（`Fig. 1a-c`）。该平台整合了确定性非神经疾病防火墙与多智能体临床停研审计器（`Fig. 1a`）、跨器官动态特异性漂移与生物物理屏障定律（`Fig. 1b`），以及经剖面似然可辨识性核验的多细胞朗之万随机微分方程（SDE）引擎（`Fig. 1c`）。我们证明，以动态重塑漂移与基质屏障门控取代静态健康 $\tau_0$，不仅显著提升了 II 期向 III/IV 期临床成熟转化的预测精度，实现了从原代人肝活检向肺与肾纤维化 IIb/III 期临床试验结局的零样本跨器官预测，还能严格基于治疗前肿瘤基线活检盲法预测患者对抗 PD-1 免疫检查点阻断的临床应答。

## 结果

### 41,396项非中枢神经系统临床试验的多智能体解卷积剥离行政损耗与生物学失败

为建立无学科交叉污染的全身及实体器官宏观临床基准，我们通过 ClinicalTrials.gov API v2 系统检索了12个非 CNS 疾病轴，涵盖实体肿瘤、代谢内分泌疾病、血液系统恶性肿瘤、肺纤维化与急性呼吸窘迫综合征（IPF/ARDS）、心血管心力衰竭与动脉粥样硬化、自身免疫性炎症性肠病与风湿病、慢性肾脏病与急性肾损伤（CKD/AKI）、皮肤免疫病、肌肉骨骼疾病、眼科视网膜疾病、肝 MASH 纤维化及感染性脓毒症（`Supplementary Table S1`）。执行自动化神经系统排斥防火墙（Firewall-2）拦截了587项伴发继发性神经系统适应证的记录，最终锁定41,396项经 SHA-256 哈希固化的唯一 II-IV 期干预性临床试验，包含23,536项 II 期、11,936项 III 期及5,924项 IV 期研究（`Fig. 2a` 与 `Extended Data Fig. 1a`）。

将干预药物与 Open Targets 药物作用机制库及 ChEMBL 分子层级交叉映射后，共有29,791项临床试验（72.0%）被精确解析至1,155个人类唯一靶基因（含586个唯一主靶基因；`Extended Data Fig. 1b,c` 与 `Supplementary Table S2`）。实体肿瘤（78.0%，797个靶基因）、血液肿瘤（83.3%，505个靶基因）、代谢内分泌疾病（72.6%，446个靶基因）与自身免疫病（71.1%，411个靶基因）的宿主靶点解析覆盖率均超过70%，而感染性脓毒症因包含大量直接针对病原体的抗微生物药物，宿主靶点解析率为26.2%（161个靶基因；`Fig. 2d`）。

随后，我们利用 OmniNiche 多智能体临床审计器对全部6,595项提前终止的靶点解析试验（`TERMINATED`、`WITHDRAWN` 或 `SUSPENDED`）进行停研归因分层。行政管理、商业管线调整及入组缓慢占全部停研试验的41.0%（2,706项），未报告或未分类停研占47.3%（3,118项），而明确的生物学失败包含466项安全性或器官毒性停研（7.1%）及305项疗效无效停研（4.6%）（`Fig. 2b`）。剔除行政与未报告停研后，共获得23,967项确定性非行政结局试验（含11,072项已完成 III/IV 期、12,124项已完成 II 期及771项明确科学失败）。按药物模态与临床分期分层可见，抗体与生物大分子在 II 期（95.3% 对比 96.4%）及 III 期（96.7% 对比 97.3%）的科学完成率均低于小分子药物（`Fig. 2c`），提示必须在靶点表达动态之外显式刻画大分子组织穿透受限。

### 跨器官单细胞特异性漂移与生物物理屏障门控显著超越健康静态特异性

为检验剔除行政噪声后健康基线细胞类型特异性能否预测临床转化，我们基于 HPA 单细胞图谱计算了20,151个人类蛋白编码基因在154种健康细胞类型中的 $\tau_0$ 与重要脏器旁观风险（`Supplementary Table S3`）。在对比 III/IV 期完成药物（`n = 11,072`）与明确科学无效或毒性失败（`n = 771`）的11,843项去混杂科学结局试验中，静态健康 $\tau_0$ 仅表现出近随机折外区分度（`AUROC = 0.5283`，Brier 评分 `= 0.0608`，每标准差比值比 `OR = 1.09`，95% 置信区间 `1.01-1.17`，`P = 0.0200`；`Fig. 3a,b` 与 `Table 1`）。在直接配对原代人类单细胞数据的6个核心疾病轴（`n = 5,815` 项试验）中，静态健康 $\tau_0$ 甚至未达统计学显著性（`AUROC = 0.5252`，`OR = 1.08`，95% 置信区间 `0.99-1.18`，`P = 0.0894`）。

为定量刻画人类活动性疾病如何重塑靶点特异性，我们统一流式处理了6个非 CNS 组织重塑轴的原代人类单细胞与单核转录组队列：肝 MASH 与肝硬化（`GSE136103` 与 `GSE212837`，17例供体）[12, 13]、肺纤维化 IPF（`GSE135893`，30例供体，114,396个细胞）[14]、慢性肾脏病与糖尿病肾病（`GSE271188`，17例肾活检）、非小细胞肺癌（`GSE207422`，约92,000个细胞）[15]、溃疡性结肠炎黏膜活检（`GSE250487`，93,900个结肠细胞）以及人类细菌性脓毒症纵向队列（`GSE167363`，12份样本）[16]（`Extended Data Fig. 2a`）。将全部1,155个临床靶基因映射至实质损伤室、髓系炎症室与肌成纤维基质室后，发现了广泛的疾病诱导谱系室切换（`Extended Data Fig. 2b,c`）：从健康或早期病变演进至晚期重塑时，肝 MASH 纤维化中27.8%的靶基因发生主导表达细胞室切换，肺 IPF 为21.5%，肾 CKD 为30.2%，实体瘤为23.2%，溃疡性结肠炎为28.2%，细菌性脓毒症高达37.1%（`Fig. 3d`）。

将单细胞动态特异性漂移（$\Delta\tau = \tau_{\text{late}} - \tau_{\text{early}}$ 及异位基质重分布）与药物模态依赖的生物物理屏障相容性（$B_{\text{modality}}$，对靶向致密基质包裹隔室的高分子量生物制剂施以空间位阻衰减）相耦合，在全部三个评测层级均取得大幅性能跃升（`Fig. 3a-c` 与 `Table 1`）。在12轴去混杂科学结局队列（`n = 11,843`）中，OmniNiche 将折外 AUROC 提升 `+0.1445` 至 `0.6728`（Brier `= 0.0594`，`OR = 1.86`，95% 置信区间 `1.72-2.01`，`P = 5.70e-55`）；在6个核心单细胞疾病轴（`n = 5,815`）中同样提升 `+0.1468` 至 `0.6720`（`OR = 1.80`，`P = 3.75e-32`）。在涵盖23,967项确定性非行政试验的 II 期向 III/IV 期临床成熟基准中，OmniNiche 的 AUROC 达 `0.6747`（静态 $\tau_0$ 为 `0.5405`，Brier `0.2279` 对比 `0.2474`，`OR = 1.87`，95% 置信区间 `1.82-1.93`，`P < 1e-300`），将最高与最低四分位（`Q4/Q1`）的 III/IV 期成熟比由静态 $\tau_0$ 的 `1.20 倍`（`53.4%` 对比 `44.6%`）扩展至 `2.56 倍`（`66.2%` 对比 `25.9%`；`Fig. 3c`）。各疾病分轴效应量一致放大，其中自身免疫病（`AUROC = 0.6993`，`OR = 2.07`，`P = 2.26e-33`）、肺 IPF/ARDS（`AUROC = 0.6104`，`OR = 1.49`，`P = 9.31e-17`）与实体肿瘤（`OR = 1.37`，`P = 1.16e-39`）尤为显著（`Fig. 3b`）。

### 零样本跨器官朗之万动力学破解 III 期抗纤维化临床失败与 F4 期协同机制

尽管动态漂移与屏障相容性解释了宏观临床成熟规律，晚期抗纤维化临床试验即便在小分子达成充分靶点结合时仍频繁折戟，提示存在非线性多细胞正反馈锁定。我们首先剖析了一个典型的高 $\tau_0$ 临床失败悖论：双靶点 CCR2/CCR5 拮抗剂 cenicriviroc（`CCR2` 健康 $\tau_0 = 0.9513`）在纳入1,778例 MASH 患者的 III 期 AURORA 试验（`NCT03028740`）中未能逆转肝纤维化 [17]。追踪 `GSE136103` 与 `GSE212837` 原代人类髓系细胞轨迹（`n = 19` 份供体剖面；`Supplementary Table S4`）可见，`CCR2` 在 NAFLD 外周血单核细胞中表达达峰值（`6.19 log2 CPM`），在健康驻留 Kupffer 巨噬细胞中降至 `4.27 log2 CPM`，而在 F2/F3 期 MASH 瘢痕相关巨噬细胞中骤降4.2倍至 `1.47 log2 CPM`（F4 肝硬化中亦仅 `4.40 log2 CPM`）；与此同时，`SPP1`（由 `0.00` 升至 `8.31 log2 CPM`）与 `TREM2`（由 `0.39` 升至 `5.77 log2 CPM`）在纤维化瘢痕内从头涌现（`Fig. 4a`）。由此可见，当人肝内微环境已跃迁至自我维持的 `SPP1+ TREM2+` 巨噬细胞-肌成纤维细胞回路后，静态健康 $\tau_0$ 所推崇的上游募集受体已丧失组织内靶向基础。

为在不引入跨器官过拟合的前提下定量刻画该多细胞回路，我们仅利用17例原代人肝供体标定了耦合实质损伤室、`SPP1+` 瘢痕巨噬细胞室与 `COL1A1+/ACTA2+` 肌成纤维细胞室的三室 Itô-朗之万 SDE 系统，并将锁定的动力学参数零样本（`0` 参数重调）迁移至原代人肺 IPF（`n = 30` 例供体）与肾 CKD/DKD（`n = 17` 例供体），盲法评测14项标志性 IIb/III 期抗纤维化临床试验（`Table 2` 与 `Supplementary Table S4`）。该基准涵盖 resmetirom（`MAESTRO-NASH`，`NCT03900429`）[18]、semaglutide（`ESSENCE`）、lanifibranor（`NATIVE`）、nintedanib（`INPULSIS-1`）[19]、pirfenidone（`ASCEND`）与 finerenone（`FIDELIO-DKD`）[20] 等阳性标杆，以及 cenicriviroc（`AURORA`）[17]、selonsertib（`STELLAR-3` 与 `STELLAR-4`）[21]、simtuzumab（肝 F3、肝 F4 及肺 IPF `RAINIER`）、elafibranor（`RESOLVE-IT`）与 pamrevlumab（`ZEPHYRUS-1`）[22] 等阴性试验，并以 bardoxolone（`BEACON`）[23] 作为4期肾毒性动力学参照。

在这一跨器官零样本基准（`人肝 -> 人肺 IPF 与人肾 CKD/DKD`）中，OmniNiche 在所有指标上均位列第一（`AUROC = 0.9583`，精确率-召回率曲线下面积 `AUPRC = 0.9583`，Brier 评分 `= 0.1257`），大幅领先静态健康 $\tau_0$（`B1`：`AUROC = 0.5000`，Brier `= 0.3181`）、多靶点计数（`B2`：`AUROC = 0.6667`，Brier `= 0.2502`）、单室基因调控网络保持度（`B3`：`AUROC = 0.5000`，Brier `= 0.2587`）、静态配体-受体中心度（`B4`：`AUROC = 0.9167`，Brier `= 0.1375`）及单纯模态先验（`B5`：`AUROC = 0.7500`，Brier `= 0.2364`）（`Fig. 4b` 与 `Table 2`）。针对4条有向跨室耦合权重的剖面似然 $\Delta\chi^2$ 扫描 [24] 证实，3条核心致病边（`Parenchyma -> SPP1_SAMac`、`SPP1_SAMac -> Myofibroblast` 与 `Myofibroblast -> Parenchyma`）均远超95%置信阈值（$\Delta\chi^2 > 3.841$，实测最小双侧 $\Delta\chi^2 = 81.39$ 至 `348.63`），而反向 `Myofibroblast -> SPP1_SAMac` 放大边呈现单侧实践松弛性（`Extended Data Fig. 3a,b`）。

将标定后的 Fokker-Planck 稳态密度投影至二维 Waddington 准势能景观 $U(\mathbf{S}) = -\ln(P_{\text{ss}}(\mathbf{S}) + \varepsilon)$ [25]，清晰揭示出分隔组织消退盆地 `(0.12, 0.12)` 与深层 F4 肝硬化锁定盆地 `(0.72, 0.72)` 的双稳态迟滞鞍点（`Fig. 4d`）。模拟延迟启动上游实质单药干预表明，病程演化至 `t50 = 17.43 h` 等效轨迹时间处存在陡峭的时间窗悬崖，越过该阈值后单药缓解率迅速衰减，并在进入 F4 锁定吸引子后坍塌至 `0.9%`（`0.0091`；`Fig. 4c`）。在锁定的 F4 肝硬化吸引子内开展系统性双药虚拟筛选表明，将低剂量肌成纤维细胞 `TGFBR1/PDGFRB` 小分子破障剂与 `SPP1/TREM2` 巨噬细胞调节剂联用可压低迟滞鞍点，恢复 `100.0%` 的模拟盆地逃逸（`Bliss 超额协同 = +0.9819`，最高单药 `HSA 超额 = +0.9909`）；反之，联用受空间位阻排斥的 `LOXL2` 抗体 simtuzumab 与靶点塌缩的 `CCR2` 拮抗剂 cenicriviroc 仅产生 `0.0%` 挽救率（`Bliss 超额 = 0.0000`；`Fig. 4c` 与 `Supplementary Table S4`）。

### 空间胶原-肌成纤维细胞排斥效应主导实体癌治疗前抗 PD-1 盲法应答预测与联合破障挽救

为将 OmniNiche 空间屏障体系推广至实体瘤免疫治疗，即 TGF-beta 活化的肌成纤维样肿瘤相关成纤维细胞（`myCAF`）在物理空间上将 `CD8+` 细胞毒性 T 细胞排斥于恶性上皮巢之外的场景 [26, 27]，我们解析了12张原代人类癌组织 10x Visium 空间转录组切片（共 `21,796` 个空间位点），涵盖膀胱癌（`GSE171351`，4张切片）[28] 以及乳腺癌、卵巢癌、肝细胞癌与胃肠间质瘤（`GSE203612`，8张切片）[29]（`Supplementary Table S5`）。跨全部12张人源切片的空间径向基邻近图谱表明，免疫抑制性 `SPP1+/CD163+` 肿瘤相关巨噬细胞（`TAM`）与 `ACTA2+/COL1A1+/TGFB1+` `myCAF` 屏障环紧密共定位（共识空间邻近权重 `= 0.6661`），而恶性肿瘤核心与 `CD8A+/GZMB+/PRF1+` 细胞毒性 T 细胞位点呈显著物理隔离（共识邻近权重 `= 0.4710`，配对 Wilcoxon 符号秩检验 `P = 0.0002`；`Fig. 5a,b`）。

基于这一实测12切片 Visium 空间邻近张量、Hill 型 `myCAF` 大分子空间排斥门控及 Holling II 型效靶比（`E:T`）肿瘤负荷饱和项，我们构建了四室空间肿瘤 SDE（`Tumor`、`TAM`、`myCAF` 与 `CD8_T`），并在人源抗 PD-1 队列中开展零信息泄露的盲法应答基准评测（`Fig. 5c`、`Table 2` 与 `Supplementary Table S5`）。在皮肤基底细胞癌队列（`GSE123813`，`n = 8` 例严格治疗前全微环境单细胞活检）[30] 中，不依赖流式分选比例的肿瘤反应性效应 T 细胞加权有效消除了独立 `CD45+ CD3+` 富集管引入的静息旁观 T 细胞虚高，取得 `AUROC = 0.8750`（`AUPRC = 0.8750`），全面超越健康 $\tau_0$（`B2 = 0.5625`）、单室 T 细胞基因调控网络评分（`B3 = 0.5000`）与静态 `PD-L1`（`CD274`）表达（`B1 = 0.4375`）。对 `GSE123813` 中15例配对治疗前后纵向患者的追踪进一步证实，临床应答者在阻断治疗后呈现特异性的细胞毒性扩增（`Extended Data Fig. 3c`）。

我们随后在独立的非小细胞肺癌（`NSCLC`）新辅助免疫联合化疗队列（`GSE207422`）[15] 中开展验证，该队列同时包含高深度微环境单细胞转录组（`n = 14` 例患者，约 `92,000` 个细胞）与严格治疗前基线诊断性穿刺转录组活检（`n = 24` 例患者）。在 `GSE207422` 单细胞队列中，OmniNiche 在全部三项指标上均居首位（`AUROC = 0.8500`，`AUPRC = 0.6083`，Brier `= 0.1751`）；在24例严格治疗前 NSCLC 诊断活检中，OmniNiche 取得了全部6种对比模型中最低的概率校准误差（`Brier = 0.2250`，`AUROC = 0.7407`，`AUPRC = 0.6025`）。将两个癌种全部32例严格治疗前基线患者活检汇总评测（`0%` 术后切除标本混入；`n = 13` 例应答者，`n = 19` 例无应答者），OmniNiche 在全部三项指标上均位列第一（`AUROC = 0.7551`，`AUPRC = 0.6385`，Brier `= 0.2366`），显著优于静态健康 $\tau_0$（`AUROC = 0.6194`，Brier `= 0.3153`）、单室 T 细胞调控网络（`AUROC = 0.6113`，Brier `= 0.2890`）及静态 `PD-L1`（`AUROC = 0.5344`，Brier `= 0.2896`）（`Fig. 5c` 与 `Table 2`）。

最后，针对19例抗 PD-1 耐药治疗前活检计算的边缘 Waddington 准势能曲线 $U(S_{\text{Tumor}})$ 揭示了免疫检查点单药在基质排斥型肿瘤中的失效根源：在抗 PD-1 单药作用下，`myCAF`-`TAM` 护盾将系统轨迹囚禁于高肿瘤负荷吸引子内（`S_Tumor ~ 0.74`，模拟主要病理缓解挽救率仅为 `16.3%`；`Fig. 5d`）。联合给予小分子 `myCAF/TGFBR1` 基质屏障重塑剂可瓦解促纤维增生排斥屏障，恢复 `CD8+` T 细胞对肿瘤核心的接触并将全局势能极小值推移至主要病理缓解阈值（`S_Tumor <= 0.36`）以下，使耐药患者的模拟缓解挽救率由 `16.3%` 跃升至 `99.2%`（`Fig. 5d`）。

## 讨论

本研究通过对41,396项非中枢神经系统临床试验开展多智能体停研解卷积，并将其与原代人类单细胞重塑轨迹及空间受限朗之万动力学相融合，在分子靶点遴选与多细胞临床疗效之间建立起定量桥梁。宏观临床审计表明，逾五分之二（`41.0%`）提前终止的靶点解析试验系因行政管理、商业策略或入组缓慢而关闭，与靶点药理学无关。在11,843项科学结局试验与23,967项确定性临床成熟试验中剥离上述行政混杂后，健康基线细胞类型特异性（$\tau_0$）仅具有近随机的临床区分能力（`AUROC = 0.5283` 至 `0.5405`）。这一局限直接源于疾病的动态组织重塑：在覆盖纤维化、肿瘤、自身免疫与脓毒症的六大原代人源单细胞队列中，高达 `21.5%` 至 `37.1%` 的临床靶基因在病程演进中发生了主导细胞室切换。

以单细胞动态特异性漂移与模态特异生物物理屏障相容性的乘积取代静态靶点指标，可提炼出三条具有直接转化价值的规律。第一，当疾病演进导致靶谱系转录沉默时，极高的健康特异性反而会成为临床陷阱。人 MASH 进展过程中从外周循环单核细胞向肝内 `SPP1+ TREM2+` 瘢痕巨噬细胞转化时发生的4.2倍 `CCR2` 表达塌缩，圆满解释了尽管 `CCR2` 拥有顶格健康 $\tau_0$（`0.9513`），cenicriviroc 仍在 III 期临床折戟的原因；相反，在重塑进程中保持或获得实质与星状细胞特异性的 `THRB` 与 `GLP1R` 则展现出稳健的晚期临床疗效。

第二，大分子药物模态与纤维化及促纤维增生基质结构之间存在强烈的非线性交互。细胞外胶原交联与致密 `myCAF` 基质护盾通过空间位阻削弱单克隆抗体穿透与 `CD8+` 淋巴细胞浸润，解释了靶向基质的生物大分子（抗 `LOXL2` 抗体 simtuzumab 与抗 `CCN2` 抗体 pamrevlumab）相较于小分子激酶与受体调节剂（`nintedanib`、`pirfenidone` 与 `finerenone`）在 IIb/III 期试验中反复受挫的现象。将上述生物物理约束嵌入仅由人肝数据标定的朗之万 SDE 后，无需针对器官重调任何动力学参数即可零样本迁移至人肺 IPF 与人肾 CKD/DKD 临床试验（`AUROC = 0.9583`）。

第三，Waddington 准势能分析表明，晚期器官纤维化（F4 肝硬化）与促纤维增生性实体癌均处于被高耸势能鞍点与缓解盆地隔绝的双稳态迟滞吸引子中。一旦越过 `17.43 h` 等效病程时间窗悬崖，上游单靶点阻断已无法瓦解自我维持的巨噬细胞-肌成纤维细胞正反馈环（F4 肝硬化单药模拟缓解率仅 `0.9%`；抗 PD-1 耐药肿瘤单药挽救率仅 `16.3%`）。跨越迟滞鞍点必须依靠正交双室联合方案，同步解除基质屏障正反馈（`TGFBR1/PDGFRB` 抑制）并重编程髓系或细胞毒性效应室（`TREM2/SPP1` 调节或 `PD-1` 阻断）。

本研究仍存在若干局限性并为前瞻性临床验证指明了方向。尽管宏观图谱解析了12个非 CNS 疾病轴中的29,791项临床试验，ClinicalTrials.gov 上仍有47.3%的停研记录缺乏细粒度公开停研说明，凸显出推行标准化无效与入组停研强制披露的必要性。此外，剖面似然审计证实核心纤维化4条反馈边中有3条具备双侧实践可辨识性，但肌成纤维细胞对巨噬细胞的反向放大边在横断面活检采样下仍呈部分参数松弛性，提示未来 II 期篮式试验需纳入纵向空间穿刺采样。最后，尽管汇总的32例免疫治疗盲法评测严格限定于治疗前基线活检，未来仍需开展将治疗前空间 Visium 活检分层与 `TGFBR1` 基质破障联合方案相配对的多中心前瞻性临床试验，以最终确证患者的长期生存获益。

## 方法

### 非中枢神经系统临床试验普查、Firewall-2 防火墙与多智能体停研解卷积

通过 ClinicalTrials.gov API v2 检索12个预设非神经疾病轴（`Solid_Oncology`、`Hepatic_MASH_Fibrosis`、`Pulmonary_IPF_ARDS`、`Renal_CKD_AKI`、`Autoimmune_IBD_Rheum`、`Cardiovascular_HF_Athero`、`Metabolic_Endocrine`、`Infectious_Sepsis`、`Hematologic_Malignancy`、`Dermatologic_Immune`、`Musculoskeletal_Bone` 与 `Ophthalmic_Retinal`）中的 II、III、IV 期干预性临床试验，并对每条原始记录生成不可篡改的 SHA-256 哈希锁。为确保与中枢神经系统研究保持零选题及零数据重叠（`Firewall-2`），自动化词法与医学主题词排斥过滤器对 `conditions`、`brief_title` 与 `keywords` 字段进行全量扫描，拦截包含脑、脊髓、卒中、阿尔茨海默病、帕金森病、胶质瘤、多发性硬化、创伤性脑损伤、癫痫及神经病理性疼痛等词条的587项记录，最终保留41,396项唯一非 CNS 临床试验。

利用 Open Targets Platform（`24.09` Parquet 版本 `molecule` 与 `mechanismOfAction` 数据表）及 ChEMBL 靶点层级映射干预药物名称与同义词，赋予人类靶基因符号（`HGNC`）及治疗模态（`Small molecule`、`Antibody`、`Protein`、`Oligo`、`Gene` 或 `Cell`）。针对全部6,595项提前终止的靶点解析试验（`TERMINATED`、`WITHDRAWN` 或 `SUSPENDED`），采用确定性多智能体词法-语义分类器联合审计 `why_stopped` 文本、实际入组人数与结果表，将其划分为四类互斥结局：`Administrative_Commercial_Enrollment`（入组缓慢、资金撤回、商业管线调整或申办方变更）、`Scientific_Efficacy_Futility`（期中无效分析停研或主终点阴性）、`Scientific_Safety_Toxicity`（不良事件、剂量限制性器官毒性或安全委员会叫停）以及 `Unreported / Other`。

### 健康基线特异性与六轴原代人类单细胞动态漂移计算

基于 Human Protein Atlas 单细胞转录组共识数据集（20,151个人类蛋白编码基因 × 154种健康单细胞类型），按 Yanai 特异性指数计算健康基线细胞类型特异性 $\tau_0$：

$$\tau = \frac{\sum_{i=1}^{N} \left(1 - \frac{x_i}{\max_{1 \le j \le N} x_j}\right)}{N - 1}$$

其中 $x_i = \log_2(\text{nTPM}_{i} + 1)$ 为细胞类型 $i$（$N = 154$）的标准化表达量。重要脏器旁观风险定义为非同源适应证中心肌细胞、肝细胞、肾近端小管上皮细胞与肺泡 II 型上皮细胞的归一化表达占比。

在严格内存上限（峰值 RAM `< 7 GB`）流式架构下处理 NCBI GEO 六大非 CNS 原代人类单细胞及单核转录组队列：肝硬化与 MASH（`GSE136103`、`GSE212837`）、肺 IPF（`GSE135893`）、肾 CKD/DKD（`GSE271188`）、NSCLC 实体瘤（`GSE207422`）、溃疡性结肠炎（`GSE250487`）与细菌性脓毒症（`GSE167363`）。将各供体细胞聚合为标准化三室或四室组织伪批量剖面（`Parenchyma_Injury`、`SPP1_SAMac` / `Myeloid`、`Myofibroblast` / `myCAF`，肿瘤队列另含 `CD8_T`）并转换为 $\log_2(\text{CPM} + 1)$。对全部1,155个临床靶基因计算早期与晚期疾病特异性（$\tau_{\text{early}}$、$\tau_{\text{late}}$）、动态特异性漂移（$\Delta\tau = \tau_{\text{late}} - \tau_{\text{early}}$）、主导谱系室切换及异位基质占比。模态依赖的生物物理屏障相容性（$B_{\text{modality}} \in [0, 1]$）刻画分子量依赖的纤维化及促纤维增生基质空间位阻衰减。采用5折分层交叉验证 L2 正则化逻辑回归生成试验级折外预测概率，报告折外 AUROC、Brier 评分、每标准差增量 Wald OR 及四分位（`Q4/Q1`）成熟比。

### 多细胞朗之万 SDE 引擎、剖面似然可辨识性审计与 Waddington 准势能构建

利用耦合 Itô-朗之万随机微分方程刻画归一化组织隔室状态 $\mathbf{S}(t) \in [0, 1]^K$ 的动态演化，并以 Euler-Maruyama 格式离散积分（步长 $\Delta t = 0.05\text{ h}$，演化时域 $T = 48\text{ h}$，每供体或给药方案模拟 `220-400` 条蒙特卡洛轨迹）：

$$dS_i = \left[ \sum_{j} v_{ji} \frac{S_j^n}{k_{ji}^n + S_j^n} - d_i S_i - \eta_i(t, \text{drug}, B_{\text{modality}}) S_i \right] dt + \sigma_i \sqrt{\max(S_i(1 - S_i), 0)} \, dW_i(t)$$

其中 $v_{ji}$ 与 $k_{ji}$ 为 Hill 激活权重与半饱和常数（$n = 2$），$d_i$ 为隔室内禀周转率，$\eta_i$ 为经基质屏障门控的药效学靶点扰动强度，$\sigma_i dW_i(t)$ 为状态有界的微环境维纳噪声。在泛器官纤维化（`K = 3`）中，动力学参数仅针对17例原代人肝供体（`GSE136103` 与 `GSE212837`）标定，随后零样本（`0` 参数重拟合）迁移至人肺 IPF（`GSE135893`，`n = 30`）与人肾 CKD/DKD（`GSE271188`，`n = 17`）以预测14项 IIb/III 期抗纤维化试验结局。在 `[0.0, 1.25]` 区间扫描各条有向耦合边权重 $v_{ji}$ 并重优化其余参数开展剖面似然审计，双侧 $\Delta\chi^2 > 3.841$（$\chi^2_1$）判定为95%置信水平实践可辨识。基于终末状态经验密度直方图按 $U(\mathbf{S}) = -\ln(P_{\text{ss}}(\mathbf{S}) + 10^{-6})$ 重建二维及边缘 Waddington 准势能面，并以 Bliss 独立超额协同度与 HSA 超额协同度量化 F4 肝硬化联合用药协同效应。

### 实体瘤 10x Visium 空间邻近参数化与治疗前抗 PD-1 盲法基准评测

对12张原代人类癌组织 10x Visium 空间切片（`GSE171351` 膀胱癌4张切片；`GSE203612` 乳腺癌、卵巢癌、肝癌及胃肠间质瘤8张切片；共 `21,796` 个空间位点）逐位点计算 `Tumor`、`SPP1+/CD163+` `TAM`、`ACTA2+/COL1A1+/TGFB1+` `myCAF` 与 `CD8A+/GZMB+/PRF1+` `CD8_T` 特征评分。基于 Visium 阵列坐标按高斯径向基衰减函数（特征尺度 $\ell = 2.5$ 倍点距）计算跨隔室空间邻近权重，并采用双侧 Wilcoxon 符号秩检验比较12张切片的 `TAM <-> myCAF` 与 `Tumor <-> CD8_T` 邻近度。

在四室实体瘤朗之万 SDE 中，肿瘤自主 Allee 自更新（`v = 0.295, k = 0.38`）与受 Hill 型 `myCAF` 位阻门控及 Holling II 型效靶比负荷饱和双重约束的 `CD8_T` 杀伤相拮抗，在主要病理缓解阈值 $S_{\text{Tumor}}(T) \le 0.36$ 处通过平滑观测似然评估缓解概率。采用跨平台中位数绝对偏差（`MAD`）逻辑斯蒂变换将患者基线状态无监督映射至 `[0, 1]^4`（零临床应答标签参与）。在 `GSE123813`（`n = 8` 例严格治疗前全微环境基底细胞癌单细胞活检，采用流式不敏感肿瘤反应性 T 细胞加权校正独立 `CD45+ CD3+` 分选管效应，另含 `n = 15` 例配对治疗前后轨迹）、`GSE207422` NSCLC 微环境单细胞队列（`n = 14`）、`GSE207422` 严格治疗前 NSCLC 诊断活检队列（`n = 24`）及汇总严格治疗前活检队列（`n = 32`）中完成盲法应答评测。

### 统计分析

采用精确 Mann-Whitney 法计算 AUROC 与精确率-召回率 AUPRC，采用 Brier 均方误差评估概率校准度。逻辑回归效应量报告为每增加1个标准差的 Wald 比值比（OR）、双侧95%置信区间及精确双侧 P 值。Visium 切片配对空间邻近比较采用双侧 Wilcoxon 符号秩检验。全部计算管线均通过27项自动化单元测试验证数学不变量、零标签泄露及 `Firewall-2` 神经系统隔离。

## 数据与代码可用性

本研究全量临床试验数据获取自 ClinicalTrials.gov API v2（`https://clinicaltrials.gov/api/v2/studies`）。药物靶点与作用机制映射数据获取自 Open Targets Platform（`https://platform.opentargets.org`）。Data from the Open Targets Platform (`https://platform.opentargets.org`) is freely available for academic and commercial use under the Open Targets licence (`https://platform-docs.opentargets.org/licence`). 健康人体单细胞参考图谱获取自 Human Protein Atlas（`https://www.proteinatlas.org`）。原代人类单细胞、单核、穿刺批量转录组及 10x Visium 空间转录组数据集均公开存储于 NCBI GEO 数据库，检索号为 `GSE136103`、`GSE212837`、`GSE135893`、`GSE271188`、`GSE207422`、`GSE250487`、`GSE167363`、`GSE171351`、`GSE203612` 与 `GSE123813`。OmniNiche v2.0 完整源代码、SHA-256 锁定汇总数据表、制图脚本及27项单元测试套件均归档于规范项目仓库（`F:/OmniNiche_Universal`）。

## 表格

### Table 1 | 29,791项靶点解析非中枢神经系统临床试验的跨病种区分度与 III/IV 期临床成熟预测表现

| 基准评测队列 | 临床试验数 (`n`) | 阳性结局数 (`n`) | 静态健康 $\tau_0$ AUROC | 疾病早期 $\tau$ AUROC | 纯漂移（无屏障）AUROC | OmniNiche 漂移 $\times$ 屏障 AUROC | 相对 $\tau_0$ AUROC 增益 | 静态 $\tau_0$ 每标准差 OR (95% CI; *P*) | OmniNiche 每标准差 OR (95% CI; *P*) | Q4/Q1 成熟比 ($\tau_0$ 对比 OmniNiche) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 去混杂科学结局（12个非 CNS 轴：III/IV 期完成对比明确科学失败） | 11,843 | 11,072 | 0.5283 | 0.6081 | 0.6073 | **0.6728** | +0.1445 | 1.09 (1.01-1.17; 2.00e-2) | **1.86** (1.72-2.01; 5.70e-55) | 1.02x vs **1.10x** |
| 去混杂科学结局（6个直接配对原代人类单细胞核心轴） | 5,815 | 5,276 | 0.5252 | 0.6434 | 0.6430 | **0.6720** | +0.1468 | 1.08 (0.99-1.18; 8.94e-2) | **1.80** (1.63-1.98; 3.75e-32) | 1.02x vs **1.16x** |
| 确定性非行政 II 期向 III/IV 期临床成熟转化（全部12个非 CNS 轴） | 23,967 | 11,072 | 0.5405 | 0.6133 | 0.6129 | **0.6747** | +0.1342 | 1.14 (1.11-1.17; 6.37e-25) | **1.87** (1.82-1.93; <1e-300) | 1.20x vs **2.56x** |

### Table 2 | 泛器官抗纤维化 IIb/III 期临床试验零样本跨器官基准与实体瘤治疗前抗 PD-1 盲法应答基准

| 转化基准评测场景 | 队列与迁移设置 | 样本量 (`n`) | OmniNiche v2.0 AUROC / AUPRC (Brier) | 静态健康 $\tau_0$ AUROC (Brier) | 单室 GRN 保持度 AUROC (Brier) | 静态 L-R / PD-L1 基线 AUROC (Brier) | 最优机制性联合挽救方案 |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| 泛器官纤维化 IIb/III 期临床试验 | 零样本 `人肝 (n=17) -> 人肺 IPF (n=30) 与人肾 CKD/DKD (n=17)` | 14项 IIb/III 期试验 | **0.9583** / **0.9583** (**0.1257**) | 0.5000 (0.3181) | 0.5000 (0.2587) | 0.9167 (0.1375) [L-R 中心度] | F4 肝硬化：`TGFBR1/PDGFRB + TREM2/SPP1`（`100%` 逃逸，Bliss `+0.9819`） |
| 皮肤基底细胞癌抗 PD-1 治疗前单细胞活检 (`GSE123813`) | 严格治疗前全微环境活检（流式去混杂） | 8例患者 (4 R, 4 NR) | **0.8750** / **0.8750** (0.2716) | 0.5625 (0.2939) | 0.5000 (0.3569) | 0.4375 (0.3892) [PD-L1] | Visium `myCAF` 屏障破除恢复 `CD8+` 核心接触 |
| NSCLC 新辅助免疫治疗微环境单细胞队列 (`GSE207422`) | 高深度人 NSCLC 约92,000细胞验证队列 | 14例患者 (4 R, 10 NR) | **0.8500** / **0.6083** (**0.1751**) | 0.7500 (0.2327) | 0.7250 (0.2280) | 0.3250 (0.3312) [PD-L1] | `抗 PD-1 + myCAF/TGFBR1` 重塑终末 $U(S_{\text{Tumor}})$ |
| 汇总严格治疗前基线肿瘤活检 (`GSE123813 + GSE207422`) | 跨 BCC 与 NSCLC 零泄露治疗前诊断性活检 | 32例患者 (13 R, 19 NR) | **0.7551** / **0.6385** (**0.2366**) | 0.6194 (0.3153) | 0.6113 (0.2890) | 0.5344 (0.2896) [PD-L1] | 耐药患者挽救 (`n=19`)：单药 `16.3%` 跃升至联合 **99.2%** |

## 图注

**Figure 1 | OmniNiche 跨病种临床解卷积与空间多细胞动力学平台总览。** **a**，跨12个非 CNS 器官轴的多智能体临床审计将临床试验档案解析为人类靶基因与药物模态（小分子对比单克隆抗体），并将提前停研分层为行政或入组停研对比科学无效与安全性失败。**b**，耦合单细胞动态特异性漂移（$\Delta\tau$，捕捉健康基线微环境向慢性组织重塑演进时的主导谱系室切换）与生物物理基质屏障门控（$B_{\text{modality}}$，刻画围绕 `SPP1+` 巨噬细胞的 `myCAF` 胶原基质护盾对 `IgG4` 抗体与 `CD8+` 细胞毒性 T 细胞的空间排斥）。**c**，随机 Itô-朗之万微环境动力学重建 Waddington 准势能景观 $U(\mathbf{S})$，识别晚期重塑中的单药迟滞陷阱并提名压低迟滞鞍点的破障联合疗法。

**Figure 2 | 41,396项非 CNS II-IV 期临床试验的宏观普查与多智能体停研解卷积。** **a**，按临床分期（II 期 `23,536` 项、III 期 `11,936` 项、IV 期 `5,924` 项）堆叠展示12个非 CNS 疾病轴共41,396项经 SHA-256 锁定的唯一干预性试验。**b**，对全部6,595项提前终止的靶点解析试验进行多智能体停研解卷积，剥离行政与入组停研（`2,706` 项，`41.0%`）及未报告或其他停研（`3,118` 项，`47.3%`）与明确的科学安全性/毒性（`466` 项，`7.1%`）及疗效无效（`305` 项，`4.6%`）失败。**c**，按药物模态类别与临床分期分层的确定性科学完成率矩阵（`n = 23,967` 项确定性非行政试验）。**d**，12个非 CNS 轴的 Open Targets 靶点解析覆盖率（总体 `72.0%`，`29,791` 项试验）与唯一人类靶基因多样性（共 `1,155` 个基因）气泡图，圆面积正比于各轴总试验数。

**Figure 3 | 跨器官动态特异性漂移与生物物理屏障相容性主导临床试验成熟转化。** **a**，四个嵌套模型（`Static Healthy tau0`、`Disease Early tau`、`Drift (No Barrier)` 与 `OmniNiche Drift x Barrier`）在12轴去混杂科学结局（`n = 11,843`）、6个核心单细胞轴（`n = 5,815`）及确定性 II 期向 III/IV 期临床成熟（`n = 23,967`）中的折外 AUROC 对比。**b**，静态健康 $\tau_0$ 与 OmniNiche 在全局队列及6个原代人类单细胞疾病轴中的标准化效应量（每标准差增量 OR）森林图。**c**，23,967项确定性试验在预测评分四分位（最低 `Q1` 至最高 `Q4`）上的实测 II 期向 III/IV 期成熟率曲线，OmniNiche 实现 `2.56 倍` `Q4/Q1` 富集（`66.2%` 对比 `25.9%`），而静态健康 $\tau_0$ 仅为 `1.20 倍`。**d**，人肝 MASH 重塑特异性漂移图谱展示健康 HPA $\tau_0$ 与疾病特异性漂移（$\Delta\tau$），标示 `27.8%` 的谱系室切换及代表性临床靶点（`GLP1R`、`CCR2`、`ACTA2`、`PDGFRB`、`THRB` 与 `ABCC8`）。

**Figure 4 | 原代人源 MASH CCR2 塌缩悖论、零样本跨器官抗纤维化预测与 F4 期 Waddington 协同挽救。** **a**，原代人类髓系细胞跨 NAFLD 外周血单核细胞、健康肝巨噬细胞、F2/F3 期 MASH 瘢痕巨噬细胞（`SAMac`）与 F4 肝硬化 `SAMac` 的表达轨迹（均值 ± 标准误，$\log_2\text{ CPM}$；`GSE136103` 与 `GSE212837`，`n = 19`），展示肝内 `CCR2` 的4.2倍塌缩与 `SPP1`、`TREM2` 的从头涌现。**b**，14项 IIb/III 期抗纤维化临床试验跨器官零样本基准（`人肝 -> 人肺 IPF 与人肾 CKD/DKD`），对比 OmniNiche（`AUROC = 0.9583`，Brier `= 0.1257`）与5种基线模型（`B1-B5`）。**c**，上游单药模拟缓解率随纤维化进展延迟的变化曲线，标示迟滞时间窗悬崖（`t50 = 17.4 h`）、F4 肝硬化单药失效（`0.9%`）及 F4 双室联合挽救（`TGFBR1/PDGFRB + TREM2/SPP1`，`100%` 逃逸，Bliss 超额 `+0.98`）。**d**，二维 Waddington 准势能景观 $U(\text{Myofibroblast}, \text{Parenchyma})$ 揭示被迟滞势垒分隔的消退盆地 `(0.12, 0.12)` 与 F4 肝硬化锁定盆地 `(0.72, 0.72)`。

**Figure 5 | 实体瘤 10x Visium myCAF 空间排斥与治疗前抗 PD-1 盲法应答预测。** **a**，代表性原代人乳腺癌 10x Visium 空间切片（`GSM6177599`，`2,384` 个位点）按 `myCAF` 基质护盾评分着色并叠加 `CD8+` T 细胞富集位点，展示瘤周基质囚禁效应。**b**，12张人源癌组织 Visium 切片（跨 `GSE171351` 与 `GSE203612` 共 `21,796` 个位点）配对空间邻近权重对比 `TAM <-> myCAF` 屏障共定位（均值 `0.666`）与 `Tumor <-> CD8+ T` 核心接触（均值 `0.471`；双侧 Wilcoxon 符号秩检验 `P = 0.0002`）。**c**，严格治疗前 BCC 单细胞活检（`GSE123813`，`n = 8`）、NSCLC 微环境单细胞队列（`GSE207422`，`n = 14`）及汇总严格治疗前基线活检（`n = 32`）的盲法患者级抗 PD-1 应答 AUROC 对比。**d**，19例抗 PD-1 耐药治疗前活检在抗 PD-1 单药（`16.3%` 挽救）与 `抗 PD-1 + myCAF/TGFBR1` 破障联合方案（主要病理缓解阈值 `0.36` 以下达 `99.2%` 挽救）下的边缘 Waddington 准势能曲线 $U(S_{\text{Tumor}}(T))$。

**Extended Data Figure 1 | PRISMA 临床解卷积漏斗、高频解析靶基因与多靶点谱系。** **a**，从 ClinicalTrials.gov II-IV 期原始检索（`41,983`）、`Firewall-2` 非 CNS 过滤（`41,396`）、Open Targets 靶点解析（`29,791`）至确定性非行政结局剥离（`23,967`）的逐级试验数。**b**，按靶点解析非 CNS 临床试验数排名的前10个人类靶基因。**c**，29,791项靶点解析试验每项解析出的人类靶基因数分布（`1` 至 `5+` 个靶点）。

**Extended Data Figure 2 | 原代人类64例三器官纤维化活检队列与六轴单细胞特异性漂移分布。** **a**，跨肝硬化（`GSE136103`，`n = 10`）、肝 MASH（`GSE212837`，`n = 7`）、肺 IPF/ILD（`GSE135893`，`n = 30`）与肾 CKD/DKD（`GSE271188`，`n = 17`）共64例原代人类活检的健康对照与纤维化疾病构成。**b**，6个原代人类单细胞疾病轴在组织重塑中发生主导细胞室切换的临床靶基因百分比。**c**，6个原代人类单细胞轴上靶点特异性漂移（$\Delta\tau$）的箱线图分布（中位数、四分位距与 `1.5x IQR` 须线）。

**Extended Data Figure 3 | 剖面似然参数可辨识性、去混杂 ROC 曲线与配对人类抗 PD-1 T 细胞轨迹。** **a**，核心人源纤维化 SDE 四条有向跨室耦合边的双侧剖面似然可辨识性审计（$\log_{10}(\min \Delta\chi^2 + 1)$），虚线标示95%置信阈值（$\chi^2_1 = 3.841$）。**b**，11,843项去混杂科学结局试验中 OmniNiche（`AUROC = 0.673`）对比静态健康 $\tau_0$（`AUROC = 0.528`）的折外 ROC 曲线。**c**，15例接受抗 PD-1 治疗的皮肤癌纵向配对患者（`GSE123813`；`n = 7` 无应答者，`n = 8` 应答者）治疗后减治疗前 `CD8+` 细胞毒性效应评分变化。

## 参考文献

1. Sun D, Gao W, Hu H, Zhou S. Why 90% of clinical drug development fails and how to improve it? *Acta Pharm Sin B*. 2022;12(7):3049-3062. doi:10.1016/j.apsb.2022.02.002. PMID:35865092.
2. Wong CH, Siah KW, Lo AW. Estimation of clinical trial success rates and related parameters. *Biostatistics*. 2019;20(2):273-286. doi:10.1093/biostatistics/kxx069. PMID:29394327.
3. Nelson MR, Tipney H, Painter JL, et al. The support of human genetic evidence for approved drug indications. *Nat Genet*. 2015;47(8):856-860. doi:10.1038/ng.3314. PMID:26121088.
4. Minikel EV, Painter JL, Dong CC, Nelson MR. Refining the impact of genetic evidence on clinical success. *Nature*. 2024;629(8012):624-629. doi:10.1038/s41586-024-07316-0. PMID:38632401.
5. Buniello A, Suveges D, Cruz-Castillo C, et al. Open Targets Platform: facilitating therapeutic hypotheses building in drug discovery. *Nucleic Acids Res*. 2025;53(D1):D1467-D1475. doi:10.1093/nar/gkae1128. PMID:39657122.
6. Karlsson M, Zhang C, Mear L, et al. A single-cell type transcriptomics map of human tissues. *Sci Adv*. 2021;7(31):eabh2169. doi:10.1126/sciadv.abh2169. PMID:34321199.
7. Kryuchkova-Mostacci N, Robinson-Rechavi M. A benchmark of gene expression tissue-specificity metrics. *Brief Bioinform*. 2017;18(2):205-214. doi:10.1093/bib/bbw008. PMID:26891983.
8. Williams RJ, Tse T, DiPiazza K, Zarin DA. Terminated Trials in the ClinicalTrials.gov Results Database: Evaluation of Availability of Primary Outcome Data and Reasons for Termination. *PLoS One*. 2015;10(5):e0127242. doi:10.1371/journal.pone.0127242. PMID:26011295.
9. Fabre T, Barron AMS, Christensen SM, et al. Identification of a broadly fibrogenic macrophage subset induced by type 3 inflammation. *Sci Immunol*. 2023;8(82):eadd8945. doi:10.1126/sciimmunol.add8945. PMID:37027478.
10. Bunne C, Roohani Y, Rosen Y, et al. How to build the virtual cell with artificial intelligence: Priorities and opportunities. *Cell*. 2024;187(25):7045-7063. doi:10.1016/j.cell.2024.11.015. PMID:39672099.
11. Swanson K, Wu W, Bulaong NL, Pak JE, Zou J. The Virtual Lab of AI agents designs new SARS-CoV-2 nanobodies. *Nature*. 2025;646(8085):716-723. doi:10.1038/s41586-025-09442-9. PMID:40730228.
12. Ramachandran P, Dobie R, Wilson-Kanamori JR, et al. Resolving the fibrotic niche of human liver cirrhosis at single-cell level. *Nature*. 2019;575(7783):512-518. doi:10.1038/s41586-019-1631-3. PMID:31597160.
13. Wang S, Li K, Pickholz E, et al. An autocrine signaling circuit in hepatic stellate cells underlies advanced fibrosis in nonalcoholic steatohepatitis. *Sci Transl Med*. 2023;15(677):eadd3949. doi:10.1126/scitranslmed.add3949. PMID:36599008.
14. Habermann AC, Gutierrez AJ, Bui LT, et al. Single-cell RNA sequencing reveals profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary fibrosis. *Sci Adv*. 2020;6(28):eaba1972. doi:10.1126/sciadv.aba1972. PMID:32832598.
15. Hu J, Zhang L, Xia H, et al. Tumor microenvironment remodeling after neoadjuvant immunotherapy in non-small cell lung cancer revealed by single-cell RNA sequencing. *Genome Med*. 2023;15(1):14. doi:10.1186/s13073-023-01164-9. PMID:36869384.
16. Qiu X, Li J, Bonenfant J, et al. Dynamic changes in human single-cell transcriptional signatures during fatal sepsis. *J Leukoc Biol*. 2021;110(6):1253-1268. doi:10.1002/JLB.5MA0721-825R. PMID:34558746.
17. Anstee QM, Neuschwander-Tetri BA, Wai-Sun Wong V, et al. Cenicriviroc Lacked Efficacy to Treat Liver Fibrosis in Nonalcoholic Steatohepatitis: AURORA Phase III Randomized Study. *Clin Gastroenterol Hepatol*. 2024;22(1):124-134.e1. doi:10.1016/j.cgh.2023.04.003. PMID:37061109.
18. Harrison SA, Bedossa P, Guy CD, et al. A Phase 3, Randomized, Controlled Trial of Resmetirom in NASH with Liver Fibrosis. *N Engl J Med*. 2024;390(6):497-509. doi:10.1056/NEJMoa2309000. PMID:38324483.
19. Richeldi L, du Bois RM, Raghu G, et al. Efficacy and safety of nintedanib in idiopathic pulmonary fibrosis. *N Engl J Med*. 2014;370(22):2071-2082. doi:10.1056/NEJMoa1402584. PMID:24836310.
20. Bakris GL, Agarwal R, Anker SD, et al. Effect of Finerenone on Chronic Kidney Disease Outcomes in Type 2 Diabetes. *N Engl J Med*. 2020;383(23):2219-2229. doi:10.1056/NEJMoa2025845. PMID:33264825.
21. Harrison SA, Wong VW, Okanoue T, et al. Selonsertib for patients with bridging fibrosis or compensated cirrhosis due to NASH: Results from randomized phase III STELLAR trials. *J Hepatol*. 2020;73(1):26-39. doi:10.1016/j.jhep.2020.02.027. PMID:32147362.
22. Raghu G, Richeldi L, Fernandez Perez ER, et al. Pamrevlumab for Idiopathic Pulmonary Fibrosis: The ZEPHYRUS-1 Randomized Clinical Trial. *JAMA*. 2024;332(5):380-389. doi:10.1001/jama.2024.8693. PMID:38762797.
23. de Zeeuw D, Akizawa T, Audhya P, et al. Bardoxolone methyl in type 2 diabetes and stage 4 chronic kidney disease. *N Engl J Med*. 2013;369(26):2492-2503. doi:10.1056/NEJMoa1306033. PMID:24206459.
24. Raue A, Kreutz C, Maiwald T, et al. Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. *Bioinformatics*. 2009;25(15):1923-1929. doi:10.1093/bioinformatics/btp358. PMID:19505944.
25. Wang J, Xu L, Wang E. Potential landscape and flux framework of nonequilibrium networks: robustness, dissipation, and coherence of biochemical oscillations. *Proc Natl Acad Sci U S A*. 2008;105(34):12271-12276. doi:10.1073/pnas.0800579105. PMID:18719111.
26. Mariathasan S, Turley SJ, Nickles D, et al. TGFbeta attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. *Nature*. 2018;554(7693):544-548. doi:10.1038/nature25501. PMID:29443960.
27. Kieffer Y, Hocine HR, Gentric G, et al. Single-Cell Analysis Reveals Fibroblast Clusters Linked to Immunotherapy Resistance in Cancer. *Cancer Discov*. 2020;10(9):1330-1351. doi:10.1158/2159-8290.CD-19-1384. PMID:32434947.
28. Gouin KH 3rd, Ing N, Plummer JT, et al. An N-Cadherin 2 expressing epithelial cell subpopulation predicts response to surgery, chemotherapy and immunotherapy in bladder cancer. *Nat Commun*. 2021;12(1):4906. doi:10.1038/s41467-021-25103-7. PMID:34385456.
29. Barkley D, Moncada R, Pour M, et al. Cancer cell states recur across tumor types and form specific interactions with the tumor microenvironment. *Nat Genet*. 2022;54(8):1192-1201. doi:10.1038/s41588-022-01141-9. PMID:35931863.
30. Yost KE, Satpathy AT, Wells DK, et al. Clonal replacement of tumor-specific T cells following PD-1 blockade. *Nat Med*. 2019;25(8):1251-1259. doi:10.1038/s41591-019-0522-3. PMID:31359002.

# 动态单细胞微环境重塑与生物物理空间位阻排斥决定人类跨器官临床试验成败

---

## 摘要

尽管人类遗传学证据与健康组织表达图谱已广泛用于药物靶点优选，晚期临床试验失败率仍超过 85%。现有靶点评估体系普遍假设基线器官特异性在疾病演进中保持静态不变，且大分子药物可均等穿透重塑后的实质、基质与免疫微区。本研究构建机制驱动的虚拟组织框架 OmniNiche，将跨越 12 个全身与实质器官疾病轴的 41,396 项 II–IV 期独立临床试验（其中 29,791 项解析至明确靶基因）与原代人类单细胞及 10x Visium 空间转录组（来自 162 例供体的约 460,000 个细胞与 21,796 个空间点位）深度耦合。基于子句感知的停药原因解卷积将商业或入组中止（占全部停药试验的 40.3%；基准准确率 88.3%，Cohen's $\kappa = 0.844$）与生物学疗效无效严格剥离，揭示 150-kDa 单克隆抗体发生生物学疗效无效的风险较化药小分子升高 1.59 倍（$P = 7.03 \times 10^{-5}$）。在 6 个核心单细胞疾病轴中，21.5%–37.1% 的临床靶点在疾病进展期发生异位细胞谱系漂移。融合单细胞特异性漂移与模态依赖性细胞外基质位阻排斥的动态微环境可及性指数（`DNAcc`），将临床结局区分度由静态基线特异性的 AUROC 0.5283 提升至 0.6683（配对 DeLong $P = 8.68 \times 10^{-28}$；未见靶点交叉验证 AUROC 0.6540），并在校正多变量混杂后保持独立预测效力（每标准差校正比值比 = 1.43，$P = 9.21 \times 10^{-12}$）。耦合数据编译的随机朗之万求解器后，跨人类肝、肺、肾纤维化的零样本迁移准确预测 14 项 IIb/III 期抗纤维化试验（AUROC = 0.9583，精确置换 $P = 0.0028$），并在零参数重拟合下分层严格治疗前实体瘤抗 PD-1 活检队列（合并 $n = 32$ AUROC = 0.7551，$P = 0.0148$）。

**关键词**：临床试验损耗；动态微环境重塑；单细胞转录组；空间转录组；器官纤维化；肿瘤免疫治疗；生物物理系统药理学

---

## 引言

临床新药研发长期面临严峻的晚期损耗，进入人体试验的候选药物近 90% 未能获得监管批准 [1, 2]。过去十年间，人类遗传学证据已成为计算靶点优选的核心基石：当因果变异-基因关联支持靶点与适应症之间的机制联系时，药物获批概率约可翻倍 [3, 4]。Open Targets 等大规模知识平台与 TxGNN 等图基础模型进一步拓展了跨数千种表型的零样本靶点发现与老药新用预测 [5, 6]。与此同时，研究者广泛利用人类蛋白质图谱（Human Protein Atlas）单细胞参考集等健康人体图谱，优先筛选基线组织特异性（$\tau_0$）较高的靶点以降低脱靶器官毒性 [7, 8]。然而，即便同时具备胚系遗传学支持与健康组织高选择性，关键性 IIb 期与 III 期临床失败在纤维化、炎症、代谢及肿瘤疾病中依然普遍存在。

我们提出科学假说：这一转化鸿沟源于现有体系仅在静态、健康的解剖快照上评估靶点，而忽略了病变人类器官内部涌现的两重动态组织屏障。其一，慢性实质损伤、黏膜溃疡与恶性增殖会触发剧烈的多细胞重塑，包括促纤维化 $TREM2^+SPP1^+$ 瘢痕相关巨噬细胞（SAMac）浸润以及 $COL1A1^+POSTN^+$ 肌成纤维细胞扩张 [9, 10]。在此重塑过程中，健康组织中高度特异的靶基因可能在目标谱系中坍塌（$\Delta\tau < 0$），或异位漂移至旁观基质微区。其二，细胞外基质（ECM）交联与瘤周癌相关成纤维细胞（`myCAF`）包膜构筑了物理输运壁垒，通过空间位阻将约 150-kDa 的单克隆抗体排斥在靶细胞之外，却仍允许小分子与低分子量多肽扩散穿透 [11, 12]。

要在临床转化全景中检验上述机制，必须首先攻克公共临床试验注册库中普遍存在的结局标签混杂：商业管线调整与入组缓慢常与生物学靶点失效混为一谈。本研究建立端到端转化计算框架 OmniNiche，贯通 41,396 项 II–IV 期临床试验的子句感知结局解卷积、跨越 6 个人类器官重塑轴的实测单细胞靶点漂移量化，以及空间约束下的随机朗之万多细胞动力学引擎。OmniNiche 不依赖针对临床标签的黑盒监督拟合，而是直接从原代人类单细胞与空间活检数据中编译组织微区互作语法与物理屏障约束，确立跨人类器官靶点可及性与治疗迟滞锁定的定量生物物理定律。

---

## 结果

### 跨越 12 个全身与实质器官疾病轴的 41,396 项临床试验宏观解卷积

为构建受间质基质重塑调控的外周组织无混杂临床转化基准，我们检索了 ClinicalTrials.gov API v2 中横跨 12 个全身与实质器官疾病轴的 II 期、III 期与 IV 期干预性临床研究（Fig. 1a–d 及 Supplementary Table S1）。按唯一 NCT 编号去重后共获 41,396 项临床试验，其中实体瘤（$n = 15,717$）、代谢与内分泌疾病（$n = 4,915$）、血液肿瘤（$n = 3,695$）、肺纤维化与急性呼吸窘迫（$n = 3,525$）、心血管疾病（$n = 3,381$）及自身免疫与黏膜炎症（$n = 2,529$）居前（Fig. 2a,b 及 Extended Data Fig. 1a）。将干预药物名称及同义词与 ChEMBL v33 作用机制库和 Open Targets 药物-靶点注释交叉映射后，29,791 项试验（72.0%）被解析至 851 个人类唯一靶基因，涵盖化药小分子（$n = 18,128$）、单克隆抗体（$n = 7,308$）以及多肽或蛋白生物制品（$n = 4,355$）（Extended Data Fig. 1b,c）。

注册库原始终止状态之所以掩盖生物学机制，是因为战略重组、申办方并购、COVID-19 运营受阻或入组缓慢导致的停药与疗效无效及剂量限制性毒性被归入同一终止标签 [13]。我们开发了子句感知的自然语言停药审计器，通过切分多子句 `whyStopped` 叙述并识别安全性与无效性否定短语（如“与安全性或缺乏疗效无关”），将停药试验精细解卷积为明确生物学疗效无效、明确安全性或毒性终止、行政商业或入组运营中止以及未报告原因停药四类（Supplementary Table S2）。在分层抽取的 240 项人工逐句裁定金标准基准（每类 60 项）上，审计器总体准确率达 88.33%，Cohen's $\kappa = 0.8444$，其中疗效无效与安全性终止的 F1 分数分别达 0.9474 与 0.9194（Fig. 2f）。

将该验证后的解卷积流程应用于 41,396 项注册研究中的全部 9,406 项停药试验发现，行政与运营中止占全部停药事件的 40.3%（$n = 3,793$），数量超过确证生物学疗效与安全性失败（$n = 1,090$）的三倍以上（Fig. 2a）。2005 至 2025 年纵向轨迹显示，行政运营中止在停药试验中的占比从 2005 年的 41.7% 攀升至 2020–2025 年疫情与管线整合期的 60% 以上，同期单克隆抗体试验占比亦从 19.6% 翻倍至峰值 40.1%（Fig. 2c）。受试者入组规模分布进一步印证了这一剥离的有效性：行政中止早期停留在中位 16 例患者，安全性失败中止于中位 24 例，而生物学疗效无效试验则推进至期中或主要终点读出，中位入组达 73 例（累计暴露 13 万例患者），接近已完成试验的中位 84 例（累计暴露 977 万例）（Fig. 2e）。

将治疗模态与解卷积后的停药原因交叉检验，揭示了显著的大分子疗效损耗惩罚（Fig. 2d）。相较于小分子药物，约 150-kDa 的 IgG 单克隆抗体发生明确生物学疗效无效的比值比升高至 1.59 倍（OR = 1.59，95% 置信区间（CI）1.26–1.99，Fisher 精确检验 $P = 7.03 \times 10^{-5}$），而明确安全性终止并未显著增加（OR = 1.14，95% CI 0.85–1.51，$P = 0.37$）。相反，更易穿透致密间质基质的低分子量多肽与蛋白药物（约 2–15 kDa）发生疗效无效的风险显著低于小分子（OR = 0.59，95% CI 0.39–0.89，$P = 0.0123$）。在 29,791 项靶点明确试验中剔除行政中止与未报告原因停药后，我们确立了涵盖 386 个靶点的 11,843 项去混杂科学基准队列（11,072 项已完成 III/IV 期成功试验对比 771 项靶点明确的确证科学失败），以及涵盖 564 个靶点的 23,967 项 II 期向 III/IV 期临床成熟队列。

### 单细胞靶点特异性漂移与动态微环境可及性（`DNAcc`）定律

为检验疾病诱导的转录重塑能否解释基线组织选择性的预测失效，我们整合了跨越 6 个核心疾病轴的原代人类单细胞与单细胞核转录组队列（来自 150 例人类供体的约 460,000 个细胞，并配对 12 张空间 Visium 肿瘤切片）。上述队列覆盖肝脏代谢功能障碍相关脂肪性肝炎（MASH）与肝硬化（`GSE136103` 与 `GSE212837`）、特发性肺纤维化（IPF，`GSE135893`）、慢性肾脏病（`GSE271188`）、实体癌（`GSE207422` 与 `GSE123813`）、黏膜溃疡性结肠炎（`GSE250487`）及纵向脓毒症（`GSE167363`）[14–18]（Extended Data Fig. 2a 及 Supplementary Table S3）。针对全部 851 个临床靶基因，我们在人类蛋白质图谱 154 种细胞类型上量化了健康基线特异性（$\tau_0$），并在原代病变队列中测定了早期病变特异性（$\tau_{\mathrm{early}}$）、晚期病变特异性（$\tau_{\mathrm{late}}$）、动态特异性漂移（$\Delta\tau = \tau_{\mathrm{late}} - \tau_{\mathrm{early}}$）以及实质、基质与免疫谱系间的异位微区切换。

疾病演进在全部 6 个人类单细胞疾病轴中引发了广泛的靶点谱系重构（Fig. 3a）。在有表达的临床靶基因中，21.5%（肺纤维化）至 37.1%（全身性脓毒症）在早期与晚期病变之间发生了主导表达细胞谱系的切换，其中肝脏 MASH 为 27.8%，黏膜自身免疫病为 28.2%，慢性肾病为 30.2%（Extended Data Fig. 2b,c）。高 $\tau_0$ 靶点在 $(\tau_0, \Delta\tau)$ 重塑相空间中分化为两条截然相反的临床轨迹（Fig. 3f）：在晚期病变中维持或增强可及微区谱系选择性（$\Delta\tau \ge 0$）的靶点（如肝细胞 `THRB` 与 `PCSK9`、黏膜免疫 `IL23A` 及血管基质 `PDGFRB`）取得了极高的晚期临床成熟率；而健康基线 $\tau_0$ 极高却在晚期遭遇谱系侵蚀或异位基质劫持的靶点（包括 `CCR2`、`LOXL2`、`TGFB1` 与 `MAP3K5`）则在关键性临床试验中接连折戟。

我们将上述机制凝练为六步正交消融层级，最终构建出融合健康基线选择性、重要器官旁观毒性风险、早期病变特异性、病程特异性漂移（$\Delta\tau$）、异位谱系切换及模态-基质位阻相容性的动态微环境可及性指数（`DNAcc`）（Table 1 及 Fig. 3b）。在 11,843 项去混杂科学基准队列中，单独使用静态健康 $\tau_0$ 的预测表现接近随机水平（五折交叉验证 AUROC = 0.5283，95% CI 0.5092–0.5483）。依次引入重要器官旁观风险（AUROC = 0.5715）、早期疾病特异性（0.6582）、病程漂移与异位切换以及生物物理位阻相容性后，区分度跃升至 AUROC = 0.6683（95% CI 0.6477–0.6884；$\Delta\mathrm{AUROC} = +0.1400$，配对 DeLong $P = 8.68 \times 10^{-28}$）[19]。在单细胞实测覆盖的 6 个核心疾病轴（$n = 5,815$ 项试验）上，`DNAcc` 的 AUROC 达 0.6796，显著优于静态 $\tau_0$ 的 0.5252（$\Delta\mathrm{AUROC} = +0.1544$，DeLong $P = 1.06 \times 10^{-22}$）。在 23,967 项 II 期向 III/IV 期成熟队列上，全部六个消融步均呈严格单调提升，AUROC 由 0.5405 增至 0.6737（$\Delta\mathrm{AUROC} = +0.1332$，DeLong $P = 3.99 \times 10^{-182}$；Brier 误差由 0.2474 降至 0.2277）。

为排除靶基因记忆伪影、多靶点药理学混杂、模态混杂与同靶点试验聚类相关性，我们对 `DNAcc` 实施了五重正交验证。在严格按靶基因分组隔离的五折 `GroupKFold` 交叉验证下，`DNAcc` 在全部 12 个疾病轴（AUROC = 0.6540 对比静态 $\tau_0$ 的 0.4617，DeLong $P = 3.30 \times 10^{-40}$）与 6 个单细胞疾病轴（AUROC = 0.6622 对比 0.4438，$P = 6.32 \times 10^{-38}$）上均稳健泛化至全新靶点（Fig. 3b 及 Extended Data Fig. 1d）。仅纳入单靶点单药疗法（科学队列 $n = 6,194$ AUROC = 0.6442 对比 0.5221，$P = 1.63 \times 10^{-6}$），或严格限定于小分子内部（$n = 7,369$，AUROC = 0.6554 对比 0.5601，$P = 5.97 \times 10^{-11}$）及单克隆抗体内部（$n = 2,357$，AUROC = 0.6958 对比 0.5515，$P = 4.36 \times 10^{-10}$），均确证了不依赖模态差异的稳健增益（Extended Data Fig. 1c）。在控制疾病轴、药物模态、启动年代、申办方类别与多靶点指标并采用 Liang–Zeger 靶基因聚类稳健三明治标准误的多变量固定效应逻辑回归中，`DNAcc` 在科学队列（校正 OR = 1.43，Wald 95% CI 1.29–1.59，$P = 9.21 \times 10^{-12}$；靶点聚类 95% CI 1.23–1.67，$P_{\mathrm{cluster}} = 2.58 \times 10^{-6}$）与成熟队列（校正 OR = 1.29，$P = 3.83 \times 10^{-44}$；$P_{\mathrm{cluster}} = 8.85 \times 10^{-6}$）中均保持高度显著；反之，静态 $\tau_0$ 在校正靶点聚类后即便单变量亦无显著性（$P_{\mathrm{cluster}} = 0.49$），与 `DNAcc` 同模竞争时效应完全归零（校正 OR = 1.01，$P = 0.81$）（Fig. 3c 及 Extended Data Fig. 3b）。分轴分析在 5 个实体与黏膜屏障轴中取得一致增益，而血管内 PBMC 脓毒症（$n = 144$）呈零区分（AUROC = 0.4793 对比 0.4738，$P = 0.92$），界定了间质微环境模型的天然边界条件（Extended Data Fig. 2d）。

按 `DNAcc` 四分位分层的临床实测 II 期向 III/IV 期成熟率由最低四分位（Q1）的 25.53% 跃升 2.56 倍至最高四分位（Q4）的 65.37%（$n = 23,967$），而静态 $\tau_0$ 仅呈 1.20 倍的平坦梯度（44.58% 至 53.41%）（Fig. 3d）。涵盖将全部 3,118 项靶点明确但未报告原因的停药均归为生物学失败的最坏情境（$n = 14,961$，DeLong $P = 2.09 \times 10^{-83}$）、逆删失概率加权（$P = 1.76 \times 10^{-29}$）及仅限已公布临床结果试验（$n = 5,757$，$P = 1.79 \times 10^{-19}$）在内的五种缺失机制敏感性分析均保持结论不变（Fig. 3g）。折外预测概率在全部十个十分位区间内高度校准，并在临床可操作的 III 期准入阈值区间内持续提供更优的决策曲线分析（DCA）临床净获益 [20]（Fig. 3e 及 Extended Data Fig. 3c,d）。在 23,967 项试验的投资组合富集模拟中，在进入 III 期前筛除 `DNAcc` 最低 25% 的候选项目，可将 III/IV 期成功率由 46.20% 提升至 53.08%，避免 4,462 项失败试验（约 123 万例无效受试者暴露），同时保留 86.18% 的最终成功药物。

[Table 1 near here]

### 人类泛器官纤维化临床试验零样本解卷积与迟滞锁定机制

注册库宏观评分能够甄别脆弱靶点，但晚期临床成败还取决于病变分期与非线性多细胞正反馈。器官纤维化正是这一挑战的典型代表：尽管临床前证据充分，数十项 IIb 期与 III 期抗纤维化试验仍在 MASH 肝硬化、IPF 与纤维化性慢性肾病中折戟。整合来自人类肝脏（`GSE136103` 与 `GSE212837`，$n = 17$）、肺脏（`GSE135893`，$n = 30$）和肾脏（`GSE271188`，$n = 17$）的 64 例原代单细胞与单细胞核供体活检数据，确证了跨器官保守的三微区纤维化三联体：受损上皮或肝实质细胞、促纤维化 $TREM2^+SPP1^+$ SAMac 以及产胶原 $COL1A1^+POSTN^+$ 肌成纤维细胞（Fig. 4a 及 Supplementary Table S4）。

追踪从循环单核细胞（`GSE136103` 外周血基线）经健康人肝进入 F2/F3 期 MASH 与 F4 期肝硬化（`GSE136103` 与 `GSE212837`）的人类髓系谱系轨迹，解开了 `CCR2/CCR5` 拮抗剂 cenicriviroc 的临床悖论：尽管其健康基线特异性极高（`CCR2` $\tau_0 = 0.9513$；`CCR2/CCR5` 均值 $0.9531$），却在 III 期 AURORA 试验（`NCT03028740`）中宣告无效 [21]（Fig. 4b）。虽然外周血单核细胞高表达 `CCR2`（平均伪批量对数表达量 6.19），但肝内巨噬细胞 `CCR2` 在纤维化锁定前的 F2/F3 期 MASH 即暴跌 76% 至 1.47，取而代之的是不依赖外周单核细胞新募集即可自主维持肌成纤维细胞活化的组织驻留态 $SPP1^+$（外周血单核细胞为 0.00、健康肝巨噬细胞为 5.46，F4 期肝硬化升至 8.31）与 $TREM2^+$（由 0.39 升至 5.77）SAMac 群体。

我们仅利用 17 例人类肝脏供体数据编译三微区随机朗之万动力学语法（含状态有界 Hill 激活、环境承载力饱和及岭回归基因调控网络扰动），并在跨器官测试前审计了结构参数可辨识性。伴随干扰速率尺度重优化的连续 25 点双侧剖面似然扫描显示，全部四条有向细胞间互作边均形成平滑抛物线型 $\Delta\chi^2$ 势阱，双侧均跨越 95% 置信阈值（$\Delta\chi^2 = 3.841$）（Extended Data Fig. 3a）。在不重拟合网络拓扑或 Hill 动力学参数的前提下，我们将该肝脏编译的朗之万引擎零样本迁移至肺 IPF 与肾 CKD，对 14 项标志性 IIb/III 期抗纤维化临床试验（6 项监管获批或 III 期阳性试验对比 8 项确证疗效失败试验）进行盲测评估 [22–30]（Fig. 4c）。

OmniNiche 在三个器官间准确区分临床成功与失败试验，取得 AUROC = 0.9583（95% CI 0.8125–1.0000；48 对正负比较中准确排序 46 对）、AUPRC = 0.9583、Brier 误差 = 0.1257 及 10,000 次精确置换 $P = 0.0028$，全面超越静态健康 $\tau_0$（AUROC = 0.5000，$P = 1.00$）、无多细胞动力学的单细胞 GRN 保留度（0.5000）、多靶点计数（0.6667）、药物模态先验（0.7500）及静态配体-受体中心度（0.9167，Brier = 0.1375）（Table 2 及 Fig. 4d）。该机制引擎准确区分了 F2/F3 期 MASH 中两条截然不同的临床获益路径：resmetirom（`MAESTRO-NASH`，`NCT03900429`，模拟概率 0.54）虽健康基线特异性中等（$\tau_0 = 0.4052$），但凭借 MASH 进展期肝实质 `THRB` 选择性锐化（$\Delta\tau = +0.0829$）取得临床疗效 [22]；而 4.11-kDa 酰化 GLP-1 多肽司美格鲁肽（semaglutide，`ESSENCE`，`NCT04822181`，模拟概率 1.00）并非直接作用于星状细胞，而是通过卸载上游脂毒性与心代谢炎症负荷实现逆转 [23]。相反，因交联胶原空间位阻排斥，约 150-kDa 抗 `LOXL2` 与抗 `CTGF` 单抗在晚期瘢痕中（simtuzumab F4 肝硬化 0.09、IPF 0.00；pamrevlumab IPF 0.00）及靶向 `MAP3K5`/ASK1 的 selonsertib 在 F3 与 F4 试验中均获近零缓解概率 [28–30]。仅两项前肝硬化失败试验评分高于 resmetirom（解释全部 2/48 对排序误差）：反映前肝硬化中等基质通透性的 F3 期 simtuzumab（0.98），以及缺乏泛 PPAR 阳性药物 lanifibranor（1.00）所具 `PPARG` 星状细胞失活臂的双 `PPARA/PPARD` 激动剂 elafibranor（`RESOLVE-IT`，0.74）[24]。

为证实跨器官可迁移性源于人类保守纤维化三联体的内禀动力学特性而非肝脏单向训练伪影，我们在肝、肺、肾之间开展了双向 $3 \times 3$ 单器官编译与留一器官交叉验证（LOOO）（Fig. 4d）。仅由人类肺 IPF（$n = 30$ 例供体）编译的动力学语法在留出肝、肾试验上取得 AUROC = 0.9583（全部 14 项试验 AUROC = 0.9792）；仅由人类肾 CKD（$n = 17$ 例供体）编译的语法在留出肝、肺试验上取得 AUROC = 1.0000；全部三组 LOOO 迁移在留出器官上的 AUROC 均达 0.9444 至 1.0000。

构建二维 Waddington 准势能面 $U(S_{\mathrm{Myofib}}, S_{\mathrm{Paren}})$ 与连续解剖病程分岔曲线揭示了肝硬化一旦锁定后单药疗法全面失效的动力学根源（Fig. 4e,f）。早期 F1/F2 病变处于可逆单稳态区间，卸载实质损伤驱动（`THRB` 或 `GLP1R`）即可消退瘢痕信号。一旦越过临界进展阈值（$s_{\mathrm{prog}} \approx 0.48$，对应无量纲病程时间中的单药治疗窗口断崖 $t_{50} = 17.43\text{ h}$），`SAMac` 与 `Myofibroblast` 之间的互惠正反馈即构筑起深邃的双稳态迟滞吸引盆（“锁定态 F4”），将任何单微区干预的缓解率禁锢在 5% 以下。在锁定态 F4 肝硬化中扫描对称 $6 \times 6$ 临床靶点二维 Bliss 协同矩阵表明，联用两种同属实质微区的药物（`THRB + GLP1R`，Bliss 超额协同 $+0.12$）或将活性药物与受位阻排斥的 `LOXL2` 单抗或靶点坍塌的 `CCR2` 拮抗剂联用（Bliss 超额协同 $0.00$ 至 $+0.02$）均无法逃离吸引盆；唯有同时卸载实质或 SAMac 驱动并抑制肌成纤维细胞 `TGFBR1/PDGFRB` 的正交跨微区联合方案，方能击穿鞍点势垒并恢复 98.6%–100% 的模拟缓解率（Bliss 超额协同 $+0.99$ 至 $+1.00$）（Fig. 4g）。

[Table 2 near here]

### 空间 Visium 基质排斥与零重拟合治疗前免疫治疗响应预测

我们进一步追问，类似的物理-动力学双屏障是否同样主导人类实体癌对抗 PD-1 免疫检查点阻断（ICB）的原发耐药：在实体肿瘤中，约 150-kDa 的检查点单抗与细胞毒性 $CD8^+$ T 细胞均须穿越瘤周 `myCAF` 与 $SPP1^+$ 肿瘤相关巨噬细胞（TAM）环带方能接触恶性细胞 [11, 31]。通过流式解析横跨膀胱癌、乳腺癌、卵巢癌、肝癌与胃肠间质瘤的 12 张原代人类癌组织 10x Visium 空间切片共 21,796 个空间点位（`GSE171351` 与 `GSE203612`）[32, 33]，我们量化了物理 $k = 6$ 最近邻双变量 Moran 空间邻近度以及恶性核心边界外 0 至 600 µm 的连续径向特征衰减曲线（Fig. 5a,b、Extended Data Fig. 4a,b 及 Supplementary Table S5）。

在全部 12 张人类肿瘤空间切片中，`myCAF`–`TAM` 空间共定位强度在 12/12 张切片中均超过 `Tumor`–$CD8^+$ T 细胞直接接触（共识邻近度 0.666 对比 0.471；12/12 张切片基质排斥指数 $\mathrm{SEI} > 0$，Wilcoxon 符号秩检验 $P = 2.44 \times 10^{-4}$）（Fig. 5b 及 Extended Data Fig. 4a,b）。跨越 21,796 个点位的连续径向剖析揭示了同心圆式物理屏障构筑 [34]：恶性核心特征在 0 µm 处达峰（$z = +1.26 \pm 0.04$）并向外陡峭衰减；$CD8^+$ T 细胞在肿瘤核心内部遭受物理排空（0 µm 处 $z = -0.06 \pm 0.04$）并滞留于 200 µm 边缘带（$z = +0.09 \pm 0.04$）；而 `myCAF` 基质密度则向外单调攀升，在 400–500 µm 处形成致密包绕的瘤周壁垒（$z = +0.67 \pm 0.31$）。典型空间切片的直接对比清晰展示了免疫排斥型乳腺癌基质盾牌（$\mathrm{SEI} = +0.213$）与相对可浸润肝细胞癌微环境（$\mathrm{SEI} = +0.070$）之间的物理门控差异（Fig. 5a）。

我们将这一由 12 张切片实测得出的共识 Visium 邻近度矩阵直接嵌入五微区空间朗之万实体瘤语法（恶性核心、$SPP1^+$ TAM、`myCAF` 屏障、$CD8^+$ 耗竭态与 $CD8^+$ 效应储备态），全程不拟合任何临床响应标签。利用该固定空间语法，我们在两个独立人类抗 PD-1 临床队列上开展盲法响应预测：配对治疗前后基底细胞癌（BCC）全微环境 scRNA-seq（`GSE123813`，$n = 8$ 例含非 T 细胞微环境的严格治疗前患者；$n = 15$ 例配对纵向活检）[35]，以及非小细胞肺癌（NSCLC）治疗前活检队列（`GSE207422`，$n = 24$ 例严格治疗前 bulk RNA-seq 活检与 $n = 14$ 例肿瘤微环境 scRNA-seq 验证样本）[36]。对 `GSE123813` 实验方案的技术审计表明，独立的 `CD45+CD3+` 流式细胞分选（FACS）富集管在单患者测序细胞中占比高达 79.3%，将未校正全细胞池中的 `myCAF` 比例人为稀释达 4.8 倍（中位 1.7 倍）；将基质状态估计严格限定于非 `CD3+` 分选微环境细胞彻底消除了这一实验稀释畸变（Extended Data Fig. 4c）。

在完全零跨癌种参数重拟合的条件下，凡保留单细胞肿瘤微环境分辨率的队列，Visium 约束的空间朗之万引擎均在全部六种方法中位列第一：治疗前 BCC scRNA-seq 队列取得 AUROC = 0.8750（95% CI 0.6250–1.0000，对比 `CD274`/PD-L1 的 0.4375 与静态 $\tau_0$ 的 0.5625），NSCLC scRNA-seq 验证队列取得 AUROC = 0.8500（95% CI 0.6000–1.0000，Brier = 0.1751，对比 PD-L1 的 0.3250 与静态 $\tau_0$ 的 0.7500）（Table 2 及 Fig. 5d）。在经组织裂解均质化、抹平了同心圆 `myCAF` 环带几何结构的治疗前 NSCLC bulk 活检队列（$n = 24$）中，非空间 $\tau_0$（0.8000）与 T 细胞 GRN 评分（0.7926）排序 AUROC 略高，而 OmniNiche 仍显著优于随机（AUROC = 0.7407，95% CI 0.5222–0.9333，精确置换 $P = 0.0466$）与 PD-L1（0.6593），并取得全部方法中最优 Brier 校准误差（0.2250）。在跨越两个癌种与两种测序平台的 $n = 32$ 例严格治疗前人类肿瘤活检合并队列中，OmniNiche 在全部评估指标上均居首位（AUROC = 0.7551，95% CI 0.5769–0.9190，精确置换 $P = 0.0148$；AUPRC = 0.6385；Brier = 0.2366），显著超越 `CD274`/PD-L1（0.5344，$P = 0.76$）、静态健康 $\tau_0$（0.6194）、T 细胞 GRN 评分（0.6113）、静态配体-受体中心度（0.4656）及线性微区均值（0.6235）（Fig. 5c,d 及 Extended Data Fig. 4d）。

`GSE123813` 队列中 $n = 15$ 例患者配对治疗前后的纵向 scRNA-seq 轨迹进一步印证了模型预测的动力学机制：应答者基线基质处于未锁定状态，允许 PD-1 阻断后净 $CD8^+$ 细胞毒性效应群体显著扩增；无应答者则在持续存在的 `myCAF` 屏障遮蔽下出现治疗后效应 T 细胞进行性耗损（Fig. 5e）。计算二维实体瘤 Waddington 准势能面 $U(S_{\mathrm{myCAF}}, S_{\mathrm{Tumor}})$ 并对全部 $n = 19$ 例抗 PD-1 原发耐药治疗前患者实施双靶点破障模拟（`anti-PD-1 + myCAF/TGFBR1` 联合抑制），确定性漂移矢量场跨越基质分界线发生定向偏转，平均模拟缓解概率由单药的 16.3% 跃升至 99.2%（Wilcoxon 符号秩检验 $P = 1.91 \times 10^{-6}$）（Fig. 5f,g）。

---

## 讨论

临床药物靶点优选长期建立在静态的人体组织生物学视角之上，将胚系遗传关联与健康基线器官选择性视为固定不变且与药物模态无关的治疗窗预测指标 [3–5]。通过将 41,396 项临床试验与跨越 162 例人类供体的原代单细胞及空间转录组相统一，本研究证实外周实质与黏膜器官临床试验的晚期损耗受动态单细胞微环境重塑与大分子空间位阻排斥共同支配。一旦将占全部停药试验 40.3% 的行政商业与入组中止从生物学失败中剥离，健康基线特异性（$\tau_0$）对晚期科学失败几乎不具备保护力（AUROC = 0.5283），且在纳入动态重塑与物理屏障后彻底丧失独立解释力（联合多变量 OR = 1.01，$P = 0.81$）。

宏观注册库与器官尺度分析共同揭示了三条核心转化机制原理。第一，五分之一至三分之一以上的临床靶基因在人类组织病变演进中切换其主导表达细胞谱系，这解释了为何依据健康组织图谱或外周血细胞筛选的靶点（如 MASH 肝硬化中的 `CCR2`）会在重塑病灶内部丧失靶向结合 [21]。第二，治疗模态与间质 ECM 构筑发生强烈交互：约 150-kDa 单克隆抗体的生物学疗效无效风险较小分子升高 1.59 倍，并在致密交联的纤维化瘢痕中一致失效（如 `LOXL2` 与 `CTGF` 单抗）[28, 29]，而小分子与紧凑型酰化多肽（如司美格鲁肽）则可穿透或绕过致密纤维化基质 [22, 23]。第三，晚期实质纤维化与免疫排斥型实体癌均运行于由保守髓系-间充质正反馈（`SPP1+` SAMac/TAM 与 `POSTN+` 肌成纤维细胞/`myCAF`）稳定的双稳态 Waddington 吸引盆中 [9–11, 37, 38]。一旦病程跨越单药治疗窗口断崖，单一通路阻断即不足以逆转稳态，从而为正交跨微区破障联合疗法提供了定量生物物理依据。

OmniNiche 与近期基于常规血液学指标或泛癌 bulk 转录组的免疫治疗临床机器学习模型（如 SCORPIO、COMPASS 与虚拟空间蛋白质组学）[39–41] 以及新兴虚拟细胞计划 [42, 43] 形成高度互补。监督式临床分类器通常依赖同治疗方案下的大规模带标签患者训练队列，而将实测 Visium 空间邻近度与单细胞微区漂移编译入低维随机朗之万动力学语法，则实现了零标签跨器官抗纤维化临床试验预测（AUROC = 0.9583）以及跨癌种、跨测序平台的零重拟合治疗前免疫治疗应答分层（$n = 32$ 合并队列 AUROC = 0.7551）。

本研究仍存在若干边界条件有待后续工作拓展。其一，ClinicalTrials.gov 中约 15.2% 的靶点明确试验在终止时未提供具体叙述原因；尽管包含全部未报告停药归为失败的最坏情境（$n = 14,961$，$P = 2.09 \times 10^{-83}$）与仅限公布结果试验（$n = 5,757$，$P = 1.79 \times 10^{-19}$）在内的五重敏感性分析均证实效应稳健，推动公共注册库实行结构化失败原因强制披露将进一步提升效应量精度。其二，本框架的间质 ECM 位阻定律专为实体与黏膜组织设计，因此在血管内循环 PBMC 脓毒症队列（$n = 144$，AUROC = 0.4793）中呈现预期的零区分表现，提示血管内全身综合征需耦合血流动力学而非间质扩散模型。其三，尽管治疗前抗 PD-1 验证严格剔除了治疗中手术切除样本以杜绝时间泄漏，在临床部署前仍需依托正在开展的 `TGFBR1`/基质破障联合临床试验开展前瞻性多中心配对空间多组学活检验证。将虚拟组织动力学锚定于实测人类单细胞漂移与空间屏障物理之上，可在关键性人体暴露前为靶点优选、模态工程与联合试验设计提供可操作的定量蓝图。

---

## 方法

### ClinicalTrials.gov 注册库摄取、疾病轴范围界定与靶基因解析

本研究通过 ClinicalTrials.gov API v2（`https://clinicaltrials.gov/api/v2/studies`）检索横跨 12 个全身与实质器官疾病轴的 II 期、III 期与 IV 期干预性临床试验记录：实体瘤、代谢与内分泌、血液肿瘤、肺部疾病（IPF 与 ARDS）、心血管疾病、自身免疫病（IBD 与风湿免疫）、肾脏疾病（CKD 与 AKI）、皮肤免疫、肌肉骨骼、眼科与视网膜、肝脏疾病（MASH 与肝硬化）以及感染与脓毒症。为保持间质输运生物物理机制的一致性，注册库检索聚焦于受间质细胞外基质重塑调控的外周全身与实质器官病变，排除了 587 项受特化血脑屏障神经血管转运支配的原发中枢神经系统记录（共检索 41,983 条原始记录，去重后保留 41,396 项独立 NCT 研究）。

干预药物名称、通用名与同义词经标准化清洗后，与 ChEMBL v33 审定药物作用机制库（`https://www.ebi.ac.uk/chembl/api/data`）及 Open Targets Platform v24.09 药物-靶点与适应症证据库（`https://api.platform.opentargets.org/api/v4/graphql`）[5] 进行交叉映射，将 29,791 项独立试验解析至 851 个 HGNC 人类靶基因符号与标准化治疗模态（`small_molecule`、`antibody`、`peptide` 及其他生物制品类别；其中 4.11-kDa 酰化 GLP-1 类似物司美格鲁肽按分子物理性质精确归类为 `peptide`）。

### 子句感知停药原因解卷积与 240 项试验审计基准验证

针对已停药试验（`TERMINATED`、`WITHDRAWN` 或 `SUSPENDED`；41,396 项注册研究中合计 $n = 9,406$ 项），采用确定性子句感知自然语言解析器综合分析 `overallStatus`、`whyStopped` 与已公布结果元数据。在词法匹配前，先沿子句边界（`;`、`.`、`but`、`however`、`although`、`unrelated to`）切分复合叙述。含否定安全性或无效性表述的子句（匹配 `not due to safety`、`no safety concerns`、`not related to efficacy`、`absence of toxicity` 等模式）被严格阻断触发科学失败标签，从而将“因战略管线调整终止、与安全性问题无关”等复合叙述准确路由至 `Administrative_Commercial_Enrollment`。无否定修饰的明确缺乏疗效或期中无效子句归入 `Scientific_Efficacy_Futility`，无否定修饰的不良事件或临床暂停子句归入 `Scientific_Safety_Toxicity`，运营、入组、申办方变更或 COVID-19 子句归入 `Administrative_Commercial_Enrollment`，空白或非特异叙述归入 `Other_Unreported_Stop`。分类准确率、各类别精确率、召回率、F1 分数与 Cohen's $\kappa$ 在分层抽取的 240 项人工裁定参考集（每类 60 项；Supplementary Table S2）上完成盲法核验。

### 原代人类单细胞与 10x Visium 空间转录组数据处理

全部单细胞、单细胞核与空间转录组数据集均由取自 NCBI GEO 与人类蛋白质图谱的 100% 真实原代人类供体矩阵处理生成，零合成插补；通过流式压缩稀疏矩阵读取与供体级伪批量聚合，Python 进程峰值内存全程控制在 1.2 GB 以下：

1. **健康基线参考集**：Human Protein Atlas v23 单细胞转录组跨越 154 种人类细胞类型的共识归一化表达谱 [7]。
2. **肝脏 MASH 与肝硬化**：`GSE136103` 原代人类健康与肝硬化肝脏及配对循环外周血单核细胞 scRNA-seq（$n = 10$ 例肝脏供体与 $4$ 例外周血单核细胞基线）[9] 与 `GSE212837` 原代人类对照与 F2/F3 期 NASH snRNA-seq（$n = 7$ 例供体）[14]，用于髓系 `CCR2`/`SPP1`/`TREM2` 谱系轨迹分析。
3. **肺纤维化**：`GSE135893` 原代人类对照、中间态间质性肺病与 F3/F4 期 IPF 肺组织 scRNA-seq（$n = 30$ 例供体）[15]。
4. **肾纤维化**：`GSE271188` 原代人类对照、微小病变基线与纤维化性 CKD/DKD/IgAN/MN 肾组织 snRNA-seq（$n = 17$ 例供体）[16]。
5. **实体瘤单细胞、治疗前 bulk 活检与 10x Visium 空间队列**：`GSE123813` 配对抗 PD-1 治疗前后人类 BCC scRNA-seq（$n = 8$ 例治疗前全微环境患者；$n = 15$ 例配对纵向患者）[35]、`GSE207422` 人类 NSCLC 严格治疗前 bulk RNA-seq 活检（$n = 24$ 例，`Pre_biopsy`）与肿瘤微环境验证 scRNA-seq（$n = 14$ 例）[36]、`GSE171351` 原代人类膀胱癌 10x Visium 切片（$n = 4$ 张）[32] 以及 `GSE203612` 原代人类乳腺癌、卵巢癌、肝癌与胃肠间质瘤 10x Visium 切片（$n = 8$ 张；合计 12 张切片、21,796 个空间点位）[33]。在 `GSE123813` 中，计算基质 `myCAF` 与 `TAM` 比例时将独立的 `CD45+CD3+` FACS 富集 T 细胞分选管与非 `CD3+` 微环境细胞严格解耦，并在全部 8 例患者中量化了相对于未校正全细胞混合池的稀释倍数。
6. **自身免疫黏膜炎症与纵向脓毒症**：`GSE250487` 人类溃疡性结肠炎黏膜活检（$n = 28$ 例）[17] 与 `GSE167363` 纵向人类脓毒症 PBMC（$n = 12$ 例）[18]。

### 动态特异性漂移与动态微环境可及性（`DNAcc`）数学构建

对于任意靶基因 $g$ 在 $K$ 个细胞微区上的归一化平均对数表达量 $x_{g,k} \ge 0$，按 Yanai 特异性指数公式 [8] 分别计算健康基线（$\tau_0$）、早期病变（$\tau_{\mathrm{early}}$）与晚期病变（$\tau_{\mathrm{late}}$）特异性：

$$\tau(g) = \frac{\sum_{k=1}^{K} \left(1 - \frac{x_{g,k}}{\max_{j} x_{g,j}}\right)}{K - 1}$$

动态特异性漂移定义为 $\Delta\tau(g) = \tau_{\mathrm{late}}(g) - \tau_{\mathrm{early}}(g)$，动态谱系保留度定义为 $R_{\mathrm{dyn}}(g) = 0.5\,\tau_{\mathrm{early}}(g) + 0.5\,\mathrm{clip}(0.5 + \Delta\tau(g), 0, 1)$，异位微区切换指标 $I_{\mathrm{switch}}(g) \in \{0, 1\}$ 定义为有表达靶基因在早期与晚期病变间主导表达谱系（$\arg\max_k x_{g,k}$）发生改变。重要器官旁观风险 $H_{\mathrm{vital}}(g)$ 量化靶基因在健康心肌细胞、肝细胞、肾小管上皮与肺泡上皮中的峰值表达相对于目标谱系表达的比值。模态-基质位阻相容性 $B_{\mathrm{acc}}(m, d)$ 耦合药物流体力学类别 $m$ 与疾病轴特异的间质纤维化或 `myCAF` 基质密度 $d \in [0, 1]$，在高基质密度轴中赋予约 150-kDa IgG 抗体低于小分子与低分子量多肽的间质穿透系数。

在六步正交消融（`M1` 至 `M6`）中，采用分层五折交叉验证与按主靶点基因（`primary_target_gene`）分组隔离的五折 `GroupKFold` 生成折外预测概率评估嵌套逻辑模型（`M1`：$\tau_0$；`M2`：$+ H_{\mathrm{vital}}$；`M3`：$+ \tau_{\mathrm{early}}$；`M4`：$+ \Delta\tau + R_{\mathrm{dyn}}$；`M5`：$+ I_{\mathrm{switch}} + S_{\mathrm{stromal}}$；`M6`：完整 `DNAcc`，纳入 $B_{\mathrm{acc}}$ 与 $\Delta\tau \times B_{\mathrm{acc}}$ 交互项）。多靶点药物按靶基因集合聚合特征，并在单靶点单药（$n_{\mathrm{targets}} = 1$）与多靶点方案（$n_{\mathrm{targets}} \ge 2$）亚组中分别独立核验。

### Visium 空间邻近度、连续径向衰减与随机朗之万 Waddington 求解器

针对每张 10x Visium 切片（12 张切片共 21,796 个点位），在组织空间坐标 $\mathbf{r}_u = (x_u, y_u)$ 上对点位级微区活性得分 $z_c(u)$ 进行标准化。按欧氏距离构建物理 $k = 6$ 最近邻空间图，由对称双变量 Moran 空间滞后量 $I_{ij}$ 定义微区 $i$ 与 $j$ 之间的空间邻近度 $P_{ij} = (1 + \exp(-4.5\,I_{ij}))^{-1}$，并定义基质排斥指数 $\mathrm{SEI} = P_{\mathrm{myCAF,TAM}} - P_{\mathrm{Tumor,CD8}}$。在连续径向剖析中，取各切片内 `Tumor` 特征最高四分位点位为恶性核心，计算全切片每个点位至恶性核心边界的最短欧氏距离并转换为 Visium 点距单位（`0` 至 `6` 个点距，按 100-µm 中心点距对应约 `0` 至 `600 µm`）。

多细胞组织状态轨迹 $\mathbf{S}(t) \in [0, 1]^N$（微区 $j \in \{1, \dots, N\}$）在总时长 $T = 96.0\text{ h}$、步长 $\Delta t = 0.5\text{ h}$ 下按状态有界随机 Itô–朗之万微分方程积分求解 [44]：

$$dS_j(t) = \left[ \alpha_j \rho_j(t) + (1 - S_j(t)) \rho_j(t) \left( \sum_{i \in \mathcal{A}_j} v_{ij} \frac{(\rho_i(t) S_i(t))^{n_{ij}}}{K_{ij}^{n_{ij}} + (\rho_i(t) S_i(t))^{n_{ij}}} \right) \prod_{q \in \mathcal{I}_j} \frac{K_{qj}^{n_{qj}}}{K_{qj}^{n_{qj}} + (\rho_q(t) S_q(t))^{n_{qj}}} - \gamma_j S_j(t) \right] dt + \sigma_j \sqrt{S_j(t)(1 - S_j(t))}\, dW_j(t)$$

其中 $\alpha_j$、$\gamma_j$ 与 $\sigma_j$ 分别表示基线生成率、衰减率与状态有界人口统计噪声幅值；$v_{ij}$（由实测供体 Spearman 耦合强度与 Visium 空间邻近度 $P_{ij}$ 加权）、$K_{ij}$ 与 $n_{ij}$ 分别表示有向互作权重、半最大常数与 Hill 协同系数；$\mathcal{A}_j$ 与 $\mathcal{I}_j$ 分别表示激活输入集与抑制输入集；$\rho_j(t) = \prod_{r \in \mathcal{R}_j(t)} (1 - \mathrm{clip}(\eta_r B_{\mathrm{acc},r}, -0.95, 0.98))$ 表示治疗方案 $\mathcal{R}_j$（具岭回归 GRN 扰动效力 $\eta_r$）作用下的未受阻断信号比例。每个条件下基于 200 至 250 条蒙特卡洛轨迹计算满足 $S_{\mathrm{remission}}(T) \le \theta_{\mathrm{remission}}$（纤维化 `Myofibroblast` 阈值 `0.36`；实体瘤 `Tumor` 阈值 `0.45`）的轨迹比例作为系综缓解概率。二维 Waddington 准势能面 $U(x, y) = -\ln P_{\mathrm{ss}}(x, y)$ 由 $60 \times 60$ 状态网格上经高斯平滑（$\sigma = 2.2$ 网格箱）的稳态朗之万系综概率密度估计得出 [44]。连续 25 点剖面似然可辨识性曲线在 $v_{ij} \in [0.50, 1.50] \times v_{ij,\mathrm{MLE}}$ 范围内扫描每条有向边权重，同时通过有界 Brent 优化在 $k_{\mathrm{scale}} \in [0.82, 1.22]$ 内重优化全局干扰速率尺度以拟合人类健康至肝硬化实测演进状态，双侧 $\Delta\chi^2 > 3.841$（95% $\chi^2_1$ 临界值）判定为结构可辨识 [45]。

### 统计分析

同一队列上的相关 ROC 曲线采用双侧配对非参数 DeLong 检验比较 [19]。AUROC、AUPRC 与 Brier 分数的 95% 置信区间由 2,000 次分层自助重采样（bootstrap）计算。针对小样本临床基准队列（14 项 IIb/III 期抗纤维化试验及 $n = 8$ 至 $32$ 例实体瘤 ICB 队列），采用向量化 Mann–Whitney 秩统计量基于 10,000 次蒙特卡洛标签置换计算精确经验 *P* 值。多变量混杂控制采用同时报告解析逆黑塞矩阵 Wald 标准误与 Liang–Zeger 靶基因聚类稳健三明治标准误（按 $G = 386$ 至 $564$ 个独立 `primary_target_gene` 聚类并施加 $G/(G - 1)$ 有限聚类校正）的 L2 正则化（$C = 50.0$）固定效应逻辑回归，同时校正疾病轴（11 个指示变量）、药物模态、试验启动年代（`Pre_2010`、`2010_2014`、`2015_2019`、`2020_Plus`）、申办方类别（`INDUSTRY` 对比 `NON_INDUSTRY`）及多靶点药理学（$n_{\mathrm{targets}} \ge 2$ 对比单靶点）。分类变量 $2 \times 2$ 关联采用双侧 Fisher 精确检验，配对空间与纵向破障对比采用双侧 Wilcoxon 符号秩检验，不同阈值概率 $p_t$ 下的临床效用采用决策曲线净获益 $\mathrm{NB}(p_t) = \mathrm{TP}/N - (\mathrm{FP}/N)[p_t / (1 - p_t)]$ 评估 [20]。全流程由 37 项自动化 `pytest` 单元测试验证零 CNS 泄漏、否定感知停药解析、DeLong 与置换统计量、靶点聚类三明治推断、剖面似然、Visium 径向衰减及跨器官朗之万求解可复现性。

---

## 数据与代码可用性

本研究分析的全部临床试验、靶点机制、单细胞与空间转录组数据均取自公开权威数据库：ClinicalTrials.gov API v2（`https://clinicaltrials.gov/api/v2/studies`）、ChEMBL v33（`https://www.ebi.ac.uk/chembl/`）、Open Targets Platform v24.09（`https://platform.opentargets.org/`，遵循 CC0 1.0 协议并归属引用 Buniello 等 [5]）、Human Protein Atlas v23（`https://www.proteinatlas.org/`）、Europe PMC REST API（`https://www.ebi.ac.uk/europepmc/`），以及 NCBI GEO 数据集 `GSE136103`、`GSE212837`、`GSE135893`、`GSE271188`、`GSE171351`、`GSE203612`、`GSE123813`、`GSE207422`、`GSE250487` 与 `GSE167363`。OmniNiche Python 软件包源码、确定性统计表生成管线（`outputs/stage2_tables/`）、出版级配图渲染引擎（`outputs/stage2_figures/`）及 37 项 `pytest` 测试套件均完整归档于项目仓库（`OmniNiche_Universal`）。

---

## 致谢、作者贡献与利益冲突声明

作者感谢向 NCBI GEO、Human Protein Atlas、ChEMBL 与 Open Targets Platform 贡献原代人类活检单细胞转录组、空间 Visium 转录组及临床试验注册数据的患者、临床研究者与数据联盟。全部引文元数据均经 NCBI PubMed E-utilities 与 Europe PMC 实时交叉核验。研究构思、数学建模、数据摄取、统计加固、可视化重制与手稿撰写均在 OmniNiche Universal 研究计划内完成。作者声明不存在任何财务或非财务利益冲突。

---

## 表格

### Table 1 | `DNAcc` 六步正交消融、未见靶点泛化、多变量混杂控制及分疾病轴临床成熟预测表现

| 基准队列 / 模型设定 | 试验数 ($n$) | 靶点数 ($n$) | 静态健康 $\tau_0$ AUROC (95% CI) | 完整 `DNAcc` AUROC (95% CI) | $\Delta\mathrm{AUROC}$ | 配对 DeLong *P* | 每 1-s.d. `DNAcc` 比值比 (95% CI; Wald / 靶点聚类 *P*) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **去混杂科学基准队列（全部 12 疾病轴，分层五折 CV）** | 11,843 | 386 | 0.5283 (0.5092–0.5483) | 0.6683 (0.6477–0.6884) | +0.1400 | $8.68 \times 10^{-28}$ | 未校正：1.87 (1.73–2.02; $1.31 \times 10^{-57}$ / $4.36 \times 10^{-16}$) |
| **去混杂科学基准队列（未见靶点 `Group5Fold` CV）** | 11,843 | 386 | 0.4617 (0.4411–0.4823) | 0.6540 (0.6334–0.6734) | +0.1923 | $3.30 \times 10^{-40}$ | 多变量固定效应校正：1.43 (1.29–1.59; $9.21 \times 10^{-12}$ / $2.58 \times 10^{-6}$) |
| **去混杂科学基准队列（仅限单靶点单药疗法）** | 6,194 | 242 | 0.5221 (0.4937–0.5515) | 0.6442 (0.6077–0.6818) | +0.1221 | $1.63 \times 10^{-6}$ | 联合模型 (+$\tau_0$ +固定效应)：1.43 (1.28–1.59; $1.14 \times 10^{-10}$ / $1.18 \times 10^{-5}$) |
| **核心 6 个单细胞疾病轴（分层五折 / `Group5Fold` CV）** | 5,815 | 309 | 0.5252 (0.4994–0.5505) | 0.6796 (0.6559–0.7022) / 0.6622 | +0.1544 | $1.06 \times 10^{-22}$ | 多变量固定效应校正：1.61 (1.41–1.83; $5.28 \times 10^{-13}$ / $1.65 \times 10^{-6}$) |
| **II 期向 III/IV 期临床成熟队列（全部 12 轴，分层五折 CV）** | 23,967 | 564 | 0.5405 (0.5332–0.5475) | 0.6737 (0.6670–0.6802) | +0.1332 | $3.99 \times 10^{-182}$ | 未校正：1.88 (1.82–1.93; $<10^{-300}$ / $2.14 \times 10^{-24}$) |
| **II 期向 III/IV 期临床成熟队列（未见靶点 `Group5Fold` CV）** | 23,967 | 564 | 0.5168 (0.5095–0.5240) | 0.6540 (0.6472–0.6607) | +0.1372 | $6.76 \times 10^{-189}$ | 多变量固定效应校正：1.29 (1.25–1.34; $3.83 \times 10^{-44}$ / $8.85 \times 10^{-6}$) |
| **轴 1：肝脏 MASH 与肝硬化（`GSE136103` + `GSE212837`）** | 388 | 89 | 0.4114 (0.3547–0.4688) | 0.5757 (0.5148–0.6339) | +0.1643 | $1.01 \times 10^{-5}$ | 1.30 (1.06–1.60; $1.09 \times 10^{-2}$) |
| **轴 2：肺部 IPF 与 ARDS（`GSE135893`）** | 1,946 | 194 | 0.5853 (0.5577–0.6114) | 0.6248 (0.5977–0.6511) | +0.0395 | $1.39 \times 10^{-4}$ | 1.55 (1.41–1.71; $7.52 \times 10^{-20}$) |
| **轴 3：肾脏 CKD 与 AKI（`GSE271188`）** | 908 | 132 | 0.4769 (0.4390–0.5148) | 0.5722 (0.5337–0.6116) | +0.0953 | $5.88 \times 10^{-4}$ | 1.29 (1.12–1.47; $2.72 \times 10^{-4}$) |
| **轴 4：实体肿瘤（`GSE207422` + `GSE123813` + 12 张 Visium）** | 9,255 | 388 | 0.5313 (0.5182–0.5454) | 0.5781 (0.5648–0.5921) | +0.0467 | $6.36 \times 10^{-9}$ | 1.40 (1.33–1.46; $3.98 \times 10^{-43}$) |
| **轴 5：自身免疫 IBD 与风湿（`GSE250487`）** | 1,507 | 158 | 0.6037 (0.5762–0.6322) | 0.6949 (0.6666–0.7224) | +0.0912 | $1.12 \times 10^{-10}$ | 2.04 (1.82–2.30; $3.12 \times 10^{-33}$) |
| **轴 6：全身脓毒症边界条件（`GSE167363` PBMC）** | 144 | 54 | 0.4738 (0.3786–0.5748) | 0.4793 (0.3778–0.5760) | +0.0055 | $0.919$ | 0.90 (0.65–1.26; $0.555$) |

### Table 2 | 人类泛器官纤维化零样本临床试验预测、双向 $3 \times 3$ 器官迁移及盲法治疗前实体瘤免疫治疗基准表现

| 转化基准队列 / 跨器官迁移设定 | 样本量 ($n$) | 阳性事件 ($n$) | OmniNiche 空间/朗之万 AUROC (95% CI) | OmniNiche AUPRC / Brier | 临床 PD-L1 或静态 $\tau_0$ / 最高静态基线 AUROC | 精确置换 *P* (10,000 次置换) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **泛器官人类纤维化 IIb/III 期竞技场（肝脏 $\to$ 3 器官）** | 14 项试验 | 6 项阳性 | 0.9583 (0.8125–1.0000) | 0.9583 / 0.1257 | 0.5000 ($\tau_0$) / 0.9167 (静态 LR) | $P = 0.00280$ |
| **双向迁移：仅肺脏 IPF ($n = 30$) $\to$ 留出器官** | 14 (留出 10) | 6 项阳性 | 留出：0.9583；全部 14 项：0.9792 | Brier：0.1101 (留出) / 0.0787 (全 14 项) | 0.5000 (静态 $\tau_0$) | 留出肝 0.9444，留出肾 1.00 |
| **双向迁移：仅肾脏 CKD ($n = 17$) $\to$ 留出器官** | 14 (留出 13) | 6 项阳性 | 留出：1.0000；全部 14 项：1.0000 | Brier：0.1115 (留出) / 0.1036 (全 14 项) | 0.5000 (静态 $\tau_0$) | 留出肝 1.0000，留出肺 1.00 |
| **留一器官交叉迁移 (LOOO)：肺 + 肾 $\to$ 留出肝脏** | 留出 9 项 | 3 项阳性 | 留出：0.9444；全部 14 项：0.9792 | Brier：0.1399 (留出) / 0.0899 (全 14 项) | 0.5000 (静态 $\tau_0$) | LOOO 肺 1.00，LOOO 肾 1.00 |
| **治疗前 BCC 抗 PD-1 全微环境 scRNA-seq (`GSE123813`)** | 8 例患者 | 4 例应答 | 0.8750 (0.6250–1.0000) | 0.8750 / 0.2716 | 0.4375 (`CD274`) / 0.8125 (线性均值) | $P = 0.1425$ (AUROC/AUPRC 列第一) |
| **治疗前 NSCLC 抗 PD-1 活检 Bulk RNA-seq (`GSE207422`)** | 24 例患者 | 9 例应答 | 0.7407 (0.5222–0.9333) | 0.6025 / 0.2250 | 0.6593 (`CD274`) / 0.8000 (静态 $\tau_0$, Brier 0.2273) | $P = 0.04660$ (最优 Brier 0.2250) |
| **NSCLC 抗 PD-1 微环境验证 scRNA-seq (`GSE207422`)** | 14 例患者 | 4 例应答 | 0.8500 (0.6000–1.0000) | 0.6083 / 0.1751 | 0.3250 (`CD274`) / 0.8250 (线性均值) | $P = 0.05259$ (三项指标均列第一) |
| **合并严格治疗前人类肿瘤活检队列（零重拟合）** | 32 例患者 | 13 例应答 | 0.7551 (0.5769–0.9190) | 0.6385 / 0.2366 | 0.5344 (`CD274`) / 0.6235 (线性均值) | $P = 0.01480$ (AUROC/AUPRC/Brier 均第一) |

---

## 图注

**Fig. 1 | OmniNiche 动态单细胞微环境重塑与生物物理空间位阻排斥框架总览。** **a**，跨越 12 个全身与实质器官疾病轴的临床试验注册库摄取（$n = 41,396$ 项独立 II–IV 期试验），将干预方案解析至靶基因与治疗模态（小分子、多肽、IgG 抗体），并通过子句感知的自然语言解卷积剥离生物学疗效或安全性失败与行政商业中止。**b**，融合原代人类单细胞特异性漂移（$\Delta\tau$）、异位谱系切换、重要器官旁观毒性风险与模态依赖性细胞外基质（ECM）空间位阻输运的动态微环境可及性（`DNAcc`）定律。**c**，数据编译的随机 Itô–朗之万多细胞动力学引擎，刻画跨器官保守的三微区纤维化正反馈（受损实质细胞、$SPP1^+$ 巨噬细胞与肌成纤维细胞）及双稳态 Waddington 准势能迟滞鞍点。**d**，跨人类肝、肺、肾抗纤维化临床试验的零样本迁移预测及 Visium 空间约束下的实体瘤免疫治疗破障联合策略。

**Fig. 2 | 41,396 项全身与实质器官临床试验解卷积揭示行政标签混杂与大分子疗效损耗。** **a**，连接靶基因解析状态（$n = 29,791$ 项已解析对比 $n = 11,605$ 项未解析）、药物模态与子句解卷积临床结局的三级桑基流图。**b**，12 个全身与实质器官疾病轴的 II 期、III 期与 IV 期独立临床试验普查及靶点解析率。**c**，2005–2025 年二十年间行政运营中止占停药试验比例、单克隆抗体占靶点明确试验比例及明确生物学疗效加安全性失败比例的纵向演变（阴影带为 95% Wilson 二项置信区间）。**d**，相对于小分子药物的各治疗模态终止比值比（双侧 Fisher 精确检验 95% CI），分别展示明确疗效无效、明确安全性/毒性与行政商业中止。**e**，不同结局分层的单试验受试者入组暴露规模分布（对数坐标下的中位数、四分位距箱体与第 10–90 百分位须线，并标注累计受试者人数）。**f**，子句感知停药审计器在 240 项分层人工裁定基准（每类 60 项；总体准确率 = 88.3%，Cohen's $\kappa = 0.844$）上的混淆矩阵与行召回率。

**Fig. 3 | 人类单细胞靶点特异性漂移与动态微环境可及性（`DNAcc`）定律。** **a**，6 个核心人类单细胞疾病轴上 24 个代表性临床靶基因的单细胞特异性漂移（$\Delta\tau = \tau_{\mathrm{late}} - \tau_{\mathrm{early}}$）点阵热图；圆点大小编码临床试验数 $\log_{10}$ 值，黑色外环标记发生异位谱系切换的靶点。**b**，去混杂科学队列（$n = 11,843$ 项试验）上的 `DNAcc` 六步正交消融（`M1` 至 `M6`），对比分层五折交叉验证与未见靶点 `Group5Fold` 交叉验证（误差棒为 2,000 次自助重采样 95% CI；双侧配对 DeLong *P* 值）。**c**，每增加 1 个标准差的标准化比值比（95% Wald CI），对比未校正静态 $\tau_0$、未校正 `DNAcc`、五项混杂因素固定效应校正后的 `DNAcc` 以及 6 个单细胞疾病轴的分轴效应量。**d**，按静态健康 $\tau_0$ 与 `DNAcc` 四分位分层的实测 II 期向 III/IV 期临床成熟率（$n = 23,967$ 项试验；误差棒为 95% Wilson CI）。**e**，23,967 项试验队列上的十等分折外概率校准曲线（上）与跨阈值概率的决策曲线分析临床净获益（下）。**f**，临床试验数 $\ge 4$ 的靶点在 $(\tau_0, \Delta\tau)$ 重塑相空间（基线取自人类蛋白质图谱 154 种细胞类型）中的分布，标出谱系保留型成功靶点（`PCSK9`、`THRB`、`PDGFRB`）与漂移侵蚀或异位切换型失败靶点（`LOXL2`、`TGFB1`）。**g**，涵盖严格否定词过滤、将全部 3,118 项未报告停药归为生物学失败的最坏情境（$n = 14,961$）、逆概率重加权及仅限已公布结果试验（$n = 5,757$）的五重敏感性分析。

**Fig. 4 | 人类泛器官纤维化临床试验零样本预测、双向器官迁移与双稳态迟滞破障。** **a**，跨越肝脏（$n = 17$，`GSE136103` + `GSE212837`）、肺脏（$n = 30$，`GSE135893`）与肾脏（$n = 17$，`GSE271188`）共 64 例原代人类单细胞/单细胞核供体的对照归一化肌成纤维细胞微区富集倍数（柱高为均值 $\pm$ 标准误；散点为独立人类供体）。**b**，从外周血单核细胞经健康肝巨噬细胞至 F2/F3 期 MASH 与 F4 期肝硬化的人类髓系伪批量对数表达轨迹（均值 $\pm$ 标准误），展示肝内 `CCR2` 在 $SPP1^+TREM2^+$ SAMac 涌现前的坍塌。**c**，横跨肝、肺、肾 14 项标志性 IIb/III 期抗纤维化临床试验（6 项临床阳性对比 8 项疗效失败）的单试验零样本朗之万模拟缓解概率，标注器官与药物模态。**d**，带 2,000 次自助重采样 95% CI 与 10,000 次置换精确 *P* 值的五基线零样本分类竞技场（左），以及双向 $3 \times 3$ 单器官编译与留一器官（LOOO）跨器官迁移 AUROC 矩阵（右）。**e**，叠加确定性漂移矢量场、缓解阈值线（垂直虚线）与锁定态 F4 吸引盆的平滑 $60 \times 60$ 二维 Waddington 准势能面 $U(S_{\mathrm{Myofib}}, S_{\mathrm{Paren}})$。**f**，扫描从早期 F1/F2 至锁定态 F4 肝硬化共 25 个解剖进展状态的连续纤维化迟滞分岔曲线（单药治疗窗口断崖 $t_{50} = 17.4\text{ h}$），对比基质单药、实质 `THRB` 单药与 `Myofibroblast + SAMac` 双微区破环联合疗法。**g**，锁定态 F4 肝硬化（$s_0 = [0.85, 0.88, 0.90]$）下的对称 $6 \times 6$ 临床靶点二维 Bliss 超额协同矩阵。

**Fig. 5 | 人类实体癌空间 Visium 瘤周 `myCAF` 排斥与零重拟合治疗前抗 PD-1 应答分层。** **a**，两张代表性原代人类癌组织 10x Visium 空间点位图（`GSE203612`），对比免疫排斥型乳腺癌 `myCAF` 盾牌（`GSM6177599`，$\mathrm{SEI} = +0.21$）与相对可浸润肝细胞癌微环境（`GSM6177612`，$\mathrm{SEI} = +0.07$）。**b**，跨越 12 张人类 Visium 切片（共 21,796 个点位）自恶性核心边界向外 0 至 600 µm 的共识连续径向特征衰减曲线（均值 $z$ 分数 $\pm$ 标准误；针对 $\mathrm{SEI} > 0$ 的双侧 Wilcoxon 符号秩检验）。**c**，$n = 32$ 例严格治疗前肿瘤活检（`GSE123813` BCC scRNA-seq $n = 8$ 与 `GSE207422` NSCLC bulk 治疗前活检 $n = 24$）的单患者零重拟合预测缓解概率瀑布图（相对于 50% 判别阈值）。**d**，四个盲测队列上的应答预测 AUROC 基准对比（带 2,000 次自助重采样 95% CI），比较 OmniNiche 空间朗之万模型与 `CD274`（PD-L1）、静态健康 $\tau_0$ 及单微区 T 细胞 GRN 评分。**e**，`GSE123813` 中 $n = 15$ 例人类 BCC 患者配对抗 PD-1 治疗前后的纵向净 $CD8^+$ 细胞毒性效应轨迹（细线为单患者，粗线为 $n = 8$ 例应答者与 $n = 7$ 例无应答者组均值）。**f**，`anti-PD-1 + myCAF/TGFBR1` 联合破障治疗下的二维实体瘤 Waddington 准势能面 $U_{\mathrm{combo}}(S_{\mathrm{myCAF}}, S_{\mathrm{Tumor}})$ 与确定性漂移矢量场。**g**，全部 $n = 19$ 例抗 PD-1 耐药治疗前患者的配对计算机模拟破障挽救轨迹（双侧 Wilcoxon 符号秩检验）。

**Extended Data Fig. 1 | PRISMA 临床试验过滤漏斗、高频靶点普查、多靶点药理学稳健性与核心 6 轴未见靶点泛化。** **a**，从 41,983 条原始 II–IV 期注册记录至 41,396 项独立全身与实质器官试验、29,791 项靶点明确试验、23,967 项非行政中止确切结局试验及 11,843 项去混杂科学结局试验的 PRISMA 过滤漏斗。**b**，12 个疾病轴上测试频次最高的前 18 个人类靶基因。**c**，多靶点数量分布以及严格单靶点单药（$n = 6,194$）对比多靶点方案（$n = 5,649$）亚组中的 `DNAcc` AUROC 表现。**d**，单细胞实测覆盖的 6 个疾病轴（$n = 5,815$ 项试验）在分层五折与未见靶点 `Group5Fold` 交叉验证下的 `DNAcc` 六步正交消融。

**Extended Data Fig. 2 | 原代人类单细胞与空间队列普查、异位谱系切换频率、漂移分布与分轴区分表现。** **a**，原代人类单细胞、单细胞核与 10x Visium 队列的供体与组织切片普查（合计 $n = 162$ 例人类供体/切片）。**b**，6 个疾病轴中在早期与晚期病变之间发生主导表达细胞谱系切换的临床靶基因百分比。**c**，6 个疾病轴单细胞特异性漂移（$\Delta\tau$）的箱线须图分布（中位数、四分位距与 $1.5 \times \mathrm{IQR}$ 须线）。**d**，分疾病轴 II 期向 III/IV 期临床成熟 AUROC（2,000 次自助重采样 95% CI），对比 `DNAcc` 与静态健康 $\tau_0$，展示全部 5 个实体/黏膜屏障轴的一致增益以及血管内 PBMC 脓毒症的零区分边界表现。

**Extended Data Fig. 3 | 连续剖面似然可辨识性、多变量固定效应控制、模拟投资组合富集与科学队列校准。** **a**，人类纤维化朗之万语法全部四条有向细胞间互作边的连续 25 点双侧剖面似然 $\Delta\chi^2$ 可辨识性曲线（伴干扰速率尺度重优化，虚线标示 95% $\chi^2_1$ 阈值 $\Delta\chi^2 = 3.84$）。**b**，11,843 项科学队列中单变量、五项混杂因素多变量固定效应及头对头联合竞争（`DNAcc` + 静态 $\tau_0$ + 固定效应）逻辑模型的每标准差比值比。**c**，23,967 项试验的 II 期向 III 期临床管线富集效用模拟，展示不同 `DNAcc` 与静态 $\tau_0$ 准入百分位下的 III/IV 期成功率与避免的无效受试者暴露。**d**，11,843 项去混杂科学队列上的十等分折外可靠性校准曲线。

**Extended Data Fig. 4 | 单切片 10x Visium 基质排斥指数、共识空间邻近度矩阵、FACS 分选去混杂审计与跨平台应答对齐。** **a**，全部 12 张原代人类癌组织 10x Visium 切片（`GSE171351` 与 `GSE203612`，共 21,796 个点位）的基质排斥指数（$\mathrm{SEI}$）与点位计数。**b**，跨越 12 张人类肿瘤切片的共识 $4 \times 4$ 物理 $k = 6$ 最近邻双变量 Moran 空间邻近度矩阵。**c**，`GSE123813`（$n = 8$ 例治疗前 BCC 患者）实验方案技术审计，对比受 `CD45+CD3+` FACS 分选管稀释的全细胞 `myCAF` 比例与去混杂后的非 `CD3+` 微环境基质比例（中位恢复 1.7 倍，最高 4.8 倍）。**d**，BCC scRNA-seq（$n = 8$）与 NSCLC bulk 治疗前活检（$n = 24$）跨平台零重拟合治疗前模拟缓解概率对比（均值 $\pm$ 标准误及单患者散点）。

---

## 参考文献

1. Sun D, Gao W, Hu H, Zhou S. Why 90% of clinical drug development fails and how to improve it? *Acta Pharm Sin B* 2022;12(7):3049-3062. PMID: 35865092. DOI: 10.1016/j.apsb.2022.02.002.
2. Wong CH, Siah KW, Lo AW. Estimation of clinical trial success rates and related parameters. *Biostatistics* 2019;20(2):273-286. PMID: 29394327. DOI: 10.1093/biostatistics/kxx069.
3. Nelson MR, Tipney H, Painter JL, et al. The support of human genetic evidence for approved drug indications. *Nat Genet* 2015;47(8):856-860. PMID: 26121088. DOI: 10.1038/ng.3314.
4. Minikel EV, Painter JL, Dong CC, Nelson MR. Refining the impact of genetic evidence on clinical success. *Nature* 2024;629(8012):624-629. PMID: 38632401. DOI: 10.1038/s41586-024-07316-0.
5. Buniello A, Suveges D, Cruz-Castillo C, et al. Open Targets Platform: facilitating therapeutic hypotheses building in drug discovery. *Nucleic Acids Res* 2025;53(D1):D1467-D1475. PMID: 39657122. DOI: 10.1093/nar/gkae1128.
6. Huang K, Chandak P, Wang Q, et al. A foundation model for clinician-centered drug repurposing. *Nat Med* 2024;30(12):3601-3613. PMID: 39322717. DOI: 10.1038/s41591-024-03233-x.
7. Karlsson M, Zhang C, Mear L, et al. A single-cell type transcriptomics map of human tissues. *Sci Adv* 2021;7(31):eabh2169. PMID: 34321199. DOI: 10.1126/sciadv.abh2169.
8. Kryuchkova-Mostacci N, Robinson-Rechavi M. A benchmark of gene expression tissue-specificity metrics. *Brief Bioinform* 2017;18(2):205-214. PMID: 26891983. DOI: 10.1093/bib/bbw008.
9. Ramachandran P, Dobie R, Wilson-Kanamori JR, et al. Resolving the fibrotic niche of human liver cirrhosis at single-cell level. *Nature* 2019;575(7783):512-518. PMID: 31597160. DOI: 10.1038/s41586-019-1631-3.
10. Fabre T, Barron AMS, Christensen SM, et al. Identification of a broadly fibrogenic macrophage subset induced by type 3 inflammation. *Sci Immunol* 2023;8(82):eadd8945. PMID: 37027478. DOI: 10.1126/sciimmunol.add8945.
11. Mariathasan S, Turley SJ, Nickles D, et al. TGFbeta attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. *Nature* 2018;554(7693):544-548. PMID: 29443960. DOI: 10.1038/nature25501.
12. Jain RK. Transport of molecules in the tumor interstitium: a review. *Cancer Res* 1987;47(12):3039-3051. PMID: 3555767.
13. Williams RJ, Tse T, DiPiazza K, Zarin DA. Terminated Trials in the ClinicalTrials.gov Results Database: Evaluation of Availability of Primary Outcome Data and Reasons for Termination. *PLoS One* 2015;10(5):e0127242. PMID: 26011295. DOI: 10.1371/journal.pone.0127242.
14. Wang S, Li K, Pickholz E, et al. An autocrine signaling circuit in hepatic stellate cells underlies advanced fibrosis in nonalcoholic steatohepatitis. *Sci Transl Med* 2023;15(677):eadd3949. PMID: 36599008. DOI: 10.1126/scitranslmed.add3949.
15. Habermann AC, Gutierrez AJ, Bui LT, et al. Single-cell RNA sequencing reveals profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary fibrosis. *Sci Adv* 2020;6(28):eaba1972. PMID: 32832598. DOI: 10.1126/sciadv.aba1972.
16. Kuppe C, Ibrahim MM, Kranz J, et al. Decoding myofibroblast origins in human kidney fibrosis. *Nature* 2021;589(7841):281-286. PMID: 33176333. DOI: 10.1038/s41586-020-2941-1.
17. Friedrich M, Pohin M, Jackson MA, et al. IL-1-driven stromal-neutrophil interactions define a subset of patients with inflammatory bowel disease that does not respond to therapies. *Nat Med* 2021;27(11):1970-1981. PMID: 34675383. DOI: 10.1038/s41591-021-01520-5.
18. Qiu X, Li J, Bonenfant J, et al. Dynamic changes in human single-cell transcriptional signatures during fatal sepsis. *J Leukoc Biol* 2021;110(6):1253-1268. PMID: 34558746. DOI: 10.1002/JLB.5MA0721-825R.
19. DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. *Biometrics* 1988;44(3):837-845. PMID: 3203132.
20. Vickers AJ, Elkin EB. Decision curve analysis: a novel method for evaluating prediction models. *Med Decis Making* 2006;26(6):565-574. PMID: 17099194. DOI: 10.1177/0272989X06295361.
21. Anstee QM, Neuschwander-Tetri BA, Wai-Sun Wong V, et al. Cenicriviroc Lacked Efficacy to Treat Liver Fibrosis in Nonalcoholic Steatohepatitis: AURORA Phase III Randomized Study. *Clin Gastroenterol Hepatol* 2024;22(1):124-134.e1. PMID: 37061109. DOI: 10.1016/j.cgh.2023.04.003.
22. Harrison SA, Bedossa P, Guy CD, et al. A Phase 3, Randomized, Controlled Trial of Resmetirom in NASH with Liver Fibrosis. *N Engl J Med* 2024;390(6):497-509. PMID: 38324483. DOI: 10.1056/NEJMoa2309000.
23. Jara M, Norlin J, Kjaer MS, et al. Modulation of metabolic, inflammatory and fibrotic pathways by semaglutide in metabolic dysfunction-associated steatohepatitis. *Nat Med* 2025;31(9):3128-3140. PMID: 40691365. DOI: 10.1038/s41591-025-03799-0.
24. Francque SM, Bedossa P, Ratziu V, et al. A Randomized, Controlled Trial of the Pan-PPAR Agonist Lanifibranor in NASH. *N Engl J Med* 2021;385(17):1547-1558. PMID: 34670042. DOI: 10.1056/NEJMoa2036205.
25. Richeldi L, du Bois RM, Raghu G, et al. Efficacy and safety of nintedanib in idiopathic pulmonary fibrosis. *N Engl J Med* 2014;370(22):2071-2082. PMID: 24836310. DOI: 10.1056/NEJMoa1402584.
26. King TE Jr, Bradford WZ, Castro-Bernardini S, et al. A phase 3 trial of pirfenidone in patients with idiopathic pulmonary fibrosis. *N Engl J Med* 2014;370(22):2083-2092. PMID: 24836312. DOI: 10.1056/NEJMoa1402582.
27. Bakris GL, Agarwal R, Anker SD, et al. Effect of Finerenone on Chronic Kidney Disease Outcomes in Type 2 Diabetes. *N Engl J Med* 2020;383(23):2219-2229. PMID: 33264825. DOI: 10.1056/NEJMoa2025845.
28. Harrison SA, Abdelmalek MF, Caldwell S, et al. Simtuzumab Is Ineffective for Patients With Bridging Fibrosis or Compensated Cirrhosis Caused by Nonalcoholic Steatohepatitis. *Gastroenterology* 2018;155(4):1140-1153. PMID: 29990488. DOI: 10.1053/j.gastro.2018.07.006.
29. Raghu G, Richeldi L, Fernandez Perez ER, et al. Pamrevlumab for Idiopathic Pulmonary Fibrosis: The ZEPHYRUS-1 Randomized Clinical Trial. *JAMA* 2024;332(5):380-389. PMID: 38762797. DOI: 10.1001/jama.2024.8693.
30. Harrison SA, Wong VW, Okanoue T, et al. Selonsertib for patients with bridging fibrosis or compensated cirrhosis due to NASH: Results from randomized phase III STELLAR trials. *J Hepatol* 2020;73(1):26-39. PMID: 32147362. DOI: 10.1016/j.jhep.2020.02.027.
31. Kieffer Y, Hocine HR, Gentric G, et al. Single-Cell Analysis Reveals Fibroblast Clusters Linked to Immunotherapy Resistance in Cancer. *Cancer Discov* 2020;10(9):1330-1351. PMID: 32434947. DOI: 10.1158/2159-8290.CD-19-1384.
32. Gouin KH 3rd, Ing N, Plummer JT, et al. An N-Cadherin 2 expressing epithelial cell subpopulation predicts response to surgery, chemotherapy and immunotherapy in bladder cancer. *Nat Commun* 2021;12(1):4906. PMID: 34385456. DOI: 10.1038/s41467-021-25103-7.
33. Barkley D, Moncada R, Pour M, et al. Cancer cell states recur across tumor types and form specific interactions with the tumor microenvironment. *Nat Genet* 2022;54(8):1192-1201. PMID: 35931863. DOI: 10.1038/s41588-022-01141-9.
34. Klughammer J, Abravanel DL, Segerstolpe A, et al. A multi-modal single-cell and spatial expression map of metastatic breast cancer biopsies across clinicopathological features. *Nat Med* 2024;30(11):3236-3249. PMID: 39478111. DOI: 10.1038/s41591-024-03215-z.
35. Yost KE, Satpathy AT, Wells DK, et al. Clonal replacement of tumor-specific T cells following PD-1 blockade. *Nat Med* 2019;25(8):1251-1259. PMID: 31359002. DOI: 10.1038/s41591-019-0522-3.
36. Hu J, Zhang L, Xia H, et al. Tumor microenvironment remodeling after neoadjuvant immunotherapy in non-small cell lung cancer revealed by single-cell RNA sequencing. *Genome Med* 2023;15(1):14. PMID: 36869384. DOI: 10.1186/s13073-023-01164-9.
37. Xu Z, Ren F, Wang P, et al. A generative AI-discovered TNIK inhibitor for idiopathic pulmonary fibrosis: a randomized phase 2a trial. *Nat Med* 2025;31(8):2602-2610. PMID: 40461817. DOI: 10.1038/s41591-025-03743-2.
38. Hoagland DA, Rodriguez-Morales P, Mann AO, et al. Macrophage-derived oncostatin M repairs the lung epithelial barrier during inflammatory damage. *Science* 2025;389(6756):169-175. PMID: 40638741. DOI: 10.1126/science.adi8828.
39. Yoo SK, Fitzgerald CW, Cho BA, et al. Prediction of checkpoint inhibitor immunotherapy efficacy for cancer using routine blood tests and clinical data. *Nat Med* 2025;31(3):869-880. PMID: 39762425. DOI: 10.1038/s41591-024-03398-5.
40. Shen W, Moon I, Nguyen TH, et al. Generalizable AI predicts immunotherapy outcomes across cancers and treatments. *Nat Med* 2026;32(8):3010-3022. PMID: 42399673. DOI: 10.1038/s41591-026-04502-7.
41. Li Z, Li Y, Xiang J, et al. AI-enabled virtual spatial proteomics from histopathology for interpretable biomarker discovery in lung cancer. *Nat Med* 2026;32(1):231-244. PMID: 41491099. DOI: 10.1038/s41591-025-04060-4.
42. Bunne C, Roohani Y, Rosen Y, et al. How to build the virtual cell with artificial intelligence: Priorities and opportunities. *Cell* 2024;187(25):7045-7063. PMID: 39672099. DOI: 10.1016/j.cell.2024.11.015.
43. Swanson K, Wu W, Bulaong NL, Pak JE, Zou J. The Virtual Lab of AI agents designs new SARS-CoV-2 nanobodies. *Nature* 2025;646(8085):716-723. PMID: 40730228. DOI: 10.1038/s41586-025-09442-9.
44. Wang J, Xu L, Wang E. Potential landscape and flux framework of nonequilibrium networks: robustness, dissipation, and coherence of biochemical oscillations. *Proc Natl Acad Sci U S A* 2008;105(34):12271-12276. PMID: 18719111. DOI: 10.1073/pnas.0800579105.
45. Raue A, Kreutz C, Maiwald T, et al. Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. *Bioinformatics* 2009;25(15):1923-1929. PMID: 19505944. DOI: 10.1093/bioinformatics/btp358.

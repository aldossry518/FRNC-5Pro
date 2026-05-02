# FRNC-5Pro
加热炉传热计算
● 1) FRNC-5PC 核心计算公式（按手册整理）

   - 总热平衡 / 炉膛热负荷与热效率
    - QABS = QR + QC
    - QABS = QF + HFUEL + HAIR - HeatLosses
    - η(热效率) = QABS / (QF + HFUEL + HAIR)（工程上常近似 QABS/QF）
    - mFUEL = QF / HCOMBUSTION
    - TFLAME = Tref + (QF + HFUEL + HAIR)/(mGAS
    
    * CpGAS)
   - 辐射段传热（Radiant）
    - 总体：QR = mGAS
    
    * CpGAS * (TFLAME - TBW) - FireboxLosses
    - 经验桥墙关系：TRAD - TBW = 40
    
    * (FireboxHeight / FireboxWidth)
    - 管元 i：
     - Qi = σ
     
    
    * Acp * F * (TRAD^4 - Touter,i^4) + ho * AT * (TRAD - Touter,i)
     - Qi = Uo,i
     
    
    * Ai * (Touter,i - Tbulk,i)
     - Qi = mOIL
     
    
    * CpOIL * (Tout,i - Tin,i)
    - 总传热系数：
     - Uo,i = 1 / (rFoul,out + rWall + rFoul,in + (1/hi)*(Ao/Ai))
   - 对流段传热（Convection）
    - 总体：QC = ΣQi
    - 管元 i：
     - Qi = Uo,i
     
    
    * Ai * LMTD
     - Qi = mOIL
     
    
    * CpOIL * (Tout,i - Tin,i)
     - Qi = mGAS
     
    
    * CpGAS * (Tg,in,i - Tg,out,i) - ConvectionLosses
   - 炉管内外壁温计算（与 API 530 关联）
    - 基本热阻网络：
     - ΔT = q
     
    
    * (Rout_foul + Rwall + Rin_foul + Rin_conv)
     - Twall,inner / Twall,outer 由内外侧热阻分配求得
    - 管内换热系数（API 530 单相常用关联）：
     - 液相：hi =
     
    
    0.023*(k/D)*Re^0.8*Pr^0.33*(μB/μW)^0.14
     - 气相：hi =
     
    
    0.021*(k/D)*Re^0.8*Pr^0.4*(TB/TW)^0.5
     - 过渡区：hi = hi(2000)*(1-F) + hi(10000)*F, F=0.000125*Re-0.25
   - 烟气温度分布计算
    - 分段能量平衡（辐射段→桥墙→对流段→烟囱）：
     - ΔTg = Qi / (mGAS
     
    
    * CpGAS)
    - 逐段迭代：由 TFLAME、TBW/TRAD、各段 Qi 收敛得到 Tstack
   - 介质压力降计算
    - 单相：
     - PT = PF + PG
     - PF 采用 Fanning/Moody 形式（L, D, f, G, ρ）
     - PG 为重力项（立管上升为正压降，下降为负）
    - 两相：
     - PT = PF + PG + PM
     - PL, PV（液相单独/气相单独压降）
     - Xtt = (PL/PV)^0.5 ≈ {(1-X)/X}^0.9*(ρV/ρL)^0.5*(μL/μV)^0.1
     - PF = PL
     
    
    * Φ^2（Lockhart-Martinelli 形式）
    - 对流段烟气侧抽力损失：
     - PF = C
     
    
    * f * G^2 * (Number of Rows)
   - 燃料燃烧计算
    - mFUEL = QF / HCOMBUSTION
    - mAIR = mFUEL
    
    * (理论空燃比) * (过量空气系数)
    - mFLUE = mFUEL + mAIR (+ 生成物修正)
    - 配合 CpGAS 与热平衡求 TFLAME/TBW/TSTACK

  -------------------------------------------------------------------------------------------------------------------------------------------------------

  2) 相关 API 530 / API 560 关键条文（用于上述计算）

   - API 530（炉管壁厚/应力）
    - 第4章：设计基础（弹性区 vs 蠕变破坏区、设计温度/设计寿命、薄壁管适用性）
    - 第5.2：平均直径应力方程（Equation 1）
    - 第5.3：弹性设计厚度（Eq.2）与最小厚度 tmin = t + CA（Eq.3）
    - 第5.4：蠕变设计厚度（Eq.4）与 tmin = t + fcorr*CA（Eq.5）
    - 第5.5：中间温区取弹性/蠕变两者中较大厚度
    - 第5.6：新管最小允许厚度下限（Table 1）
    - 第5.7：最小壁厚与平均壁厚换算按 ASTM 厚度公差
    - 第5.8：等效金属温度 Teq（Eq.6）及参数 V, N（Eq.7, Eq.8）
    - 附录A（信息性）：允许管壁温度、退役厚度、剩余寿命估算方法
   - API 560（加热炉整体设计）
    - 规定加热炉工艺/热设计边界：热负荷、效率、烟囱温度、过量空气、热损失
    - 要求辐射段/对流段热分配与温度校核（含桥墙温度、烟气侧校核）
    - 要求工艺侧与烟气侧压降限值及其核算
    - 炉管金属温度控制与材质选用；炉管强度/厚度计算通常与 API 530 联用
    - 机械设计、制造与试验要求与相关规范（如 ASME 体系）衔接执行

  -------------------------------------------------------------------------------------------------------------------------------------------------------

  3) FRNC-5PC 软件界面主要模块与布局

   - 主界面布局
    - 上部：项目/工况摘要（Case、单位、收敛控制）
    - 左侧：输入数据模块导航（分数据表/分段）
    - 中央：参数输入表单（可编辑字段、方法选择）
    - 下部/右侧：计算结果窗口（分表格与曲线）
   - 主要功能模块
    - 工况与流程配置（流股、入口出口条件、流量）
    - 炉膛/辐射段模块（炉膛几何、火焰/桥墙参数、辐射计算）
    - 对流段模块（管排、翅片/钉头/光管、对流换热与烟气侧损失）
    - 燃料与燃烧模块（燃料组成、空燃比、过量空气、烟气性质）
    - 管内换热与压降模块（单相/两相方法切换）
    - 计算选项模块（分段数、收敛容差、关联式选择）
    - 结果模块（热平衡、温度分布、壁温、压降、效率、报警/校核）
    - 报告模块（计算书导出、工况对比）

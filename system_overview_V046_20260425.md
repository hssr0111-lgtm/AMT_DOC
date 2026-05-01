# AMeDAS気温空間分布解析システム — 全体詳細解説
**VER.0.4.6.000　2026/04/25**

---

## 1. ファイル構成

| ファイル | 行数 | 役割 |
|---|---|---|
| `amt.cpp` | 2,696 | 解析ロジック・実行制御の全体（設定読込・初期化・実行モード） |
| `AnalAMeDAS.h` | 1,160 | 全構造体・定数・インライン関数の定義 |
| `AnalAMeDAS.cpp` | 7,052 | 初期化・解放・データ読込・出力関数群 |
| `Interpolation.cpp` | 2,480 | 補間計算・近傍点テーブル構築・エリア統計集計 |
| `OpticalFlow.cpp` | 616 | 勾配フロー・Horn-Schunck実装（AVX-512） |
| `HsCommon.h` | 1,373 | 汎用クラス（CHsString / CHsIniFile / CHsLog） |
| `WriteGeoTIFF.h` | 600 | GeoTIFF出力（GDAL）|
| `WriteGRD.h` | 1,493 | 独自バイナリ(GRD)出力・読込、BLOSC2圧縮 |
| `ConvMatrix.h` | 404 | float/int16変換・カラーマップ適用・PNG行生成 |
| `TMD.h` | 705 | 汎用時系列データロガー形式の読み書き |
| `Nifsq8.cpp/.h` | 1,060 | 補助モジュール（データ圧縮補助） |
| `toml.hpp` | — | TOML設定パーサ（ヘッダオンリー、外部ライブラリ） |
| `stb_truetype.h` | — | TrueTypeフォントレンダリング（ヘッダオンリー、外部ライブラリ） |
| `*.cfg` | — | 実行設定ファイル（INI形式） |

---

## 2. 構造体体系（AnalAMeDAS.h）

### 引数順の統一規則

```
DataDEM → DataAMeDAS → DataAnal → DataFlow
```

原則として構造体は全て参照渡し（const付き or 非const）。

---

### 2-1. AMeDASバイナリ入力系（`#pragma pack(1)`）

#### `AMeDAS_STDT`（観測所マスタ）
固定長バイト列で観測所番号・名称（漢字/カタカナ）・緯度・経度・標高を保持。`parse_fixed_int()` / `parse_fixed_char()` で数値化する。

#### `AMeDAS_TMDT`（時別値）
```c
int16_t  Rain;      // 降水量 [mm]   32767=欠測
uint8_t  WndDirect; // 風向 [16方位] 255=欠測
uint8_t  WndSpeed;  // 風速 [m/s]    255=欠測
int16_t  Sun;       // 日照時間 [0.1h] 32767=欠測
int16_t  Temp;      // 気温 [0.1℃]  32767=欠測
```

#### `AMeDAS_DYDT`（日別値）
日合計雨量・平均/最高/最低気温・平均風速など、リマークコード（0=正常/1=推定/9=欠測）付き。

#### `AMeDAS_DAY`（1日分レコード）
256バイト固定。1観測所×1日の `AMeDAS_TMDT[24]` + `AMeDAS_DYDT` のコンテナ。

---

### 2-2. 空間データ系

#### `DataDEM`
| フィールド | 内容 |
|---|---|
| `dLonOrigin`, `dLatOrigin` | 左上ピクセルの経緯度 |
| `dResX`, `dResY` | 1px当たりの経度・緯度幅（dResYは負値） |
| `nWidth` / `nHeight` | グリッドサイズ（日本全国250m: 7329×12057） |
| `nAlignedWidth` | 64bit境界アライメント後の1行要素数 |
| `total_allocated_pix` | 実確保ピクセル数（nAlignedWidth × nHeight） |
| `pElevation` | 標高配列 [float] |
| `pSeaMask` | 陸地マスク（1=陸, 0=海）[uint8_t] |
| `pCoastDist` | 全ピクセルの海岸距離 [m]（Meijster法EDT） |

#### `DataAMeDAS`
1ヶ月分の観測所データを保持。月替わり時に `FreeAMeDAS` → `LoadAMeDAS` で再ロードされる。

```c
int      count;          // 観測所数（約1300）
int      alignedNum;     // AVX-512アライメント後の観測所数
int32_t* numbers;        // 観測所番号
float*   lats, *lons, *alts;  // 緯度・経度・標高
int32_t* pixX, *pixY;         // DEM上のピクセル座標
AMeDAS_DAY* DT;          // [count × Days] の観測データ
float*   fTemps;         // 補間計算用気温（NaN=欠測）
float*   fRains;         // 補間計算用降水量（NaN=欠測）
// 複数時刻並列補間用バッファ
float*   pfTempsPool;    // 全スロット連続確保バッファ
float*   pfTempsArr[];   // スロット別ポインタ（[0]=fTemps）
```

#### `InterpConfig`（補間設定パラメータ）
```c
int   method;            // INTERP_* 定数
float idw_power;         // IDW/Shepard べき乗値（デフォルト 1.5）CFGキー: IDW_POWER
float rbf_sigma;         // RBF バンド幅 [m]（デフォルト 80000）  CFGキー: RBF_SIGMA
float barnes_kappa;      // Barnes 収束パラメータ（CalcInterpParamsで自動計算）
float barnes_kappa_fac;  // Barnes: kappa = d_nn² × この値（デフォルト 2.0、文献値 5.052）CFGキー: BARNES_KAPPA_FAC
float barnes_gamma;      // Barnes 2Pass スキャン間収縮率（デフォルト 0.3）CFGキー: BARNES_GAMMA
float coast_L;           // 海岸距離補正スケール [m]（0で補正なし）CFGキー: COAST_L
float tps_lambda;        // TPS 正則化パラメータ
```

`barnes_kappa` は `CalcInterpParams()` が全観測点ペアの平均最近傍距離 d_nn から `kappa = d_nn² × barnes_kappa_fac` として自動計算する。`rbf_sigma` / `idw_power` / `barnes_gamma` は自動計算されず、CFG値またはデフォルト値をそのまま使用する。

#### `DataAnal`（補間設定・作業バッファ・結果バッファ）

近傍点探索パラメータ（`PixCtx_*`）・気温減率テーブル・補間設定・TPS係数を保持するほか、以下の主要バッファを含む。

**近傍点テーブル（SoA構造、月1回BuildNbrTable()で構築）**

インデックス規則: `base = landIdx * NBR_MAX`（NBR_MAX=16）

| フィールド | 内容 | 型 |
|---|---|---|
| `pNbrN` | 有効近傍点数 [nLandCount] | uint8_t |
| `pNbrIdx` | 観測点インデックス [nLandCount × 16] | uint16_t |
| `pNbrD2` | 距離² [m²] [nLandCount × 16] | float |
| `pNbrLapse` | 高度補正係数 [°C/m] [nLandCount] | float |
| `pNbrAltDiff` / `pNbrAltDiffF` | 標高差 (int16 or float) | — |
| `pNbrCoastPx` / `pNbrCoastPxF` | ピクセル海岸距離 (int16 or float) | — |
| `pNbrCoastStn` / `pNbrCoastStnF` | 観測点海岸距離 (int16 or float) | — |
| `pNbrW` / `pNbrWF` | 重み (uint16 or float) | — |
| `pNbrAltW` / `pNbrAltWF` | lapse_rate×altDiff (int16 or float) | — |

`nNbrFloatMode` でint/float版を切り替え（0=全整数, 1=SRCのみfloat, 2=Wのみfloat, 3=全float）。

**気温補間結果（リングバッファ）**

```c
float*    pAnalTemp;              // 気温補間データ float（オプション確保）
int16_t*  mShtTemp[SLOTS];       // int16_t リングバッファ（HSS_GRD_FAC_TEMP倍）
int       idx_order[SLOTS];      // idx_order[0]=書込先（最新）、idx_order[1]=1ステップ前(t0)
int       nParallelSlices;       // 補間並列時刻数
int       nShtTempSlots;         // スロット数（= nParallelSlices + 1）
```

**等温線境界マスク（ビットパック）**

```c
uint16_t* mBoundPacked[SLOTS];   // bit c = 閾値c の境界マスク（HSS_GRD_TYPE_B16として出力）
int       nCntThresholds;        // 閾値数（最大16）
int16_t   nvThresholds[16];      // 閾値リスト（設定値℃にHSS_GRD_FAC_TEMPを乗じた値）
```

**エリア別統計集計**

```c
int          nAreas;                   // エリア数
AreaDef      areaDefs[AREA_DEF_MAX];   // エリア定義配列
uint32_t*    pAreaMask;                // エリアマスク（bit[a]=1→エリアa所属）
AreaMonthStat vAreaStat[OBTAIN_AREA];  // 月次統計バッファ（[0]=全領域）
```

**確率マップ**

```c
uint32_t* pCountsA;    // 全期間カウント(lo) [total_allocated_pix × nThresholds]
uint16_t* pCountsM;    // 月次カウント(lo)
uint32_t* pCountsA_hi; // 全期間カウント(hi) DIR=2時
uint16_t* pCountsM_hi; // 月次カウント(hi)
```

**気温統計（時別/日別/月別/全期間）**

各集計期間で `fTempMin/fTempMax/dTempSum/nTempCnt[OBTAIN_AREA]` を保持。

**気温ヒストグラム**

```c
uint32_t* pTempHistBuf;   // 時別ヒストグラム 0.1℃精度（中央値導出用）
uint32_t* pTHistHourBuf;  // 時別ヒストグラム 0.5℃精度（出力用）
uint32_t* pTHistDayBuf;   // 日別累積ヒストグラム
uint32_t* pTHistMthBuf;   // 月別累積ヒストグラム
uint32_t* pTHistAllBuf;   // 全期間累積ヒストグラム
```

#### `DataFlow`（ベクトル解析結果）
```c
float* pFlowU;      // 経度方向速度 [pixel/step]
float* pFlowV;      // 緯度方向速度 [pixel/step]
float* pMagnitude;  // |∇T|
float* pDiverg;     // 発散: ∂u/∂x + ∂v/∂y
float* pVorticity;  // 渦度: ∂v/∂x - ∂u/∂y
```

---

### 2-3. 評価・制御系

#### `DateTime`
年月日時を保持し、`DateTimeAdvance()` で1時間ずつ繰り上がり処理（月末・年末対応）。`DateTimeYearMonthEqual()` でロードの要否を判定。

#### `ConfigCore`（INI読込可能な実行設定）
旧 `AppConfig` に相当。`AppConfigDefault()` でデフォルト値セット後、`LoadConfig()` でINIから上書き。主要セクション：[MAIN] [PATH] [RUN] [EXT_ANAL] [INTERP] [PARAM] [OUTPUT] [OUTPUT_VECTOR] [OUTPUT_GEOTIFF] [OUTPUT_PNG] [OUTPUT_GRD] [TEMP_STAT] [TEMP_HIST] [THRESHOLD_ANAL] [THRESHOLD] [AREA] [POI]

#### `ConfigRun`（実行動作定義）
旧 `RunConfig` に相当。開始日時・ステップ数・年範囲・対象月リスト等。

#### `ConfigRain`（雨量相関解析用設定）
開始日時・ステップ数・収束判定閾値・降水判定閾値。

#### `AreaDef`（エリア定義）
矩形指定（lat_s/lat_n/lon_w/lon_e）またはSHPポリゴン指定（szShpPath/szFieldName/szFieldValue）。

#### `AreaMonthStat`（エリア別月次統計）
閾値別の存在時間数・面積積算・最大面積・出現/消滅回数をT<閾値（lo）とT≥閾値（hi）の両方向で管理。

#### `DivergRainResult`（1ステップ評価）
収束域での降水予測結果（TP/FP/FN/TN）+ 相関係数 + 降水あり/なし観測点での発散統計。

#### `DivergRainAccum`（累積評価）
複数ステップの結果を集計。HR/FAR/CSI/HSS算出に使用。

#### `FlowSnapshot`（ParamScan用スナップショット）
```c
float*   pDiverg;  // 発散場（全ピクセル）
float*   fRains;   // 観測点降水量（alignedNum点）
int      stCount;  // 有効観測点数
DateTime dt;       // この時刻
```
Stage1で保存しStage2で閾値を変えて再評価する2段階構造。

---

## 3. ユーティリティクラス（HsCommon.h）

### `CHsString`
`std::string` ベースのVC6 `CString` 互換ラッパ。`Format()`, `Mid()`, `Find()`, `TrimLeft/Right()` 等を提供。Clang on Windowsでの互換性維持目的。

### `CHsIniFile`
Windowsスタイルのセクション/キー/値形式のINIファイルを `std::map` にキャッシュし、型付きGetter（`GetInt32`, `GetFloat`, `GetString`, `GetDouble`）を提供。コメント行（`;`, `#`）はスキップ。

### `CHsLog`（グローバルインスタンス `g_LF`）
標準出力とログファイルへの同時出力を担うクラス。
- `OpenMakPath()` で実行ファイル名+日時のログファイルを自動生成（`log/`サブディレクトリ）
- 16KBの内部バッファ（`_IOFBF`）でIO負荷を軽減
- `PRINTF` マクロで `g_LF.Print()` を呼び出す

---

## 4. 処理フロー（main関数）

```
起動
│
├─ CHsLog::OpenMakPath()          # ログ開始
├─ GDAL / Blosc2 初期化
├─ AppConfigDefault()             # デフォルト設定
├─ LoadConfig(argv[1], CFG)       # 第1引数のcfg読み込み
│
├─ InitDEM(CFG.szDemPath)         # GeoTIFF DEM読み込み
└─ CalcCoastDist(DEM)             # 海岸距離計算（Meijster法）
    │
    └─ 設定ファイルループ（argv[1..N]）
        │
        ├─ DEMパスが変わった場合: FreeDEM → InitDEM → CalcCoastDist
        ├─ InitAnal(DEM, CFG, ANAL)
        ├─ InitFlow(DEM, FL)
        │
        └─ switch(CFG.nRunMode)
            ├─ CONTINUOUS(0):     RunContinuous()    # 通常連続実行
            ├─ MULTI_YEAR(1):     RunMultiYear()     # 同月複数年
            ├─ PARAM_SCAN(2):     RunParamScan()     # 閾値パラメータスキャン
            ├─ COMPARE(3):        ComparInterpolationMethod()  # 補間法比較
            ├─ VALUE_SURVEY(4):   RunValueSurvey()   # 値調査
            ├─ EVAL_INTERP(5):    RunEvalInterp()    # 補間法精度評価（LOOCV）
            ├─ OPTIM_PARAMS(6):   RunOptimParams()   # 補間パラメータ最適化スキャン
            └─ AREA_CHAR(7):      RunAreaChar()      # 解析エリア特性評価
```

**複数cfgファイルの連続実行**: DEMが共通であれば再初期化なしで続けて実行。DEMパスが変わった場合のみ再ロード。

---

## 5. 主要関数群の詳細

### 5-1. 初期化・解放系

#### `InitDEM()` / `FreeDEM()`
GDALでGeoTIFFを読み込み、`pElevation`・`pSeaMask`・`pCoastDist` を確保。`nAlignedWidth` はAVX-512の64バイト境界（16float）に合わせて切り上げ。

#### `InitAnal()` / `FreeAnal()`
`DataAnal` の全ポインタをAVX-512アライメント確保（`_mm_malloc`）。`mShtTemp[nShtTempSlots]` はリングバッファとして動作。CFG.nAreaStatMode・nProveMode等の設定に応じてオプションバッファ（pCountsM系・pTempHistBuf系等）の確保を制御。

#### `InitFlow()` / `FreeFlow()`
`DataFlow` の5配列（pFlowU/V/Magnitude/Diverg/Vorticity）を確保。

#### `InitTPS()`
`DataAnal` の `tps_w`, `tps_stLon`, `tps_stLat` を確保（TPS使用時のみ）。

#### `FreeAMeDAS()`
月替わり時に `_mm_free` で観測データを全解放。

---

### 5-2. データ読み込み系

#### `LoadAMeDAS()`
```
1. STNファイル（観測所マスタ）読み込み → 緯度・経度・標高をfloatに変換
2. IDXファイル（観測所番号リスト）と照合
3. 各観測所のAMDファイル（日別データ）を読み込み
4. ValidateStationElevation() でDEM標高と観測所標高を比較・補正
5. pixX/pixY にDEMグリッド上の対応ピクセル座標を設定
```

#### `SetCalcTemp()` / `SetCalcRain()`
指定日時の `AMeDAS_TMDT` から `DT.fTemps[]` / `DT.fRains[]` を設定。欠測（32767）は `NAN` に変換。

---

### 5-3. 前処理系

#### `CalcCoastDist(DataDEM&)`
**Meijster法**による正確なユークリッド距離変換（EDT）。
```
STEP1: 海ピクセル=0、陸ピクセル=INF で初期化
STEP2: x方向パス（左→右・右→左）で1D距離²を計算
STEP3: y方向Meijster法（放物線エンベロープ構築 → 最小距離²確定）
STEP4: sqrt → pCoastDist[m] に格納
```
STEP2/3/4を `#pragma omp parallel for` で並列化。

#### `CalcInterpParams()`
観測点の実配置からBarnes法のkappa値とRBF sigma値を自動計算。
全観測点ペアO(N²)で最近傍距離を求め、平均最近傍距離 d_nn を算出。

#### `BuildNbrTable()`
全陸地ピクセルに対して近い観測所最大16点を選択し、SoA形式の近傍点テーブルを構築。月替わり時に再計算。

#### `BuildWeightTable()`
近傍点テーブルから事前計算済みの重み（`pNbrW`/`pNbrAltW`）を更新。補間手法変更時・月替わり時に実行。

---

### 5-4. 補間手法

`InterpConfig.method` で切り替え可能。

| 定数 | 値 | 手法 |
|---|---|---|
| `INTERP_IDW` | 0 | 逆距離加重（p=2固定） |
| `INTERP_RBF_GAUSS` | 1 | ガウシアンRBF |
| `INTERP_SHEPARD` | 2 | Shepard修正IDW（べき乗可変） |
| `INTERP_KRIGING_S` | 3 | 簡易クリギング（球形バリオグラム） |
| `INTERP_BARNES` | 4 | Barnes解析（推奨） |
| `INTERP_BARNES_H` | 5 | Barnes高精度版（Schraudolph近似exp） |
| `INTERP_BARNES2PASS` | 6 | Barnes 2Passスキャン |
| `INTERP_TPS` | 7 | Thin Plate Spline（計算コスト大） |

**海岸距離補正**（`GetCoastFactor()`）: ピクセルと観測点の海岸距離差に基づいてIDW重みを減衰。`coast_L=0` で補正なし。

**気温減率補正**（`GetLapseRate()`）: `lapse_rate_table[緯度帯5][季節4]` から減率を取得し、標高差による気温差を補正。

**Schraudolph近似exp**（`INTERP_BARNES_H`使用時）:
```c
union { float f; int i; } u;
u.i = (int)(12102203.0f * x + 1065353216.0f);
```
float標準expfより高速な整数演算近似。

---

### 5-5. 補間実行系

#### `TemparatureInterpolation()`
時間ループの1ステップとして以下を実行。
```
1. SetCalcTemp(DT, day, hour)
2. CalcTPSCoeff()（TPS使用時のみ）
3. InterpolationMain() または InterpolationMainMulti()
   - 全陸地ピクセルを並列処理（#pragma omp parallel for）
   - 近傍点テーブルから重みを取得し各手法で補間
   - 結果を mShtTemp[idx_order[0]] に格納（int16_t）
```

#### `RotateTimeSlice()`
`idx_order[]` を1つずつシフト（リングバッファローテーション）。`idx_order[0]` が常に最新書込先。

---

### 5-6. ベクトル解析系（OpticalFlow.cpp）

#### `CalcGradientFlow()`（2点局所法）
```
∇T を各ピクセルのSobel勾配で計算
T_t0 = mShtTemp[idx_order[1]]（1ステップ前）
T_t1 = mShtTemp[idx_order[0]]（最新）
dT/dt = (T_t1 - T_t0) / dt

FlowU = -dT/dt × Gx / (Gx² + Gy² + ε)
FlowV = -dT/dt × Gy / (Gx² + Gy² + ε)
```
AVX-512 SIMD（`_mm512_fmadd_ps` 等）で加速。

#### `Calc3PointGradientFlow()`
中心差分法（TIME_SLICES > 2 の場合）。現在は無効化（評価上2点法が優勢）。

#### `CalcHornSchunckFlow()`
Horn-Schunckオプティカルフロー（HS_ALPHA_SQ=0.1, HS_ITER=32）。Sobel勾配 + 反復ラプラシアン平滑化。現在は無効化。

#### 発散・渦度の計算
```
∂u/∂x, ∂v/∂y の中心差分 → pDiverg（発散）
∂v/∂x - ∂u/∂y の中心差分 → pVorticity（渦度）
```

---

### 5-7. 評価系

#### `EvalInterp()`
**LOOCV（Leave-One-Out Cross Validation）**: 各観測点を1つ除外して補間し、実測値との誤差を評価。
- `loocv_RMSE`, `loocv_MAE`, `loocv_maxErr`
- 緯度帯別（南/<33°N、中/33-40°N、北/>40°N）集計
- ラプラシアン（空間的滑らかさ指標）も計算

#### `EvalDivergRain()`
観測点位置の発散値と降水量を照合。
```
発散 < fDivergThreshold → 収束（降水予測あり）
降水量 >= fRainThreshold → 降水あり
→ TP/FP/FN/TN を DivergRainResult に格納
```

#### HSS計算（`PrintDivergRainAccum()`）
```
Heidke Skill Score (HSS):
  E = (TP+FN)(TP+FP)/n² + (TN+FP)(TN+FN)/n²
  H = (TP+TN)/n
  HSS = (H-E)/(1-E)
```
HSS > 0 でランダム予測より優れる。

#### `CalcBoundaryMaskAndCollect()`
気温場から各閾値の境界マスクを生成し、エリア別統計集計を同時実行。x方向・y方向いずれかで閾値をまたぐピクセルを境界(1)とする。

---

### 5-8. 実行モード系

#### RunOneStep（時間ループの中核）
```
for step in range(nSteps):
    1. 月替わり判定 → LoadAMeDAS / BuildNbrTable / BuildWeightTable
    2. RotateTimeSlice
    3. TemparatureInterpolation（補間 + int16化）
    4. CalcGradientFlow（step≥1から開始）
    5. CalcBoundaryMaskAndCollect（エリア統計集計）
    6. OutputBoundCentroid（重心計算、nCentroidMode>0時）
    7. [GeoTIFF出力] step % nThinningGeoTIFF == 0
    8. [GRD出力]    step % nThinningGRD == 0
    9. [PNG出力]    step % nThinningGeoTIFF == 0
    10. DateTimeAdvance
```

#### `RunContinuous()`
CFGの開始日時から `nRunCount` ステップを1回実行。

#### `RunMultiYear()`
`nYearStart`〜`nYearEnd` の同月（`nTargetMonthList[]`）をループし、年別統計を1行ずつ出力。終了後に全年集計を出力。

#### `RunParamScan()`
**2段階方式**:
- **Stage1**: 全ステップを実行し、各ステップの発散値（観測点位置のみ）と降水量をスナップショット保存
- **Stage2**: 発散閾値 × 降水閾値 × 時間ラグの全組み合わせでスナップショットを参照し高速再評価

#### `RunEvalInterp()`（MODE=EVAL_INTERP）
複数の補間手法を指定月・年に対してLOOCV評価し、精度比較結果をCSVとして出力する。
- `LoadAMeDAS / BuildNbrTable` の後、手法を切り替えながら `EvalInterp()` を各手法で呼び出す
- RMSE / MAE / 最大誤差 / 緯度帯別誤差 / ラプラシアン（空間的滑らかさ）を集計
- 出力: `*_eval_summary.txt`（手法×指標のTSV）、`*_eval_optimize.txt`（各手法のパラメータ最適値一覧）

#### `RunOptimParams()`（MODE=OPTIM_PARAMS）
Barnes / Barnes2Pass / RBF-Gauss / Shepard の各補間パラメータをグリッドスキャンし、LOOCV最適値を導出する。
- **スキャン対象パラメータ**:
  - Barnes/Barnes2Pass: `kappa_fac` ∈ {1.0, 2.0, 3.0, 5.052, 8.0, 12.0, 15.0, 20.0}
  - Barnes2Pass: `gamma` ∈ {0.1, 0.2, 0.3, 0.5, 0.7}
  - RBF-Gauss: `rbf_sigma` ∈ {10000, 20000, 30000, 50000, 80000, 100000, 150000}（m）
  - Shepard: `idw_power` ∈ {1.0, 1.5, 2.0, 2.5, 3.0}
- 各組み合わせでLOOCV RMSEを計算し、手法ごとにベストパラメータを記録
- **境界警告**: 最優パラメータがスキャン範囲の上下端に一致した場合、コンソールとファイルに警告を出力（範囲拡大の要否を通知）
- 出力: `*_optim_result.txt`（手法別最適パラメータ・LOOCV RMSE）

#### `RunAreaChar()`（MODE=AREA_CHAR）
全陸地ピクセルの補間難易度・地形特性を定量化し、GeoTIFFおよびCSVとして出力する。単発調査用モード。
- `LoadAMeDAS / CalcInterpParams / BuildNbrTable` を実行（`BuildWeightTable` は不要）
- 各陸地ピクセルに対して以下の指標を並列（OpenMP）計算:

| 指標 | 説明 |
|---|---|
| `d_nn` | 最近傍観測点距離 √(pNbrD2[base+0]) [m] |
| `mean_d` | 近傍NBR_MAX点の平均距離 [m] |
| `cv_d` | 近傍距離の変動係数（std/mean）（観測点分布の均一性指標） |
| `alt_std` | 近傍観測点との標高差の標準偏差 [m]（地形不均一性・補間難易度指標） |
| `alt_range` | 近傍観測点との標高差の最大−最小 [m] |
| `nbr_n` | 有効近傍点数（pNbrN） |

- **GeoTIFF出力**: 指標ごとに `*_area_char_d_nn.tif`, `*_area_char_mean_d.tif`, `*_area_char_cv_d.tif`, `*_area_char_alt_std.tif`, `*_area_char_alt_range.tif`, `*_area_char_nbr_n.tif`（float型, NODATA=-9999）
- **サマリーCSV出力** (`*_area_char_summary.csv`): 標高帯4区分（全体/低地<200m/中高地200-800m/高山≥800m）× 6指標 × 統計値（N/mean/std/min/max/p10/p50/p90）
- **観測点別CSV出力** (`*_area_char_stations.csv`): 各観測点ピクセル座標でサンプリングした6指標値（補間困難点の特定・OPTIM_PARAMSとの相互参照用）

**解析上の意義**: `alt_std` が高いピクセルは平地-山岳遷移帯（群馬～宇都宮、飛騨、福島市付近、阿蘇周辺等）に集中し、近傍点の地形多様性が高く補間誤差が大きい傾向がある。OPTIM_PARAMSで得られた最適`kappa_fac`とこのエリア特性の関係を把握することで、解析エリア選定の根拠を定量化できる。

---

## 6. 出力ファイル系

### 6-1. GeoTIFF（WriteGeoTIFF.h）
- `nOutputGeoTIFF2=1`: 気温・フロー結果（TEMP, FLOWU, FLOWV, MAGNITUDE, DIVERG）
- `nOutputGeoTIFF3=1`: 気温閾値確率分析結果
- `nOutputGeoTIFF4=1`: 気温境界面マスク（月次代表）
- `nThinningGeoTIFF` で時間間引き

### 6-2. GRD バイナリ（WriteGRD.h）
独自フォーマット。BLOSC2（ZSTD lv1、8スレッド）で圧縮。

**ファイル構造**:
```
HSS_GRD_HEADER（3519 byte）
  + HSS_GRD_CHINFO × nChNum（80 byte × CH数）
  + 各チャンネルの圧縮データ
```

**主要チャンネル構成**:
| タグ名 | 内容 | 型 | スケール |
|---|---|---|---|
| TEMP | 気温 | int16 | HSS_GRD_FAC_TEMP（デフォルト100=0.01℃精度） |
| FLOWU | 経度方向速度 | float | — |
| FLOWV | 緯度方向速度 | float | — |
| MAGNITUDE | 速度スカラー | float | — |
| DIVERG | 発散 | float | — |
| VORTICITY | 渦度 | float | — |
| BOUND | 等温線境界マスク | B16（ビットパック） | bit c = 閾値c |

`nThinningGRD` で出力頻度を制御。

### 6-3. PNG出力（AnalAMeDAS.cpp）
- `nOutputPNG2=1`: 気温分布画像 + 等温線別PNG（`nOutputPNG4=1` 時）
- `nOutputPNG2=2`: マージ出力（気温分布+等温線合成）
- `nOutputPNG3=1`: 閾値確率分析結果画像

### 6-4. テキスト統計出力（TSV形式）

| ファイル名パターン | 内容 |
|---|---|
| `*_threshold_stat_hourly.txt` | 閾値別エリア別ピクセル数・比率（時別） |
| `*_threshold_lo_stat_monthly.txt` | 閾値以下の月別統計（面積・出現/消滅回数等） |
| `*_threshold_hi_stat_monthly.txt` | 閾値以上の月別統計 |
| `*_temp_stat_hourly.txt` | エリア別気温統計（平均/最低/最高/中央値、時別） |
| `*_temp_stat_daily.txt` | エリア別気温統計（日別） |
| `*_temp_stat_monthly.txt` | エリア別気温統計（月別） |
| `*_temp_stat_all.txt` | エリア別気温統計（全期間） |
| `*_temp_hist_hourly.txt` | エリア別気温ヒストグラム（時別） |
| `*_temp_hist_monthly.txt` | エリア別気温ヒストグラム（月別） |
| `*_temp_hist_all.txt` | エリア別気温ヒストグラム（全期間） |
| `*_bound_cen.txt` | 等温線境界面重心（時別） |
| `*_poi_centroid.txt` | POI別最近傍境界面重心追跡 |
| `*_area_land_pixels.txt` | エリア別陸地ピクセル数 |
| `*_eval_summary.txt` | 補間手法別LOOCV精度比較（RMSE/MAE/最大誤差/ラプラシアン等）[MODE=EVAL_INTERP] |
| `*_eval_optimize.txt` | 各補間手法の最適パラメータ一覧 [MODE=EVAL_INTERP] |
| `*_optim_result.txt` | パラメータスキャン結果・手法別最適値・境界警告 [MODE=OPTIM_PARAMS] |
| `*_area_char_d_nn.tif` | 最近傍観測点距離分布GeoTIFF [MODE=AREA_CHAR] |
| `*_area_char_mean_d.tif` | 近傍点平均距離分布GeoTIFF [MODE=AREA_CHAR] |
| `*_area_char_cv_d.tif` | 近傍距離変動係数分布GeoTIFF [MODE=AREA_CHAR] |
| `*_area_char_alt_std.tif` | 近傍点標高差標準偏差分布GeoTIFF（補間難易度指標）[MODE=AREA_CHAR] |
| `*_area_char_alt_range.tif` | 近傍点標高差レンジ分布GeoTIFF [MODE=AREA_CHAR] |
| `*_area_char_nbr_n.tif` | 有効近傍点数分布GeoTIFF [MODE=AREA_CHAR] |
| `*_area_char_summary.csv` | 標高帯別×指標別統計サマリー [MODE=AREA_CHAR] |
| `*_area_char_stations.txt` | 観測点ピクセル位置でサンプリングした各指標値 [MODE=AREA_CHAR] |

---

## 7. 設定ファイル（CFGファイル）仕様

```
ini
[MAIN]
TYPE         = main                          # 固定(設定ファイルの適合性確認用)

[PATH]
DEM_PATH     = U:\DEM\japan_250m.tif         # DEM GeoTIFFファイルパス(末尾に「\」が必要)
RESULT_PATH  = U:\PRG3\AMT_RES\20260425A\    # 結果出力先ディレクトリ
FILE_KEY     = japan                         # 画像出力ファイルの先頭文字列
FILE_KEY_ST  = RES_JPN                       # 統計結果テキストファイルの先頭文字列
AMEDAS_PATH  = U:\AMeDAS\                    # AMeDASデータディレクトリ（STN/IDX/AMDファイルが存在するパス）

[RUN]
# MODE : CONTINUOUS / MULTI_YEAR / PARAM_SCAN / COMPARE / VALUE_SURVEY / EVAL_INTERP / OPTIM_PARAMS / AREA_CHAR
MODE         = MULTI_YEAR
CONT_START   = 2000/2/15/0                   # CONTINUOUSモードの実行開始日時（年/月/日/時）
CONT_DUR     = 2/1                           # CONTINUOUSモードの実行期間（時/日）、日設定を優先
MULTI_YEAR   = 1999/2001                     # MULTI_YEARモードの対象年範囲（開始年/終了年）
MULTI_MONTH  = 1/2/3/12                      # MULTI_YEARモードの対象月リスト（1-12のスラッシュ区切り）
PARALLEL_SLICES = 4                          #複数時刻同時計算
USE_NBR_FLOAT   = 3                          #近傍点関連をfloatで処理する（0:すべて整数、1:SRCのみfloat、2:Wのみfloat、3：全てfloat）

[EXT_ANAL]
PROB_MODE    = 3                             # 閾値以上/以下の確率分析(0:なし、1：以下、2：以上、3:両方)
CENTROID     = 1                             # 重心分析(0:なし、1:実施)
AREA_STAT    = 1                             # 統計分析(0:なし、1:実施)
AREA_HIST    = 1                             # ヒストグラム分析(0:なし、1:実施)
PROB_ANNUAL  = 0                             # 全期間の確率分析を行う（0:行わない、1:行う）
OPT_FLOW     = 0                             # 時系列オプティカルフロー解析（0:無効、1:勾配法ベクトル解、2:Horn-Schunck法）

[INTERP]
METHOD          = 4                          # 気温空間補間の補間方法 (0=IDW, 1=RBF_GAUSS, 2=SHEPARD, 3=KRIGING_S, 4=BARNES, 5=BARNES_H, 6=BARNES2PASS, 7=TPS)
COAST_L         = 20000.0                   # 海岸距離補正スケール [m]（0で補正なし）
IDW_POWER       = 1.5                        # IDW/Shepard べき乗値（デフォルト 1.5）
RBF_SIGMA       = 80000.0                    # RBF-Gauss バンド幅 [m]（デフォルト 80000）
BARNES_KAPPA_FAC = 2.0                       # Barnes: kappa = d_nn² × この値（デフォルト 2.0、文献値 5.052）
BARNES_GAMMA    = 0.3                        # Barnes 2Pass スキャン間収縮率（デフォルト 0.3）

[PARAM]
DIVERG_RAIN      = 0
DIVERG_THRESHOLD = -2.200
RAIN_THRESHOLD   = 1.5

[OUTPUT]
#ログ出力関連
DEBUG        = 0
COMMENT      = 1
ANAL0        = 0
ANAL1        = 0
ANAL2        = 0

[TEMP_STAT]
HOUR_OUT     = 0                             # 気温統計解析の時別値を出力（0:なし、1:あり）
DAY_OUT      = 1                             # 気温統計解析の日別値を出力（0:なし、1:あり）
MONTH_OUT    = 1                             # 気温統計解析の月別値を出力（0:なし、1:あり）
TOTAL_OUT    = 0                             # 気温統計解析の全体を出力（0:なし、1:あり）
NUM_DEC      = 1                             # 気温統計解析の小数点以下桁数（1 or 2）
THINNING     = 6                             # 気温統計解析の時別値の出力間隔

[TEMP_HIST]
# MIN=-20、MAX=20、BIN=1とした場合は
# -20：T<-19.5
# -19:-19.5≦T<-18.5
#   ...
# 19:18.5≦T<19.5
# 20:19.5≦T
BIN          = 1                             # 気温ヒストグラム解析のBIN幅 ℃
MIN          = -35                           # 気温ヒストグラム解析の最低気温 ℃
MAX          = 20                            # 気温ヒストグラム解析の最高気温 ℃
HOUR_OUT     = 1                             # 気温ヒストグラム解析の時別値を出力（0:なし、1:あり）
DAY_OUT      = 1                             # 気温ヒストグラム解析の日別値を出力（0:なし、1:あり）
MONTH_OUT    = 1                             # 気温ヒストグラム解析の月別値を出力（0:なし、1:あり）
TOTAL_OUT    = 0                             # 気温ヒストグラム解析の全体を出力（0:なし、1:あり）
NUM_DEC      = 1                             # 気温ヒストグラム解析の小数点以下桁数（1 or 2）
FG_CNT       = 0                             # 気温ヒストグラム解析で該当メッシュ数を出力（0:なし、1:あり）
THINNING     = 6                             # 気温ヒストグラム解析の時別値の出力間隔

[THRESHOLD_ANAL]
HOUR_OUT     = 1                             # 気温閾値解析の時別値（閾値以下）を出力（0:なし、1:あり）
LO_OUT       = 1                             # 気温閾値解析の閾値以下統計値を出力（0:なし、1:あり）
HI_OUT       = 1                             # 気温閾値解析の閾値以上統計値を出力（0:なし、1:あり）
NUM_DEC      = 1                             # 気温閾値解析の小数点以下桁数（1 or 2）
FG_CNT       = 0                             # 気温閾値解析で該当メッシュ数を出力（0:なし、1:あり）
THINNING     = 6                             # 気温閾値解析の時別値の出力間隔

[OUTPUT_VECTOR]
FLOW         = 1                             # オプティカルフロー解析のベクトルデータをGPKGファイルとして出力（0:なし、1:あり）
THINNING     = 1                             # オプティカルフロー解析のベクトルデータの出力時間間隔
THINNING_XY = 20                             # オプティカルフロー解析のベクトルデータの空間間引き間隔（20⇒20メッシュ毎）

[OUTPUT_GEOTIFF]
GEOTIFF1     = 0                             # 陸マスクなど確認用GeoTIFF出力（0:無効、1:有効）
TEMP         = 0                             # 気温解析結果GeoTIFF出力（0:無効、1:有効）
FLOW         = 0                             # フロー結果GeoTIFF出力（0:無効、1:有効）
GEOTIFF3     = 0                             # 気温閾値確率分析結果GeoTIFF出力（0:無効、1:有効）
GEOTIFF4     = 0                             # 気温境界面マスクGeoTIFF出力（0:無効、1:有効）
THINNING     = 1                             # GeoTIFF出力の出力時間間隔
THINNING_XY  = 1                             # GeoTIFF出力の空間間引き間隔

[OUTPUT_PNG]
PNG2         = 0                             # 気温補間PNG出力（0:無効、1:有効、2：3:マージ出力）、PNG=2の場合はPNG4=0となる
TEMP_MIN     = -30                           # PNG2の気温範囲(下限)
TEMP_MAX     = 20                            # PNG2の気温範囲(上限)
PNG3         = 0                             # 気温閾値確率分析結果画像出力（0:無効、1:有効）
PNG4         = 0                             # 気温境界面PNG出力（全閾値重ね合わせ）（0:無効、1:有効）

[OUTPUT_GRD]
GRD          = 0                             # GRDバイナリ出力（0:無効、1:有効）
THINNING     = 1                             # バイナリ出力時間間隔
CH_NUM  = 1                                  # 出力チャンネル
CH_TAG1 = TEMP                               # 番号は1からNUM

[PNG_OPT]                                    # PNG画像出力の年月日時刻文字
FONT_PATH    = C:/Windows/Fonts/arial.ttf    # フォントファイルパス
FONT_HEIGHT  = 64                            # 文字高さ(ピクセル)
IMPOSE_X     = 10                            # 描画位置(左上を基準ピクセル)
IMPOSE_Y     = 10                            # 描画位置(左上を基準ピクセル)

[LEGEND_TP]                                  # 気温分布PNG画像出力の凡例文字設定
FONT_PATH    = C:/Windows/Fonts/arial.ttf    # フォントファイルパス
FONT_HEIGHT  = 48                            # 文字高さ(ピクセル)
IMPOSE_X     = 1200                          # 描画位置(左上を基準ピクセル)
IMPOSE_Y     = 10                            # 描画位置(左上を基準ピクセル)
BAR_HIGHT    = 40                            # バーの高さ(ピクセル)
BAR_BIN_WIDTH  = 2                           # バーの1INDEXの幅(ピクセル)

[LEGEND_PR]                                  # 閾値確率分析PNG画像出力の凡例文字設定
FONT_PATH    = C:/Windows/Fonts/arial.ttf    # フォントファイルパス
FONT_HEIGHT  = 48                            # 文字高さ(ピクセル)
IMPOSE_X     = 1200                          # 描画位置(左上を基準ピクセル)
IMPOSE_Y     = 10                            # 描画位置(左上を基準ピクセル)
BAR_HIGHT    = 40                            # バーの高さ(ピクセル)
BAR_BIN_WIDTH  = 2                           # バーの1INDEXの幅(ピクセル)

[LEGEND_TC]                                  # 等温線PNG画像出力の凡例文字設定
FONT_PATH    = C:/Windows/Fonts/arial.ttf    # フォントファイルパス
FONT_HEIGHT  = 40                            # 文字高さ(ピクセル)
IMPOSE_X     = 3500                          # 描画位置(左上を基準ピクセル)
IMPOSE_Y     = 10                            # 描画位置(左上を基準ピクセル)

[THRESHOLD]                                  # 閾値設定
VALUE    = -30,-25,-20,-15,-10,-5,0,5,10,15

[AREA]                                       # エリア解析条件
NUM    = 32                                  # エリア解析設定数(最大32)
#番号は0からNUM-1
AREA_0 = 札幌周辺,42.8,43.4,141.0,141.8                                   # 矩形領域の設定
AREA_1 = 札幌市,shp:U:\GIS\AdmDiv\N03-20250101.shp,N03_007,01101-01110    # ポリゴン設定

[POI_CENTROID]                               # 指定点重心解析設定
COUNT        = 3                             # 設定数
SEARCH_R     = 50                            # 検索範囲メッシュ数
POI_1        = 札幌,43.060000,141.326667
POI_2        = 名寄,44.370000,142.455000
POI_3        = 陸別,43.468333,143.738333

[LAPSE]
# 緯度ゾーン境界（LAT_0～LAT_3の昇順4値で5ゾーン分割）
# デフォルト: 33/40/999/999 → 実質3ゾーン
LAT   = 33.0/40.0/999.0/999.0
# 季節開始月（デフォルト: 春=3月、夏=6月、秋=9月、冬=12月）
SPRING_MONTH = 3
SUMMER_MONTH = 6
AUTUMN_MONTH = 9
# 高度補正値 [℃/m]（負値=高度上昇で低温、全ゼロ=高度補正なし）
# キー: WINTER/SPRING/SUMMER/AUTUMN + _1～5（_1=最南ゾーン）
# ゾーン1: lat < 33.0（南帯）
VALUE1   = -0.0055/-0.0060/-0.0070/-0.0062
# ゾーン2: 33.0 <= lat < 40.0（中帯）
VALUE2   = -0.0050/-0.0060/-0.0065/-0.0060
# ゾーン3: 40.0 <= lat（北帯・北海道）
VALUE3   = -0.0045/-0.0055/-0.0060/-0.0055
# ゾーン4,5: LAT_2,LAT_3が999.0のため未使用（将来の拡張用）
VALUE4   = -0.0045/-0.0055/-0.0060/-0.0055
VALUE5   = -0.0045/-0.0055/-0.0060/-0.0055
```

---

## 8. パフォーマンス特性
### 評価DEM条件
エリアとメッシュ
|  |  | 北海道 | 日本 | 群馬～長野～富山 |
| :--- | :---: | ---: | ---: | ---: |
| メッシュサイズ | m | 250 | 250 | 50 |
| 東西サイズ | px | 2,208 | 7,329 | 3,999 |
| 南北サイズ | px | 2,075 | 12,057 | 2,887 |
| 総ピクセル数 | px | 4,581,600 | 88,365,753 | 11,545,113 |
| 陸地ピクセル数 | px | 1,238,181 | 5,955,323 | 11,545,113 |
| 陸地比率 | % | 27.0 | 6.7 | 100.0 |

### 実行環境①
AMD 7940HS
メモリ:16GB（DDR5-4800 x 1）
Greekom A7max(オリジナルまま)
　※上記条件はスワップの発生なく処理可能
　　　※Windows初期設定でメモリ消費を抑えていれば

### 処理時間@実行環境①
|  |  | 日本全国 | 北海道 | 北海道 | 群馬～長野～富山 | 群馬～長野～富山 | 群馬～長野～富山 | 九州 | 九州 | 九州 |
|:---|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| メッシュサイズ | m | 250 | 250 | 100 | 250 | 100 | 50 | 250 | 100 | 50 |
| 東西サイズ | px | 8,295 | 2,208 | 5,521 | 800 | 2,000 | 3,999 | 1,340 | 3,350 | 6,700 |
| 南北サイズ | px | 9,763 | 2,075 | 5,185 | 577 | 1,443 | 2,887 | 1,360 | 3,401 | 6,801 |
| 総ピクセル数 | px | 80,984,085 | 4,581,600 | 28,626,385 | 461,600 | 2,886,000 | 11,545,113 | 1,822,400 | 11,393,350 | 45,566,700 |
| 陸地ピクセル数 | px | 6,193,139 | 1,250,946 | 7,784,422 | 451,827 | 2,824,156 | 11,296,741 | 653,257 | 4,013,483 | 15,948,564 |
| 陸地比率 | % | 7.6 | 27.3 | 27.2 | 97.9 | 97.9 | 97.8 | 35.8 | 35.2 | 35.0 |
| 条件 |  | 計算のみ | 計算のみ | 計算のみ | 計算のみ | 計算のみ | 計算のみ | 計算のみ | 計算のみ | 計算のみ |
| USE_NBR_FLOAT  |  | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 |
| メモリ消費（主要、参考） |  |  |  |  |  |  |  |  |  |  |
| pNbr* (mode=3) | MB | 2108 | 425 | 2650 | 153 | 961 | 3846 | 222 | 1366 | 5429 |
| ShtTemp+Bound x5 | MB | 1546 | 87 | 547 | 8 | 55 | 220 | 34 | 217 | 869 |
| ProbMap(N=16,A=1) | MB | 3092 | 174 | 1094 | 17 | 110 | 440 | 62 | 392 | 1565 |
| Total | MB | 8082 | 764 | 4785 | 185 | 1182 | 4736 | 348 | 2174 | 8661 |
| 3か月実行時間 | sec. | 127.1 | 30.4 | 121.9 | 23.1 | 52.1 | 159.8 | 27.3 | 73.7 | 246.1 |
| エリア条件読込 | sec. | 19.2 | 11.1 | 13.5 | 16.6 | 17.4 | 22.9 | 17.5 | 20.3 | 32.6 |
| 平均処理時間 | sec./月 | 36.0 | 6.5 | 36.2 | 2.2 | 11.6 | 45.6 | 3.3 | 17.8 | 71.2 |
| 近傍点計算 | sec./月 | 2.3 | 0.3 | 1.6 | 0.2 | 1.0 | 4.1 | 0.4 | 2.2 | 9.0 |
| Analyze1Step | msec. | 40.77 | 7.87 | 43.97 | 2.53 | 13.61 | 53.58 | 3.64 | 20.01 | 77.69 |
| TemparatureInterpolation | msec. | 10.49 | 2.10 | 13.05 | 0.77 | 4.59 | 18.57 | 1.11 | 6.65 | 26.19 |
| CalcBoundaryMaskAndCollect | msec. | 20.72 | 3.97 | 24.15 | 1.35 | 7.74 | 30.30 | 1.95 | 10.79 | 42.22 |
| AreaStatics | msec. | 2.83 | 0.62 | 2.97 | 0.21 | 0.89 | 3.58 | 0.31 | 1.66 | 5.44 |
| OutputThrStatHour | msec. | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| CalcProb | msec. | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| CalcCentroid | msec. | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| OutputTempHistHour | msec. | 0.04 | 0.07 | 0.07 | 0.04 | 0.05 | 0.06 | 0.03 | 0.04 | 0.04 |
| OutputTempHistDay | msec. | 0.01 | 0.01 | 0.02 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 |
| OutputBoundCentroid | msec. | 6.65 | 1.06 | 3.67 | 0.13 | 0.30 | 1.02 | 0.21 | 0.84 | 3.75 |

### 処理時間まとめ@実行環境①
16GBメモリ搭載PCではスワップなしで実行できる限界の条件
|  |  | 日本 | 北海道 | 群馬～長野～富山 | 九州 |
| :--- | :---: | ---: | ---: | ---: | ---: |
| メッシュサイズ | m | 250 | 250 | 50 | 100 |
| 東西サイズ | px | 7,329 | 2,208 | 3,999 | 3,350 |
| 南北サイズ | px | 12,057 | 2,075 | 2,887 | 3,401 |
| 総ピクセル数 | px | 88,365,753 | 4,581,600 | 11,545,113 | 11,393,350 |
| 陸地ピクセル数 | px | 5,955,323 | 1,250,946 | 11,296,741 | 4,013,483 |
| 陸地比率 | % | 6.7 | 27.3 | 97.8 | 35.2 |
| 平均処理時間  | sec./月 | 36.0 | 6.5 | 45.6 | 17.8 |
| 気温の空間補間計算  | msec./h | 10.5 | 2.1 | 18.6 | 6.7 |

### 参考消費メモリサイズ
近接点情報floatモード
|  |  | 日本 | 北海道 | 群馬～長野～富山 | 九州 |
|:---|:---:|---:|---:|---:|---:|
| メッシュサイズ | m | 250 | 250 | 50 | 100 |
| 東西サイズ | px | 7,329 | 2,208 | 3,999 | 3,350 |
| 南北サイズ | px | 12,057 | 2,075 | 2,887 | 3,401 |
| 総ピクセル数 | px | 88,365,753 | 4,581,600 | 11,545,113 | 11,393,350 |
| 陸地ピクセル数 | px | 5,955,323 | 1,250,946 | 11,296,741 | 4,013,483 |
| 陸地比率 | % | 6.7 | 27.3 | 97.8 | 35.2 |
| DEM | MB | 695 | 39 | 99 | 98 |
| DataAnal fixed | MB | 332 | 22 | 87 | 58 |
| pNbr* (mode=3) | MB | 2108 | 425 | 3846 | 1366 |
| ShtTemp+Bound x5 | MB | 1546 | 87 | 220 | 217 |
| pAreaMask | MB | 309 | 17 | 44 | 43 |
| ProbMap(N=16,A=1) | MB | 3092 | 174 | 440 | 392 |
| Total | MB | 8082 | 764 | 4736 | 2174 |

### ファイル出力関連実行時間
ストレージのID律速になるため環境により異なる
| 出力対象 | 時間 |
:--- | ---: |
| PNG出力 | |
| GeoTIFF出力(float) | 134～137msec.
| GeoTIFF出力(int16) | 97msec.
| GPKG出力 | 20msec. |

### 実行環境②
AMD 7940HS
メモリ:64GB（DDR5-5600 x 2）
Greekom A7max
メモリ:Crucial CT2K32G56C46S5

### 処理時間@実行環境②
|  |  | 日本全国 | 北海道 | 北海道 | 群馬～長野～富山 | 群馬～長野～富山 | 群馬～長野～富山 | 九州 | 九州 | 九州 |
|:---|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| メッシュサイズ | m | 250 | 250 | 100 | 250 | 100 | 50 | 250 | 100 | 50 |
| 東西サイズ | px | 8,295 | 2,208 | 5,521 | 800 | 2,000 | 3,999 | 1,340 | 3,350 | 6,700 |
| 南北サイズ | px | 9,763 | 2,075 | 5,185 | 577 | 1,443 | 2,887 | 1,360 | 3,401 | 6,801 |
| 総ピクセル数 | px | 80,984,085 | 4,581,600 | 28,626,385 | 461,600 | 2,886,000 | 11,545,113 | 1,822,400 | 11,393,350 | 45,566,700 |
| 陸地ピクセル数 | px | 6,193,139 | 1,250,946 | 7,784,422 | 451,827 | 2,824,156 | 11,296,741 | 653,257 | 4,013,483 | 15,948,564 |
| 陸地比率 | % | 7.6 | 27.3 | 27.2 | 97.9 | 97.9 | 97.8 | 35.8 | 35.2 | 35.0 |
| 条件 |  | 計算のみ | 計算のみ | 計算のみ | 計算のみ | 計算のみ | 計算のみ | 計算のみ | 計算のみ | 計算のみ |
| メモリ消費（主要、参考） |  | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 |
| pNbr* (mode=3) | MB | 2108 | 425 | 2650 | 153 | 961 | 3846 | 222 | 1366 | 5429 |
| ShtTemp+Bound x5 | MB | 1546 | 87 | 547 | 8 | 55 | 220 | 34 | 217 | 869 |
| ProbMap(N=16,A=1) | MB | 3092 | 174 | 1094 | 17 | 110 | 440 | 62 | 392 | 1565 |
| Total | MB | 8082 | 764 | 4785 | 185 | 1182 | 4736 | 348 | 2174 | 8661 |
| 3か月実行時間 | sec. | 102.3 | 26.1 | 108.9 | 23.2 | 49.6 | 145.0 | 26.8 | 72.0 | 214.0 |
| エリア条件読込 | sec. | 19.2 | 10.2 | 16.6 | 17.2 | 18.5 | 23.0 | 18.3 | 22.8 | 36.6 |
| 平均処理時間 | sec./月 | 27.7 | 5.3 | 30.8 | 2.0 | 10.4 | 40.7 | 2.8 | 16.4 | 59.1 |
| 近傍点計算 | sec./月 | 2.2 | 0.3 | 1.6 | 0.2 | 1.0 | 3.9 | 0.3 | 2.2 | 8.4 |
| 陸地ピクセル⇒補間計算 | 1.64389057944836E-06 | 9.5 | 1.9 | 12.87 | 0.81 | 4.71 | 17.3 | 1.15 | 6.2 | 26.29 |
|  | 0.0711925627857486 | 1.0782 | 1.1891 | 1.1860 | 1.2920 | 1.2503 | 1.1018 | 1.2866 | 1.0892 | 1.1726 |
| Analyze1Step | msec. | 31.57 | 6.39 | 37.56 | 2.29 | 12.07 | 47.72 | 3.12 | 18.20 | 65.30 |
| TemparatureInterpolation | msec. | 8.81 | 1.61 | 10.85 | 0.63 | 3.77 | 15.73 | 0.89 | 5.65 | 22.42 |
| CalcBoundaryMaskAndCollect | msec. | 17.22 | 3.33 | 21.83 | 1.28 | 7.23 | 28.94 | 1.72 | 10.59 | 37.63 |
| AreaStatics | msec. | 1.75 | 0.39 | 1.96 | 0.15 | 0.61 | 2.17 | 0.22 | 1.13 | 3.06 |
| OutputThrStatHour | msec. | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| CalcProb | msec. | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| CalcCentroid | msec. | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| OutputTempHistHour | msec. | 0.06 | 0.04 | 0.05 | 0.04 | 0.06 | 0.06 | 0.04 | 0.05 | 0.04 |
| OutputTempHistDay | msec. | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 |
| OutputBoundCentroid | msec. | 3.68 | 0.98 | 2.81 | 0.14 | 0.36 | 0.78 | 0.21 | 0.74 | 2.11 |
| メモリ影響 |  |  |  |  |  |  |  |  |  |  |
| Analyze1Step改善 | % | -22.6 | -18.8 | -14.6 | -9.5 | -11.3 | -10.9 | -14.3 | -9.0 | -15.9 |
| TemparatureInterpolation改善 | % | -16.0 | -23.3 | -16.9 | -18.2 | -17.9 | -15.3 | -19.8 | -15.0 | -14.4 |
| CalcBoundaryMaskAndCollect改善 | % | -16.9 | -16.1 | -9.6 | -5.2 | -6.6 | -4.5 | -11.8 | -1.9 | -10.9 |
| AreaStatics改善 | % | -38.2 | -37.1 | -34.0 | -28.6 | -31.5 | -39.4 | -29.0 | -31.9 | -43.8 |

### 処理時間まとめ@実行環境②
|  |  | 日本 | 北海道 | 群馬～長野～富山 | 九州 |
| :--- | :---: | ---: | ---: | ---: | ---: |
| メッシュサイズ | m | 250 | 250 | 50 | 100 |
| 東西サイズ | px | 7,329 | 2,208 | 3,999 | 3,350 |
| 南北サイズ | px | 12,057 | 2,075 | 2,887 | 3,401 |
| 総ピクセル数 | px | 88,365,753 | 4,581,600 | 11,545,113 | 11,393,350 |
| 陸地ピクセル数 | px | 5,955,323 | 1,250,946 | 11,296,741 | 4,013,483 |
| 陸地比率 | % | 6.7 | 27.3 | 97.8 | 35.2 |
| 平均処理時間  | sec./月 | 27.7 | 5.3 | 40.7 | 16.4 |
| 気温の空間補間計算  | msec./h | 8.8 | 1.6 | 15.7 | 5.7 |
| 月処理時改善率 | % | -22.9 | -17.7 | -10.8 | -7.8 |
| 温度補間処理改善率 | % | -16.0 | -23.3 | -15.3 | -15.0 |
```
デュアルチャンネル化でメモリ帯域が2倍になっているが、
SIMD処理化、並列処理の最適化、L1～L3キャッシュの利用最適化などを行っているため2割程度の改善となっている
　※メモリ帯域律速になるべくならないようにしているため
　※AreaStatics処理はSINMD化、並列処理化が困難なため改善率が高い
```

---

## 9. 関数一覧

各ファイルで定義される関数の一覧。引数の `DataDEM/DataAMeDAS/DataAnal/DataFlow` は全て参照渡し。

---

### 9-1. amt.cpp

#### `InitAppConfig`
- **処理概要**: AppConfig 構造体にデフォルト値を設定
- **引数**: `AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `PrintCFG`
- **処理概要**: AppConfig の設定内容をログ出力
- **引数**: `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `LoadConfig`
- **処理概要**: CFGファイルを読み込み AppConfig に反映
- **引数**: `const char* pszPath`, `AppConfig& CFG`
- **戻り値**: bool — 読み込み成功/失敗
- **呼び出し元**: main()

#### `EvalInterp`
- **処理概要**: 複数の補間手法を比較検証し精度評価を実施
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `DataAnal& AN`, `InterpEvalResult& RES`
- **戻り値**: void
- **呼び出し元**: main()

#### `PrintEvalResult`
- **処理概要**: 補間精度評価結果をログ出力
- **引数**: `const char* label`, `const InterpEvalResult& r`, `const DataAnal* pAN`
- **戻り値**: void
- **呼び出し元**: EvalInterp()

#### `ValidateStationElevation`
- **処理概要**: AMeDAS観測点の標高をDEMと照合・修正
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`
- **戻り値**: void
- **呼び出し元**: EvalInterp()

#### `CheckPixelContext`
- **処理概要**: PixelContext（近傍観測点情報）の整合性検証
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`
- **戻り値**: void
- **呼び出し元**: パラメータ検証処理

#### `ScanFactorTPS`
- **処理概要**: TPS補間の正則化パラメータ λ を最適化スキャン
- **引数**: `const DataDEM& DEM`, `DataAnal& ANAL`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()（パラメータスキャンモード）

#### `OutputAnalyzeImage`
- **処理概要**: 気温分布・フロー解析結果を GeoTIFF/PNG で出力
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `const ConfigCore& CFG`, `DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `OutputAnalyzeBinary`
- **処理概要**: 解析結果を GRD バイナリ形式で出力
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `const AppConfig& CFG`, `DateTime& dtWr`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `BuildPixelContext`
- **処理概要**: 全陸地ピクセルの近傍16観測点情報（PixelContext）を構築
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `DataAnal& AN`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main(), RunParamScan(), RunCompare()

#### `Eval_Coast_L_2Time`
- **処理概要**: 海岸距離補正係数 coast_L を2時間データで評価
- **引数**: `const DataDEM& DEM`, `DataAnal& ANAL`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: EvalInterp()

#### `Eval_Coast_L`
- **処理概要**: 海岸距離補正係数 coast_L を1か月データで最適化スキャン
- **引数**: `const DataDEM& DEM`, `DataAnal& ANAL`, `DataAMeDAS& DT`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: EvalInterp()

#### `ComparInterpolationMethod`
- **処理概要**: 複数補間手法を同一データで比較
- **引数**: `const DataDEM& DEM`, `DataAnal& ANAL`, `DataAMeDAS& DT`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: RunCompare()

#### `RunOneBlock`
- **処理概要**: 指定期間の全解析処理（補間→境界マスク→フロー→出力）を実行するメインループ
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `DataAMeDAS& DT`, `const AppConfig& CFG`, `const ConfigRain& ITR_CFG`, `DivergRainAccum* pAccum`
- **戻り値**: void
- **呼び出し元**: RunContinuous(), RunMultiYear()

#### `RunContinuous`
- **処理概要**: 指定開始日時から TIME_SLICES ステップを連続実行（CONTINUOUS モード）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `DataAMeDAS& DT`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `RunMultiYear`
- **処理概要**: 同月・複数年を順次処理（MULTI_YEAR モード）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `DataAMeDAS& DT`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `RunParamScan`
- **処理概要**: 2段階パラメータスキャン（Stage1:スナップショット保存 / Stage2:閾値スキャン）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `DataAMeDAS& DT`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `RunCompare`
- **処理概要**: 補間手法比較モード（ComparInterpolationMethod を呼び出し）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `DataAMeDAS& DT`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `RunValueSurvey`
- **処理概要**: pNbr* 配列の値域を調査して CSV 出力（整数化設計用）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataAMeDAS& DT`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `main`
- **処理概要**: エントリポイント。CFG読み込み・DEM/AMeDAS初期化・実行モード振り分け
- **引数**: `int argc`, `char* argv[]`
- **戻り値**: int — 0:正常終了, -1:異常終了
- **呼び出し元**: OS

---

### 9-2. AnalAMeDAS.cpp

#### `CHsTextRenderer::Init`
- **処理概要**: STB TrueType フォントを読み込みグリフキャッシュを準備
- **引数**: `const std::string& fontPath`, `int h` — フォント高さ(px)
- **戻り値**: bool — 成功/失敗
- **呼び出し元**: CHsImageLegend::Init

#### `CHsImageLegend::get_text_width`
- **処理概要**: TrueType 文字列のピクセル幅を計算（カーニング対応）
- **引数**: `const stbtt_fontinfo* font`, `const char* text`, `float scale`
- **戻り値**: int — 描画幅(px)
- **呼び出し元**: CHsImageLegend::Init

#### `CHsImageLegend::Init`
- **処理概要**: 凡例画像を生成（気温/確率/等温線マップ用、mode 1-5）
- **引数**: `const ConfigCore& CFG`, `int mode`
- **戻り値**: bool — 成功/失敗
- **呼び出し元**: main() 初期化処理

#### `CHsTextRenderer::cacheGlyph`
- **処理概要**: 指定文字のグリフビットマップをキャッシュに登録
- **引数**: `char c`
- **戻り値**: void
- **呼び出し元**: CHsTextRenderer::Init

#### `CHsTextRenderer::renderTextRow`
- **処理概要**: テキスト1行を画像バッファに描画（アルファブレンディング）
- **引数**: `uint8_t* rowBuffer`, `int rowY`, `int imgWidth`, `const std::string& text`, `int posX`, `int posY`, `uint8_t textValue`
- **戻り値**: void
- **呼び出し元**: SavePNG8wSIMD

#### `NanInitFloatAVX512`
- **処理概要**: AVX-512 で float 配列を NaN で高速初期化
- **引数**: `float* p`, `size_t n`
- **戻り値**: void
- **呼び出し元**: InitAnal(), InitFlow()

#### `InitConfigCore`
- **処理概要**: ConfigCore 構造体にデフォルト値を設定
- **引数**: `ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `SanitizeFileKey`
- **処理概要**: Windows 禁止文字をアンダースコアに置換（ファイルキー用）
- **引数**: `char* sz`, `const char* szKeyName`
- **戻り値**: void
- **呼び出し元**: LoadConfigCore()

#### `LoadConfigCore`
- **処理概要**: INIファイルから ConfigCore を読み込み（未記載キーはデフォルト維持）
- **引数**: `const char* pszPath`, `ConfigCore& CFG`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: main()

#### `InitConfigRun`
- **処理概要**: ConfigRun 構造体にデフォルト値を設定
- **引数**: `ConfigRun& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `LoadConfigRun`
- **処理概要**: INIファイルから実行動作設定（実行モード・年月・ループ設定等）を読み込み
- **引数**: `const char* pszPath`, `ConfigRun& CFG`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: main()

#### `PathExistCheck`
- **処理概要**: パスの存在確認、存在しない場合はディレクトリ生成
- **引数**: `const char* szPath`
- **戻り値**: void
- **呼び出し元**: main()

#### `PrintDivergRainResult`
- **処理概要**: 発散場・降水量照合の評価結果（HSS/Hit Rate/CSI 等）をコンソール出力
- **引数**: `const DivergRainResult& r`, `int year`, `int month`, `int day`, `int hour`
- **戻り値**: void
- **呼び出し元**: main() or OpticalFlow.cpp

#### `PrintDivergRainAccum`
- **処理概要**: 発散場・降水量照合の月間累積統計を出力
- **引数**: `const DivergRainAccum& a`
- **戻り値**: void
- **呼び出し元**: main()

#### `FreeDEM`
- **処理概要**: DataDEM 構造体の全メモリを解放
- **引数**: `DataDEM& DEM`
- **戻り値**: void
- **呼び出し元**: main()

#### `InitDEM`
- **処理概要**: GeoTIFF から DEM（標高・陸マスク）を読み込み、海岸距離(CalcCoastDist)を計算
- **引数**: `DataDEM& DEM`, `const char* pPath`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: main()

#### `CalcCoastDist`
- **処理概要**: 全陸地ピクセルから最近接海岸までのユークリッド距離を Meijster EDT で計算
- **引数**: `DataDEM& DEM`
- **戻り値**: void
- **呼び出し元**: InitDEM()

#### `FreeAMeDAS`
- **処理概要**: DataAMeDAS 構造体の全メモリを解放
- **引数**: `DataAMeDAS& DT`
- **戻り値**: void
- **呼び出し元**: main()

#### `LoadAMeDAS`
- **処理概要**: 年月指定で AMeDAS バイナリを読み込み、グリッド座標へ変換
- **引数**: `const char* path`, `int year`, `int month`, `const DataDEM& DEM`, `DataAMeDAS& DT`, `const ConfigCore& CFG`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: main()

#### `SetCalcTemp`
- **処理概要**: 指定日時の気温データを fTemps バッファにセット（0.1℃整数→float変換）
- **引数**: `DataAMeDAS& DT`, `int day`, `int hour`
- **戻り値**: void
- **呼び出し元**: TemparatureInterpolation()

#### `SetCalcTempToSlot`
- **処理概要**: 指定日時の気温を複数時刻並列用スロット pfTempsArr[slot] にセット
- **引数**: `DataAMeDAS& DT`, `int day`, `int hour`, `int slot`
- **戻り値**: void
- **呼び出し元**: TemparatureInterpolationMulti()

#### `SetCalcRain`
- **処理概要**: 指定日時の降水量を fRains バッファにセット（メモリ未確保時は確保）
- **引数**: `DataAMeDAS& DT`, `int day`, `int hour`
- **戻り値**: int — 0:成功, -1:メモリ確保失敗
- **呼び出し元**: RunOneBlock()

#### `FreeAnal`
- **処理概要**: DataAnal 構造体の全メモリを解放
- **引数**: `DataAnal& AN`
- **戻り値**: void
- **呼び出し元**: InitAnal()（失敗時）, main()（終了時）

#### `InitAnal`
- **処理概要**: DataAnal を初期化（補間パラメータ設定・各種バッファ確保・近傍テーブル構築）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `ConfigCore& CFG`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: main()

#### `MatchFieldValue`
- **処理概要**: OGR フィールド値がフィルタ条件（単値/範囲/複数値）と一致するか判定
- **引数**: `const char* szFieldVal`, `const char* szSpec`
- **戻り値**: bool — 一致時 true
- **呼び出し元**: BuildAreaMask()

#### `BuildAreaMask`
- **処理概要**: SHP ポリゴンから GDAL OGR 経由でエリアマスクをグリッドに焼き込み
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `CalcAreaStat`
- **処理概要**: 気温閾値エリア別統計（発生数・消滅数等）を時別集計
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const DateTime& dtCurDT`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `GetRunMonthString`
- **処理概要**: 実行月設定から月文字列（例: "1-3"）を生成
- **引数**: `const ConfigRun& RCFG`
- **戻り値**: CHsString
- **呼び出し元**: 各種初期化・出力関数

#### `InitOutputTempStat`
- **処理概要**: 気温統計ファイルのヘッダを書き込む共通ヘルパー（内部用）
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`, `int nToutIdx`, `int nOutMode`
- **戻り値**: void
- **呼び出し元**: InitOutputTempStatHour/Day/Mth/All

#### `InitOutputTempStatHour`
- **処理概要**: 時別気温統計ファイルを初期化（ヘッダ書込）
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputTempStatHour`
- **処理概要**: 時別気温統計（平均/最低/最高/中央値）を計算して TSV 出力
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`, `int mode`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `InitOutputTempStatDay`
- **処理概要**: 日別気温統計ファイルを初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputTempStatDay`
- **処理概要**: 日別気温統計を TSV 出力してバッファリセット
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `InitOutputTempStatMth`
- **処理概要**: 月別気温統計ファイルを初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputTempStatMth`
- **処理概要**: 月別気温統計を TSV 出力してバッファリセット
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `InitOutputTempStatAll`
- **処理概要**: 全期間気温統計ファイルを初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputTempStatAll`
- **処理概要**: 全期間気温統計を TSV 出力（走行完了時に1回）
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: main()（走行完了後）

#### `InitOutputAreaLandPx`
- **処理概要**: エリア別陸地ピクセル数出力ファイルを初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`, `const DataAnal& AN`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `InitOutputThrStatHour`
- **処理概要**: 時別閾値エリア統計ファイルを初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputThrStatHour`
- **処理概要**: 時別閾値エリア統計を TSV 出力してバッファリセット
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `WriteThrStatMonthHeader`
- **処理概要**: 月次統計ファイルの共通ヘッダを書き込む（lo/hi 共用ヘルパー）
- **引数**: `int nTout`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: InitOutputThrStatMonthLo/Hi()

#### `InitOutputThrStatMonthLo`
- **処理概要**: 月次統計ファイル（T < 閾値側）を初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `InitOutputThrStatMonthHi`
- **処理概要**: 月次統計ファイル（T ≥ 閾値側）を初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputThrStatMonth`
- **処理概要**: 月次閾値エリア統計を TSV 出力してバッファリセット
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `FlushProvMonthToAccum`
- **処理概要**: 月次確率バッファを全期間バッファに加算してリセット（AVX-512）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()（月末）

#### `InitOutputProbText`
- **処理概要**: 確率マップ出力用 TSV ファイル（4種）をオープンしてヘッダを書き込み
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `CalcProbMonthStat`
- **処理概要**: 月次/全期間の確率統計（平均/最大/ヒストグラム）を計算
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `int mode` — 月次(0)/全期間(1), `int dir` — lo(0)/hi(1)
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `InitOutputBoundCentroid`
- **処理概要**: 境界面重心出力ファイルを初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputBoundCentroid`
- **処理概要**: POI 別境界面重心を BFS（8近傍）で計算して時別 TSV 出力
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `SavePNG16wSIMD`
- **処理概要**: float 配列から 16bit PNG を AVX-512 対応で高速出力
- **引数**: `const char* path`, `ConvMatrixParam& DAT`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: OutputAnalyzeImage()

#### `SavePNG8wSIMD`
- **処理概要**: float 配列から 8bit PNG を出力（パレット/凡例インポーズ対応）
- **引数**: `const char* path`, `ConvMatrixParam& DAT`, `const uint8_t* palette`, `const CHsImageLegend* pIL`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: OutputAnalyzeImage()

#### `SavePNG16woSIMD`
- **処理概要**: float 配列から 16bit PNG を出力（スカラー処理版）
- **引数**: `const char* path`, `ConvMatrixParam& DAT`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: OutputAnalyzeImage()

#### `SavePNG8woSIMD`
- **処理概要**: float 配列から 8bit PNG を出力（スカラー処理版）
- **引数**: `const char* path`, `ConvMatrixParam& DAT`, `const uint8_t* palette`, `const CHsImageLegend* pIL`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: OutputAnalyzeImage()

#### `InitProcTimeMeas`
- **処理概要**: 処理時間計測用 TMD ファイルを準備
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `FlushAllTextOut`
- **処理概要**: 全テキスト出力バッファをフラッシュ
- **引数**: なし
- **戻り値**: void
- **呼び出し元**: main()（終了時）

#### `InitAllTextOut`
- **処理概要**: 設定に応じた全テキスト出力ファイルを一括初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `OutputTempHistCore`
- **処理概要**: 気温ヒストグラムバッファを集約して TSV 出力する共通ヘルパー
- **引数**: `const DataAnal& AN`, `const ConfigCore& CFG`, `int mode`, `const uint32_t* pHistBuf`, `int nToutIdx`, `const char* szTimeKey`, `int nOutMode`, `const DateTime* pDT`
- **戻り値**: void
- **呼び出し元**: OutputTempHistHour/Day/Mth/All()

#### `InitOutputTempHist`
- **処理概要**: 気温ヒストグラムファイルのヘッダを書き込む共通ヘルパー
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`, `int nToutIdx`, `int nOutMode`
- **戻り値**: void
- **呼び出し元**: 各 InitOutputTempHist*()

#### `InitOutputTempHistHour`
- **処理概要**: 時別気温ヒストグラムファイルを初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputTempHistHour`
- **処理概要**: 時別気温ヒストグラム生データを TSV 出力
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`, `int mode`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `InitOutputTempHistDay`
- **処理概要**: 日別気温ヒストグラムファイルを初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputTempHistDay`
- **処理概要**: 日別気温ヒストグラム生データを TSV 出力
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `InitOutputTempHistMth`
- **処理概要**: 月別気温ヒストグラムファイルを初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputTempHistMth`
- **処理概要**: 月別気温ヒストグラム生データを TSV 出力
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `InitOutputTempHistAll`
- **処理概要**: 全期間気温ヒストグラムファイルを初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputTempHistAll`
- **処理概要**: 全期間気温ヒストグラム生データを TSV 出力（走行完了時に1回）
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: main()（走行完了後）

---

### 9-3. Interpolation.cpp

#### `CalcInterpParams`
- **処理概要**: 観測点配置から補間パラメータ（Barnes kappa, RBF sigma 等）を自動計算
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `InterpConfig& ITR_CFG`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: main()（初期化時）, RunParamScan(), RunCompare()

#### `BuildNbrTable`
- **処理概要**: 全陸地ピクセルの近傍16観測点を探索し SoA 配列（pNbr*）を構築（AVX-512最適化）
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `DataAnal& AN`
- **戻り値**: void
- **呼び出し元**: InitAnal(), BuildPixelContext(), RunValueSurvey()

#### `BuildWeightTable`
- **処理概要**: 補間手法別重み・高度補正係数を事前計算（pNbrW, pNbrAltW）
- **引数**: `const DataDEM& DEM`, `const DataAMeDAS& DT`, `DataAnal& AN`
- **戻り値**: void
- **呼び出し元**: TemparatureInterpolation(), ChgInterpMethod(), TemparatureInterpolationMulti()

#### `ChgInterpMethod`
- **処理概要**: 補間手法を切り替え（Barnes2Pass 用 pNbrLapse 確保を含む）
- **引数**: `const DataDEM& DEM`, `const DataAMeDAS& DT`, `DataAnal& AN`, `int nNewMode`
- **戻り値**: void
- **呼び出し元**: main()（パラメータ検証時）

#### `InterpolationMain`
- **処理概要**: IDW/Shepard/RBF/Barnes 共通補間ループ（AVX-512 SIMD 化）
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `DataAnal& AN`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: TemparatureInterpolation()

#### `InterpolationMainMulti`
- **処理概要**: N 時刻分の補間を共有ループで一括処理（キャッシュ効率向上）
- **引数**: `const DataDEM& DEM`, `float* const* ppFTemps`, `int16_t** ppDst`, `int nBatch`, `DataAnal& AN`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: TemparatureInterpolationMulti()

#### `CollectAreaCntFromShtTemp`
- **処理概要**: mShtTemp から閾値別エリアカウント（areaCnt）を収集して pAreaCntBuf を更新
- **引数**: `const DataAnal& AN`
- **戻り値**: void
- **呼び出し元**: TemparatureInterpolationMulti()

#### `TemparatureInterpolationMulti`
- **処理概要**: N 時刻分の気温補間を並列実行（InterpolationMainMulti を呼び出し）
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `DataAnal& AN`, `const ConfigCore& CFG`, `const int* pSlots`, `int nBatch`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `InitTPS`
- **処理概要**: TPS 補間用作業バッファ（tps_w, tps_stLon, tps_stLat 等）を確保・初期化
- **引数**: `const DataAMeDAS& DT`, `DataAnal& AN`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: main()（初期化時）

#### `CalcBoundaryMaskCoreInt16`
- **処理概要**: int16 気温場から各閾値の境界マスクを生成する内部ヘルパー（AVX-512）
- **引数**: `const DataDEM& DEM`, `const DataAnal& AN`, `const int16_t* pSht`, `uint16_t* pPacked`
- **戻り値**: void
- **呼び出し元**: CalcBoundaryMaskAndCollect()

#### `CalcBoundaryMaskAndCollect`
- **処理概要**: 境界マスク生成とエリア統計集計を同時実行
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`
- **戻り値**: void
- **呼び出し元**: TemparatureInterpolation(), RunOneBlock()

#### `CalcBoundaryMaskAndCollectMulti`
- **処理概要**: N 時刻の境界マスク生成・エリア統計集計を一括処理
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const int* pBatchSlots`, `int nBatch`
- **戻り値**: void
- **呼び出し元**: RunOneBlock()

#### `TemparatureSetInt16`
- **処理概要**: float 気温補間結果を int16 にスケーリングしてパック
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `DataAnal& AN`, `const ConfigCore& CFG`, `int idx`
- **戻り値**: void
- **呼び出し元**: main()（補間処理内）

#### `CalcTPSCoeff`
- **処理概要**: TPS 補間の係数をガウス消去法で計算
- **引数**: `DataAMeDAS& DT`, `DataAnal& AN`, `const float* pTempOverride` — nullptr 時は DT.fTemps 使用
- **戻り値**: void
- **呼び出し元**: InterpolationTPS()

#### `InterpolationTPS`
- **処理概要**: Thin Plate Spline 法による気温補間
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `DataAnal& AN`, `int day`, `int hour`
- **戻り値**: void
- **呼び出し元**: TemparatureInterpolation()

#### `InterpolationBarnes2Pass`
- **処理概要**: Barnes 法 2 パス補間（第1パス: 基本補間, 第2パス: 残差補間）
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `DataAnal& AN`, `int day`, `int hour`
- **戻り値**: void
- **呼び出し元**: TemparatureInterpolation()

#### `TemparatureInterpolation`
- **処理概要**: 補間計算のエントリ関数。`AN.interp.method` で手法を切り替え（IDW/RBF/Barnes/TPS/Barnes2Pass）
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `DataAnal& AN`, `const ConfigCore& CFG`, `int day`, `int hour`
- **戻り値**: void
- **呼び出し元**: RunOneBlock(), main()（パラメータスキャン内）

---

### 9-4. OpticalFlow.cpp

#### `FreeFlow`
- **処理概要**: DataFlow 構造体の全メモリを解放
- **引数**: `DataFlow& FL`
- **戻り値**: void
- **呼び出し元**: main()（終了時・初期化失敗時）

#### `InitFlow`
- **処理概要**: DataFlow を確保・初期化（pFlowU/V/Magnitude/Diverg/Vorticity を NaN 初期化）
- **引数**: `const DataDEM& DEM`, `DataFlow& FL`, `const ConfigCore& CFG`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: main()

#### `CalcGradientFlow`
- **処理概要**: 勾配法オプティカルフロー（2時刻の気温場から空間勾配+時間差分でU/Vを算出）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `float alpha_sq` — 正則化パラメータ
- **戻り値**: void
- **呼び出し元**: RunOneBlock(), RunParamScan()

#### `Calc3PointGradientFlow`
- **処理概要**: 3時刻の気温場を使った勾配法（2次精度中心差分による時間微分）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `float alpha_sq`
- **戻り値**: void
- **呼び出し元**: RunOneBlock(), RunParamScan()（TIME_SLICES ≥ 3 時）

#### `CalcHornSchunckFlow`
- **処理概要**: Horn-Schunck 法オプティカルフロー（反復法で滑らかなフロー場を生成）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `float alpha_sq`, `int nIter` — 反復回数
- **戻り値**: void
- **呼び出し元**: RunOneBlock(), RunParamScan()

---

*以上、VER.0.4.6.000時点のシステム全体解説。*

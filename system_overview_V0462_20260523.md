# AMeDAS気温空間分布解析システム — 全体詳細解説
**VER.0.4.6.140　2026/05/11**

本書は `system_overview_V046_20260425.md` を元に、2026/05/11時点の現行コード（`amt.cpp` 先頭の `VER.0.4.6.140`）へ合わせて修正・追記したもの。

---

## 0. 設計思想・コーディング方針

本プロジェクトのコードを読み書きする際（AIエージェント・人間とも）に前提とすべき設計思想。

### 0-1. C言語スタイルの計算コアを核とする
- **解析・補間・統計の計算部分は C 言語スタイル**で記述する（`typedef struct` によるデータ構造、ポインタ直接アクセス、関数ベースの処理）。
- これは速度のための意図的な設計。AVX-512 SIMD・OpenMP 並列・メモリ階層（L3キャッシュ）を意識した最適化は、ポインタを直接制御できる C/C++ だからこそ達成できる。**おそらく他の言語ではこの速度は出せない**。本プロジェクトの計算部分はこの思想で結果を出している。
- C でも構造化は十分可能（zlib/libpng/libjpeg が C で構造化されている前例）。**ポインタを過度に隠蔽するのは本末転倒**で、C/C++ の強み（直接メモリ制御）を捨てることになる。それなら別言語を使うべき。

### 0-2. C++ はユーティリティクラスに限定する
- C++ クラスでのカプセル化は **「中の構造を細かく見る必要がないユーティリティ」に限定**する（`CHsString`/`CHsIniFile`/`CHsLog`/`CHsTextRenderer`/NetCDF出力クラス等）。
- **`std` namespace の抽象（`std::string_view` 等の C++17 機能）を計算部分に持ち込まない**。計算側を C 言語シンプル記述にしている意味がなくなり、C++ 言語のハードルを不必要に上げる。

### 0-3. `CHsString`（CString互換クラス）の使い方
- `CHsString` は VC6 の `CString` 互換＋発展。`operator const char*()` を持ち、**通常の関数引数（`const char*` を取る）には暗黙変換が効く**ため明示変換は不要。
- **例外: 可変長引数（`printf`/`fprintf` 等の `...`）には C++ 仕様上、暗黙変換が効かない**（`va_arg` は型を見ない）。可変長引数に渡すときのみ `(const char*)str` または `str.GetBuffer()` で明示する。CHsString を裸で `%s` に渡すとクラッシュする（`-Wformat` で検出可能）。
- `GetBuffer()` への全面統一は不要。「関数引数は暗黙変換に任せ、可変長引数のみ明示」が CString 互換として自然。

### 0-4. amt.cpp の位置づけ（責務境界）
- **当初の「解析コアを完全分離しライブラリ化（別アプリ解析パス）」構想は撤回（2026-05-21、review 5-B-1 撤回）**。理由は処理速度が DISK-IO 律速まで到達し「GRD生成→別アプリ解析」より複数回走行が有利になったこと、RUN_MODE 整理で異機能組込が容易になり別アプリ化の必要が薄れたこと。
- 現行方針は **「amt.cpp＝全体コントロール（main/設定/実行制御＋モード別 Run* と評価コア）／AnalAMeDAS.cpp 等＝ライブラリ（複数モード共通の汎用機能）」の責務境界明確化**。汎用出力（GeoTIFF/GPKG/AnalyzeBinary）・補間準備（BuildPixelContext）はライブラリ側へ移設済み、`Run*` は月単位純化（RunOneMonth）＋共通ヘルパー（SetupForMonth/FlushMonth/ResetAnalBuffers/ZoneClassifier）化済み（Phase 7 完了、VER.0.4.6.179）。`EvalInterp`/`ComparInterpolationMethod` 等の評価コアは特定モード専用のため amt.cpp に同居（ユーザー方針）。

---

## 1. ファイル構成

| ファイル | 行数 | 役割 |
|:---|---:|:---|
|`amt.cpp`|4233|解析ロジック・実行制御の全体（設定読込・初期化・実行モード）|
|`AnalAMeDAS.h`|1311|全構造体・定数・インライン関数の定義|
|`AnalAMeDAS.cpp`|8504|初期化・解放・データ読込・出力関数群|
|`Interpolation.cpp`|2769|補間計算・近傍点テーブル構築・エリア統計集計|
|`OpticalFlow.cpp`|616|勾配フロー・Horn-Schunck実装（AVX-512）|
|`LogStDatNetCDF.h`|479|NetCDF出力クラス定義、観測点計算結果時別データ出力|
|`HourlyStatNetCDF.h`|1061|NetCDF出力クラス定義、時別データ出力|
|`HsCommon.h`|1394|汎用クラス（CHsString/CHsIniFile/CHsLog）|
|`WriteGeoTIFF.h`|600|GeoTIFF出力（GDAL）|
|`WriteGRD.h`|1493|独自バイナリ(GRD)出力・読込、BLOSC2圧縮|
|`ConvMatrix.h`|404|float/int16変換・カラーマップ適用・PNG行生成|
|`TMD.h`|705|汎用時系列データロガー形式の読み書き|
|`Nifsq8.cpp/.h`|1059|補助モジュール（データ圧縮補助）|
|`palette1.h`|259|8bitPNG用パレット定義|
|`palette2.h`|259|8bitPNG用パレット定義|
|`stb_truetype.h`|5079|TrueTypeフォントレンダリング（ヘッダオンリー、外部ライブラリ）|
|`*.cfg`|—|実行設定ファイル（INI形式）|


---

## 1-1. VER.0.4.6.140時点の主な更新点

- `HSS_GRD_FAC_TEMP=500.0f` に変更済み。気温GRD/int16温度面は0.002℃精度。
- `pAnalTemp` は通常Barnes経路では確保しない。Barnes2Pass/TPSの一時バッファとしてのみ確保する。
- 補間手法定数は `0=IDW, 1=RBF_GAUSS, 2=SHEPARD, 3=KRIGING_S未完成, 4=BARNES, 5=BARNES2PASS, 6=TPS`。旧資料の `INTERP_BARNES_H` は現行コードに存在しない。
- `PARALLEL_SLICES` により複数時刻同時補間を行う。日境界・lapse時刻帯・地形補正時刻帯をまたぐ場合はバッチを分割する。
- `USE_NBR_FLOAT` により近傍点テーブルのSRC/Wをfloat/intへ切り替える。`ALTW_CACHE` は時刻帯別の高度/地形補正項をキャッシュする。
- `AREA_CONFIG.TERRAIN_CORR` により地形特徴量を使った補正を追加。`TERRAIN_DEM_PATH` が必要。
- `MODE=ST_NBR_EVAL` を追加。観測点近傍、地形特徴量、補正項の検証用に通常出力を抑止して限定日時を処理する。
- 気温統計・気温ヒストグラム・閾値統計の一部はNetCDF出力にも対応している。

## 1-2. VER.0.4.6.142→200 リファクタリングの主な更新点（2026-05-23 クローズ）

`system_review_V0461_20260516.md` に基づくリファクタリングを VER.0.4.6.200 で完了。本概要書も以下の構造変更を反映している:

- **`DataAnal` を5サブ構造体に分割**（Phase6）: `nbr`(NbrTable) / `stats`(TempStats) / `prob`(ProbMap) / `area`(AreaStats) / `multi`(MultiPool)。参照は `AN.area.vAreaStat` のようにサブ構造体経由。
- **設定構造の整理**（Phase5・3-3）: `[INTERP]` 5float と `[OUTPUT_GRD]` チャンネル配列は `ConfigCore` に統合済み。`AppConfig` は `{ConfigCore C; ConfigRun R}` の薄いラッパで、**ConfigCore の既定値・読込は `InitConfigCore`/`LoadConfigCore` に一元化**（`InitAppConfig`/`LoadConfig` は C と R の Init/Load を束ねるだけ）。引数規約: **Run* エントリ＝AppConfig／リーフ計算・出力＝ConfigCore**。
- **時間ループ関数の再編**（柱C、VER.0.4.6.178）: 旧 `RunOneBlock` を **`RunOneMonth`（月単位純化、月境界判定を内部に持たない）** へ。月初準備＝`SetupForMonth`、月末出力＝`FlushMonth`、バッファ初期化＝`ResetAnalBuffers` に抽出。**CONTINUOUS=期間動作**（月チャンク分割ループ＋TOTAL_OUT、MONTH_OUT廃止）／**MULTI_YEAR=月動作**（月別独立リセット＋月別ProbImage）。
- **責務境界の明確化**（柱A）: 汎用出力（`OutputTempGeoTIFF`/`OutputFlowGeoTIFF`/`OutputFlowGPKG`/`OutputAnalyzeBinary`）を AnalAMeDAS.cpp へ、補間準備 `BuildPixelContext` を Interpolation.cpp へ移設。**amt.cpp 完全解体・ライブラリ化は撤回**（§0-4）。
- **テーブル駆動化**: 補間手法名（`InterpMethodInfo`/`GetInterpMethodName`）、`ComparInterpolationMethod` の手法ループ（8-8）、出力定義（`g_TempTextDefs`）。**LOOCV ゾーン分類は `ZoneClassifier` 構造体に集約**（5-B-8）。
- **出力ファイル名の統一**（C3d）: 統計TXT/NC・TMD・LOG が共通ステム `GetRunStemName()` を共有。MULTI_YEAR=`KEY_YYYY-YYYY_MM`／その他（CONTINUOUS等）=`KEY_YYYYMMDD_HH`。
- **RunMode別固定値の単一窓口化**: `ApplyRunModeFixedConfig()` を `LoadConfig` 末尾で適用（CONTINUOUS の MONTH_OUT 無効化、ST_NBR_EVAL の出力抑止等）。
- **検証規約**: GRD(.tva) は Blosc2 並列圧縮で**バイト非決定的**のためリグレッション比較に使えない。TXT/NC の値・SHA256 で判定する。

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
| `pTerrainElevation` / `pTerrainMask` | 地形補正用DEMと有効マスク（`AREA_CONFIG.TERRAIN_CORR=1`時に使用） |

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
float*   pTerrainFeature;// 観測点側の地形補正特徴量
// 複数時刻並列補間用バッファ
float*   pfTempsPool;    // 全スロット連続確保バッファ
float*   pfTempsArr[];   // スロット別ポインタ（[0]=fTemps）
char     (*szKNames)[16];// 観測点漢字名（Shift-JIS）
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

> **【Phase6 分割】** `DataAnal` は機能別の5サブ構造体に分割済み。以下の各バッファは対応するサブ構造体経由でアクセスする（コード上は `AN.area.vAreaStat` 等）:
> | サブ構造体 | 主な内容 |
> |---|---|
> | `nbr`（NbrTable） | 近傍点テーブル（pNbrN/pNbrIdx/pNbrD2/pNbrW/pNbrAltW… ＋ nPreMethod 等） |
> | `stats`（TempStats） | 気温統計（時別/日別/月別/全期間）＋気温ヒストグラム各バッファ |
> | `prob`（ProbMap） | 確率マップ（pCountsA/pCountsM/pCountsA_hi/pCountsM_hi/pValid*） |
> | `area`（AreaStats） | エリア別統計（vAreaStat/pAreaMask/nAreaLandPx/pHourlyBuf/nHourlyCapacity 等） |
> | `multi`（MultiPool） | 複数時刻並列バッチ状態（nBatchNext/batchSlots/nBatchN ＋状態遷移メソッド） |
>
> 以下の表・コードは論理的内容を示すもので、実アクセスは上記サブ構造体経由（例: `pNbrN`→`AN.nbr.pNbrN`、`vAreaStat`→`AN.area.vAreaStat`、`pCountsM`→`AN.prob.pCountsM`、`pTHistHourBuf`→`AN.stats.pTHistHourBuf`）。

**近傍点テーブル（SoA構造、月1回BuildNbrTable()で構築、`AN.nbr`）**

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
| `pNbrAltW` / `pNbrAltWF` | lapse_rate×altDiff + 地形補正項 (int16 or float) | — |
| `pNbrAltWCache` / `pNbrAltWFCache` | 時刻帯別の高度/地形補正キャッシュ（`ALTW_CACHE=1`時） | — |

`nNbrFloatMode` でint/float版を切り替え（0=全整数, 1=SRCのみfloat, 2=Wのみfloat, 3=全float）。

**気温補間結果（リングバッファ）**

```c
float*    pAnalTemp;              // Barnes2Pass/TPS用のfloat一時バッファ（通常Barnesでは未確保）
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
AreaMonthStat vAreaStat[OBTAIN_AREA];  // 月次統計バッファ（[AREA_DEF_MAX]=全域）
uint32_t      nAreaLandPx[OBTAIN_AREA];// 陸地ピクセル数（[AREA_DEF_MAX]=全域）
```

**確率マップ**

```c
uint32_t* pCountsA;    // 全期間カウント(lo) [pixel * nThresholds + threshold]、PROB_ANNUAL=1時
uint16_t* pCountsM;    // 月次カウント(lo)。PROB_MODE=2時はhi側カウントを兼用格納
uint32_t* pCountsA_hi; // 全期間カウント(hi)。PROB_MODE=3かつPROB_ANNUAL=1時のみ
uint16_t* pCountsM_hi; // 月次カウント(hi)。PROB_MODE=3時のみ
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

#### `ConfigCore`（INI読込可能な解析設定）
解析共通設定。`InitConfigCore()` でデフォルト値セット後、`LoadConfigCore()` でINIから上書き。主要セクション：[MAIN] [PATH] [RUN] [EXT_ANAL] [INTERP] [PARAM] [OUTPUT] [OUTPUT_VECTOR] [OUTPUT_GEOTIFF] [OUTPUT_PNG] [PNG_OPT] [LEGEND_*] [OUTPUT_GRD] [TEMP_STAT] [TEMP_HIST] [THRESHOLD_ANAL] [THRESHOLD] [AREA] [POI_CENTROID] [LAPSE] [AREA_CONFIG]。**[INTERP] 5float（COAST_L/IDW_POWER/RBF_SIGMA/BARNES_KAPPA_FAC/BARNES_GAMMA）と [OUTPUT_GRD] チャンネルタグ（CH_NUM/CH_TAG）も ConfigCore に統合済み**で、既定値は `InitConfigCore`、読込は `LoadConfigCore` に一元化（3-3）。

#### `ConfigRun`（実行動作定義）
旧 `RunConfig` に相当。開始日時・ステップ数・年範囲・対象月リスト等を保持する。現行の `nRunMode` は 0=CONTINUOUS, 1=MULTI_YEAR, 2=PARAM_SCAN, 3=COMPARE, 4=VALUE_SURVEY, 5=EVAL_INTERP, 6=OPTIM_PARAMS, 7=AREA_CHAR, 8=ST_NBR_EVAL。

#### `AppConfig`（`ConfigCore C` + `ConfigRun R` の薄いラッパ）
`InitAppConfig()`＝`InitConfigCore`+`InitConfigRun`、`LoadConfig()`＝`LoadConfigCore`+`LoadConfigRun`+type検証+`ApplyRunModeFixedConfig()`（RunMode別固定値の単一窓口）。**引数規約**: Run* エントリ関数と main 直下の設定処理は `AppConfig`、補間/出力などのリーフ計算・出力関数は `ConfigCore`（`CFG.C`）のみを受け取る（3-3 で明文化）。

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
├─ CHsLog::Open()                 # ログ開始
├─ GDAL / Blosc2 初期化
├─ InitAppConfig()                # デフォルト設定
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
            ├─ COMPARE(3):        RunCompare()       # 補間法比較
            ├─ VALUE_SURVEY(4):   RunValueSurvey()   # 値調査
            ├─ EVAL_INTERP(5):    RunEvalInterp()    # 補間法精度評価（LOOCV）
            ├─ OPTIM_PARAMS(6):   RunOptimParams()   # 補間パラメータ最適化スキャン
            ├─ AREA_CHAR(7):      RunAreaChar()      # 解析エリア特性評価
            └─ ST_NBR_EVAL(8):    RunStNbrEval()     # 観測点近傍・補正項評価用出力
```

**複数cfgファイルの連続実行**: DEMが共通であれば再初期化なしで続けて実行。DEMパスが変わった場合のみ再ロード。

---

## 5. 主要関数群の詳細

### 5-1. 初期化・解放系

#### `InitDEM()` / `FreeDEM()`
GDALでGeoTIFFを読み込み、`pElevation`・`pSeaMask`・`pCoastDist` を確保。`nAlignedWidth` はAVX-512の64バイト境界（16float）に合わせて切り上げ。

#### `InitTerrainDEM()`
`AREA_CONFIG.TERRAIN_CORR=1` の場合に、`PATH.TERRAIN_DEM_PATH` で指定した地形補正用DEMを読み込む。読み込み失敗、メモリ確保失敗、RasterIO失敗時は警告を出して `TERRAIN_CORR=0` 扱いへ落とす。

#### `InitAnal()` / `FreeAnal()`
`DataAnal` の主要ポインタをAVX-512アライメント確保（`_mm_malloc`）。`mShtTemp[nShtTempSlots]` はリングバッファとして動作。`pAnalTemp` は通常Barnesでは確保せず、Barnes2Pass/TPS時のみ確保する。CFG.nAreaStatMode・nAreaHistMode・nProveMode等の設定に応じてオプションバッファ（pCountsM系・pTempHistBuf系・NetCDF系等）の確保を制御。

#### `InitFlow()` / `FreeFlow()`
`DataFlow` の5配列（pFlowU/V/Magnitude/Diverg/Vorticity）を確保。

#### 補間手法別バッファの初期化（`ChgInterpMethod()` 経由、VER.0.4.6.172〜173）
補間手法に依存するバッファ（Barnes2Pass: `pNbrLapse`、TPS: `tps_*`、両者共通: `pAnalTemp`）は `ChgInterpMethod()` 内の static ヘルパー（`InitInterpLapseCache`/`InitInterpTPSWork`/`EnsureAnalTempBuffer`）で**選択手法に必要な分のみ**確保する（TPSを選ばない限り `tps_*` を確保しない＝メモリ節約）。旧 `InitTPS()`（public）は `InitInterpTPSWork`（Interpolation.cpp 内 static）に統合され、評価モード・通常運用とも `ChgInterpMethod()` 経由に一本化された（呼び出し漏れ防止）。通常運用では `SetupForMonth()` 内（月初ロード直後）で `ChgInterpMethod()` を呼ぶことで初期化を保証する。

#### `FreeAMeDAS()`
月替わり時に `_mm_free` で観測データを全解放。

---

### 5-2. データ読み込み系

#### `LoadAMeDAS()`
```
1. STNファイル（観測所マスタ）読み込み → 緯度・経度・標高をfloatに変換
2. IDXファイル（観測所番号リスト）と照合
3. 各観測所のAMDファイル（日別データ）を読み込み
4. pixX/pixY にDEMグリッド上の対応ピクセル座標を設定（旧 `ValidateStationElevation` はデッドコードのため削除済み、柱A A1）
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
| `INTERP_KRIGING_S` | 3 | 簡易クリギング（未完成、通常使用しない） |
| `INTERP_BARNES` | 4 | Barnes解析（推奨） |
| `INTERP_BARNES2PASS` | 5 | Barnes 2Passスキャン（観測点密度不足で逆効果の評価） |
| `INTERP_TPS` | 6 | Thin Plate Spline（計算コスト大） |

**海岸距離補正**（`GetCoastFactor()`）: ピクセルと観測点の海岸距離差に基づいてIDW重みを減衰。`coast_L=0` で補正なし。

**気温減率補正**（`GetLapseRate()`）: `lapse_rate_table[緯度帯5][季節4]` から減率を取得し、標高差による気温差を補正。

旧資料にあった `INTERP_BARNES_H` / Barnes高精度版は現行コードでは定数として存在しない。Barnes系の通常運用は `INTERP_BARNES`、評価・特殊用途として `INTERP_BARNES2PASS` を使用する。

**補間手法と複数時刻並列（PARALLEL_SLICES）の関係**（VER.0.4.6.173）:
- `INTERP_IDW`/`INTERP_RBF_GAUSS`/`INTERP_SHEPARD`/`INTERP_BARNES` は逐次パス・並列パス（`TemparatureInterpolationMulti`）の両方に対応する。
- `INTERP_BARNES2PASS`/`INTERP_TPS` は**逐次パス専用**。並列パスは Barnes系1Passのカーネル（`InterpolationMainMulti`）に特化しているため、これらの手法は並列パスでは正しく補間できない。
- そのため `LoadConfigCore` で **Barnes2Pass/TPS かつ PARALLEL_SLICES>1 の場合は警告を出力して `nParallelSlices=1`（逐次）へ自動矯正し、動作を継続**する。Barnes2Pass/TPS は元々低速（特にTPSは数百ms/ステップ）で並列の費用対効果も低く、これらを使うのは評価・研究目的で速度より正しさが優先されるための設計判断（方針: 逐次フォールバック）。
- 各手法に必要なバッファ（Barnes2Pass: `pNbrLapse`、TPS: `tps_*`、両者共通: `pAnalTemp`）は、通常運用でも `SetupForMonth()`（月初ロード直後）が `ChgInterpMethod()` を経由して確保する（評価モードと同経路に統一、VER.0.4.6.172）。

---

### 5-5. 補間実行系

#### `TemparatureInterpolation()`
時間ループの1ステップとして以下を実行。
```
1. SetCalcTemp(DT, day, hour)
2. 手法/時刻帯変更時に BuildWeightTable() を更新
3. InterpolationMain() / InterpolationBarnes2Pass() / InterpolationTPS()
   - 全陸地ピクセルを並列処理（#pragma omp parallel for）
   - 近傍点テーブルから重みを取得し各手法で補間
   - 結果を mShtTemp[idx_order[0]] に格納（int16_t）
   - 通常Barnes系は `pAnalTemp` を使わず `mShtTemp` に直接格納
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

|指標|正式名称|概要|
|---|---|---|
|RMSE|Root Mean Squared Error（二乗平均平方根誤差）|誤差を二乗して平均し、その平方根をとったもの。大きな誤差に敏感（ペナルティが大きい）。|
|MAE|Mean Absolute Error（平均絶対誤差）誤差の絶対値の平均。|外れ値の影響を受けにくく、直感的に理解しやすい。|
|BIAS|Bias（バイアス）|予測値と実測値の差の平均。モデルが全体として過大評価（＋）か過小評価（－）かを示す。|
|NSE|Nash-Sutcliffe Efficiency（ナッシュ・サトクリフ効率係数）|主に水文学で使用される。1に近いほど精度が高く、0以下は「平均値を使うより精度が低い」ことを意味する。|
|R|Correlation Coefficient（相関係数）|予測値と実測値の連動性（相関）の強さを示す。通常、1に近いほど正の相関が強い。|

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

#### `RunOneMonth()`（月単位処理の中核、VER.0.4.6.178 で RunOneBlock から純化）
月境界判定を内部に持たず、呼出側が「月内に収まる nSteps」を保証する純粋関数。リセット/flush は呼出側責務。
```
RunOneMonth(dtStart, nSteps, nGlobalStepBase, dtPrev[in/out]):
    0. SetupForMonth（LoadAMeDAS / CalcInterpParams / BuildPixelContext / ChgInterpMethod / 月初初期化）※ループ前に1回
       AN.multi.Reset()（並列バッチ状態は月単位リセット）
    for s in range(nSteps):
        step = nGlobalStepBase + s   # グローバル: フロー有無/間引き位相/GRD判定用。s は月ローカル: 並列バッチ判定用
        1. RotateTimeSlice
        2. （0時かつ day>1 のとき）前日分の日別統計出力
        3. TemparatureInterpolation / TemparatureInterpolationMulti（補間 + int16化）
        4. CalcBoundaryMaskOnly または CalcBoundaryMaskAndCollectMulti（境界マスク・統計）
        5. 観測点別ログ、時別気温統計、時別ヒストグラム出力
        6. PNG出力（気温、等温線、マージ。実行判断は各関数先頭）
        7. CalcAreaStat / OutputBoundCentroid / OutputBoundMaskGeoTIFF
        8. CalcGradientFlow / CalcHornSchunckFlow（OPT_FLOW有効かつ step>=1）
        9. GeoTIFF / GPKG / GRD / DivergRain評価
        10. dtPrev=t0 / DateTimeAdvance / 処理時間ログ更新
```
共通ヘルパー: `SetupForMonth`（月初ロード+補間準備）、`FlushMonth`（月末/期間末の統計・確率出力）、`ResetAnalBuffers`（idx_order/エリアバッファ初期化）。

#### `RunContinuous()`（期間動作）
CFGの開始日時から `nRunCount` ステップを**月チャンクに分割**し `RunOneMonth` を反復。`ResetAnalBuffers` は期間先頭で1回（月跨ぎで idx_order/フロー連続性維持）、`nGlobalStepBase`/`dtPrev` を月跨ぎで引継ぎ。各チャンク末で `FlushMonth(bProbImage=false)`、期間末で `FlushAllTempStat`（TOTAL）。MONTH_OUT は `ApplyRunModeFixedConfig` で無効化（期間統計は TOTAL_OUT）。

#### `RunMultiYear()`（月動作）
`nYearStart`〜`nYearEnd` の同月（`nTargetMonthList[]`）をループ。各月で `ResetAnalBuffers`（月別独立）→ `RunOneMonth(nGlobalStepBase=0)` → `FlushMonth(bProbImage=true)`。年別統計を1行ずつ出力し、終了後に全期間確率マップ＋全期間統計を出力。

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

#### `RunStNbrEval()`（MODE=ST_NBR_EVAL）
観測点近傍、地形特徴量、補正項確認用の補助モード。`SetStNbrEvalFixedConfig()` によりPNG/GeoTIFF/GRD/確率/フロー等の通常出力を抑止し、指定年・月・日・時刻だけを処理する。結果は観測点ログNetCDFとmanifestに記録され、地形補正や近傍点構成の検証に使う。

---

## 6. 出力ファイル系

### 6-1. GeoTIFF（WriteGeoTIFF.h）
- `nOutputTempGeoTIFF=1`: 気温解析結果
- `nOutputFlowGeoTIFF=1`: フロー結果（FLOWU, FLOWV, MAGNITUDE, DIVERG, VORTICITY）
- `nOutputGeoTIFF3=1`: 気温閾値確率分析結果
- `nOutputGeoTIFF4=1`: 気温境界面マスク
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
| TEMP | 気温 | int16 | HSS_GRD_FAC_TEMP=500（0.002℃精度） |
| FLOWU | 経度方向速度 | int16 | 1.0e-4、offset=3.2767 |
| FLOWV | 緯度方向速度 | int16 | 1.0e-4、offset=3.2767 |
| MAGNITUDE | 速度スカラー | int16 | 1.0e-4、offset=3.2767 |
| DIVERG | 発散 | int16 | 1.0e-3、offset=32.767 |
| VORTICITY | 渦度 | int16 | 1.0e-3、offset=32.767 |
| BOUND | 等温線境界マスク | B16（ビットパック） | bit c = 閾値c |

`nThinningGRD` で出力頻度を制御。

> **【検証上の注意】GRD(.tva) はバイト非決定的**: BLOSC2 の8スレッド圧縮はブロックのスレッド分配が実行ごとに変わるため、**同一入力・同一バイナリでもバイト列が一致しない**（解凍後の値は常に同一）。リグレッション検証では GRD のバイト/SHA256 比較は使えず、**TXT/NC の値・SHA256 で判定**すること（決定的にしたい場合は圧縮スレッド数を1にする選択肢がある）。

> **【出力ファイル名の共通ステム】** 統計TXT/NC・TMD・LOG は `GetRunStemName()` で生成する共通ステムを共有する（VER.0.4.6.179、C3d）。`{RESULT_PATH}{FILE_KEY_ST}` に続けて、MULTI_YEAR=`_YYYY-YYYY_MM`／その他（CONTINUOUS等）=`_YYYYMMDD_HH` を付す。GeoTIFF/GRD/PNG 等の地理データは `FILE_KEY`＋日時ベース。

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

```ini
[MAIN]
TYPE         = main                          # 固定(設定ファイルの適合性確認用)

[PATH]
DEM_PATH     = U:\DEM\japan_250m.tif         # DEM GeoTIFFファイルパス(末尾に「\」が必要)
TERRAIN_DEM_PATH = U:\DEM\terrain_50m.tif    # 地形補正用DEM（AREA_CONFIG.TERRAIN_CORR=1時）
RESULT_PATH  = U:\PRG3\AMT_RES\20260425A\    # 結果出力先ディレクトリ
FILE_KEY     = japan                         # 画像出力ファイルの先頭文字列
FILE_KEY_ST  = RES_JPN                       # 統計結果テキストファイルの先頭文字列
AMEDAS_PATH  = U:\AMeDAS\                    # AMeDASデータディレクトリ（STN/IDX/AMDファイルが存在するパス）

[RUN]
# MODE : CONTINUOUS / MULTI_YEAR / PARAM_SCAN / COMPARE / VALUE_SURVEY / EVAL_INTERP / OPTIM_PARAMS / AREA_CHAR / ST_NBR_EVAL
MODE         = MULTI_YEAR
CONT_START   = 2000/2/15/0                   # CONTINUOUSモードの実行開始日時（年/月/日/時）
CONT_DUR     = 2/1                           # CONTINUOUSモードの実行期間（時/日）、日設定を優先
MULTI_YEAR   = 1999/2001                     # MULTI_YEARモードの対象年範囲（開始年/終了年）
MULTI_MONTH  = 1/2/3/12                      # MULTI_YEARモードの対象月リスト（1-12のスラッシュ区切り）
PARALLEL_SLICES = 4                          #複数時刻同時計算
USE_NBR_FLOAT   = 3                          #近傍点関連をfloatで処理する（0:すべて整数、1:SRCのみfloat、2:Wのみfloat、3：全てfloat）
ALTW_CACHE      = 0                          #時刻帯別pNbrAltWキャッシュ（0:無効、1:有効）

[EXT_ANAL]
PROB_MODE    = 3                             # 閾値以上/以下の確率分析(0:なし、1：以下、2：以上、3:両方)
CENTROID     = 1                             # 重心分析(0:なし、1:実施)
AREA_STAT    = 1                             # 統計分析(0:なし、1:実施)
AREA_HIST    = 1                             # ヒストグラム分析(0:なし、1:実施)
PROB_ANNUAL  = 0                             # 全期間の確率分析を行う（0:行わない、1:行う）
OPT_FLOW     = 0                             # 時系列オプティカルフロー解析（0:無効、1:勾配法ベクトル解、2:Horn-Schunck法）

[INTERP]
METHOD          = 4                          # 気温空間補間の補間方法 (0=IDW, 1=RBF_GAUSS, 2=SHEPARD, 3=KRIGING_S未完成, 4=BARNES, 5=BARNES2PASS, 6=TPS)
COAST_L         = 20000.0                   # 海岸距離補正スケール [m]（0で補正なし）
IDW_POWER       = 1.5                        # IDW/Shepard べき乗値（デフォルト 1.5）
RBF_SIGMA       = 80000.0                    # RBF-Gauss バンド幅 [m]（デフォルト 80000）
BARNES_KAPPA_FAC = 2.0                       # Barnes: kappa = d_nn² × この値（デフォルト 2.0、文献値 5.052）
BARNES_GAMMA    = 0.3                        # Barnes 2Pass スキャン間収縮率（デフォルト 0.3）

[AREA_CONFIG]
TERRAIN_CORR           = 0                    # 地形補正（0:無効、1:有効）
TERRAIN_CORR_CLIP      = 0.0                  # 地形補正項クリップ幅[℃]（0以下で無効）
TERRAIN_CORR_COAST_MIN = 0.0                  # 適用最小海岸距離[m]（0以下で無効）
WINTER_CORR            = 0/0/0/0              # 季節別・時刻帯別の地形補正係数
SPRING_CORR            = 0/0/0/0
SUMMER_CORR            = 0/0/0/0
AUTUMN_CORR            = 0/0/0/0

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
# 季節開始月（春/夏/秋、冬は12月扱い）
SEASON = 3 / 6 / 9
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

各ファイルで定義される主要関数の一覧。全ての `static` helper までは列挙せず、処理フロー理解・保守・デバッグで参照頻度が高い関数を中心に記載する。引数の `DataDEM/DataAMeDAS/DataAnal/DataFlow` は原則として参照渡し。

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

> **削除済み（柱A A1、VER.0.4.6.175）**: `ValidateStationElevation` / `CheckPixelContext` は呼び出しのないデッドコードのため削除（計82行）。

#### `ScanFactorTPS`
- **処理概要**: TPS補間の正則化パラメータ λ を最適化スキャン
- **引数**: `const DataDEM& DEM`, `DataAnal& ANAL`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()（パラメータスキャンモード）

> **移設（柱A A4、VER.0.4.6.175）**: 以下の汎用出力4関数（`OutputTempGeoTIFF`/`OutputFlowGeoTIFF`/`OutputFlowGPKG`/`OutputAnalyzeBinary`）は **amt.cpp → AnalAMeDAS.cpp** へ移動済み（複数モード共通の汎用機能をライブラリ側へ集約）。`OutputAnalyzeBinary` は引数を `AppConfig`→`ConfigCore` に変更。

#### `OutputTempGeoTIFF`（→ AnalAMeDAS.cpp）
- **処理概要**: `mShtTemp` の気温面をGeoTIFFとして出力
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const ConfigCore& CFG`, `DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `OutputFlowGeoTIFF`
- **処理概要**: FlowU/FlowV/Magnitude/Diverg/Vorticity をGeoTIFFとして出力
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `const ConfigCore& CFG`, `DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `OutputFlowGPKG`
- **処理概要**: フロー解析結果をベクトルGPKGとして間引き出力
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `const ConfigCore& CFG`, `DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `OutputAnalyzeBinary`（→ AnalAMeDAS.cpp）
- **処理概要**: 解析結果を GRD バイナリ形式で出力
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `const ConfigCore& CFG`, `DateTime& dtWr`（AppConfig→ConfigCore に変更）
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `BuildPixelContext`（→ Interpolation.cpp、柱A A3）
- **処理概要**: 全陸地ピクセルの近傍16観測点情報（PixelContext）を構築。**amt.cpp → Interpolation.cpp** へ移動（3モード共通の補間準備をライブラリ側へ）、引数を `AppConfig`→`ConfigCore` に変更
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& AMD`, `DataAnal& AN`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: SetupForMonth(), RunParamScan(), RunCompare()

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

#### `RunOneMonth`（旧 RunOneBlock、static、VER.0.4.6.178 で月単位純化）
- **処理概要**: 1か月（または期間内の月チャンク）を処理する純粋関数。先頭で `SetupForMonth` を1回呼び、月内 nSteps ステップをループ。月境界判定・月次flush・バッファリセットは内部に持たず呼出側責務。月ローカル `s`（並列バッチ判定）とグローバル `step = nGlobalStepBase + s`（フロー有無/間引き位相/GRD判定）を使い分ける。
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `DataAMeDAS& AMD`, `const AppConfig& CFG`, `const ConfigRain& ITR_CFG`, `DivergRainAccum* pAccum`, `int nGlobalStepBase`, `DateTime& dtPrev`, `double& init_time_msec`
- **戻り値**: int — 0:成功 / -1:月初ロード失敗
- **呼び出し元**: RunContinuous(), RunMultiYear(), RunStNbrEval()

#### `SetupForMonth`（static、VER.0.4.6.177）
- **処理概要**: 月初の AMeDAS ロード＋補間準備（CalcInterpParams/BuildPixelContext/ChgInterpMethod）＋月初初期化（観測点ログ/エリアマスク/時別統計NC/解析概要）を一括実行
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataAMeDAS& AMD`, `const AppConfig& CFG`, `const DateTime& dt`, `double& init_time_msec`
- **戻り値**: int — 0:成功 / -1:LoadAMeDAS失敗
- **呼び出し元**: RunOneMonth()

#### `FlushMonth`（static、VER.0.4.6.177）
- **処理概要**: 月末/期間末の統計・確率出力を一括実行（閾値統計月別/時別、日別・月別の気温統計/ヒストグラム、確率マップ月別統計＋テキスト、`bProbImage` 時は確率マップ画像、TMD flush）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const AppConfig& CFG`, `const DateTime& dtMonth`, `const DateTime& dtPrev`, `bool bProbImage`
- **戻り値**: void
- **呼び出し元**: RunContinuous()（bProbImage=false）, RunMultiYear()（bProbImage=true）

#### `ResetAnalBuffers`（static、VER.0.4.6.178）
- **処理概要**: idx_order とエリア統計バッファ（vAreaStat/nHourlySteps/pHourlyBuf）を初期化。MULTI_YEAR は月ごと、CONTINUOUS は期間先頭で1回呼ぶ（AN.multi のリセットは RunOneMonth 先頭で別途実施）
- **引数**: `DataAnal& AN`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: RunContinuous(), RunMultiYear(), RunStNbrEval()

#### `RunContinuous`（期間動作）
- **処理概要**: `CONT_START`/`CONT_DUR` の期間を月チャンクに分割し `RunOneMonth` を反復。先頭で `ResetAnalBuffers` 1回、`nGlobalStepBase`/`dtPrev` を月跨ぎ引継ぎ（フロー連続性）。各チャンク末で `FlushMonth(false)`、期間末で `FlushAllTempStat`（TOTAL）。MONTH_OUT は無効（TOTAL_OUT で期間統計）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `DataAMeDAS& AMD`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `RunMultiYear`（月動作）
- **処理概要**: 同月・複数年を順次処理。各月で `ResetAnalBuffers`（月別独立）→ `RunOneMonth(nGlobalStepBase=0)` → `FlushMonth(true)`。終了後に全期間確率マップ＋全期間統計を出力
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `DataAMeDAS& AMD`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `SetStNbrEvalFixedConfig`
- **処理概要**: ST_NBR_EVAL用に通常出力・並列設定を固定し、観測点近傍評価向けの条件にする
- **引数**: `AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main(), RunStNbrEval()

#### `RunStNbrEval`
- **処理概要**: 指定年月日・時刻だけを処理し、観測点近傍・地形特徴量・補正項確認用の出力を行う
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
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataAMeDAS& DT`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `RunEvalInterp`
- **処理概要**: 指定日時の観測点LOOCVを実行し、補間手法別の誤差・地域帯別指標を出力
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataAMeDAS& DT`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `RunOptimParams`
- **処理概要**: Barnes/RBF/Shepard/Barnes2Passの補間パラメータをグリッドスキャンし、LOOCV最適値を出力
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataAMeDAS& DT`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `RunAreaChar`
- **処理概要**: 解析エリアの近傍距離・標高差・近傍点数などをGeoTIFF/CSVで出力
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataAMeDAS& DT`, `const AppConfig& CFG`
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

#### `EvalDivergRain`
- **処理概要**: 観測点位置の発散値と降水量を照合し、TP/FP/FN/TN・相関・降水有無別統計を算出
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& AMD`, `DataAnal& AN`, `DataFlow& FL`, `DivergRainResult* pR`, `float fDivergThreshold`, `float fRainThreshold`
- **戻り値**: void
- **呼び出し元**: RunOneMonth(), RunParamScan()

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
- **処理概要**: GeoTIFF から DEM（標高・陸マスク）を読み込む。海岸距離は呼び出し側で `CalcCoastDist()` を実行する
- **引数**: `DataDEM& DEM`, `const char* pPath`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: main()

#### `FreeTerrainDEM`
- **処理概要**: 地形補正用DEMのメモリを解放
- **引数**: `DataDEM& DEM`
- **戻り値**: void
- **呼び出し元**: FreeDEM(), InitTerrainDEM()

#### `InitTerrainDEM`
- **処理概要**: `TERRAIN_CORR` 用のDEMを読み込み、失敗時は警告を出して地形補正を無効化
- **引数**: `DataDEM& DEM`, `ConfigCore& CFG`
- **戻り値**: int — 0:成功または無効化, -1:失敗
- **呼び出し元**: main()

#### `CalcTerrainFeatureAtLonLat`
- **処理概要**: 指定経緯度周辺の地形特徴量（中心標高と周辺低位標高の差）を算出
- **引数**: `const DataDEM& DEM`, `double lon`, `double lat`
- **戻り値**: float — 地形特徴量
- **呼び出し元**: InitAnal(), LoadAMeDAS()

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
- **呼び出し元**: RunOneMonth()

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
- **呼び出し元**: RunOneMonth()

#### `GetRunMonthString`
- **処理概要**: 実行月設定から月文字列（例: "1-3"）を生成
- **引数**: `const ConfigRun& RCFG`
- **戻り値**: CHsString
- **呼び出し元**: 各種初期化・出力関数

#### `GetCommonOutputFileName`
- **処理概要**: 結果パス、ファイルキー、実行月条件を使って共通出力ファイル名を生成
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`, `const char* pszName`
- **戻り値**: CHsString
- **呼び出し元**: 各種InitOutput*(), OutputAnalysisSummary()

#### `OutputAnalysisSummary`
- **処理概要**: 実行条件、DEM、補間、出力ファイル一覧、エリア情報などの概要TSVを出力
- **引数**: `const DataDEM& DEM`, `const DataAMeDAS& AMD`, `const DataAnal& AN`, `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

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
- **呼び出し元**: RunOneMonth()

#### `InitOutputTempStatDay`
- **処理概要**: 日別気温統計ファイルを初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputTempStatDay`
- **処理概要**: 日別気温統計を TSV 出力してバッファリセット
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `InitOutputTempStatMth`
- **処理概要**: 月別気温統計ファイルを初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputTempStatMth`
- **処理概要**: 月別気温統計を TSV 出力してバッファリセット
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

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
- **呼び出し元**: RunOneMonth()

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
- **呼び出し元**: RunOneMonth()

#### `FlushProvMonthToAccum`
- **処理概要**: 月次確率バッファを全期間バッファに加算してリセット（AVX-512）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()（月末）

#### `InitOutputProbText`
- **処理概要**: 確率マップ出力用 TSV ファイル（4種）をオープンしてヘッダを書き込み
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `CalcProbMonthStat`
- **処理概要**: 月次/全期間の確率統計（平均/最大/ヒストグラム）を計算
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `int mode` — 月次(0)/全期間(1), `int dir` — lo(0)/hi(1)
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `OutputProbText`
- **処理概要**: `CalcProbMonthStat()` の結果を月次/全期間TSVとして出力
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`, `int dir`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `OutputProbImage`
- **処理概要**: 月次/全期間の確率マップをGeoTIFF/PNGとして出力
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `InitOutputBoundCentroid`
- **処理概要**: 境界面重心出力を初期化。現行のTXT出力は廃止され、NetCDFは `InitHourlyStatNC()` で初期化する
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputBoundCentroid`
- **処理概要**: 全域/POI別の境界面重心を BFS（8近傍）で計算してNetCDFへ出力
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `OutputBoundMaskGeoTIFF`
- **処理概要**: `mBoundPacked` の境界マスクをGeoTIFFとして出力
- **引数**: `const DataDEM& DEM`, `const DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `OutputBoundContourPNG`
- **処理概要**: 等温線境界マスクを8bit PNGとして出力
- **引数**: `const DataDEM& DEM`, `const DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `OutputMergedTempContourPNG`
- **処理概要**: 気温分布と等温線を合成したPNGを出力
- **引数**: `const DataDEM& DEM`, `const DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `OutputTempDistPNG`
- **処理概要**: `mShtTemp` の気温分布を8bit PNGとして出力
- **引数**: `const DataDEM& DEM`, `const DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

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

#### `InitOutputStationTemprature`
- **処理概要**: 観測点別時別データNetCDFと、必要に応じて近傍ダンプ出力を初期化
- **引数**: `const DataDEM& DEM`, `const DataAMeDAS& AMD`, `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `OutputStDatHour`
- **処理概要**: 観測点別の実測/補間/近傍/補正項データを時別NetCDFへ出力
- **引数**: `const DataDEM& DEM`, `const DataAMeDAS& AMD`, `const DataAnal& AN`, `const ConfigCore& CFG`, `int y`, `int m`, `int d`, `int h`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `InitHourlyStatNC`
- **処理概要**: 時別気温統計、閾値統計、境界重心などのNetCDF出力を初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`, `const DataAnal& AN`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

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
- **呼び出し元**: RunOneMonth()

#### `InitOutputTempHistDay`
- **処理概要**: 日別気温ヒストグラムファイルを初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputTempHistDay`
- **処理概要**: 日別気温ヒストグラム生データを TSV 出力
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `InitOutputTempHistMth`
- **処理概要**: 月別気温ヒストグラムファイルを初期化
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `OutputTempHistMth`
- **処理概要**: 月別気温ヒストグラム生データを TSV 出力
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

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
- **呼び出し元**: RunOneMonth()

#### `InitInterpTPSWork`（旧 `InitTPS`、VER.0.4.6.172 で static 化）
- **処理概要**: TPS 補間用作業バッファ（tps_w, tps_stLon, tps_stLat）を確保・初期化
- **引数**: `const DataAMeDAS& DT`, `DataAnal& AN`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: `ChgInterpMethod()`（TPS選択時のみ）。旧 public `InitTPS` を Interpolation.cpp 内 static に統合し、評価モード・通常運用とも ChgInterpMethod 経由に一本化（呼び出し漏れ防止）

#### `CalcBoundaryMaskCoreInt16`
- **処理概要**: int16 気温場から各閾値の境界マスクを生成する内部ヘルパー（AVX-512）
- **引数**: `const DataDEM& DEM`, `const DataAnal& AN`, `const int16_t* pSht`, `uint16_t* pPacked`
- **戻り値**: void
- **呼び出し元**: CalcBoundaryMaskAndCollect()

#### `CalcBoundaryMaskAndCollect`
- **処理概要**: 境界マスク生成とエリア統計集計を同時実行
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`
- **戻り値**: void
- **呼び出し元**: TemparatureInterpolation(), RunOneMonth()

#### `CalcBoundaryMaskAndCollectMulti`
- **処理概要**: N 時刻の境界マスク生成・エリア統計集計を一括処理
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const int* pBatchSlots`, `int nBatch`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

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
- **呼び出し元**: RunOneMonth(), main()（パラメータスキャン内）

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
- **呼び出し元**: RunOneMonth(), RunParamScan()

#### `Calc3PointGradientFlow`
- **処理概要**: 3時刻の気温場を使った勾配法（2次精度中心差分による時間微分）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `float alpha_sq`
- **戻り値**: void
- **呼び出し元**: RunOneMonth(), RunParamScan()（TIME_SLICES ≥ 3 時）

#### `CalcHornSchunckFlow`
- **処理概要**: Horn-Schunck 法オプティカルフロー（反復法で滑らかなフロー場を生成）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `float alpha_sq`, `int nIter` — 反復回数
- **戻り値**: void
- **呼び出し元**: RunOneMonth(), RunParamScan()

---

## 10. 設定キー一覧

設定ファイルはINI形式。未記載キーは `InitConfigCore()` / `InitConfigRun()` のデフォルトを維持する。値域は読み込み時に `HsClip()` で丸めるものが多い。ここでは運用・保守で重要な主要キーを記載する。

### 10-1. 基本・パス

| セクション | キー | 型 | デフォルト/範囲 | 対応メンバ | 内容 |
|---|---|---:|---|---|---|
| MAIN | TYPE | string | main想定 | `ConfigCore::szType` | 設定ファイル種別。通常実行は `main` |
| PATH | DEM_PATH | path | 環境依存 | `szDemPath` | 解析用DEM GeoTIFF |
| PATH | TERRAIN_DEM_PATH | path | DEM_PATH | `szTerrainDemPath` | 地形補正用DEM。`TERRAIN_CORR=1`時に必須 |
| PATH | RESULT_PATH | path | 環境依存 | `szResultPath` | 出力先ディレクトリ |
| PATH | FILE_KEY | string | 環境依存 | `szFileKey` | 画像/GRD/GeoTIFF系のファイルキー |
| PATH | FILE_KEY_ST | string | 環境依存 | `szStatFileKey` | 統計/補助出力系のファイルキー |
| PATH | AMEDAS_PATH | path | 環境依存 | `szAMeDASPath` | AMeDAS年報バイナリ格納ルート |

### 10-2. 実行制御

| セクション | キー | 型 | デフォルト/範囲 | 対応メンバ | 内容 |
|---|---|---:|---|---|---|
| RUN | MODE | string | CONTINUOUS等 | `ConfigRun::nRunMode` | `CONTINUOUS` / `MULTI_YEAR` / `PARAM_SCAN` / `COMPARE` / `VALUE_SURVEY` / `EVAL_INTERP` / `OPTIM_PARAMS` / `AREA_CHAR` / `ST_NBR_EVAL` |
| RUN | CONT_START | y/m/d/h | 2000/1/1/0相当 | `nStartYear`等 | 連続実行開始日時 |
| RUN | CONT_DUR | hour/day | 設定値 | `nRunCount`, `nRunDays` | 連続実行期間。日指定を優先 |
| RUN | MULTI_YEAR | y0/y1 | 設定値 | `nYearStart`, `nYearEnd` | 複数年実行の年範囲 |
| RUN | MULTI_MONTH | list | 1-12 | `nTargetMonthList` | 対象月リスト。`1/2/3/12`形式 |
| RUN | DAYS | list | ST_NBR_EVAL用 | `nTargetDayList` | ST_NBR_EVAL対象日 |
| RUN | HOURS | list | ST_NBR_EVAL用 | `nTargetHourList` | ST_NBR_EVAL対象時刻 |
| RUN | EVAL_METHODS | list | 0-6 | `nEvalMethods` | EVAL_INTERP対象補間手法 |
| RUN | PARALLEL_SLICES | int | 1 / 1-8 | `ConfigCore::nParallelSlices` | 複数時刻同時補間数。日境界・時刻帯境界で分割 |
| RUN | USE_NBR_FLOAT | int | 3 / 0-3 | `nNbrFloatMode` | 0=全整数、1=SRCのみfloat、2=Wのみfloat、3=全float |
| RUN | ALTW_CACHE | int | 0 / 0-1 | `nAltWCache` | 時刻帯別 `pNbrAltW` キャッシュ |

### 10-3. 解析・補間

| セクション | キー | 型 | デフォルト/範囲 | 対応メンバ | 内容 |
|---|---|---:|---|---|---|
| EXT_ANAL | PROB_MODE | int | 0 / 0-3 | `nProveMode` | 0=無効、1=閾値以下、2=閾値以上、3=両方 |
| EXT_ANAL | CENTROID | int | 0 / 0-1 | `nCentroidMode` | 境界面重心解析 |
| EXT_ANAL | AREA_STAT | int | 0 | `nAreaStatMode` | 閾値エリア統計 |
| EXT_ANAL | AREA_HIST | int | 0 / 0-1 | `nAreaHistMode` | 気温ヒストグラム集計 |
| EXT_ANAL | OPT_FLOW | int | 0 / 0-2 | `nOptFlowMode` | 0=無効、1=勾配法、2=Horn-Schunck |
| EXT_ANAL | PROB_ANNUAL | int | 0 / 0-2 | `nProbAnnual` | 全期間確率バッファ確保。省メモリ時は0 |
| INTERP | METHOD | int | 4 / 0-6 | `nInterpMethod` | 0=IDW、1=RBF、2=Shepard、3=未完成、4=Barnes、5=Barnes2Pass、6=TPS |
| INTERP | COAST_L | float | 20000相当 | `AppConfig::fCoastL` | 海岸距離補正スケール[m] |
| INTERP | IDW_POWER | float | 1.5 | `fIdwPower` | Shepard/IDW系べき乗 |
| INTERP | RBF_SIGMA | float | 80000 | `fRbfSigma` | RBF sigma[m] |
| INTERP | BARNES_KAPPA_FAC | float | 2.0 | `fBarnesKappaFac` | `kappa = d_nn^2 * factor` |
| INTERP | BARNES_GAMMA | float | 0.3 | `fBarnesGamma` | Barnes2Passの第2パス係数 |
| PARAM | DIVERG_RAIN | int | 0 | `nDivergRain` | 発散・降水照合 |
| PARAM | DIVERG_THRESHOLD | float | -2.2 | `fDivergThreshold` | 収束判定閾値 |
| PARAM | RAIN_THRESHOLD | float | 1.5 | `fRainThreshold` | 降水判定閾値[mm/h] |

### 10-4. 出力制御

| セクション | キー | 型 | デフォルト/範囲 | 対応メンバ | 内容 |
|---|---|---:|---|---|---|
| OUTPUT | DEBUG | int | 0 / 0-1 | `nDebug` | デバッグ出力 |
| OUTPUT | COMMENT | int | 1 / 0-1 | `nComment` | 条件ログ出力 |
| OUTPUT | ANAL0 | int | 1 / 0-1 | `nOutputAnal0` | 実行時間表示 |
| OUTPUT | ANAL1 | int | 0 / 0-1 | `nOutputAnal1` | 数値評価表示 |
| OUTPUT | ANAL2 | int | 0 / 0-1 | `nOutputAnal2` | coast_L評価等の補助 |
| OUTPUT | ST_DAT | int | 0 / 0-1 | `nOutputStDat` | 観測点別NetCDF |
| OUTPUT | ST_NBR_DUMP | list | 任意 | `nStNbrDumpDT` | 近傍点ダンプ日時 `yyyymmddhh/...` |
| OUTPUT_VECTOR | FLOW | int | 0 / 0-1 | `nOutputFlowGPKG` | GPKGベクトル出力 |
| OUTPUT_VECTOR | THINNING | int | 6 / 1-1000 | `nThinningGPKG` | GPKG時間間引き |
| OUTPUT_VECTOR | THINNING_XY | int | 20 / 1-1000 | `nThinningGPKGXY` | GPKG空間間引き |
| OUTPUT_GEOTIFF | TEMP | int | 0 / 0-1 | `nOutputTempGeoTIFF` | 気温GeoTIFF |
| OUTPUT_GEOTIFF | FLOW | int | 0 / 0-1 | `nOutputFlowGeoTIFF` | フローGeoTIFF |
| OUTPUT_GEOTIFF | GEOTIFF3 | int | 0 / 0-1 | `nOutputGeoTIFF3` | 確率GeoTIFF |
| OUTPUT_GEOTIFF | GEOTIFF4 | int | 0 / 0-1 | `nOutputGeoTIFF4` | 境界マスクGeoTIFF |
| OUTPUT_GEOTIFF | THINNING | int | 6 / 1-1000 | `nThinningGeoTIFF` | GeoTIFF時間間引き |
| OUTPUT_GEOTIFF | THINNING_XY | int | 1 / 1-1000 | `nThinningGeoXY` | GeoTIFF空間間引き |
| OUTPUT_PNG | PNG2 | int | 0 / 0-2 | `nOutputPNG2` | 0=無効、1=気温+等温線別、2=マージ |
| OUTPUT_PNG | TEMP_MIN/MAX | int | -40/40 | `nDistTempMin/Max` | 気温PNGレンジ |
| OUTPUT_PNG | PNG3 | int | 0 / 0-1 | `nOutputPNG3` | 確率PNG |
| OUTPUT_PNG | PNG4 | int | 0 / 0-1 | `nOutputPNG4` | 等温線PNG。`PNG2=1`時のみ有効 |
| OUTPUT_GRD | GRD | int | 1 / 0-1 | `nOutputGRD` | GRDバイナリ出力 |
| OUTPUT_GRD | THINNING | int | 1 / 1-1000 | `nThinningGRD` | GRD時間間引き |
| OUTPUT_GRD | CH_NUM / CH_TAGn | int/string | 任意 | `nGrdChNum`, `szGrdChTags` | GRD出力チャンネル指定 |

`TEMP_STAT.HOUR_OUT`、`TEMP_HIST.HOUR_OUT`、`THRESHOLD_ANAL.HOUR_OUT` は 0=無効、1=TXT、2=NetCDF、3=TXT+NetCDF。日別/月別/全期間系は現行ではTXT出力フラグとして扱う。

### 10-5. 閾値・エリア・補正

| セクション | キー | 型 | デフォルト/範囲 | 対応メンバ | 内容 |
|---|---|---:|---|---|---|
| THRESHOLD | VALUE | list | 最大16 | `nvThresholds` | 昇順の気温閾値リスト。`,`区切り |
| AREA | NUM | int | 0 / 0-32 | `nAreas` | エリア数 |
| AREA | AREA_n | string | 任意 | `areaDefs` | 矩形 `name,lat_s,lat_n,lon_w,lon_e` またはSHP `name,shp:path,field,value` |
| POI_CENTROID | COUNT | int | 0 / 0-16 | `nPOIs` | POI数 |
| POI_CENTROID | POI_n | string | 任意 | `poiDefs` | `name,lat,lon` |
| POI | SEARCH_R | int | 50 | `nPoiSearchR` | 近傍境界探索半径[px] |
| LAPSE | LAT | list | 33/40/999/999 | `lapse_lat_zone` | 緯度帯境界 |
| LAPSE | SEASON | list | 3/6/9 | `lapse_*_month` | 春/夏/秋開始月 |
| LAPSE | HOUR | list | 未指定 | `lapse_hour_th` | 4時刻帯使用時の境界 |
| LAPSE | VALUEz / VALUEz_b | list | デフォルト減率 | `lapse_rate_table` | 季節別/時刻帯別の気温減率[℃/m] |
| AREA_CONFIG | TERRAIN_CORR | int | 0 / 0-1 | `nTerrainCorr` | 地形特徴量補正 |
| AREA_CONFIG | TERRAIN_CORR_CLIP | float | 4.0 | `fTerrainCorrClip` | 地形補正項クリップ幅[℃] |
| AREA_CONFIG | TERRAIN_CORR_COAST_MIN | float | 0.0 | `fTerrainCorrCoastMin` | 適用最小海岸距離[m] |
| AREA_CONFIG | SEASON / HOUR | list | LAPSE準拠 | `terrain_*` | 地形補正用の季節/時刻帯 |
| AREA_CONFIG | WINTER_CORR等 | list | 0 | `terrain_corr_table` | 季節別/時刻帯別の地形補正係数 |

---

## 11. 出力ファイル仕様

出力ファイル名は原則として `GetCommonOutputFileName()` で生成され、`RESULT_PATH`、`FILE_KEY_ST`、実行年/月条件、末尾名を組み合わせる。画像・GeoTIFF・GRD等の一部は `FILE_KEY` を使う。

### 11-1. 共通・補助出力

| ファイル | 生成条件 | 内容 |
|---|---|---|
| `*_analysis_summary.txt` | 初回エリアマスク構築後 | 実行条件、DEM情報、補間設定、出力ファイル一覧、エリア陸地ピクセル数 |
| `*_area_land_pixels.txt` | `AREA.NUM>0` | 全域と各エリアの陸地ピクセル数。全域は `AREA_DEF_MAX` |
| `*_stdat.nc` | `OUTPUT.ST_DAT=1` またはST_NBR_EVAL | 観測点別の実測値、補間値、近傍、地形特徴量、補正項 |
| `*_manifest.txt` | ST_NBR_EVAL | ST_NBR_EVALで処理した日時と出力先の対応 |
| `*_value_survey_YYYYMM.txt` | `MODE=VALUE_SURVEY` | pNbr*等の値域調査 |
| `*_area_char_YYYYMM_*.tif` | `MODE=AREA_CHAR` | d_nn, mean_d, cv_d, alt_std, alt_range, nbr_n |
| `*_area_char_YYYYMM_summary.txt` | `MODE=AREA_CHAR` | エリア特性の統計サマリ |
| `*_area_char_YYYYMM_stations.txt` | `MODE=AREA_CHAR` | 観測点位置でのエリア特性サンプル |
| `*_eval_loocv_YYYYMMDD_HH.txt` | `MODE=EVAL_INTERP` | 観測点別LOOCV結果 |
| `*_eval_summary_YYYYMMDD_HH.txt` | `MODE=EVAL_INTERP` | 手法別・地域帯別のLOOCVサマリ |
| `*_eval_optimize_YYYYMMDD_HH.txt` | `MODE=OPTIM_PARAMS` | パラメータスキャン全結果 |
| `*_optim_result_YYYYMMDD_HH.txt` | `MODE=OPTIM_PARAMS` | 推奨パラメータ結果 |

### 11-2. 気温統計・ヒストグラム

| ファイル | 生成条件 | 内容 |
|---|---|---|
| `*_temp_stat_hourly.txt` | `TEMP_STAT.HOUR_OUT=1 or 3` | 時別の平均/最低/最高/中央値。全域+エリア別 |
| `*_temp_stat_hourly.nc` | `TEMP_STAT.HOUR_OUT=2 or 3` | 時別気温統計NetCDF |
| `*_temp_stat_daily.txt` | `TEMP_STAT.DAY_OUT=1` | 日別気温統計 |
| `*_temp_stat_monthly.txt` | `TEMP_STAT.MONTH_OUT=1` | 月別気温統計 |
| `*_temp_stat_total.txt` | `TEMP_STAT.TOTAL_OUT=1` | 全期間気温統計 |
| `*_temp_hist_hourly.txt` | `TEMP_HIST.HOUR_OUT=1 or 3` | 時別気温ヒストグラム。BIN幅/範囲は `TEMP_HIST` |
| `*_temp_hist_hourly.nc` | `TEMP_HIST.HOUR_OUT=2 or 3` | time × area × bin のNetCDF |
| `*_temp_hist_daily.txt` | `TEMP_HIST.DAY_OUT=1` | 日別累積ヒストグラム |
| `*_temp_hist_monthly.txt` | `TEMP_HIST.MONTH_OUT=1` | 月別累積ヒストグラム |
| `*_temp_hist_total.txt` | `TEMP_HIST.TOTAL_OUT=1` | 全期間累積ヒストグラム |

時別気温統計は `mShtTemp` の `int16_t` 値を `HSS_GRD_FAC_TEMP=500` で戻した値を基準に集計する。中央値は0.1℃内部BINのヒストグラムから簡易算出する。

### 11-3. 閾値・確率・重心

| ファイル | 生成条件 | 内容 |
|---|---|---|
| `*_threshold_stat_hourly.txt` | `THRESHOLD_ANAL.HOUR_OUT=1 or 3` | 時別の閾値以下ピクセル比率/ピクセル数。閾値 × 全域/エリア |
| `*_threshold_stat_hourly.nc` | `THRESHOLD_ANAL.HOUR_OUT=2 or 3` | 時別閾値統計NetCDF。NetCDFではcntも保存 |
| `*_threshold_lo_stat_monthly.txt` | `THRESHOLD_ANAL.LO_OUT=1` | 月別T<閾値統計。発生時間、平均/最大面積、消滅/出現 |
| `*_threshold_hi_stat_monthly.txt` | `THRESHOLD_ANAL.HI_OUT=1` | 月別T≥閾値統計 |
| `*_prob_monthly.txt` | `PROB_MODE=1 or 2` | 月別確率統計 |
| `*_prob_histgram.txt` | `PROB_MODE=1 or 2` | 月別確率ヒストグラム |
| `*_prob_lo_monthly.txt` / `*_prob_hi_monthly.txt` | `PROB_MODE=3` | lo/hi別の月別確率統計 |
| `*_prob_lo_histgram.txt` / `*_prob_hi_histgram.txt` | `PROB_MODE=3` | lo/hi別の確率ヒストグラム |
| `*_bound_centroid.nc` | `CENTROID=1` | 全域/POI別の境界面重心NetCDF |

確率マップ内部配列は `pixel * nCntThresholds + threshold`。`PROB_MODE=2` では `pCountsM/pCountsA` を閾値以上側として兼用し、`PROB_MODE=3` のみ `pCountsM_hi/pCountsA_hi` を別確保する。

### 11-4. 画像・GIS・バイナリ

| 出力 | 生成条件 | ファイル名例/内容 |
|---|---|---|
| 気温PNG | `OUTPUT_PNG.PNG2=1` | `TempPNG\YYYY\temp_KEY_YYYYMMDD_HH.png` |
| 等温線PNG | `OUTPUT_PNG.PNG2=1 && PNG4=1` | `BoundMask\ContourPNG\..\bound_ctr_KEY_YYYYMMDD_HH.png` |
| 気温+等温線PNG | `OUTPUT_PNG.PNG2=2` | `TempContourPNG\YYYY\tempctr_KEY_YYYYMMDD_HH.png` |
| 確率PNG | `OUTPUT_PNG.PNG3=1 && PROB_MODE>0` | `prob_KEY_YYYYMM_+nnC_lo/hi.png` |
| 気温GeoTIFF | `OUTPUT_GEOTIFF.TEMP=1` | `KEY_Temparature_YYYYMMDD_HH.tif` |
| Flow GeoTIFF | `OUTPUT_GEOTIFF.FLOW=1` | `KEY_FlowU_*.tif`, `KEY_FlowV_*.tif`, `KEY_Diverg_*.tif`, `KEY_Vorticity_*.tif` |
| 境界マスクGeoTIFF | `OUTPUT_GEOTIFF.GEOTIFF4=1` | `BoundMask\...\bound_mask_KEY_YYYYMM_DD_HH_+nnC.tif` |
| 確率GeoTIFF | `OUTPUT_GEOTIFF.GEOTIFF3=1 && PROB_MODE>0` | `prob_KEY_YYYYMM_+nnC_lo/hi.tif` |
| フローGPKG | `OUTPUT_VECTOR.FLOW=1 && OPT_FLOW>0` | `KEY_Vector_YYYYMMDD_HH.gpkg` |
| GRD | `OUTPUT_GRD.GRD=1` | `BIN_TEMP` または `BIN_DATA` 配下。TEMP/FLOW/BOUND等をBLOSC2圧縮格納 |

GRDのTEMPは `int16_t`、スケールは `HSS_GRD_FAC_TEMP=500`（0.002℃精度）。BOUNDは `HSS_GRD_TYPE_B16` で、bit c が `THRESHOLD.VALUE[c]` の境界マスクを表す。

---

*以上、VER.0.4.6.140時点のシステム全体解説。*

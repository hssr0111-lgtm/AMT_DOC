# AMT(Advanced Mesh Temperature）解析システム — 全体詳細解説
**VER.0.4.6.500　2026/06/11**

本書は `system_overview_V046_20260425.md` を元に、2026/06/11時点の現行コード（`amt.cpp` main内の `VER.0.4.6.494`）へ合わせて修正・追記したもの。VER.0.4.6.413 の LOOCV一本化・COMPARE廃止・coast_L評価のOPTIM_PARAMS統合・ANAL2廃止は **§1-9**、その後の保守リファクタ（pCFG完全廃止・引数 ConfigCore 統一・重複共通化・デッドコード削除等、〜VER.0.4.6.460）は **§1-10**、リファクタリング4B（`DataAnal` 3サブ構造体化＝pixctx/tps/ring・`nvThresholdsScaled`改名・`MergeHistU32AVX512`集約・6-1/7-1/7-4 等、VER.0.4.6.470〜475、`system_review_V0464_20260606.md` 準拠）は **§12-11**、リファクタリング6（画像名のFILE_KEY先頭統一・スペル統一・interp集約・境界マスク統一・DataAnal static化・非フローGRD位相・HS magnitude外周・int16 PNG×SIMDクラッシュ修正・SavePNG8汎用化・回帰ハーネス新設、VER.0.4.6.480〜488、`system_review_V0464_20260607.md` 準拠）は **§12-12**、PNG出力本体のクラス化・外出し（`HsOutputPNG.h`/`.cpp` 新設・`CHsOutputPNG`/`COutputPNG`、VER.0.4.6.491〜494）は **§12-13** に集約して記載する。

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
- 現行方針は **「amt.cpp＝全体コントロール（main/設定読込/実行制御＋モード別 Run* と実行モード固有評価）／AnalCore.cpp＝複数モード共通のデータ処理・解析補助／Output.cpp＝出力サブシステム」** の責務境界明確化。**amt.cpp の完全解体は行わない**。`Run*`、`LoadConfig`、`ApplyRunModeFixedConfig`、ST_NBR_EVALのmanifest、LOOCVゾーン統計の共有部品（`ZoneStat`/`ZoneStatAccum`/`ZoneStatFinal`/`ZoneClassifier`/`CalcLaplacianSmoothness`）等の特定モード専用処理は amt.cpp に残す。**※ VER.0.4.6.413 で `EvalInterp`(static)/`ComparInterpolationMethod`/`Eval_Coast_L*`/`ScanFactorTPS`/`PrintEvalResult` は削除し、LOOCV は本番補間ベースの単一系統（`RunEvalInterp`/`RunOptimParams`）へ統合した（§1-9）。**
- `SetupForMonth`/`ResetAnalBuffers`/`FlushMonth` は複数モード共通の月次ライフサイクル処理として AnalCore.cpp に配置する。`EvalDivergRain` 関連は解析本筋ではないが複数モードから利用し得る妥当性評価処理として AnalCore.cpp に残す。
- VER.0.4.6.300 で出力系を `Output.cpp` に分離。TXT/NetCDF/PNG/GeoTIFF/GPKG/GRD、出力ディレクトリ・ファイル名ヘルパ、出力系グローバルは Output.cpp 側へ集約した。さらに `AnalAMeDAS.cpp` を `AnalCore.cpp`、`AnalAMeDAS.h` を `AMTCore.h` に名称変更し、AMeDAS専用に見える名称依存を弱めた。現行処理はAMeDAS時別データを主対象とする一方、原理上はAMeDAS年報バイナリ形式に準拠した観測点・観測データにも適用できる。

### 0-5. 近傍点テーブルの int版/float版 二重フィールド維持は意図的設計（レビュー対象外）

- **「最低限のメモリで動作すること」を前提とする本ソフトウェアの基本思想に基づき、`NbrTable` の int版/float版 二重フィールド
  （`pNbrAltDiff/F`・`pNbrCoastPx/F`・`pNbrCoastStn/F`・`pNbrW/F`・`pNbrAltW/F`）と `USE_NBR_FLOAT`(`exec.nNbrFloatMode`) による切替は、
  保守コスト増加を許容したうえで意図的に継続する確定方針である。**
- 理由: 近傍点情報は解析条件によらず**最大メモリ消費となる部分**であり、整数化による省メモリ効果が大きい。家庭用PC（16GB級）での
  広域エリア・詳細メッシュ解析では、この整数化が「解析可能/不可能」を分ける現実的な手段になる。
  - 参考（README 留意点より、4000×4000全陸地解析）: 全float 約6.9GB に対し、`USE_NBR_FLOAT=2`（一部int）で約1GB、`=0`（全int）で約2.1GB のメモリ節約。
  - 精度影響も評価済みで、温度統計・ヒストグラムで実用上問題ない範囲（README「整数モードに関して」参照）。
- **したがって、本二重フィールド維持・`exec.nNbrFloatMode` 分岐の「保守負担」は今後のコードレビューで減点・指摘の対象としない。**
  （`pNbr*` 読取は `NbrGetCoastPx`/`NbrGetCoastStn`/`NbrGetAltDiff`（`AMTCore.h`）等のヘルパに集約し、分岐の散在を抑える設計とする。）
- 関連: §0-1（C言語スタイルの計算コア）。READMEの「整数モードに関して」「DEMデータ変換」節に運用上の指針を記載。

### 0-6. 地形補正用DEMの `DataDEM` インライン同居は意図的設計（レビュー対象外）

- **`DataDEM` は「DEMに由来する空間情報を一括して保持するコンテナ」と位置づける思想に基づき、地形補正用DEM
  （`terrain_dLonOrigin`/`terrain_dResX`/`terrain_dLatOrigin`/`terrain_dResY`/`terrain_nWidth`/`terrain_nHeight`/
  `terrain_nAlignedWidth`/`terrain_total_allocated_pix`/`pTerrainElevation`/`pTerrainMask`）を本体DEMと同階層にインライン同居させる。
  これは確定方針であり、`TerrainDEM` サブ構造体への切り出しは行わない。**
- 理由: 本体DEM（標高・陸海マスク・海岸距離）と地形補正用DEMは、ともに「同一解析空間に対するDEM系入力」であり、
  解析全体を通じて一括ロード（`InitDEM`/`InitTerrainDEM`）・一括解放（`FreeDEM`/`FreeTerrainDEM`）・複数 cfg 連続実行時の
  差し替え判定（`amt.cpp` の DEMパス/TERRAIN_CORR 変化検出）を**1つの `DataDEM` 単位で扱う**ことが運用上自然なため。
  サブ構造体化すると引数経路（`DEM.terrain_*` → `DEM.terrain.*`）の波及と、DEMコンテナとしての一体性低下を招く。
- 補足: 地形補正用DEMは解析DEM本体と**解像度・範囲が異なってよい**（解析エリア外近傍点の地形特徴量取得のため範囲外も含めて作成する）。
  名称プレフィックス `terrain_*` で本体DEMフィールドと明確に区別しており、混線リスクは低い。
- **したがって、`DataDEM` への地形補正DEMインライン同居（サブ構造体化の見送り）は今後のコードレビューで減点・指摘の対象としない。**
- 関連: §1-7（地形補正方式 TERRAIN_CORR=1/2）、§0-5（同様に省メモリ/コンテナ一体性を優先した確定方針）。

---

## 1. ファイル構成

| ファイル | 行数 | 役割 |
|:---|---:|:---|
|`amt.cpp`|3217|AMT実行制御（main、設定読込、Run*、実行モード固有評価/デバッグ処理）|
|`AMTCore.h`|1691|全構造体・定数・インライン関数・関数宣言の定義|
|`AnalCore.cpp`|4331|設定Core読込、AMeDAS形式データ読込、DEM/地形/共通解析ヘルパ、月次ライフサイクル、DivergRain評価|
|`Output.cpp`|4154|出力サブシステム（TXT/NetCDF/PNG/GeoTIFF/GPKG/GRD、出力名・出力ディレクトリ、出力系状態）。PNG本体は HsOutputPNG.h/.cpp ＋ 派生 COutputPNG へ移設|
|`HsOutputPNG.h`|122|PNG出力本体クラス CHsOutputPNG 宣言（AMT非依存。Save8/Save16/SaveByBuilder8/SaveIndexed8、overlayRow/buildRow8仮想）|
|`HsOutputPNG.cpp`|514|CHsOutputPNG 実装（libpng定型集約・float/int16×SIMD/スカラー行変換）|
|`Interpolation.cpp`|3180|補間計算・近傍点テーブル構築・エリア統計集計|
|`OpticalFlow.cpp`|686|勾配フロー・Horn-Schunck実装（AVX-512）|
|`LogStDatNetCDF.h`|479|NetCDF出力クラス定義、観測点計算結果時別データ出力|
|`HourlyStatNetCDF.h`|1086|NetCDF出力クラス定義、時別データ出力|
|`HsCommon.h`|1667|汎用クラス（CHsString/CHsIniFile/CHsLog/CHsCrono/CHsTextRenderer）|
|`WriteGeoTIFF.h`|626|GeoTIFF出力（GDAL）|
|`WriteGRD.h`|1515|独自バイナリ(GRD)出力・読込、BLOSC2圧縮|
|`ConvMatrix.h`|434|float/int16変換・カラーマップ適用・PNG行生成|
|`TMD.h`|705|汎用時系列データロガー形式の読み書き|
|`Nifsq8.cpp/.h`|1059|補助モジュール（データ圧縮補助）|
|`palettes.h`|534|8bitPNG用パレット定義|
|`stb_truetype.h`|5079|STB TrueTypeフォントラスタライザ（外部ライブラリ、CHsTextRenderer が使用）|
|`Constants.h`|42|全体定数の補助定義（Kriging安定化・旧評価定数削除注記等）|
|`*.cfg`|—|実行設定ファイル（INI形式）|


---

## 2. 構造体体系（AMTCore.h）

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
| `terrain_*`（座標系・サイズ） | 地形補正用DEMの起点経緯度・解像度・幅高・アライン幅・確保数（本体DEMと解像度/範囲が異なってよい） |

**地形補正用DEM（`terrain_*` / `pTerrainElevation` / `pTerrainMask`）の本体DEMインライン同居は、DEM系入力を `DataDEM` に集約する意図的設計であり、サブ構造体化は行わない（レビュー対象外、§0-6 参照）。**

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

近傍点探索パラメータ（`pixctx`）・気温減率テーブル・補間設定・TPS係数を保持するほか、以下の主要バッファを含む。`pixctx` の既定値は `InitAnal` で、最小観測点数16、初期半径2000m、最大半径250000m、拡大閾値30000m、拡大ステップ2000m/10000m、拡大時最小観測点数16に設定する。

> **【Phase6＋保守リファクタ分割／VER.0.4.6.475 で pixctx/tps/ring 追加】** `DataAnal` は機能別サブ構造体に分割済み。以下の各バッファ・状態は対応するサブ構造体経由でアクセスする（コード上は `AN.area.vAreaStat`、`AN.lapse.nCurHour`、`AN.ring.mShtTemp` 等）:
> | サブ構造体 | 主な内容 |
> |---|---|
> | `pixctx`（PixCtxParam） | 近傍点探索パラメータ（min_station/radius_init/radius_max/extent_range/extent_step1/extent_step2/extent_min_station）。旧 `PixCtx_*`、接頭辞除去（11-2-1） |
> | `interp`（InterpConfig） | 補間手法パラメータ（method/idw_power/rbf_sigma/barnes_kappa*/barnes_gamma/coast_L/tps_lambda） |
> | `tps`（TpsWork） | TPS補間の作業領域（w/a/stLon/stLat/valid/n）。旧 `tps_*`、接頭辞除去（11-2-1） |
> | `ring`（RingBuffer） | 気温int16リング＋境界マスク＋スロット（pShtTempPool/mShtTemp[]/idx_order[]/mBoundPacked[]/**nvThresholdsScaled**[]/nParallelSlices/nShtTempSlots）。旧トップレベル、`nvThresholds`は改名（11-2-1/11-2-2） |
> | `nbr`（NbrTable） | 近傍点テーブル（pNbrN/pNbrIdx/pNbrD2/pNbrW/pNbrAltW… ＋ nPreMethod 等） |
> | `stats`（TempStats） | 気温統計（時別/日別/月別/全期間）＋気温ヒストグラム各バッファ |
> | `prob`（ProbMap） | 確率マップ（pCountsA/pCountsM/pCountsA_hi/pCountsM_hi/pValid*） |
> | `area`（AreaStats） | エリア別統計（vAreaStat/pAreaMask/nAreaLandPx/pHourlyBuf/nHourlyCapacity 等） |
> | `multi`（MultiPool） | 複数時刻並列バッチ状態（nBatchNext/batchSlots/nBatchN ＋状態遷移メソッド） |
> | `lapse`（LapseState） | lapse推定・時刻帯カーソル状態（nCurHour/nCurTimeBand/fLapseEstMonth/fLapseEstHour 等） |
>
> トップレベル据え置きは `pAnalTemp`/`nLandCount`/`pLandToFull`/`pLandIdxMap`/`nFgFlow`。
> 以下の表・コードは論理的内容を示すもので、実アクセスは上記サブ構造体経由（例: `pNbrN`→`AN.nbr.pNbrN`、`mShtTemp`→`AN.ring.mShtTemp`、`idx_order`→`AN.ring.idx_order`、`nvThresholds`→`AN.ring.nvThresholdsScaled`、`vAreaStat`→`AN.area.vAreaStat`、`pCountsM`→`AN.prob.pCountsM`、`nCurHour`→`AN.lapse.nCurHour`）。

**近傍点テーブル（SoA構造、月1回BuildNbrTable()で構築、`AN.nbr`）**

インデックス規則: `base = landIdx * NBR_MAX`（NBR_MAX=16）

| フィールド | 内容 | 型 |
|---|---|---|
| `pNbrN` | 有効近傍点数 [nLandCount] | uint8_t |
| `pNbrIdx` | 観測点インデックス [nLandCount × 16] | uint16_t |
| `pNbrD2` | 距離² [m²] [nLandCount × 16] | float |
| `pLandLapse` | 全補間手法共通の陸地ピクセル別lapse [°C/m] [nLandCount] | float |
| `pNbrAltDiff` / `pNbrAltDiffF` | 標高差 (int16 or float) | — |
| `pNbrCoastPx` / `pNbrCoastPxF` | ピクセル海岸距離 (int16 or float) | — |
| `pNbrCoastStn` / `pNbrCoastStnF` | 観測点海岸距離 (int16 or float) | — |
| `pNbrW` / `pNbrWF` | 重み (uint16 or float) | — |
| `pNbrAltW` / `pNbrAltWF` | lapse_rate×altDiff + 地形補正項 (int16 or float) | — |
| `pNbrAltWCache` / `pNbrAltWFCache` | 時刻帯別の高度/地形補正キャッシュ（`ALTW_CACHE=1`時） | — |

`exec.nNbrFloatMode` でint/float版を切り替え（0=全整数, 1=SRCのみfloat, 2=Wのみfloat, 3=全float）。
**この int版/float版 二重フィールド維持は省メモリ運用を前提とした意図的設計であり、保守負担はレビュー対象外（§0-5 参照）。**
`pCurNbrAltW` / `pCurNbrAltWF` は時刻帯別キャッシュ参照用で、キャッシュ未使用時の `nullptr` は正常状態である。参照側は本体面 `pNbrAltW` / `pNbrAltWF` へfallbackする。

**気温補間結果（リングバッファ）**　※ `pAnalTemp` 以外は `ring`（RingBuffer）配下（11-2-1）

```c
float*    pAnalTemp;                  // [トップレベル] Barnes2Pass/TPS用float一時バッファ（通常Barnesでは未確保）
int16_t*  ring.mShtTemp[SLOTS];       // int16_t リングバッファ（HSS_GRD_FAC_TEMP倍）
int       ring.idx_order[SLOTS];      // idx_order[0]=書込先（最新）、idx_order[1]=1ステップ前(t0)
int       ring.nParallelSlices;       // 補間並列時刻数
int       ring.nShtTempSlots;         // スロット数（= nParallelSlices + 1）
int16_t*  ring.pShtTempPool;          // mShtTemp全スロット連続確保バッファ
```

`pAnalTemp == nullptr` は通常Barnes等では正常で、`ring.mShtTemp` へ直接格納する。Barnes2Pass/TPSでは `ChgInterpMethod()` 経由で必要時に確保する。

**等温線境界マスク（ビットパック）**　※ `ring`（RingBuffer）配下（11-2-1）。閾値リストは `nvThresholdsScaled` に改名（11-2-2）

```c
uint16_t* ring.mBoundPacked[SLOTS];      // bit c = 閾値c の境界マスク（HSS_GRD_TYPE_B16として出力）
int16_t   ring.nvThresholdsScaled[16];   // 閾値リスト（設定値℃にHSS_GRD_FAC_TEMPを乗じた値、旧 nvThresholds）
// 閾値数 nCntThresholds は ConfigCore.threshold 側（引数 CFG 経由で参照）
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

確率マップ系ポインタは設定依存で未確保が正常なものを含む。`nAnalProb=0` では `pCountsM` 系は未確保、`PROB_ANNUAL=0` では `pCountsA` 系は未確保、`nAnalProb=1/2` では `_hi` 系は未確保となる。出力・集計関数は `nullptr` を「対象出力なし」として早期returnする。

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
解析共通設定。`InitConfigCore()` でデフォルト値セット後、`LoadConfigCore()` でINIから上書き。主要セクション：[MAIN] [PATH] [RUN] [EXT_ANAL] [INTERP] [PARAM] [OUTPUT] [OUTPUT_VECTOR] [OUTPUT_GEOTIFF] [OUTPUT_PNG] [PNG_OPT] [LEGEND_*] [OUTPUT_GRD] [TEMP_STAT] [TEMP_HIST] [THRESHOLD_ANAL] [THRESHOLD] [AREA] [POI_CENTROID] [LAPSE] [AREA_CONFIG]。`ConfigCore` 本体はカテゴリ別サブ構造体を保持し、実参照は `CFG.path.szDemPath`、`CFG.exec.nParallelSlices`、`CFG.interp.fCoastL`、`CFG.threshold.nCntThresholds` のようにカテゴリ経由とする。旧フラット参照は残していない。

#### `ConfigRun`（実行動作定義）
旧 `RunConfig` に相当。開始日時・ステップ数・年範囲・対象月リスト等を保持する。現行の `nRunMode` は `RunMode` enum（`RUN_MODE_CONTINUOUS` ～ `RUN_MODE_ST_NBR_EVAL`）で扱う。`RUN_MODE_COMPARE_DEPRECATED` は欠番として保持し、`MODE=COMPARE` 指定時はロードエラーにする。

#### `AppConfig`（`ConfigCore C` + `ConfigRun R` の薄いラッパ）
`InitAppConfig()`＝`InitConfigCore`+`InitConfigRun`、`LoadConfig()`＝`LoadConfigCore`+`LoadConfigRun`+type検証+`ApplyRunModeFixedConfig()`（RunMode別固定値の単一窓口）。**引数規約**: Run* エントリ関数と main 直下の設定処理は `AppConfig`、補間/出力などのリーフ計算・出力関数は `ConfigCore`（`CFG.C`）のみを受け取る（3-3 で明文化）。

#### `RunChunkConfig`（月内チャンク実行範囲）
`RunOneMonth` に渡す月内実行範囲。`startDT` と `nSteps` のみを持ち、雨量相関設定とは分離する。

#### `ConfigRain`（雨量相関解析用設定）
雨量相関評価用の収束判定閾値・降水判定閾値のみを保持する。月内チャンク範囲は `RunChunkConfig` が保持する。

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

### `CHsTextRenderer`
STB TrueType（`stb_truetype.h`）でフォントを読み込み、PNG画像への文字描画（日時インポーズ・凡例ラベル）を行うクラス。HsCommon.h 内のインライン実装。
- `Init(fontPath, h)` でグリフキャッシュを準備、`cacheGlyph(c)` で文字単位キャッシュ
- `renderTextRow()` がテキスト1行を画像行バッファへアルファブレンディング描画。`COutputPNG::overlayRow`（Output.cpp）から行単位で呼ばれる（§12-13）

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
├─ InitDEM(CFG.path.szDemPath)         # GeoTIFF DEM読み込み
└─ CalcCoastDist(DEM)             # 海岸距離計算（Meijster法）
    │
    └─ 設定ファイルループ（argv[1..N]）
        │
        ├─ DEMパスが変わった場合: FreeDEM → InitDEM → CalcCoastDist
        ├─ InitAnal(DEM, ANAL, CFG)
        ├─ InitFlow(DEM, FL)
        │
        └─ switch(CFG.nRunMode)
            ├─ CONTINUOUS(0):     RunContinuous()    # 通常連続実行
            ├─ MULTI_YEAR(1):     RunMultiYear()     # 同月複数年
            ├─ PARAM_SCAN(2):     RunParamScan()     # 閾値パラメータスキャン
            │  （3=旧COMPARE は VER.0.4.6.413 で廃止。EVAL_INTERP/OPTIM_PARAMS に統合）
            ├─ VALUE_SURVEY(4):   RunValueSurvey()   # 値調査
            ├─ EVAL_INTERP(5):    RunEvalInterp()    # 補間法精度評価（LOOCV、本番補間ベース）
            ├─ OPTIM_PARAMS(6):   RunOptimParams()   # 補間パラメータ最適化スキャン（coast_L含む）
            ├─ AREA_CHAR(7):      RunAreaChar()      # 解析エリア特性評価
            └─ ST_NBR_EVAL(8):    RunStNbrEval()     # 観測点近傍・補正項評価用出力
```

**複数cfgファイルの連続実行**: 各cfgは読込前に `InitAppConfig` で初期化し、未記載キーは既定値へ戻す。前cfgの設定値は継承しない。DEMが共通であればDEM本体は再初期化なしで続けて実行し、DEMパスが変わった場合のみ再ロード。

---

## 5. 主要関数群の詳細

### 5-1. 初期化・解放系

#### `InitDEM()` / `FreeDEM()`
GDALでGeoTIFFを読み込み、`pElevation`・`pSeaMask`・`pCoastDist` を確保。`nAlignedWidth` はAVX-512の64バイト境界（16float）に合わせて切り上げ。

#### `InitTerrainDEM()`
`AREA_CONFIG.TERRAIN_CORR>0` の場合に、`PATH.TERRAIN_DEM_PATH` で指定した地形特徴量抽出用DEMを読み込む。読み込み失敗、メモリ確保失敗、RasterIO失敗時は警告を出して `TERRAIN_CORR=0` 扱いへ落とす。`TERRAIN_CORR=1` は値補正、`TERRAIN_CORR=2` はBarnes重み補正で同じ地形特徴量を参照する。

#### `InitAnal()` / `FreeAnal()`
`DataAnal` の主要ポインタをAVX-512アライメント確保（`_mm_malloc`）。`mShtTemp[nShtTempSlots]` はリングバッファとして動作。`pAnalTemp` は通常Barnesでは確保せず、Barnes2Pass/TPS時のみ確保する。CFG.anal.nAnalAreaStat・nAnalAreaHist・nAnalProb等の設定に応じてオプションバッファ（pCountsM系・pTempHistBuf系・NetCDF系等）の確保を制御。`InitAnal()` は `const ConfigCore&` を参照するだけで設定値を書き換えない。POIのピクセル座標は `ConfigCore::poi.poiDefs` ではなく `DataAnal::area.poiDefs` 側へ解決して保持する。

#### `InitFlow()` / `FreeFlow()`
`DataFlow` の5配列（pFlowU/V/Magnitude/Diverg/Vorticity）を確保。

#### 補間手法別バッファの初期化（`ChgInterpMethod()` 経由、VER.0.4.6.172〜173）
補間手法に依存するバッファ（TPS: `tps_*`、Barnes2Pass/TPS共通: `pAnalTemp`）は `ChgInterpMethod()` 内の static ヘルパー（`InitInterpTPSWork`/`EnsureAnalTempBuffer`）で**選択手法に必要な分のみ**確保する（TPSを選ばない限り `tps_*` を確保しない＝メモリ節約）。旧 `InitTPS()`（public）は `InitInterpTPSWork`（Interpolation.cpp 内 static）に統合され、評価モード・通常運用とも `ChgInterpMethod()` 経由に一本化された（呼び出し漏れ防止）。`pLandLapse` は手法共通キャッシュとして `InitAnal()` で確保し、通常運用では `SetupForMonth()` 内（月初ロード直後）で `BuildPixelContext()` と `ChgInterpMethod()` を呼ぶことで初期化を保証する。

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
近傍点テーブルから事前計算済みの重み（`pNbrW`/`pNbrAltW`）を更新。補間手法変更時・月替わり時に実行。`BuildAltWeightTable()` は現在時刻条件で `pLandLapse` を更新した上で `pNbrAltW` を生成する。`TERRAIN_CORR=1` の地形値補正は `pNbrAltW` に加算し、`TERRAIN_CORR=2` の地形類似度重みはBarnes1Passの `pNbrW` に乗じる。

---

### 5-4. 補間手法

`InterpConfig.method` で切り替え可能。

| 定数 | 値 | 手法 |
|---|---|---|
| `INTERP_IDW` | 0 | 逆距離加重（p=2固定） |
| `INTERP_RBF_GAUSS` | 1 | ガウシアンRBF |
| `INTERP_SHEPARD` | 2 | Shepard修正IDW（べき乗可変） |
| `INTERP_BARNES` | 3 | Barnes解析（推奨） |
| `INTERP_BARNES2PASS` | 4 | Barnes 2Passスキャン（観測点密度不足で逆効果の評価） |
| `INTERP_TPS` | 5 | Thin Plate Spline（計算コスト大） |
| `INTERP_KRIGING_S` | 6 | 簡易クリギング（球形バリオグラム由来の共分散重み。Ordinary Krigingの行列解法ではない） |
| `INTERP_KRIGING_OK` | 7 | Ordinary Kriging近似（近傍8点固定、球形共分散、重み事前計算） |

**海岸距離補正**（`GetCoastFactor()`）: ピクセルと観測点の海岸距離差に基づいてIDW重みを減衰。`coast_L=0` で補正なし。

**気温減率補正**（`GetLapseRate()` / `pLandLapse`）: 標準運用では `lapse_rate_table[緯度帯5][季節4][時刻帯]` から減率を取得し、標高差による気温差を補正する。`LAPSE.EST_MODE=1/2` 指定時は月内観測値から推定した減率を使用する。陸地ピクセル別の現在値は `pLandLapse` に保持し、通常補間、Barnes2Pass、TPSのピクセル側高度補正で共通参照する。標準運用は `EST_MODE=0` とする。

**地形特徴量補正**（`AREA_CONFIG.TERRAIN_CORR`）: `TERRAIN_CORR=1` は地形特徴量差に基づく値補正を `pNbrAltW` に加算する。`TERRAIN_CORR=2` はBarnes1Passの距離重みに地形類似度重みを乗じる。`TERRAIN_CORR=2` の強度は `TERRAIN_CORR_L` で指定し、小さいほど地形特徴差に敏感になる。現行評価では `TERRAIN_CORR=1` がLOOCV最良だが、観測点が存在しない地形・標高レンジでは外挿リスクがあるため、長期統計では `TERRAIN_CORR=0/1` の感度分析を併用する。

旧資料にあった `INTERP_BARNES_H` / Barnes高精度版は現行コードでは定数として存在しない。Barnes系の通常運用は `INTERP_BARNES`、評価・特殊用途として `INTERP_BARNES2PASS` を使用する。

**補間手法と複数時刻並列（PARALLEL_SLICES）の関係**（VER.0.4.6.173）:
- `INTERP_IDW`/`INTERP_RBF_GAUSS`/`INTERP_SHEPARD`/`INTERP_BARNES`/`INTERP_KRIGING_S`/`INTERP_KRIGING_OK` は逐次パス・並列パス（`TemperatureInterpolationMulti`）の両方に対応する。
- `INTERP_BARNES2PASS`/`INTERP_TPS` は**逐次パス専用**。並列パスは Barnes系1Passのカーネル（`InterpolationMainMulti`）に特化しているため、これらの手法は並列パスでは正しく補間できない。
- そのため `LoadConfigCore` で **Barnes2Pass/TPS かつ PARALLEL_SLICES>1 の場合は警告を出力して `nParallelSlices=1`（逐次）へ自動矯正し、動作を継続**する。Barnes2Pass/TPS は元々低速（特にTPSは数百ms/ステップ）で並列の費用対効果も低く、これらを使うのは評価・研究目的で速度より正しさが優先されるための設計判断（方針: 逐次フォールバック）。
- `EST_MODE=2` は時刻別 `pLandLapse` を単一面で管理するため、補間手法に関係なく `PARALLEL_SLICES>1` を `1` へ自動矯正する。スロット別 `pLandLapse` は実装複雑度とメモリ増加に対して効果が限定的なため採用していない。
- 各手法に必要なバッファ（TPS: `tps_*`、Barnes2Pass/TPS共通: `pAnalTemp`）は、通常運用でも `SetupForMonth()`（月初ロード直後）が `ChgInterpMethod()` を経由して確保する（評価モードと同経路に統一、VER.0.4.6.172）。

---

### 5-5. 補間実行系

#### `TemperatureInterpolation()`
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

#### LOOCV評価方式（`RunEvalInterp`/`RunOptimParams`、VER.0.4.6.413で本番補間ベースに一本化）
**LOOCV（Leave-One-Out Cross Validation）**: 各観測点を1つ除外（観測値を欠測マスク）して**本番補間（`TemperatureInterpolation`）を再実行**し、当該ピクセルの推定値と実測値の誤差を評価する。地形補正・KrigingS/KrigingOK・海岸補正など本番の重み計算がそのまま反映される（旧 static `EvalInterp` の手計算LOOCVは廃止、§1-9）。
- 指標: RMSE / MAE / MaxErr / Bias / NSE / r（共有 `ZoneStat`/`ZoneStatAccum`/`ZoneStatFinal` で集計）
- ゾーン別集計（`ZoneClassifier`）: 全体 + 緯度3帯 + 経度3帯 + 標高3帯（DEM範囲を3等分、標高は200m/600m境界）
- 空間的滑らかさ（ラプラシアン）は `CalcLaplacianSmoothness` で算出し `eval_summary` の `LapMean`/`Lap95p` に出力

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
        1. RunOneMonthPrepareStep（RotateTimeSlice、前日分の日別統計出力、時刻帯設定）
        2. RunOneMonthInterpolateStep（逐次/複数時刻並列補間、境界マスク）
        3. RunOneMonthOutputStep
           - RunOneMonthOutputCommonStep（観測点ログ、時別統計、PNG、エリア統計、境界面出力）
           - RunOneMonthOutputFlowStep（フロー計算、フロー系GeoTIFF/GPKG/GRD、DivergRain評価）
           - RunOneMonthOutputNoFlowStep（非フロー時の気温GeoTIFF/GRD）
        4. dtPrev=t0 / DateTimeAdvance / 処理時間ログ更新
```
共通ヘルパー: `SetupForMonth`（月初ロード+補間準備）、`InitOutputAfterAreaMask`（DataAnal依存の後段出力初期化）、`FlushMonth`（月末/期間末の統計・確率出力。内部は `FlushMonthAreaThreshold` / `FlushMonthTemperature` / `FlushMonthProbability` に分割）、`ResetAnalBuffers`（idx_order/エリアバッファ初期化）。

#### `RunContinuous()`（期間動作）
CFGの開始日時から `nRunCount` ステップを**月チャンクに分割**し `RunOneMonth` を反復。`ResetAnalBuffers` は期間先頭で1回（月跨ぎで idx_order/フロー連続性維持）、`nGlobalStepBase`/`dtPrev` を月跨ぎで引継ぎ。各チャンク末で `FlushMonth(bProbImage=false)`、期間末で `FlushAllTempStat`（TOTAL）。MONTH_OUT は `ApplyRunModeFixedConfig` で無効化（期間統計は TOTAL_OUT）。

#### `RunMultiYear()`（月動作）
`nYearStart`〜`nYearEnd` の同月（`nTargetMonthList[]`）をループ。各月で `ResetAnalBuffers`（月別独立）→ `RunOneMonth(nGlobalStepBase=0)` → `FlushMonth(bProbImage=true)`。年別統計を1行ずつ出力し、終了後に全期間確率マップ＋全期間統計を出力。

#### `RunParamScan()`
**2段階方式**:
- **Stage1**: 全ステップを実行し、各ステップの発散値（観測点位置のみ）と降水量をスナップショット保存
- **Stage2**: 発散閾値 × 降水閾値 × 時間ラグの全組み合わせでスナップショットを参照し高速再評価

#### `RunEvalInterp()`（MODE=EVAL_INTERP）
複数の補間手法を指定日時に対してLOOCV評価し、精度比較結果をTSVとして出力する。
- `PrepareEvalModeData`（LoadAMeDAS / BuildNbrTable 等）の後、手法を切り替えながら、各手法で「通常補間」＋「観測点マスク＋本番補間再実行のLOOCV」を実行する（旧 static `EvalInterp` 呼び出しは廃止、本番補間ベースに統一）
- `RUN.EVAL_METHODS` で評価対象・評価順・結果出力順を指定する。例: `3/6/7` は Barnes1Pass → KrigingS → KrigingOK の順。TPSは低速のため通常評価では除外する。
- RMSE / MAE / MaxErr / Bias / NSE / r をゾーン別（全体・緯度3帯・経度3帯・標高3帯）に集計。手法別に空間的滑らかさ `LapMean`/`Lap95p` も算出（`eval_summary` 末尾2列、「全体」行のみ値）
- 出力: `*_eval_loocv_YYYYMMDD_HH.txt`（観測点別LOOCV結果）、`*_eval_summary_YYYYMMDD_HH.txt`（手法×ゾーン×指標＋平滑度のTSV）。`geotiff.nOutGeoTIFF_Temp`(TEMP) 有効時は手法別気温 GeoTIFF `*_eval_temp_{手法名}_YYYYMMDD_HH.tif` も出力

#### `RunOptimParams()`（MODE=OPTIM_PARAMS）
Barnes / Barnes2Pass / RBF-Gauss / Shepard / IDW / coast_L の補間パラメータをグリッドスキャンし、本番補間ベースLOOCVの最適値を導出する。LOOCVは共有 `ZoneStat` でゾーン別集計。
- **スキャン対象パラメータ（Phase順）**:
  - Phase1 Barnes/Barnes2Pass: `kappa_fac` ∈ {0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 5.052, 8.0, 12.0, 15.0, 20.0}
  - Phase2 Barnes2Pass: `gamma` ∈ {0.1, 0.2, 0.3, 0.5}（Phase1最適 kappa_fac 固定）
  - Phase3 RBF-Gauss: `rbf_sigma` ∈ {10000, 20000, 30000, 50000, 80000, 100000, 150000}（m）
  - Phase4 Shepard: `idw_power` ∈ {1.0, 1.5, 2.0, 2.5, 3.0}
  - Phase5 IDW: 固定 power=2.0 の単発評価（参考値、CalcWeight内固定のため）
  - **Phase6 Barnes: `coast_L` ∈ {0, 20000, 30000, 50000, 80000, 120000}（m）（VER.0.4.6.413、旧 Eval_Coast_L 統合。Phase1最適 kappa_fac 下で評価し推奨CFGの前提を一致させる）**
- 各組み合わせでLOOCV RMSEを計算し、手法ごとにベストパラメータを記録
- **境界警告**: 最優パラメータがスキャン範囲の上下端に一致した場合、コンソールとファイルに警告を出力（範囲拡大の要否を通知）
- 出力: `*_eval_optimize_YYYYMMDD_HH.txt`（全パラメータ×ゾーンの詳細）、`*_optim_result_YYYYMMDD_HH.txt`（推奨CFG: `BARNES_KAPPA_FAC`/`BARNES_GAMMA`/`RBF_SIGMA`/`IDW_POWER`/`COAST_L`＋境界警告）

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
観測点近傍、地形特徴量、補正項確認用の補助モード。`SetStNbrEvalFixedConfig()` によりPNG/GeoTIFF/GRD/確率/フロー等の通常出力を抑止し、指定年・月・日・時刻だけを処理する。結果は `StNbrDump\YYYY\{FILE_KEY_ST}_stnbr_yyyymmdd_hh.tsv` と `{FILE_KEY_ST}_eval_manifest.tsv` に記録され、地形補正や近傍点構成の検証に使う。`FILE_KEY` ではなく `FILE_KEY_ST` をファイル名プレフィックスに使う。`RunOneMonth()` 失敗時は当該日の処理をスキップして次の指定日へ進み、manifestは成功後に追記する。
**HOURS仕様**: 不連続指定（例 `HOURS = 0/23`）の場合も min～max の全時刻（この例では24時刻）を補間計算する（`RunOneMonth` が連続ステップ実行を前提とするため。ダンプ出力は指定時刻のみ）。まばらな時刻指定では中間時刻の計算が無駄になる（レビュー6-10、仕様として明文化）。

---

## 6. 出力ファイル系

### 6-1. GeoTIFF（WriteGeoTIFF.h）
- `nOutGeoTIFF_Temp=1`: 気温解析結果
- `nOutGeoTIFF_Flow=1`: フロー結果（FLOWU, FLOWV, MAGNITUDE, DIVERG, VORTICITY）
- `nOutGeoTIFF_Prob=1`: 気温閾値確率分析結果
- `nOutGeoTIFF_Cont=1`: 気温境界面マスク
- `geotiff.nCfgGeoTIFF_Thinning` で時間間引き

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

`grd.nCfgGRD_Thinning` で出力頻度を制御。

> **【検証上の注意】GRD(.tva) はバイト非決定的**: BLOSC2 の8スレッド圧縮はブロックのスレッド分配が実行ごとに変わるため、**同一入力・同一バイナリでもバイト列が一致しない**（解凍後の値は常に同一）。リグレッション検証では GRD のバイト/SHA256 比較は使えず、**TXT/NC の値・SHA256 で判定**すること（決定的にしたい場合は圧縮スレッド数を1にする選択肢がある）。

> **【出力ファイル名の共通ステム】** 統計TXT/NC・TMD・LOG は `GetRunStemName()` で生成する共通ステムを共有する（VER.0.4.6.179、C3d）。`{RESULT_PATH}{FILE_KEY_ST}` に続けて、MULTI_YEAR=`_YYYY-YYYY_MM`／その他（CONTINUOUS等）=`_YYYYMMDD_HH` を付す。GeoTIFF/GRD/PNG 等の地理データは `FILE_KEY`＋日時ベース。

### 6-3. PNG出力（Output.cpp）
- `nOutPNG_Temp=1`: 気温分布画像 + 等温線別PNG（`nOutPNG_Cont=1` 時）
- `nOutPNG_Temp=2`: マージ出力（気温分布+等温線合成）
- `nOutPNG_Prob=1`: 閾値確率分析結果画像

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
TERRAIN_DEM_PATH = U:\DEM\terrain_50m.tif    # 地形特徴量抽出用DEM（AREA_CONFIG.TERRAIN_CORR>0時）
RESULT_PATH  = U:\PRG3\AMT_RES\20260425A\    # 結果出力先ディレクトリ
FILE_KEY     = japan                         # 画像、GeoTIFF、GRD等の地理データ出力プレフィックス
FILE_KEY_ST  = RES_JPN                       # 統計、評価、ログ、ST_NBR_DUMP出力プレフィックス
AMEDAS_PATH  = U:\AMeDAS\                    # AMeDASデータディレクトリ（STN/IDX/AMDファイルが存在するパス）

[RUN]
# MODE : CONTINUOUS / MULTI_YEAR / PARAM_SCAN / VALUE_SURVEY / EVAL_INTERP / OPTIM_PARAMS / AREA_CHAR / ST_NBR_EVAL
#        （COMPARE は VER.0.4.6.413 で廃止。指定するとロードエラー）
MODE         = MULTI_YEAR
CONT_START   = 2000/2/15/0                   # CONTINUOUSモードの実行開始日時（年/月/日/時）
CONT_DUR     = 2/1                           # CONTINUOUSモードの実行期間（時/日）、日設定を優先
EVAL_METHODS = 3/6/7                         # EVAL_INTERP/OPTIM_PARAMS評価対象・評価順・結果出力順
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
METHOD          = 3                          # 気温空間補間の補間方法 (0=IDW, 1=RBF_GAUSS, 2=SHEPARD, 3=BARNES, 4=BARNES2PASS, 5=TPS, 6=KRIGING_S, 7=KRIGING_OK)
COAST_L         = 20000.0                   # 海岸距離補正スケール [m]（0で補正なし）
IDW_POWER       = 1.5                        # IDW/Shepard べき乗値（デフォルト 1.5）
RBF_SIGMA       = 80000.0                    # RBF-Gauss バンド幅 [m]（デフォルト 80000）
BARNES_KAPPA_FAC = 2.0                       # Barnes: kappa = d_nn² × この値（デフォルト 2.0、文献値 5.052）
BARNES_GAMMA    = 0.3                        # Barnes 2Pass スキャン間収縮率（デフォルト 0.3）

[AREA_CONFIG]
TERRAIN_CORR           = 0                    # 地形補正（0:無効、1:値補正、2:地形類似度重み）
TERRAIN_CORR_CLIP      = 0.0                  # 地形補正項クリップ幅[℃]（0以下で無効）
TERRAIN_CORR_COAST_MIN = 0.0                  # 適用最小海岸距離[m]（0以下で無効）
TERRAIN_CORR_L         = 200.0                # TERRAIN_CORR=2 の地形類似度重みスケール
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
# ANAL2 は VER.0.4.6.413 で廃止（旧 coast_L 評価フラグ。coast_L 評価は MODE=OPTIM_PARAMS に統合）

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
# EST_MODE: 0=既存TABLE（標準）, 1=月別推定（評価用）, 2=月別×24時刻別推定（研究用）
EST_MODE = 0
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
| TemperatureInterpolation | msec. | 10.49 | 2.10 | 13.05 | 0.77 | 4.59 | 18.57 | 1.11 | 6.65 | 26.19 |
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
| TemperatureInterpolation | msec. | 8.81 | 1.61 | 10.85 | 0.63 | 3.77 | 15.73 | 0.89 | 5.65 | 22.42 |
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
| TemperatureInterpolation改善 | % | -16.0 | -23.3 | -16.9 | -18.2 | -17.9 | -15.3 | -19.8 | -15.0 | -14.4 |
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

> **削除済み（VER.0.4.6.413、§1-9）**: `EvalInterp`(static)（手計算LOOCV、地形補正/Kriging非反映）/ `PrintEvalResult` / `InterpEvalResult` 構造体 / `ScanFactorTPS`（TPS λ スキャン、呼び出し元なし）。LOOCV は本番補間ベースの `RunEvalInterp`/`RunOptimParams` へ一本化。

> **削除済み（柱A A1、VER.0.4.6.175）**: `ValidateStationElevation` / `CheckPixelContext` は呼び出しのないデッドコードのため削除（計82行）。

> **移設（柱A A4、VER.0.4.6.175 → VER.0.4.6.300）**: 以下の汎用出力4関数（`OutputTempGeoTIFF`/`OutputFlowGeoTIFF`/`OutputFlowGPKG`/`OutputAnalyzeBinary`）は **amt.cpp → AnalCore.cpp → Output.cpp** へ移動済み。複数モード共通の出力機能として `Output.cpp` に集約する。`OutputAnalyzeBinary` は引数を `AppConfig`→`ConfigCore` に変更済み。

#### `OutputTempGeoTIFF`（→ Output.cpp）
- **処理概要**: `mShtTemp` の気温面をGeoTIFFとして出力
- **引数**: `const DataDEM& DEM`, `const DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`, `int nShtSlot`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()。フロー有効時は `dtPrev` + `idx_order[1]`、非フロー時は `dtCurDT` + `idx_order[0]` を渡す。

#### `OutputFlowGeoTIFF`
- **処理概要**: FlowU/FlowV/Magnitude/Diverg/Vorticity をGeoTIFFとして出力
- **引数**: `const DataDEM& DEM`, `const DataAnal& AN`, `const DataFlow& FL`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `OutputFlowGPKG`
- **処理概要**: フロー解析結果をベクトルGPKGとして間引き出力
- **引数**: `const DataDEM& DEM`, `const DataAnal& AN`, `const DataFlow& FL`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `OutputAnalyzeBinary`（→ Output.cpp）
- **処理概要**: 解析結果を GRD バイナリ形式で出力（タグ列構築・mode導出は WriteGRD 内で CFG から実施）
- **引数**: `const DataDEM& DEM`, `const DataAnal& AN`, `const DataFlow& FL`, `const ConfigCore& CFG`, `const DateTime& dtWr`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `BuildPixelContext`（→ Interpolation.cpp、柱A A3）
- **処理概要**: 全陸地ピクセルの近傍16観測点情報（PixelContext）を構築。**amt.cpp → Interpolation.cpp** へ移動（3モード共通の補間準備をライブラリ側へ）、引数を `AppConfig`→`ConfigCore` に変更
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& AMD`, `DataAnal& AN`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: SetupForMonth(), RunParamScan()

> **削除済み（VER.0.4.6.413、§1-9）**: `Eval_Coast_L`（旧 ANAL2 経由の coast_L 手計算スキャン）/ `Eval_Coast_L_2Time`（呼び出し元なし）/ `ComparInterpolationMethod`（COMPARE本体）。coast_L 評価は `RunOptimParams` Phase6（本番LOOCV）へ統合。

#### `RunOneMonth`（旧 RunOneBlock、static、VER.0.4.6.178 で月単位純化）
- **処理概要**: 1か月（または期間内の月チャンク）を処理する純粋関数。先頭で `SetupForMonth` を1回呼び、月内 nSteps ステップをループ。月境界判定・月次flush・バッファリセットは内部に持たず呼出側責務。月ローカル `s`（並列バッチ判定）とグローバル `step = nGlobalStepBase + s`（フロー有無/間引き位相/GRD判定）を使い分ける。
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataFlow& FL`, `DataAMeDAS& AMD`, `const AppConfig& CFG`, `const RunChunkConfig& CHUNK_CFG`, `DivergRainAccum* pAccum`, `int nGlobalStepBase`, `DateTime& dtPrev`, `double& init_time_msec`
- **戻り値**: int — 0:成功 / -1:月初ロード失敗
- **呼び出し元**: RunContinuous(), RunMultiYear(), RunStNbrEval()

> **配置注記**: `SetupForMonth`/`FlushMonth`/`ResetAnalBuffers` は複数モード共通の月次ライフサイクル処理として **AnalCore.cpp の非static関数**（§0-4 参照）。本一覧では **§9-2** に記載する。

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

> **削除済み（VER.0.4.6.413、§1-9）**: `RunCompare`（mode3、COMPARE）。手法別LOOCV比較は `RunEvalInterp` へ統合。`MODE=COMPARE` は `LoadConfigRun` がエラーで停止する。

#### `RunEvalInterp`
- **処理概要**: 指定日時の観測点LOOCV（観測点マスク＋本番補間再実行）を手法別に実行し、`eval_loocv`（観測点別）・`eval_summary`（手法×ゾーン×指標）を出力。`eval_summary` には空間的滑らかさ `LapMean`/`Lap95p`（`CalcLaplacianSmoothness`、「全体」行のみ）を含む。`geotiff.nOutGeoTIFF_Temp`(TEMP) 有効時は手法別気温 GeoTIFF（`{FILE_KEY_ST}_eval_temp_{手法名}_*.tif`）を出力（旧COMPARE固有機能を移植、VER.0.4.6.413）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataAMeDAS& DT`, `const AppConfig& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `RunOptimParams`
- **処理概要**: 補間パラメータをグリッドスキャンし LOOCV 最適値を出力。Phase1: Barnes/Barnes2P kappa_fac、Phase2: Barnes2P gamma、Phase3: RBF rbf_sigma、Phase4: Shepard idw_power、Phase5: IDW参考値、**Phase6: Barnes coast_L（VER.0.4.6.413、旧 Eval_Coast_L 統合。Phase1 最適 kappa_fac 下で評価）**。LOOCV はゾーン別共有統計（`ZoneStat`）で集計。`optim_result` に `BARNES_KAPPA_FAC`/`BARNES_GAMMA`/`RBF_SIGMA`/`IDW_POWER`/`COAST_L` と境界警告を出力
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

### 9-2. AnalCore.cpp

> VER.0.4.6.300 以降、出力関数群は `Output.cpp` に分離済み。ここでは設定Core読込、データ読込、解析共通ヘルパ、月次ライフサイクル、DivergRain評価などを `AnalCore.cpp` の主責務とする。

> **配置注記**: `CHsTextRenderer`（Init/cacheGlyph/renderTextRow）は **HsCommon.h 内のインライン実装**（§3 参照）であり AnalCore.cpp には存在しない。`renderTextRow` は `COutputPNG::overlayRow`（Output.cpp）から呼ばれる（旧 SavePNG8wSIMD 呼び出しは §12-13 で廃止）。

#### `CHsImageLegend::get_text_width`
- **処理概要**: TrueType 文字列のピクセル幅を計算（カーニング対応）
- **引数**: `const stbtt_fontinfo* font`, `const char* text`, `float scale`
- **戻り値**: int — 描画幅(px)
- **呼び出し元**: CHsImageLegend::Init

#### `CHsImageLegend::Init`
- **処理概要**: 凡例画像を生成（気温/確率/等温線マップ用、mode 1-5）
- **引数**: `const ConfigCore& CFG`, `int mode`
- **戻り値**: bool — 成功/失敗
- **呼び出し元**: InitOutputLegendsCore()（Output.cpp。LoadConfigCore 末尾の InitOutputLegends 経由）

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
- **処理概要**: INIファイルから ConfigCore を読み込み（未記載キーはデフォルト維持）。VER.0.4.6.500 でセクション群毎の static 7関数（`LoadConfigBase`/`LoadConfigAnal`/`LoadConfigStatOut`/`LoadConfigFileOut`/`LoadConfigThresholdArea`/`LoadConfigLapse`/`LoadConfigAreaCorr`）＋`FinalizeConfigCore`（パレット・凡例・FileKey検証）の順次呼出へ分割（§12-14）。読込定型は `GetIniString`（空初期化保証）＋`SplitValueList`（/,両対応リスト分割）に集約
- **引数**: `CHsIniFile& ReadINI`, `ConfigCore& CFG`（VER.0.4.6.497 で `const char*`→`CHsIniFile&` へ変更。LoadConfig が生成した単一インスタンスを LoadConfigRun と共有し、未参照キー検出を可能にする）
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: LoadConfig()

#### `InitConfigRun`
- **処理概要**: ConfigRun 構造体にデフォルト値を設定
- **引数**: `ConfigRun& CFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `LoadConfigRun`
- **処理概要**: INIファイルから実行動作設定（実行モード・年月・ループ設定等）を読み込み。`MULTI_MONTH` / `DAYS` / `HOURS` / `EVAL_METHODS` は `CHsIniFile::GetIntList` で共通読込
- **引数**: `CHsIniFile& ReadINI`, `ConfigRun& CFG`（VER.0.4.6.497 で `const char*`→`CHsIniFile&` へ変更。LoadConfigCore と同一インスタンス共有）
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: LoadConfig()

#### `HsPathExistCheck`（VER.0.4.6.500 で HsCommon.h へ外出し）
- **処理概要**: パスの存在確認、存在しない場合はディレクトリ生成。旧 AnalCore.cpp の `PathExistCheck` を AMT非依存ユーティリティとして HsCommon.h の inline 関数へ移設・改名
- **引数**: `const char* szPath`
- **戻り値**: int — 0:成功, -1:失敗（空パス・確認失敗・生成失敗）
- **呼び出し元**: main()、OutputContext のパス初期化（Output.cpp）

#### `EvalDivergRain`
- **処理概要**: 観測点位置の発散値と降水量を照合し、TP/FP/FN/TN・相関・降水有無別統計を算出（相関は `CorrAccum`/`CalcPearsonCorr` 共通ヘルパ使用）
- **引数**: `const DataDEM& DEM`, `const DataAMeDAS& AMD`, `const DataFlow& FL`, `DivergRainResult* pResult`, `const ConfigRain& RAIN_CFG`（VER.0.4.6.500 で閾値2個ばら渡しを ConfigRain 受けへ変更）
- **戻り値**: void
- **呼び出し元**: RunOneMonth(), RunParamScan()

#### `PrintDivergRainResult`
- **処理概要**: 発散場・降水量照合の評価結果（HSS/Hit Rate/CSI 等）をコンソール出力
- **引数**: `const DivergRainResult& r`, `const DateTime& dt`
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
- **引数**: `DataAMeDAS& DT`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: TemperatureInterpolation()

#### `SetCalcTempToSlot`
- **処理概要**: 指定日時の気温を複数時刻並列用スロット pfTempsArr[slot] にセット
- **引数**: `DataAMeDAS& DT`, `const DateTime& dt`, `int slot`
- **戻り値**: void
- **呼び出し元**: RunOneMonth(), TemperatureInterpolationMulti()

#### `SetCalcRain`
- **処理概要**: 指定日時の降水量を fRains バッファにセット（メモリ未確保時は確保）
- **引数**: `DataAMeDAS& DT`, `const DateTime& dt`
- **戻り値**: int — 0:成功, -1:メモリ確保失敗
- **呼び出し元**: RunOneMonth(), RunParamScan()

#### `FreeAnal`
- **処理概要**: DataAnal 構造体の全メモリを解放
- **引数**: `DataAnal& AN`
- **戻り値**: void
- **呼び出し元**: InitAnal()（失敗時）, main()（終了時）

#### `InitAnal`
- **処理概要**: DataAnal を初期化（補間パラメータ設定・各種バッファ確保・近傍テーブル構築）。VER.0.4.6.500 で機能ブロック毎の static 8関数（`InitAnalDefaults`/`InitAnalGrid`/`InitAnalNbrTable`/`InitAnalAreaStats`/`InitAnalMultiPool`/`InitAnalProb`/`InitAnalThreshold`/`InitAnalPoi`）＋`PrintInitAnalMemSummary` の順次呼出へ分割（§12-14）。確保定型は `ANAL_ALLOC` マクロに集約し、確保失敗時の `FreeAnal` は本体で一括実施
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const ConfigCore& CFG`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: main()。呼出側で `DataAnal` を0クリアしてから渡す。関数先頭で `FreeAnal()` を行うため、既存確保済みANの再初期化でも旧メモリを失わない。

#### `MatchFieldValue`
- **処理概要**: OGR フィールド値がフィルタ条件（単値/範囲/複数値）と一致するか判定
- **引数**: `const char* szFieldVal`, `const char* szSpec`
- **戻り値**: bool — 一致時 true
- **呼び出し元**: BuildAreaMask()

#### `BuildAreaMask`
- **処理概要**: SHP ポリゴンから GDAL OGR 経由でエリアマスクをグリッドに焼き込み
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: main()。内部は `LoadAreaShpPolygons` / `BuildAreaBBoxList` / `RasterizeAreaMask` / `CountAreaLandPixels` / `PrintAreaMaskSummary` に分割。SHP読込失敗・フィールド不一致・リング0件・陸地ピクセル0件は警告を出す。

#### `CalcAreaStat`
- **処理概要**: 気温閾値エリア別統計（発生数・消滅数等）を時別集計
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const DateTime& dtCurDT`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### `FlushProvMonthToAccum`
- **処理概要**: 月次確率バッファを全期間バッファに加算してリセット（AVX-512）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: FlushMonth()（月末）

#### `InitOutputAfterAreaMask`（static）
- **処理概要**: `BuildAreaMask` 後に必要な DataAnal 依存出力を初期化。`InitHourlyStatNC` と `OutputAnalysisSummary` を1回だけ実行
- **引数**: `const DataDEM& DEM`, `const DataAMeDAS& AMD`, `DataAnal& AN`, `const ConfigCore& CFG`, `const ConfigRun& RCFG`（VER.0.4.6.500 で AppConfig 一括受けを層別に分解）
- **戻り値**: void
- **呼び出し元**: SetupForMonth()

#### `SetupForMonth`（VER.0.4.6.177）
- **処理概要**: 月初の AMeDAS ロード＋補間準備（EstimateLapseRateForMonth/CalcInterpParams/BuildPixelContext/ChgInterpMethod）＋月初初期化（観測点ログ/エリアマスク/時別統計NC/解析概要）を一括実行
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataAMeDAS& AMD`, `const ConfigCore& CFG`, `const ConfigRun& RCFG`, `const DateTime& dt`, `double& init_time_msec`（VER.0.4.6.500 で AppConfig 一括受けを層別に分解）
- **戻り値**: int — 0:成功 / -1:LoadAMeDAS失敗
- **呼び出し元**: RunOneMonth()

#### `ResetAnalBuffers`（VER.0.4.6.178）
- **処理概要**: idx_order とエリア統計バッファ（vAreaStat/nHourlySteps/pHourlyBuf）を初期化。MULTI_YEAR は月ごと、CONTINUOUS は期間先頭で1回呼ぶ（AN.multi のリセットは RunOneMonth 先頭で別途実施）
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: RunContinuous(), RunMultiYear(), RunStNbrEval()

#### `FlushMonth`（VER.0.4.6.177）
- **処理概要**: 月末/期間末の統計・確率出力の順序管理を行う。内部で `FlushMonthAreaThreshold`（閾値統計月別/時別）、`FlushMonthTemperature`（日別・月別の気温統計/ヒストグラム）、`FlushMonthProbability`（確率マップ月別統計＋テキスト、`bProbImage` 時は確率マップ画像、全期間累積）を順に呼び、最後に TMD をflushする。
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dtMonth`, `const DateTime& dtPrev`, `bool bProbImage`
- **戻り値**: void
- **呼び出し元**: RunContinuous()（bProbImage=false）, RunMultiYear()（bProbImage=true）

#### `PrepareEvalModeData`
- **処理概要**: 評価系RunMode共通の AMeDASロード＋idx_order初期化＋近傍表構築（Run* プロローグ共通化、§12-10）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `DataAMeDAS& AMD`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: int — 0:成功 / -1:LoadAMeDAS失敗
- **呼び出し元**: RunEvalInterp(), RunOptimParams(), RunValueSurvey(), RunAreaChar()

### 9-3. Output.cpp

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
- **処理概要**: 実行条件、DEM、補間、出力ファイル一覧、エリア情報などの概要TSVを出力。`MODE` 行は VER.0.4.6.500 から「番号+名前」のタブ区切り併記（例 `MODE	0	CONTINUOUS`、`GetRunModeName` で逆引き。回帰ハーネスは番号部のみ比較）
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
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`, `int step`（THINNING判定は関数内、Phase7-B）
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

#### `InitOutputProbText`
- **処理概要**: 確率マップ出力用 TSV ファイル（4種）をオープンしてヘッダを書き込み
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: InitAllTextOut()

#### `CalcProbMonthStat`
- **処理概要**: 月次/全期間の確率統計（平均/最大/ヒストグラム）を計算
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const ConfigCore& CFG`, `int mode` — 月次(0)/全期間(1), `int dir` — lo(0)/hi(1)
- **戻り値**: void
- **呼び出し元**: FlushMonth(), RunMultiYear()

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

> **削除済み（VER.0.4.6.491〜494、§12-13）**: `SavePNG8`/`SavePNG8wSIMD`/`SavePNG8woSIMD`/`SavePNG16`/`SavePNG16wSIMD`/`SavePNG16woSIMD` の自由関数6本は廃止。PNG出力本体は `CHsOutputPNG`（HsOutputPNG.h/.cpp、**§9-6**）へ移設し、Output.cpp 側はファイルローカル派生 `COutputPNG`（文字・凡例の `overlayRow`／マージPNG行生成の `buildRow8`＋`SetMergedParams`）で利用する。

#### `InitProcTimeMeas`
- **処理概要**: 処理時間計測用 TMD ファイルを準備
- **引数**: `const ConfigCore& CFG`, `const ConfigRun& RCFG`
- **戻り値**: void
- **呼び出し元**: main()

#### `InitOutputStationTemperature`
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
- **呼び出し元**: InitOutputAfterAreaMask()

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
- **備考**: cfg確定直後に初期化できる出力を担当。DataAnal依存のNetCDF/解析概要は `InitOutputAfterAreaMask` が担当する

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
- **引数**: `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`, `int step`（THINNING判定は関数内、Phase7-B）
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

### 9-4. Interpolation.cpp

#### `CalcInterpParams`
- **処理概要**: 観測点配置から補間パラメータ（Barnes kappa, RBF sigma 等）を自動計算
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `InterpConfig& ITR_CFG`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: SetupForMonth(), PrepareEvalModeData(), RunParamScan(), RunValueSurvey(), RunAreaChar()

#### `BuildNbrTable`
- **処理概要**: 全陸地ピクセルの近傍16観測点を探索し SoA 配列（pNbr*）を構築（AVX-512最適化）。pLandLapse 初期値も格納（GetLapseRate に CFG が必要）
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `DataAnal& AN`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: BuildPixelContext(), PrepareEvalModeData()

#### `BuildWeightTable`
- **処理概要**: 補間手法別重み・高度補正項を事前計算（pNbrW, pNbrAltW）。pNbrAltW生成時にpLandLapseを更新・参照
- **引数**: `const DataDEM& DEM`, `const DataAMeDAS& DT`, `DataAnal& AN`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: TemperatureInterpolation(), BuildPixelContext(), TemperatureInterpolationMulti()

#### `ChgInterpMethod`
- **処理概要**: 補間手法を切り替え（TPS用tps_*、Barnes2Pass/TPS用pAnalTempなど手法固有バッファを必要時に確保）
- **引数**: `const DataDEM& DEM`, `const DataAMeDAS& DT`, `DataAnal& AN`, `InterpMethod nNewMode`
- **戻り値**: void
- **呼び出し元**: SetupForMonth(), RunEvalInterp(), RunOptimParams()

#### `InterpolationMain`
- **処理概要**: IDW/Shepard/RBF/Barnes 共通補間ループ（AVX-512 SIMD 化）
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `DataAnal& AN`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: TemperatureInterpolation()

#### `InterpolationMainMulti`
- **処理概要**: N 時刻分の補間を共有ループで一括処理（キャッシュ効率向上）
- **引数**: `const DataDEM& DEM`, `float* const* ppFTemps`, `int16_t** ppDst`, `int nBatch`, `DataAnal& AN`, `const ConfigCore& CFG`
- **戻り値**: void
- **呼び出し元**: TemperatureInterpolationMulti()

> `CollectAreaCntFromShtTemp` は VER.0.4.6.413 で削除（呼び出し0・デッドコード。areaCnt収集は InterpolationMain / CalcBoundaryMaskAndCollect[Multi] が担う）。

#### `TemperatureInterpolationMulti`
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
- **処理概要**: int16 気温場から各閾値の境界マスクを生成する内部ヘルパー（per-pixel は共通 `CalcBoundaryBits()` へ集約、AVX-512）
- **引数**: `const DataDEM& DEM`, `const DataAnal& AN`, `const ConfigCore& CFG`, `const int16_t* pSht`, `uint16_t* pPacked`
- **戻り値**: void
- **呼び出し元**: CalcBoundaryMaskOnly()

#### `CalcBoundaryMaskAndCollect`
- **処理概要**: 境界マスク生成とエリア統計集計を同時実行（境界bitは `CalcBoundaryBits()` 共通化）。`bGenMask=false` で**集計のみ・マスク非生成**（B2P/TPS逐次が利用、マスクは外側 `CalcBoundaryMaskOnly` で一律生成、VER.0.4.6.488）
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const ConfigCore& CFG`, `bool bGenMask=true`
- **戻り値**: void
- **呼び出し元**: TemperatureInterpolation(), RunOneMonth()

#### `CalcBoundaryMaskAndCollectMulti`
- **処理概要**: N 時刻の境界マスク生成・エリア統計集計を一括処理
- **引数**: `const DataDEM& DEM`, `DataAnal& AN`, `const ConfigCore& CFG`, `const int* pBatchSlots`, `int nBatch`
- **戻り値**: void
- **呼び出し元**: RunOneMonth()

#### ~~`TemperatureSetInt16`~~（VER.0.4.6.500 で削除）
- 呼出0のデッドコードのため削除（レビューP3-5）。int16化は各補間関数が `TempToI16` で直接実施する。TMD処理時間チャンネル3は欠番（再利用しない）

#### `CalcTPSCoeff`
- **処理概要**: TPS 補間の係数をガウス消去法で計算
- **引数**: `DataAMeDAS& DT`, `DataAnal& AN`, `const float* pTempOverride` — nullptr 時は DT.fTemps 使用
- **戻り値**: void
- **呼び出し元**: InterpolationTPS()

#### `InterpolationTPS`
- **処理概要**: Thin Plate Spline 法による気温補間
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: TemperatureInterpolation()

#### `InterpolationBarnes2Pass`
- **処理概要**: Barnes 法 2 パス補間（第1パス: 基本補間, 第2パス: 残差補間）
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: TemperatureInterpolation()

#### `TemperatureInterpolation`
- **処理概要**: 補間計算のエントリ関数。`AN.interp.method` で手法を切り替え（IDW/RBF/Barnes/TPS/Barnes2Pass/KrigingS/KrigingOK）
- **引数**: `const DataDEM& DEM`, `DataAMeDAS& DT`, `DataAnal& AN`, `const ConfigCore& CFG`, `const DateTime& dt`
- **戻り値**: void
- **呼び出し元**: RunOneMonth(), RunEvalInterp(), RunOptimParams()

---

### 9-5. OpticalFlow.cpp

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

### 9-6. HsOutputPNG.cpp（VER.0.4.6.491〜494 新設、§12-13）

PNG出力本体クラス `CHsOutputPNG`（AMT非依存ユーティリティ、§0-2準拠）の実装。文字・凡例インポーズは仮想 `overlayRow`、独自行生成は仮想 `buildRow8` を派生クラス（Output.cpp の `COutputPNG`）で注入する。

#### `CHsOutputPNG::Save8`
- **処理概要**: ConvMatrixParam（float/int16スケール済ソース）から8bitパレット/グレースケールPNGを出力。`nUseSIMD` でAVX-512/スカラー、`nSrcType` でfloat/int16を分岐し、int16は内部でfloat[℃]へ復元してfloatパスを共用
- **引数**: `const char* path`, `ConvMatrixParam& DAT`, `const uint8_t* palette`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: OutputProbImage(), OutputTempDistPNG()（Output.cpp、COutputPNG経由）

#### `CHsOutputPNG::Save16`
- **処理概要**: ConvMatrixParam（floatソース）から16bitグレースケールPNGを出力（オーバーレイなし）
- **引数**: `const char* path`, `ConvMatrixParam& DAT`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: DEM等の生データグレースケール出力

#### `CHsOutputPNG::SaveByBuilder8`
- **処理概要**: 行内容を仮想 `buildRow8()` で都度生成しながら8bitパレットPNGを出力（1行分のみのメモリfootprint）。気温＋等温線マージPNG（横2倍行）用
- **引数**: `const char* path`, `int w`, `int h`, `int alignedW`, `const uint8_t* palette`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: OutputMergedTempContourPNG()（Output.cpp、COutputPNG::SetMergedParams 設定後）

#### `CHsOutputPNG::SaveIndexed8`
- **処理概要**: 事前生成済みインデックスバッファから8bitパレットPNGを出力（行を作業バッファへ複製してから overlayRow を適用し元バッファを破壊しない）
- **引数**: `const char* path`, `const uint8_t* pIndex`, `int w`, `int h`, `int alignedW`, `const uint8_t* palette`
- **戻り値**: int — 0:成功, -1:失敗
- **呼び出し元**: OutputBoundContourPNG()（Output.cpp）

---

## 10. 設定キー一覧

設定ファイルはINI形式。未記載キーは `InitConfigCore()` / `InitConfigRun()` のデフォルトを維持する。値域は読み込み時に `HsClip()` で丸めるものが多い。ここでは運用・保守で重要な主要キーを記載する。

### 10-1. 基本・パス

| セクション | キー | 型 | デフォルト/範囲 | 対応メンバ | 内容 |
|---|---|---:|---|---|---|
| MAIN | TYPE | string | main想定 | `ConfigCore::meta.szType` | 設定ファイル種別。通常実行は `main` |
| PATH | DEM_PATH | path | 環境依存 | `path.szDemPath` | 解析用DEM GeoTIFF |
| PATH | TERRAIN_DEM_PATH | path | DEM_PATH | `path.szTerrainDemPath` | 地形補正用DEM。`TERRAIN_CORR>0`時に使用 |
| PATH | RESULT_PATH | path | 環境依存 | `path.szResultPath` | 出力先ディレクトリ |
| PATH | FILE_KEY | string | 環境依存 | `path.szFileKey` | 画像/GRD/GeoTIFF系のファイルキー |
| PATH | FILE_KEY_ST | string | 環境依存 | `path.szStatFileKey` | 統計/補助出力系のファイルキー |
| PATH | AMEDAS_PATH | path | 環境依存 | `path.szAMeDASPath` | AMeDAS年報バイナリ格納ルート |

### 10-2. 実行制御

| セクション | キー | 型 | デフォルト/範囲 | 対応メンバ | 内容 |
|---|---|---:|---|---|---|
| RUN | MODE | string | CONTINUOUS等 | `ConfigRun::nRunMode` | `CONTINUOUS` / `MULTI_YEAR` / `PARAM_SCAN` / `VALUE_SURVEY` / `EVAL_INTERP` / `OPTIM_PARAMS` / `AREA_CHAR` / `ST_NBR_EVAL`（`COMPARE` は VER.0.4.6.413 で廃止＝ロードエラー） |
| RUN | CONT_START | y/m/d/h | 2000/1/1/0相当 | `nStartYear`等 | 連続実行開始日時 |
| RUN | CONT_DUR | hour/day | 設定値 | `nRunCount`, `nRunDays` | 連続実行期間。日指定を優先 |
| RUN | MULTI_YEAR | y0/y1 | 設定値 | `nYearStart`, `nYearEnd` | 複数年実行の年範囲 |
| RUN | MULTI_MONTH | list | 1-12 | `nTargetMonthList` | 対象月リスト。`1/2/3/12`または`1,2,3,12`形式 |
| RUN | DAYS | list | ST_NBR_EVAL用 | `nTargetDayList` | ST_NBR_EVAL対象日。`/`と`,`区切りに対応 |
| RUN | HOURS | list | ST_NBR_EVAL用 | `nTargetHourList` | ST_NBR_EVAL対象時刻。`/`と`,`区切りに対応 |
| RUN | EVAL_METHODS | list | 0-7 | `nEvalMethods` | EVAL_INTERP/OPTIM_PARAMS評価対象・評価順・結果出力順。`/`と`,`区切りに対応 |
| RUN | PARALLEL_SLICES | int | 1 / 1-8 | `ConfigCore::exec.nParallelSlices` | 複数時刻同時補間数。日境界・時刻帯境界で分割。Barnes2Pass/TPSまたはEST_MODE=2では1へ自動矯正 |
| RUN | USE_NBR_FLOAT | int | 3 / 0-3 | `exec.nNbrFloatMode` | 0=全整数、1=SRCのみfloat、2=Wのみfloat、3=全float |
| RUN | ALTW_CACHE | int | 0 / 0-1 | `exec.nAltWCache` | 時刻帯別 `pNbrAltW` キャッシュ |

`RUN.MODE` は未記載/空値の場合のみ既定値 `CONTINUOUS` のまま扱う。非空の未知値は `LoadConfigRun` がエラーを返し、`LoadConfig` 失敗として当該cfgの実行をスキップする。

### 10-3. 解析・補間

| セクション | キー | 型 | デフォルト/範囲 | 対応メンバ | 内容 |
|---|---|---:|---|---|---|
| EXT_ANAL | PROB_MODE | int | 0 / 0-3 | `anal.nAnalProb` | 0=無効、1=閾値以下、2=閾値以上、3=両方 |
| EXT_ANAL | CENTROID | int | 0 / 0-1 | `anal.nAnalCentroid` | 境界面重心解析 |
| EXT_ANAL | AREA_STAT | int | 0 | `anal.nAnalAreaStat` | 閾値エリア統計 |
| EXT_ANAL | AREA_HIST | int | 0 / 0-1 | `anal.nAnalAreaHist` | 気温ヒストグラム集計 |
| EXT_ANAL | OPT_FLOW | int | 0 / 0-2 | `anal.nAnalOptFlow` | 0=無効、1=勾配法、2=Horn-Schunck |
| EXT_ANAL | PROB_ANNUAL | int | 0 / 0-2 | `anal.nProbAnnual` | 全期間確率バッファ確保。省メモリ時は0 |
| INTERP | METHOD | int | 3 / 0-7 | `interp.nInterpMethod` | 0=IDW、1=RBF、2=Shepard、3=Barnes、4=Barnes2Pass、5=TPS、6=Kriging_S、7=Kriging_OK |
| INTERP | COAST_L | float | 20000相当 | `AppConfig::fCoastL` | 海岸距離補正スケール[m] |
| INTERP | IDW_POWER | float | 1.5 | `interp.fIdwPower` | Shepard/IDW系べき乗 |
| INTERP | RBF_SIGMA | float | 80000 | `interp.fRbfSigma` | RBF sigma[m] |
| INTERP | BARNES_KAPPA_FAC | float | 2.0 | `interp.fBarnesKappaFac` | `kappa = d_nn^2 * factor` |
| INTERP | BARNES_GAMMA | float | 0.3 | `interp.fBarnesGamma` | Barnes2Passの第2パス係数 |
| PARAM | DIVERG_RAIN | int | 0 | `rain.nDivergRain` | 発散・降水照合 |
| PARAM | DIVERG_THRESHOLD | float | -2.2 | `rain.fDivergThreshold` | 収束判定閾値 |
| PARAM | RAIN_THRESHOLD | float | 1.5 | `rain.fRainThreshold` | 降水判定閾値[mm/h] |

### 10-4. 出力制御

| セクション | キー | 型 | デフォルト/範囲 | 対応メンバ | 内容 |
|---|---|---:|---|---|---|
| OUTPUT | DEBUG | int | 0 / 0-1 | `log.nOutLog_Debug` | デバッグ出力 |
| OUTPUT | COMMENT | int | 1 / 0-1 | `log.nOutLog_Comment` | 条件ログ出力 |
| OUTPUT | TIME | int | 1 / 0-1 | `log.nOutLog_Time` | 実行時間表示 |
| OUTPUT | KEY_TIME | int | 0 / 0-1 | `log.nOutLog_KeyTime` | キー項目実行時間表示（メンバ名は VER.0.4.6.500 で KetTime→KeyTime の typo 修正） |
| OUTPUT | ~~ANAL2~~ | — | — | （廃止） | VER.0.4.6.413で廃止。coast_L評価は MODE=OPTIM_PARAMS Phase6 へ統合（§1-9） |
| OUTPUT | ST_DAT | int | 0 / 0-1 | `log.nOutLog_StDat` | 観測点別NetCDF |
| OUTPUT | ST_NBR_DUMP | list | 任意 | `log.nStNbrDumpDT` | 近傍点ダンプ日時 `yyyymmddhh/...` |
| OUTPUT_VECTOR | FLOW | int | 0 / 0-1 | `gpkg.nOutGPKG_Flow` | GPKGベクトル出力 |
| OUTPUT_VECTOR | THINNING | int | 6 / 1-1000 | `gpkg.nCfgGPKG_Thinning` | GPKG時間間引き |
| OUTPUT_VECTOR | THINNING_XY | int | 20 / 1-1000 | `gpkg.nCfgGPKG_ThinningXY` | GPKG空間間引き |
| OUTPUT_GEOTIFF | TEMP | int | 0 / 0-1 | `geotiff.nOutGeoTIFF_Temp` | 気温GeoTIFF |
| OUTPUT_GEOTIFF | FLOW | int | 0 / 0-1 | `geotiff.nOutGeoTIFF_Flow` | フローGeoTIFF |
| OUTPUT_GEOTIFF | GEOTIFF3 | int | 0 / 0-1 | `geotiff.nOutGeoTIFF_Prob` | 確率GeoTIFF |
| OUTPUT_GEOTIFF | GEOTIFF4 | int | 0 / 0-1 | `geotiff.nOutGeoTIFF_Cont` | 境界マスクGeoTIFF |
| OUTPUT_GEOTIFF | THINNING | int | 6 / 1-1000 | `geotiff.nCfgGeoTIFF_Thinning` | GeoTIFF時間間引き |
| OUTPUT_GEOTIFF | THINNING_XY | int | 1 / 1-1000 | `geotiff.nCfgGeoTIFF_ThinningXY` | GeoTIFF空間間引き |
| OUTPUT_PNG | PNG2 | int | 0 / 0-2 | `png.nOutPNG_Temp` | 0=無効、1=気温+等温線別、2=マージ |
| OUTPUT_PNG | TEMP_MIN/MAX | int | -40/40 | `nCfgPNG_TempMin/Max` | 気温PNGレンジ |
| OUTPUT_PNG | PNG3 | int | 0 / 0-1 | `png.nOutPNG_Prob` | 確率PNG |
| OUTPUT_PNG | PNG4 | int | 0 / 0-1 | `png.nOutPNG_Cont` | 等温線PNG。`PNG2=1`時のみ有効 |
| OUTPUT_GRD | GRD | int | 1 / 0-1 | `grd.nOutGRD` | GRDバイナリ出力 |
| OUTPUT_GRD | THINNING | int | 1 / 1-1000 | `grd.nCfgGRD_Thinning` | GRD時間間引き |
| OUTPUT_GRD | CH_NUM / CH_TAGn | int/string | 任意 | `grd.nGrdChNum`, `grd.szGrdChTags` | GRD出力チャンネル指定 |

`TEMP_STAT.HOUR_OUT`、`TEMP_HIST.HOUR_OUT`、`THRESHOLD_ANAL.HOUR_OUT` は 0=無効、1=TXT、2=NetCDF、3=TXT+NetCDF。日別/月別/全期間系は現行ではTXT出力フラグとして扱う。

### 10-5. 閾値・エリア・補正

| セクション | キー | 型 | デフォルト/範囲 | 対応メンバ | 内容 |
|---|---|---:|---|---|---|
| THRESHOLD | VALUE | list | 最大16 | `threshold.nvThresholds` | 昇順の気温閾値リスト。`,`区切り。未記載時は閾値なし（`nCntThresholds=0`） |
| AREA | NUM | int | 0 / 0-32 | `area.nAreas` | エリア数 |
| AREA | AREA_n | string | 任意 | `area.areaDefs` | 矩形 `name,lat_s,lat_n,lon_w,lon_e` またはSHP `name,shp:path,field,value` |
| POI_CENTROID | COUNT | int | 0 / 0-16 | `poi.nPOIs` | POI数 |
| POI_CENTROID | POI_n | string | 任意 | `poi.poiDefs` | `name,lat,lon` |
| POI | SEARCH_R | int | 50 | `poi.nPoiSearchR` | 近傍境界探索半径[px] |
| LAPSE | EST_MODE | int | 0 / 0-2 | `lapse.nLapseEstMode` | 0=既存TABLE（標準）、1=月別推定（評価用）、2=月別×24時刻別推定（研究用、PARALLEL_SLICES=1へ矯正） |
| LAPSE | LAT | list | 33/40/999/999 | `lapse.lapse_lat_zone` | 緯度帯境界 |
| LAPSE | SEASON | list | 3/6/9 | `lapse_*_month` | 春/夏/秋開始月 |
| LAPSE | HOUR | list | 未指定 | `lapse.lapse_hour_th` | 4時刻帯使用時の境界 |
| LAPSE | VALUEz / VALUEz_b | list | デフォルト減率 | `lapse.lapse_rate_table` | 季節別/時刻帯別の気温減率[℃/m] |
| AREA_CONFIG | TERRAIN_CORR | int | 0 / 0-2 | `terrain.nTerrainCorr` | 地形特徴量補正。0=無効、1=値補正、2=Barnes地形類似度重み |
| AREA_CONFIG | TERRAIN_CORR_CLIP | float | 4.0 | `terrain.fTerrainCorrClip` | 地形補正項クリップ幅[℃] |
| AREA_CONFIG | TERRAIN_CORR_COAST_MIN | float | 0.0 | `terrain.fTerrainCorrCoastMin` | 適用最小海岸距離[m] |
| AREA_CONFIG | TERRAIN_CORR_L | float | 200.0 | `terrain.fTerrainCorrL` | `TERRAIN_CORR=2` の地形類似度重みスケール |
| AREA_CONFIG | SEASON / HOUR | list | LAPSE準拠 | `terrain_*` | 地形補正用の季節/時刻帯 |
| AREA_CONFIG | WINTER_CORR等 | list | 0 | `terrain.terrain_corr_table` | 季節別/時刻帯別の地形補正係数 |

---

## 11. 出力ファイル仕様

出力ファイル名は原則として `GetCommonOutputFileName()` で生成され、`RESULT_PATH`、`FILE_KEY_ST`、実行年/月条件、末尾名を組み合わせる。画像・GeoTIFF・GRD等の一部は `FILE_KEY` を使う。

### 11-1. 共通・補助出力

| ファイル | 生成条件 | 内容 |
|---|---|---|
| `*_analysis_summary.txt` | 初回エリアマスク構築後 | 実行条件、DEM情報、補間設定、出力ファイル一覧、エリア陸地ピクセル数 |
| `*_area_land_pixels.txt` | `AREA.NUM>0` | 全域と各エリアの陸地ピクセル数。全域は `AREA_DEF_MAX` |
| `*_stdat.nc` | `OUTPUT.ST_DAT=1` またはST_NBR_EVAL | 観測点別の実測値、補間値、近傍、地形特徴量、補正項 |
| `StNbrDump\{FILE_KEY_ST}_eval_manifest.tsv` | ST_NBR_EVAL | ST_NBR_EVALで処理した日時と出力先の対応 |
| `StNbrDump\YYYY\{FILE_KEY_ST}_stnbr_yyyymmdd_hh.tsv` | ST_NBR_EVAL | 観測点別の近傍点、距離、標高差、altW、重み、LOOCV寄与確認 |
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
| 気温PNG | `OUTPUT_PNG.PNG2=1` | `TempPNG\YYYY\KEY_temp_YYYYMMDD_HH.png` |
| 等温線PNG | `OUTPUT_PNG.PNG2=1 && PNG4=1` | `BoundMask\ContourPNG\..\KEY_bound_ctr_YYYYMMDD_HH.png` |
| 気温+等温線PNG | `OUTPUT_PNG.PNG2=2` | `TempContourPNG\YYYY\KEY_tempctr_YYYYMMDD_HH.png` |
| 確率PNG | `OUTPUT_PNG.PNG3=1 && PROB_MODE>0` | `KEY_prob_YYYYMM_+nnC_lo/hi.png` |
| 気温GeoTIFF | `OUTPUT_GEOTIFF.TEMP=1` | `KEY_Temperature_YYYYMMDD_HH.tif` |
| Flow GeoTIFF | `OUTPUT_GEOTIFF.FLOW=1` | `KEY_FlowU_*.tif`, `KEY_FlowV_*.tif`, `KEY_Diverg_*.tif`, `KEY_Vorticity_*.tif` |
| 境界マスクGeoTIFF | `OUTPUT_GEOTIFF.GEOTIFF4=1` | `BoundMask\...\KEY_bound_mask_YYYYMM_DD_HH_+nnC.tif` |
| 確率GeoTIFF | `OUTPUT_GEOTIFF.GEOTIFF3=1 && PROB_MODE>0` | `KEY_prob_YYYYMM_+nnC_lo/hi.tif` |
| フローGPKG | `OUTPUT_VECTOR.FLOW=1 && OPT_FLOW>0` | `KEY_Vector_YYYYMMDD_HH.gpkg` |
| GRD | `OUTPUT_GRD.GRD=1` | `BIN_TEMP` または `BIN_DATA` 配下。TEMP/FLOW/BOUND等をBLOSC2圧縮格納 |

GRDのTEMPは `int16_t`、スケールは `HSS_GRD_FAC_TEMP=500`（0.002℃精度）。BOUNDは `HSS_GRD_TYPE_B16` で、bit c が `THRESHOLD.VALUE[c]` の境界マスクを表す。

---

# 12. VER.0.4.6.140時点の主な更新点

## 12-1.
- `HSS_GRD_FAC_TEMP=500.0f` に変更済み。気温GRD/int16温度面は0.002℃精度。
- `pAnalTemp` は通常Barnes経路では確保しない。Barnes2Pass/TPSの一時バッファとしてのみ確保する。
- 補間手法定数は `0=IDW, 1=RBF_GAUSS, 2=SHEPARD, 3=BARNES, 4=BARNES2PASS, 5=TPS, 6=KRIGING_S, 7=KRIGING_OK`。旧資料の `INTERP_BARNES_H` は現行コードに存在しない。
- `PARALLEL_SLICES` により複数時刻同時補間を行う。日境界・lapse時刻帯・地形補正時刻帯をまたぐ場合はバッチを分割する。
- `USE_NBR_FLOAT` により近傍点テーブルのSRC/Wをfloat/intへ切り替える。`ALTW_CACHE` は時刻帯別の高度/地形補正項をキャッシュする。
- `AREA_CONFIG.TERRAIN_CORR` により地形特徴量を使った補正を追加。`TERRAIN_DEM_PATH` が必要。
- `MODE=ST_NBR_EVAL` を追加。観測点近傍、地形特徴量、補正項の検証用に通常出力を抑止して限定日時を処理する。
- 気温統計・気温ヒストグラム・閾値統計の一部はNetCDF出力にも対応している。

## 12-2. VER.0.4.6.142→200 リファクタリングの主な更新点（2026-05-23 クローズ）

`system_review_V0461_20260516.md` に基づくリファクタリングを VER.0.4.6.200 で完了。本概要書も以下の構造変更を反映している:

- **`DataAnal` を5サブ構造体＋`LapseState`に分割**（Phase6＋保守リファクタ）: `nbr`(NbrTable) / `stats`(TempStats) / `prob`(ProbMap) / `area`(AreaStats) / `multi`(MultiPool) に加え、lapse推定・時刻帯カーソル状態を `lapse`(LapseState) へ分離。参照は `AN.area.vAreaStat`、`AN.lapse.nCurHour` のようにサブ構造体経由。
- **設定構造の整理**（Phase5・3-3＋2026/06/06保守）: `[INTERP]` 5float と `[OUTPUT_GRD]` チャンネル配列は `ConfigCore` に統合済み。さらに `ConfigCore` は `meta/path/exec/anal/interp/rain/log/gpkg/geotiff/png/grd/tempStat/tempHist/thrOut/threshold/area/poi/lapse/terrain` のカテゴリ別サブ構造体へ分割済み。`AppConfig` は `{ConfigCore C; ConfigRun R}` の薄いラッパで、**ConfigCore の既定値・読込は `InitConfigCore`/`LoadConfigCore` に一元化**（`InitAppConfig`/`LoadConfig` は C と R の Init/Load を束ねるだけ）。引数規約: **Run* エントリ＝AppConfig／リーフ計算・出力＝ConfigCore**。
- **時間ループ関数の再編**（柱C、VER.0.4.6.178）: 旧 `RunOneBlock` を **`RunOneMonth`（月単位純化、月境界判定を内部に持たない）** へ。月初準備＝`SetupForMonth`、月末出力＝`FlushMonth`、バッファ初期化＝`ResetAnalBuffers` に抽出。さらに `RunOneMonthPrepareStep` / `RunOneMonthInterpolateStep` / `RunOneMonthOutputStep` と、通常出力・フロー出力・非フロー出力の小関数へ段階分割済み。**CONTINUOUS=期間動作**（月チャンク分割ループ＋TOTAL_OUT、MONTH_OUT廃止）／**MULTI_YEAR=月動作**（月別独立リセット＋月別ProbImage）。
- **責務境界の明確化**（柱A）: 汎用出力（`OutputTempGeoTIFF`/`OutputFlowGeoTIFF`/`OutputFlowGPKG`/`OutputAnalyzeBinary`）を一度 AnalCore.cpp へ、補間準備 `BuildPixelContext` を Interpolation.cpp へ移設。**amt.cpp 完全解体・ライブラリ化は撤回**（§0-4）。その後 VER.0.4.6.300 で出力系は `Output.cpp` に分離。
- **テーブル駆動化**: 補間手法名（`InterpMethodInfo`/`GetInterpMethodName`）、出力定義（`g_TempTextDefs`）。**LOOCV ゾーン分類は `ZoneClassifier` 構造体に集約**（5-B-8）。旧 `ComparInterpolationMethod` は VER.0.4.6.413 で削除済み。
- **出力ファイル名の統一**（C3d）: 統計TXT/NC・TMD・LOG が共通ステム `GetRunStemName()` を共有。MULTI_YEAR=`KEY_YYYY-YYYY_MM`／その他（CONTINUOUS等）=`KEY_YYYYMMDD_HH`。
- **RunMode別固定値の単一窓口化**: `ApplyRunModeFixedConfig()` を `LoadConfig` 末尾で適用（CONTINUOUS の MONTH_OUT 無効化、ST_NBR_EVAL の出力抑止等）。`RUN.MODE` の文字列対応は `RunMode` enum と `RunModeDef` テーブルへ集約し、旧 `COMPARE` は廃止モードとしてロードエラーを返す。
- **検証規約**: GRD(.tva) は Blosc2 並列圧縮で**バイト非決定的**のためリグレッション比較に使えない。TXT/NC の値・SHA256 で判定する。

## 12-3. VER.0.4.6.300 出力系分割（2026-05-26）

VER.0.4.6.300 では `AnalCore.cpp` から出力系を `Output.cpp` へ分離した。分割後は `AnalCore.cpp`、`Output.cpp` とも4000行台で、内包関数数も概ね同程度となった。`AnalCore.cpp` の追加分割は可能だが、現時点では分割しすぎによる保守性低下を避けるため保留する。

- **Output.cpp へ移動した範囲**: 出力共通ヘルパ（`MakeYearDir`/`MakeSubDir`/`GetRunMonthString`/`GetRunStemName`/`GetCommonOutputFileName`）、出力初期化・flush（`InitAllTextOut`/`FlushAllTextOut`/`InitHourlyStatNC`/`InitProcTimeMeas`）、気温統計TXT/NetCDF、気温ヒストグラム、閾値・確率・エリア系出力、境界・PNG・画像出力、観測点ログ・ST_NBR_DUMP出力、GeoTIFF/GPKG/GRD出力。
- **Output.cpp へ移動した出力状態**: `OutputContext` にテキスト出力、観測点ログ/時別NetCDF、画像文字描画、画像凡例を集約。`Output.cpp` 外から `g_aTOUT` / `g_TR` / `g_IL*` / NetCDF系グローバルを直接参照せず、`InitAllTextOut` / `FlushAllTextOut` / `InitHourlyStatNC` / `InitOutputStationTemperature` / `InitOutputRenderDefaults` / `LoadOutputRenderConfig` / `InitOutputLegends` 経由で操作する。`Output.cpp` 内部の主要出力処理は `*Core(OutputContext& CTX, ...)` に分離し、公開APIラッパーだけを `OCX` 境界とする。処理時間TMD（`g_TML`）とログ（`g_LF`）は引き続き全体共通状態として `AnalCore.cpp` 側に置く。
- **AnalCore.cpp に残した範囲**: `LoadConfigCore`/`InitConfigCore`、パレット設定、データ読込、DEM/AMeDAS形式データ処理、`SetupForMonth`/`ResetAnalBuffers`/`FlushMonth`（内部はエリア閾値/気温/確率flushへ分割）、`FlushProvMonthToAccum`、`EvalDivergRain` 関連。
- **amt.cpp に残した範囲**: `LoadConfig`、`ApplyRunModeFixedConfig`、全 `Run*`、`AppendStNbrEvalManifest`/`WriteStNbrEvalManifestHeader` 等の実行モード固有処理。（`ComparInterpolationMethod`/`Eval_Coast_L*`/`ScanFactorTPS`/`EvalInterp`(static) は VER.0.4.6.413 で削除、§1-9）
- **ビルド変更**: `mk.bat` のコンパイル対象に `Output.cpp` を追加。

## 12-4. VER.0.4.6.400後 KrigingS 実装・評価（2026-05-27）

VER.0.4.6.400後の作業で、未完成扱いだった `INTERP_KRIGING_S` を簡易Krigingとして実装した。手法番号は将来拡張と既存cfgの互換性を考慮し、Barnes標準を `METHOD=3`、KrigingSを `METHOD=6` とした。

- **実装内容**: `CalcWeight()` に球形バリオグラム由来の共分散重みを追加し、既存の近傍点重み補間経路に載せた。これは Ordinary Kriging の行列解法ではなく、距離共分散を重みに使う簡易Krigingである。
- **並列対応**: `INTERP_KRIGING_S` は `INTERP_BARNES` と同じ逐次/複数時刻並列パスに対応する。`PARALLEL_SLICES=8` の北海道250m評価で速度劣化は確認されていない。
- **速度評価**: 1999〜2001年1月（3か月、北海道250m）では Barnes1Pass 46.8秒、KrigingS 45.8秒。`TemperatureInterpolation` は Barnes1Pass 1.68msec/step、KrigingS 1.67msec/step で実用上同等。
- **画質評価**: PNG/GeoTIFF目視では等温境界に細かな差異はあるが、どちらが優位か判断できる明確差は確認されていない。
- **精度評価**: `EVAL_INTERP` のLOOCVでは、KrigingSはBarnes1Passを上回る結果ではない。2026-05-26評価では全体LOOCV RMSEが Barnes1Pass 1.5107、KrigingS 1.6062。
- **結論**: KrigingSは「比較評価可能な簡易Kriging」として実装完了。ただし現時点ではBarnes1Passから置き換える根拠はなく、標準手法は引き続きBarnes1Passとする。

## 12-5. Ordinary Kriging近似（KRIGING_OK）実装方針（2026-05-27）

VER.0.4.6.400後の追加検証として、`INTERP_KRIGING_OK = 7` を追加した。これは厳密な全点Ordinary Krigingではなく、既存の近傍点テーブルを利用する近似実装である。

- **実装単位**: 月初の `BuildWeightTable()` で、各陸地ピクセルの近傍上位 `KRIGING_OK_NBR=8` 点を使い、Ordinary Kriging制約 `sum(lambda)=1` の線形方程式を解く。
- **共分散モデル**: `KRIGING_S` と同じ球形共分散を使用する。レンジは `range=sqrt(barnes_kappa)` とし、`BARNES_KAPPA_FAC` の評価軸を流用する。
- **海岸距離補正**: target-station側の共分散ベクトルに既存の `GetCoastFactor()` を掛ける。station-station共分散には海岸補正を入れない。
- **重み保持**: Kriging重みは負値を取り得るため、`KRIGING_OK` 使用時は `NBR_MODE_W_FLOAT` を強制する。時別補間本体は既存の `InterpolationMain` / `InterpolationMainMulti` を流用する。
- **欠測時動作**: 時別欠測で有効近傍点が減る場合は、既存経路と同様に有効点の `sum(w*T)/sum(w)` で再正規化する。分母が極小の場合は欠測扱いとする。
- **位置づけ**: 計算コストと精度改善の評価用。観測点密度不足によりBarnes1Passを上回る保証はないため、標準手法は引き続きBarnes1Passとする。
- **安定化**: `KRIGING_OK_FALLBACK_DEN=0.05` を `Constants.h` に定義し、OK重みの有効和が小さい場合は同じ有効近傍点集合で `KRIGING_S` 重みに退避する。これはLOOCV時の自己点除外で負重み・小分母が作る極端外挿を抑制するための安全弁である。
- **評価結果**: `HOK250V46404`（北海道250m、1999-01-20 00:00、EVAL_INTERP）では `KRIGING_OK` の通常RMSEは 0.2594 と低いが、LOOCV RMSEは 2.3711、MAEは 1.7021 で Barnes1Pass（RMSE 1.5107、MAE 1.1460）を下回った。通常RMSEの良さは観測点自己再現性を示すもので、未観測点推定性能はLOOCVを優先して判断する。
- **結論**: `KRIGING_OK` は観測点自己再現では有利だが、近傍固定OK重みの負値・外挿性によりLOOCVではBarnes1Passを上回らない。標準手法はBarnes1Passを維持し、`KRIGING_OK` は比較評価・研究用手法とする。

## 12-6. 高度補正係数推定（LAPSE EST_MODE）実装・評価（2026-05-27）

標高を共変数として扱う簡易的な regression kriging 相当の検証として、既存の高度補正経路に `LAPSE.EST_MODE` を追加した。これは新しい補間手法ではなく、`lapse_rate` の決め方を切り替える機能である。HOK250V46412対応で全陸地ピクセル別の共通キャッシュ `pLandLapse` を追加し、通常補間の `pNbrAltW = pLandLapse × altDiff + terrainTerm`、Barnes2Pass、TPSのピクセル側高度補正が同じlapse参照元を使うよう整理した。

- **設定**: `EST_MODE=0` は既存TABLE、`EST_MODE=1` は月別推定、`EST_MODE=2` は月別×24時刻別推定。
- **推定方法**: 月内の観測値と観測点標高から `T = a + b × z` の単回帰を行い、傾き `b` を推定lapseとして使用する。`EST_MODE=2` は同じ処理を時刻別に行う。
- **安全策**: サンプル数・標高分散の下限を設け、推定lapseを `LAPSE_EST_MIN_RATE=-0.009` ～ `LAPSE_EST_MAX_RATE=-0.003` [℃/m] にクリップする。`EST_MODE=2` は時刻別に `pLandLapse` と `pNbrAltW` が変化するため、`PARALLEL_SLICES>1` 指定時は `LoadConfigCore()` が警告を出して `PARALLEL_SLICES=1` に矯正する。
- **評価結果1**: `HOK250V46405A～C`（クリップ前、推定下限 -0.012）では推定lapseが `-0.011031` となり、Barnes1PassのLOOCV RMSEは `EST_MODE=0: 1.5107`、`EST_MODE=1: 1.5454`、`EST_MODE=2: 1.5739` と悪化した。全域一括推定では強すぎる減率を拾い、過補正になる。
- **評価結果2**: `HOK250V46406A～C`（推定範囲 -0.009～-0.003）では、`EST_MODE=1/2` とも推定lapseが `-0.009000」にクリップされ、Barnes1PassのLOOCV RMSEは `EST_MODE=0: 1.5107`、`EST_MODE=1: 1.4941'、`EST_MODE=2: 1.4941' となった。小幅改善はあるが、時刻別推定の追加効果は確認されていない。
- **共通化後評価**: `HOK250V46412_BASE0/BASE1` と `HOK250V46412_R0/R1` の `eval_summary`/`eval_loocv` はSHA256一致し、`TERRAIN_CORR=0/1, EST_MODE=0` では互換性が保たれた。`pLandLapse` 追加により北海道250m条件のメモリ見積もりは約750MBから約755MBへ増加した。`HOK250V46412_R2`（TERRAIN_CORR=2、EST_MODE=1）ではBarnes1Pass LOOCV RMSEは `1.5596`、Barnes2Pass `1.7321`、TPS `1.8999` となり、Barnes2Pass/TPSにも推定lapseが反映された。`HOK250V46412_R3`（MULTI_YEAR、EST_MODE=2、PARALLEL_SLICES=8指定）ではログ上で `PARALLEL_SLICES=1` への矯正と `hour_valid=24/24` を確認した。
- **結論**: 標準運用は従来条件 `METHOD=3`、`EST_MODE=0` を維持する。`EST_MODE=1` は年月・地域により改善する可能性があるが、観測点配置や内陸/沿岸差を拾って悪化するリスクもあるため評価用機能とする。`EST_MODE=2` は現時点で `EST_MODE=1` を上回る根拠がなく、研究用扱いとする。実装上は全補間方法のlapse参照元を `pLandLapse` に寄せたため、制限付き機能ではなく共通高度補正として扱える。

## 12-7. 地形補正方式（TERRAIN_CORR=1/2）実装・評価（2026-05-27）

`AREA_CONFIG.TERRAIN_CORR` を `0/1/2` の3状態に整理した。`TERRAIN_DEM_PATH` は地形特徴量抽出用DEMであり、現在の北海道評価では50mメッシュを使用している。解析エリア外近傍点の地形特徴量も取得できるよう、地形DEMは解析範囲外も含めて作成する。

- **`TERRAIN_CORR=0`**: 地形補正なし。標高補正は既存の `lapse_rate × altDiff` のみ。
- **`TERRAIN_CORR=1`**: 地形特徴量に基づく値補正。`pNbrAltW = lapse_rate × altDiff + terrainTerm` として補正項を加算する。補正係数は観測点の情報と観測データから統計的に求めた値で、観測点LOOCVでは最も改善する一方、観測点が存在しない標高・地形レンジでは外挿になる。過剰補正は `TERRAIN_CORR_CLIP` で抑制するが、`CLIP=1.0` も経験的な安全弁であり物理的に一意な根拠値ではない。
- **`TERRAIN_CORR=2`**: 地形類似度重み。値補正は行わず、Barnes1Passの近傍重みに `exp(-(terrainDiff^2)/(2*TERRAIN_CORR_L^2))` を乗じる。値を直接動かさないため外挿リスクは小さいが、今回評価では改善幅は限定的。`INTERP.METHOD=3`（Barnes1Pass）専用であり、他手法では `LoadConfigCore` が警告を出して `TERRAIN_CORR=0` に強制補正する。
- **評価結果（北海道250m、1999-01-20 22:00、EVAL_INTERP）**: `TERRAIN_CORR=1` のBarnes1Pass LOOCV RMSEは `1.5107`、MAEは `1.1460`。`TERRAIN_CORR=2` は `TERRAIN_CORR_L=50` が全体RMSE最良で `1.5759`、`L=35` がMAE最良で `1.1901`。`L=25` は過強調でRMSE `1.5851` と悪化した。
- **ゾーン別評価**: `TERRAIN_CORR=2` は低地側では小幅改善するが、中高地（200-600m）では `L=50` でLOOCV RMSE `2.1471`、`L=200` で `2.1407` と、強い地形重みほど安定性が下がる傾向がある。`TERRAIN_CORR=1` は中高地でも RMSE `1.8773` と明確に良いが、これは値補正が直接効くためで、非観測域での外挿妥当性は別途確認が必要。
- **ST_NBR_EVAL確認**: `TERRAIN_CORR=2` の `L=50` と `L=200` は出力行数・近傍構成は同一で、重み分布のみ変化する。`L=50` は一部地点で自己点や地形類似点への重み集中が強まるが、平均有効近傍数はおおむね8点台で、1点支配へ極端に崩れてはいない。
- **結論**: 観測点LOOCVでは `TERRAIN_CORR=1` が最良。ただし公開・長期解析では、非観測高標高域への外挿リスクを明記し、`TERRAIN_CORR=0/1` の長期統計感度分析で主要な傾向結論が変わらないかを確認する。`TERRAIN_CORR=2` は物理的には自然だが改善幅が小さいため、比較実装・補助評価用として残す。

## 12-8. EVAL_METHODSと評価モードの整理（2026-05-27）

`RUN.EVAL_METHODS` は `EVAL_INTERP` / `OPTIM_PARAMS` の評価対象・評価順・結果出力順を制御するリストとして明確化した。例: `EVAL_METHODS=3/6/7` は Barnes1Pass、KrigingS、KrigingOK の順に評価し、`eval_loocv` の列順および `eval_summary` の出力順も同じになる。TPSは計算時間が大きいため、通常の地形補正評価では除外する。

- 未指定時は `EVAL_INTERP` では `3`（Barnes1Pass）のみをデフォルトとする。
- `RunEvalInterp()` は実行開始時に `EVAL_METHODS: 3(BARNES_1PASS), ...` の形で評価順をログ出力する。
- `ST_NBR_DUMP` の実体ファイル名は `FILE_KEY_ST` を使用するよう修正し、`ST_NBR_EVAL` のmanifestと実体TSVのパスを一致させた。`FILE_KEY` は画像、GeoTIFF、独自バイナリ等の地理データ出力用、`FILE_KEY_ST` は統計・評価・ログ系出力用である。

## 12-9. LOOCV一本化・COMPARE廃止・coast_L統合・ANAL2廃止（VER.0.4.6.413、2026-05-30）

評価系に二系統存在した LOOCV を、本番補間ベースの単一系統へ統合した。従来の static `EvalInterp`（手計算で近傍重みを再構成する簡易LOOCV）は **地形補正（TERRAIN_CORR）も KrigingS/KrigingOK も反映しない**ため、本番補間（`TemperatureInterpolation`）と数値が一致しない問題があった。これを廃止し、観測点をマスクして本番補間を再実行する方式（`RunEvalInterp`/`RunOptimParams`）へ一本化した。

- **COMPARE（mode3）廃止**: `RunCompare`/`ComparInterpolationMethod` を削除。手法別LOOCV比較は上位互換の `EVAL_INTERP`（手法可変・Kriging対応・地形補正込み・緯度/経度/標高ゾーン別集計）に統合。`MODE=COMPARE` 指定時は `LoadConfigRun` がエラーを返して停止する（mode3 は欠番、再採番しない）。
- **coast_L 評価の OPTIM_PARAMS 統合**: 旧 `Eval_Coast_L`（`ANAL2=1` 時に通常運用ループ内から呼ばれていた手計算スキャン）と未使用の `Eval_Coast_L_2Time`/`ScanFactorTPS` を削除。coast_L スキャンは `RunOptimParams` の **Phase6** として追加し、本番LOOCV（`DoLOOCV`）で評価する。`optim_result` に `COAST_L=` と境界警告を出力。Phase6 は Phase1 で求めた Barnes 最適 `kappa_fac` の下でスキャンし、推奨 `BARNES_KAPPA_FAC` と `COAST_L` の前提を一致させる。
- **ANAL2（`nOutputAnal2`）廃止**: 旧 `Eval_Coast_L` 呼び出し用フラグを `ConfigCore` から削除。cfg の `ANAL2=` 行は不活性（読み捨て）。通常運用ループから coast_L スキャン呼び出しが消え、通常パスが純化した。
- **ZoneStat 共有化**: `RunEvalInterp`/`RunOptimParams` に重複していた LOOCVゾーン統計の構造体・集計（`ZoneStat` + RMSE/MAE/Bias/NSE/r 算出）を `ZoneStat`/`ZoneStatAccum`/`ZoneStatFinal`/`AccumulateZoneStats`/`WriteZoneStatRows`（amt.cpp ファイルスコープ）へ一本化。
- **ラプラシアン平滑度・手法別GeoTIFF の EVAL_INTERP 移植**: 旧 COMPARE 固有だった空間的滑らかさ指標（`CalcLaplacianSmoothness`）を `eval_summary` の末尾2列 `LapMean`/`Lap95p`（グリッド全体指標のため「全体」行のみ値）として追加。手法別の気温 GeoTIFF は `geotiff.nOutGeoTIFF_Temp`(TEMP) 有効時に `{FILE_KEY_ST}_eval_temp_{手法名}_YYYYMMDD_HH.tif` で出力。
- **検証（HOK250、1999-01-20 00:00）**: Barnes 全体 LOOCV RMSE は `EVAL_INTERP`＝`OPTIM_PARAMS Phase1`＝`OPTIM_PARAMS Phase6(coast_L=20000)`＝**1.5107** で end-to-end 一致。KrigingOK の LOOCV RMSE 2.3711 も従来 HOK250V46404 と一致。旧二系統の乖離が解消したことを実データで確認した。
- **`LoadAMeDAS` の不具合修正**: 観測データ本体ファイルを開けない場合に NULL ポインタへ `fseek`/`fread` していた箇所を、`FreeAMeDAS`＋`return -1`（ロード失敗）へ修正（クラッシュ回避、同関数内の他失敗パスと同規約）。

## 12-10. 保守リファクタリング（VER.0.4.6.413〜440、2026-05-30〜06-05）

§1-9 完了後に実施した保守整理。**すべて byte-identical**（HOK250 各構成で統計出力が変更前と一致することを確認済み。GRD は Blosc 非決定のため TEMP 値で確認）。

### 引数規約の統一（pCFG 廃止）
- **`DataAnal.pCFG`（ConfigCore* 二重参照）を完全廃止**。従来「設定値を AN に複製しない」目的で `InitAnal` が `AN.pCFG = &CFG` を保持していたが、`AN.pCFG->` と引数 `CFG.` の二重アクセスが「どちらが正典か」を曖昧にしていた。設計を **「設定は全関数が引数 `ConfigCore` で受け取る」** に一本化。
  - 第1段: `WriteGRD`/`ReadGRD`/manifest 2関数（`AppConfig`→`ConfigCore`）。`WriteGRD` は `(…,dt,const ConfigCore& CFG,prm)` に変更し、**出力タグ列構築（`AN.nFgFlow<1` 時 TEMP/BOUND のみ）と mode 導出を内包**。これに伴い旧 `OutputAnalyzeBinary` の「`nPos` 個構築なのに `nTags=nGrdChNum` を渡す」**潜在バグ（未初期化tag参照）も解消**。`ReadGRD` は `const_cast<pAN->pCFG>` を廃し引数 `ConfigCore&` へ直接書き戻し。
  - 第2段: CFG 引数を既に持つ関数内の `AN.pCFG->` を `CFG.` へ統一（Output/AnalCore/Interpolation/NetCDF 等 計63箇所、`IsMultiBatchBoundary` は引数 `C` へ）。
  - 第3段: CFG 引数を持たなかった `CalcBoundaryMaskOnly`/`CalcBoundaryMaskCoreInt16`/`CalcBoundaryMaskAndCollect`/`CalcBoundaryMaskAndCollectMulti`/`FlushProvMonthToAccum`/`CalcProbMonthStat` に `const ConfigCore& CFG` を追加。
  - 最終段: inline ヘルパ `GetTimeBand` を `(const ConfigCore& CFG,int hour)`（AN除去）、`GetLapseRate` を `(const DataAnal& AN,const ConfigCore& CFG,float lat,int month)` に変更。これに伴い `BuildNbrTable`/`UpdateLandLapseTable`/`InterpolationTPS`/`InterpolationBarnes2Pass` の署名にも `CFG` を追加（呼出元は全て CFG 保持で波及収束）。**結果 `DataAnal` から設定ポインタが消滅**。
- **`const` 化**: `OutputTempGeoTIFF`/`OutputFlowGeoTIFF`/`OutputFlowGPKG`/`OutputAnalyzeBinary` の `DateTime&`→`const DateTime&`、`EvalDivergRain` の読取専用引数（DT/AN/FL）を const。

### 重複コード共通化
- **`CalcDivergSkillScores(TP,FP,FN,TN)`**（AMTCore.h）: HR/AR/CSI/HSS 算出を `RunMultiYear`/`RunParamScan`/`PrintDivergRainAccum` の3重複から集約（HSS は `E<1.0` ガード統一）。
- **`CalcBoundaryBits()`**（Interpolation.cpp）: 等温線境界ビット算出（中心と右/下隣の `T<閾値` XOR）を `CalcBoundaryMaskCoreInt16`/`CalcBoundaryMaskAndCollect`/`…Multi` の3コピーから集約。
- **`ZoneStat`/`ZoneStatAccum`/`ZoneStatFinal`/`AccumulateZoneStats`/`WriteZoneStatRows`**（§1-9・amt.cpp）: LOOCV ゾーン統計の累積・確定・TSV行出力を `RunEvalInterp`/`RunOptimParams` から共有。
- **`EvaluateOneStationLOOCV`**（amt.cpp）: 観測点1点の欠測化、本番補間再実行、推定値取得、観測値復元を集約し、`RunEvalInterp`/`RunOptimParams` のLOOCVマスク処理重複を排除。
- **`RunModeDef`**（AnalCore.cpp）: `RUN.MODE` 文字列、内部 `RunMode` 定数、廃止フラグをテーブル化し、`LoadConfigRun` の `strcmp` 連鎖を排除。
- **Run* プロローグ共通化**: `RunValueSurvey`/`RunAreaChar` のロード〜近傍構築を `PrepareEvalModeData` に集約。

### 削除（デッドコード等）
- `CollectAreaCntFromShtTemp`（呼び出し0。areaCnt収集は `InterpolationMain`/`CalcBoundaryMaskAndCollect[Multi]` が担う）。
- デッドフィールド `NbrTable.pBarnes2PassLapse`（常に nullptr）。
- 未使用定数 `LATBAND_*`/`ELVBAND_*`（旧 static EvalInterp 専用、Constants.h）。

### 不具合・コメント修正
- **7-2**: `CalcHornSchunckFlow`（`OPT_FLOW=2`）で `pMagnitude` を作業バッファ（Ix）に流用後 `|∇T|` へ復元していなかった点を修正（`CalcDivergVorticitySIMD` 直前に `sqrt(Ix²+Iy²)` を確定）。
- **`value_survey` 出力名**: 統計系のため `FILE_KEY`→`FILE_KEY_ST`（`path.szStatFileKey`）に修正。
- **設定読込の明確化**: `RUN.MODE` は未記載/空値のみ既定値 `CONTINUOUS` を維持し、非空未知値はロードエラーとする。複数cfg連続実行では `argv[2]` 以降も読込前に `InitAppConfig` を行い、未記載キーの前cfg継承を廃止する。`THRESHOLD.VALUE` 未記載時は `nCntThresholds=0` とし、直前文字列の誤読を防止する。RUN系の整数リストキーは `CHsIniFile::GetIntList` へ共通化済み。
- **気温GeoTIFFの時刻・スロット明示化**: `OutputTempGeoTIFF` に出力対象スロット引数を追加。フロー有効時は `dtPrev` + `idx_order[1]`、非フロー時はPNGと同じ `dtCurDT` + `idx_order[0]` を出力する。
- **スロット参照・出力const・InitAnal契約整理**: `GetCurrentTempSlot` / `GetPrevTempSlot` を追加し、主要な現在/前時刻参照をヘルパ経由へ整理。出力4関数は読み取り主体の `DataAnal` / `DataFlow` を `const` 参照化。`InitAnal` は呼出側0クリア済み前提で先に `FreeAnal` し、既存確保済みANの再初期化時にメモリを失わない契約へ修正。
- **`BuildAreaMask` の段階分割と警告強化**: SHP読込、BBOX準備、マスク生成、陸地数集計、サマリ出力を静的ヘルパへ分割。SHPファイルを開けない、レイヤ0がない、指定フィールドがない、条件に一致するリングが0件、構築後のエリア陸地ピクセルが0件の場合に `WARNING` を出す。指定フィールドなしの場合は全フィーチャ誤採用を避けるため当該SHPエリアを空扱いにする。
- **`RunChunkConfig` / `ConfigRain` の責務分離**: 月内チャンク実行範囲（`startDT`/`nSteps`）を `RunChunkConfig` として新設し、`RunOneMonth` はこれを受け取る。`ConfigRain` は雨量相関閾値（`rain.fDivergThreshold`/`rain.fRainThreshold`）専用に戻し、月内実行範囲との混在を解消。
- `exec.nParallelSlices` 補助モード固定値を 0→1 に統一、`AMeDAS_TMDT.Rain` 単位コメント（mm→0.1mm）、`PROB_ANNUAL` のコメント/clip 整合（0..2→0..1、挙動不変）、`PixCtx_*` デフォルト値コメント、`CalcWeight(KRIGING_S)` 遠方重み床の意図コメント 等。

---

## 12-11. リファクタリング4B（VER.0.4.6.470〜475、2026-06-06）

`system_review_V0464_20260606.md` のレビューに基づく保守整理。**6-1 と 7-4 を除き byte-identical**（HOK 1999-2001 正常系で変更前後のデータTXTが完全一致を確認）。6-1 は失敗月のみの挙動差、7-4 は中央値の丸め統一（本体C++＋変換スクリプト両側を同時修正し一致維持）。

### DataAnal 機能別サブ構造体化（11-2-1）＋ nvThresholds 改名（11-2-2）
- **`DataAnal` トップレベルの残フィールドを3サブ構造体へ分離**（既存の `nbr`/`stats`/`prob`/`area`/`multi`/`lapse`/`interp` に追加）:
  - **`pixctx`（`PixCtxParam`）**: 近傍点探索パラメータ。旧 `PixCtx_*` を接頭辞除去し `min_station`/`radius_init`/`radius_max`/`extent_range`/`extent_step1`/`extent_step2`/`extent_min_station`。アクセスは `AN.pixctx.radius_init` 等（`InitAnal` 設定 / `BuildNbrTable` 参照）。
  - **`tps`（`TpsWork`）**: TPS補間の作業領域。旧 `tps_w`/`tps_a`/`tps_stLon`/`tps_stLat`/`tps_valid`/`tps_n` を接頭辞除去し `w`/`a`/`stLon`/`stLat`/`valid`/`n`。アクセスは `AN.tps.w` 等。※正則化パラメータ `tps_lambda` は `InterpConfig` 側（`AN.interp.tps_lambda`）で本構造体に含めない。
  - **`ring`（`RingBuffer`）**: 気温int16リングバッファ＋等温線境界マスク＋スロット管理。`pShtTempPool`/`mShtTemp[]`/`idx_order[]`/`mBoundPacked[]`/`nvThresholdsScaled[]`/`nParallelSlices`/`nShtTempSlots`。アクセスは `AN.ring.mShtTemp`/`AN.ring.idx_order`/`pAN->ring.idx_order` 等。
- **`nvThresholds`→`ring.nvThresholdsScaled` 改名（11-2-2）**: スケール済（×`HSS_GRD_FAC_TEMP`）の SIMD 用閾値。ConfigCore 側の生℃ `threshold.nvThresholds` との同名混同を解消。
- `pAnalTemp`/`nLandCount`/`pLandToFull`/`pLandIdxMap`/`nFgFlow` はトップレベル据え置き。波及は本体7ファイル（amt/AnalCore/Interpolation/Output/OpticalFlow/AMTCore.h/WriteGRD.h）。`test/recalc.cpp`（mk.bat 非対象）は未更新。

### 重複共通化・整理（純粋リファクタ、出力不変）
- **`MergeHistU32AVX512(dst,src,nBins)`**（Interpolation.cpp）: スレッドローカル uint32 ヒストグラムの AVX-512 集約マージ（16幅で総和0ブロックskip＋スカラー末端）を `InterpolationMain`/`CalcBoundaryMaskAndCollect`/`CalcBoundaryMaskAndCollectMulti` の計6箇所から集約。※5fine→1coarse ダウンサンプルは3変種でターゲットが異なるため非対象（見送り継続）。
- **`BuildAltWeightTable` を `static` 化**（Interpolation.cpp、file-local 関数の linkage 統一）。
- **`nRunMode` マジックナンバー→ 名前付き enum 定数**（Output.cpp）: `GetRunStemName`/`OutputAnalysisSummaryCore` の `==1`→`RUN_MODE_MULTI_YEAR`、`InitHourlyStatNC`/`InitAllTextOut` の `==4..7`→`RUN_MODE_VALUE_SURVEY` 他、`InitOutputStationTemperature` の `>1`→`>RUN_MODE_MULTI_YEAR`。
- **未使用 `nPix` 削除**（`InterpolationMain`/`CalcBoundaryMaskAndCollect`）。
- ※ 重み再計算一本化 `EnsureWeightTableForStep`、多段ステップ引数集約 `RunStepContext`、`PrepareParamScanMonthData`、`WriteZoneStatRows` は先行リファクタ（〜VER.0.4.6.472）で実装済み。

### 不具合・防御・丸め
- **6-1**: `RunMultiYear` で月初ロード失敗時に `FlushMonth`/DivergRain 出力を行わず `continue`（`RunContinuous` の `break` と整合、コメント「Skipping to next month」と挙動一致）。失敗月のゼロ埋め確率テキスト・空 ProbImage を抑止（正常系は不変）。
- **7-1**: `LoadAMeDAS` データ本体ループの防御強化。観測点番号が index 未掲載のレコードを skip（`nIdxST=-1` 表現で station0 誤上書き防止）、`DAY.Day` 範囲外を skip（範囲外書込防止）。WARNING は各1回に抑制。正規データでは出力不変。
- **7-4**: 気温統計の**中央値**を avg/min/max と同じ四捨五入 `HsRoundNoNegZero(fMed, num_dec+1)` に統一（`OutputTempStat{Hour,Day,Mth,All}Core` 計4箇所）。変換スクリプト `tools/cnv_tstat_nc2tsv.py` の `fmt_median` も `hs_round(num_dec+1)` 経由へ揃え、本体TXT＝NC変換TSV の一致を維持（HOK 1999 dec1/dec2・hourly/monthly で byte 一致を確認）。`hs_round` は C++ `HsRound` を float32・half-away-from-zero で正確再現するため .X5 半端値でも完全一致。NC には従来どおり素の fMed を格納し読取側で丸める（avg/min/max と同様）。

---

## 12-12. リファクタリング6（VER.0.4.6.480〜488、2026-06-07）

`system_review_V0464_20260607.md` のレビューに基づく整理。**出力不変のリファクタ／意図的な挙動変更／既存不具合の修正**を含む。検証は新設の回帰ハーネス（下記）で行う方針（ビルド・実行確認はユーザー側責務）。

### 回帰ハーネス新設（`dbg_tools/`）
- リファクタ前後（または別モード）の**結果ディレクトリを総当たり比較**する開発者ツール群を新設。
  - `regress_dir.py`（入口）＋ `regress.bat`：`regress.bat <RESULT_DIR> <PREFIX_before> <PREFIX_after>`。**両PREFIXの積集合**のみ比較し、片側のみは `ONLY_BEFORE/ONLY_AFTER`（エラー扱いしない）。事前検証で PERIOD複数/ MODE・PERIOD不一致はエラー終了（`analysis_summary.txt` から取得、実行cfgは非参照）。
  - 拡張子別: **TXT/TSV=SHA256**（`analysis_summary`/`log`除外）、**NC=データレイヤ値**（`compare_nc_values.py`）、**GeoTIFF/PNG=画素**（`compare_geotiff.py`）、**GPKG=`sqlite3`でフィーチャテーブル**（`compare_gpkg.py`）、**GRD(.tva)=スキップ**（Blosc非決定）。
  - 依存は numpy/netCDF4/rasterio（既存）＋標準 `sqlite3`/`hashlib`（追加インストールなし）。

### 画像出力ファイル名の FILE_KEY 先頭統一（出力名変更）
- PNG/prob/bound 系（旧「種別先頭」`temp_KEY…`等）を **`{FILE_KEY}_{種別}_{日時}`（KEY先頭）**へ統一。GeoTIFF/GPKG/GRD は元からKEY先頭で不変。KEYと日時の無区切り連結（`{KEY}{日時}`）も `_` 挿入で解消。→ §11-4 の出力名表に反映済。回帰ハーネスの先頭プレフィックス探索が全画像に効くようになる。

### スペル統一 `Temparature`/`Temprature` → `Temperature`（識別子＝出力不変／出力名＝変更）
- コア識別子（`TemperatureInterpolation`/`…Multi`/`TemperatureGeoTIFF`/`TemperatureSetInt16`/`PROC_TM_Temperature*`/`InitOutputStationTemperature`）を一括統一。出力面では GeoTIFFサブDir `TemperatureGeoTIFF`・ファイル名 `_Temperature_`・TMD処理時間ラベルも変更。`test/recalc.cpp`（mk.bat非対象・公開対象外）は旧綴り残置。

### 補間設定の `InitAnal` 集約（出力不変）
- `InitAnal` が `CFG.interp.*`（`idw_power`/`rbf_sigma`/`barnes_kappa_fac`/`barnes_gamma`/`coast_L`）を直接設定。main 側の後付け代入5行を廃止。`InitAnal` の唯一の呼出元は main、かつ `InitConfigCore` デフォルト＝旧ハードコード値一致のため出力不変。`barnes_kappa`（`CalcInterpParams` 上書き）/`tps_lambda`（1.0f）は据置。

### 境界マスク生成の統一（構造リファクタ・出力不変）
- `CalcBoundaryMaskAndCollect(..., bool bGenMask=true)` を追加。`bGenMask=false` で**集計のみ・マスク非生成**。`TemperatureInterpolation` の Barnes2Pass/TPS 逐次分岐を `bGenMask=false` とし、**境界マスクは外側 `CalcBoundaryMaskOnly` で全手法一律生成**へ統一（逐次 B2P/TPS の二重マスク計算を解消）。並列パス（`CalcBoundaryMaskAndCollectMulti`／フォールバックの既定 `bGenMask=true`）は不変。

### `DataAnal` の static 化（出力不変）
- `main` の `DataAnal ANAL;` → `static DataAnal ANAL;`（~190KB をスタック→.bss）。直後の `memset` で cfg ループ毎に0クリア＝各cfg独立・出力不変。将来の `OBTAIN_AREA`/`TEMP_HIST_BINS` 拡張時のスタック逼迫を回避。

### 非フローGRD出力位相の温度位相化（**意図的な挙動変更**）
- `RunOneMonthOutputNoFlowStep` の GRD 条件 `(step-1)%N` → **`step%N`**。非フローGRD＝温度データのため**温度位相(step 0,N,2N…)**に統一し、同経路の気温GeoTIFF(`step%N`)と出力時刻が一致。**N=1は不変／N>1で出力step集合が `{1,N+1,…}`→`{0,N,2N,…}` に変化**。フロー時GRDは `(step-1)%N`（t0基準）を維持。

### HS magnitude 外周1pxのNaN統一（不具合修正）
- `CalcHornSchunckFlow` は `pMagnitude` を Ix 作業バッファに流用するため、`InitFlow` の事前NaN初期化が外周でも消えていた。magnitude復元（内部 [1,H-2]×[1,nW-2]）後に**外周1px（x=0/nW-1, y=0/H-1）を NAN 充填**し、U/V/Diverg/Vorticity と同じ外周NaN規約に統一。OPT_FLOW=2(HS) のフロー面 magnitude の外周1pxのみ変化（内部不変）。

### int16 PNG × USE_SIMD クラッシュ修正（既存バグ）
- `SavePNG8wSIMD` は float ソース(`fSrcData`)専用だったが、温度PNG(`OUTPUT_PNG.TEMP=1`)は int16 ソース(`CM_SRC_INT16_SCALED`/`pShtSrcData`)を渡す。VER.0.4.6.170 の int16 型追加時に woSIMD だけ `nSrcType` 分岐対応し SIMD版が取り残され、int16時に未設定 `fSrcData`(NULL) を読み**クラッシュ（0バイトPNG）**していた（通常運用が TEMP=2＝別経路 `WriteMergedPNG` のため潜在）。**SIMD版に int16 経路を正規実装**：`_mm512_cvtepi16_epi32`→`cvtepi32_ps`×`fShtFacTemp` で float[℃]へ1行 SIMD 展開(`row_f`)し、欠損 `INT16_MIN` は `fInv` へ写像、以降は既存 float パイプラインを共用。

### `SavePNG8` の汎用関数化（リファクタ・出力不変）
- `SavePNG8*Core(OutputContext& CTX)` ＋ OCX注入ラッパを廃止し、`SavePNG8`/`SavePNG8wSIMD`/`SavePNG8woSIMD` を **`CHsTextRenderer& TR` を明示引数で取る汎用関数**へ一本化（唯一の OutputContext 依存 `CTX.TR` を除去）。`SavePNG16` と同格のユーティリティ化＝将来のヘッダ外出しが容易に。`AMTCore.h` に `class CHsTextRenderer;` 前方宣言を追加、宣言から `pIL=NULL` デフォルトを廃止（省略呼出0件）。

### 検証・コミット
- 出力不変項目（interp集約／境界マスク統一／DataAnal static化／SavePNG8汎用化／スペル識別子）は byte/値一致で確認。出力変化項目（画像名統一・スペル出力名・非フローGRD位相・HS magnitude外周・int16 PNG修正）は変化点を明示記録。VER.0.4.6.488 としてコミット済。
- 公開準備（公開用最小テストの同梱／ビルド警告棚卸し／README_EN）は同梱データ不可・ROI等の理由で本リファクタでは見送り（`system_review_V0464_20260607.md` §10-2 参照）。

---

## 12-13. PNG出力本体のクラス化・外出し（VER.0.4.6.491〜494、2026-06-08）

`SavePNG*` 一群を **AMTプロジェクト非依存の汎用PNG出力クラス**へ外出しし、`Output.cpp` 内に散在していた libpng 定型コードの重複を排除した。

### 新設 `HsOutputPNG.h` / `HsOutputPNG.cpp`（コア・AMT非依存）
- **基底クラス `CHsOutputPNG`**：libpngライフサイクル（fopen/create/setjmp/IHDR/PLTE/write_info/filter/write_row/end/destroy）を集約し、`ConvMatrix.h` 以外のAMT固有型（`CHsImageLegend`/`CHsTextRenderer`/`ConfigCore`等）に一切依存しない。
  - `Save8` / `Save16`：`ConvMatrixParam`（float/int16スケール済ソース、`nUseSIMD`でAVX-512/スカラー）から8bit/16bit PNGを出力。旧 `SavePNG8/16(w/woSIMD)` の行変換を内部1本へ集約。
  - `SaveByBuilder8`：仮想 `buildRow8(row,y)` で行を都度生成（1行メモリfootprint維持）。横2倍マージPNG用。
  - `SaveIndexed8`：事前生成済みインデックスバッファから出力。等温線PNG用。
  - 旧 `PNG_MODE` マクロはコンストラクタ引数化。`HsClip` 依存はコア内 `hs_clipf` で内包。
- **インポーズの仮想化**：文字・凡例の行オーバーレイは仮想 `overlayRow(row,y,width)` に分離（基底は何もしない）。画像出力はユーティリティのためC++クラス化を許容（§0-2準拠）。
- **setjmp/longjmp対策**：各Saveメソッド内で fp/png/info をローカル保持し手動解放。RAII・例外巻き戻しに後始末を依存させない。

### `Output.cpp` 側（AMT従属を閉じ込め）
- **派生クラス `COutputPNG`**（ファイルローカル）：`CHsTextRenderer*` ＋ 凡例最大2枚を保持し `overlayRow` で文字＋凡例（`BlitLegendRow`）を適用。マージPNGは `SetMergedParams()`＋`buildRow8` override で横2倍行（左:気温0-236/右:等温線237-252）を生成。
- 置換：`OutputProbImage`/`OutputTempDistPNG`→`Save8`、`OutputBoundContourPNG`→`SaveIndexed8`、`OutputMergedTempContourPNG`（旧 `WriteMergedPNG` ラムダ廃止）→`SaveByBuilder8`。
- 旧 `SavePNG8/16(w/woSIMD・ディスパッチャ)` 計6自由関数（約640行）と `AMTCore.h` の同宣言・`class CHsTextRenderer;` 前方宣言を削除。`mk.bat` に `HsOutputPNG.cpp` を追加。
- 既存呼び出しの後方互換は維持しない（全面クラス化）。`test/test1.cpp`・`test2.cpp`（mk.bat非対象）は旧 `SavePNG16*` 呼び出しを残置（公開対象外）。

### 意図的な挙動変更（出力変化）
- **等温線PNG（`OutputBoundContourPNG`）にIL3凡例が出現**：旧実装は凡例memcpyの書込先が `&pBuf[nImgPosX]`（行0絶対位置）かつソース行無オフセットで、行出力 `&pBuf[y*nWA]` に反映されない**実質no-opの潜在バグ**だった。統一により `overlayRow` で正しく per-row 描画されるため、等温線PNGに等温線凡例が表示されるようになる（温度/確率/マージPNGの凡例実装と整合）。
- **凡例X位置の画像幅自動フィット（`BlitLegendRow`、VER.0.4.6.491A）**：`LEGEND_TC` の `IMPOSE_X` はマージPNG（横2倍幅）と単独ContourPNG（横1倍幅）で共有のため、マージ向けに大きい値（例 3500）を設定すると単独ContourPNG（例:北海道250m=幅2177）で画面外になり凡例が出ない。`IMPOSE_X + 凡例幅 > 画像幅` の場合は**右寄せクランプ**し、狭い画像でも必ず表示する。温度/確率/マージPNGは凡例が画像内に収まるためクランプ不発＝出力不変。
- **等温線凡例のパレットインデックス基点修正（`CHsImageLegend::Init`、VER.0.4.6.491A）**：等温線凡例(IL3)のラベル色は出力先パレットの閾値色インデックスに一致させる必要がある。単独ContourPNGは `g_IsothermPAL` の閾値色=**0～15**、マージPNGは `g_MergedPAL` の閾値色=**237～252**。旧 `nOffIdx = (mode==3 ? 0 : 237)` は単独(mode4)も237を使い、`g_IsothermPAL[237..]`（未初期化≒黒）を参照して**色化け**していた（凡例no-opバグで従来は露呈せず）。`nOffIdx = (mode==5 ? 237 : 0)` に修正し、単独ContourPNGは0基点でcontour線と同色に、マージ(mode5)は237基点で従来どおり。
- 温度/確率/マージPNGの可視出力は不変（emitted byte一致。マージは横2倍・左右生成・文字/凡例とも旧 `WriteMergedPNG` と同一）。
- **VER.0.4.6.494 としてコミット済**（amt.cpp main内のバージョン文字列も `VER.0.4.6.494`）。

---

## 12-14. レビューP1/P2/P3・残課題修正（VER.0.4.6.495〜500、2026/06/11〜13）

`system_review_V0464_20260611.md` のP1（バグ修正6件）・§6-11（PARAM_SCAN＋OPT_FLOW=0クラッシュ）・P2（防御強化5件）・P3（出力不変リファクタ8件）を実施した。詳細と検証結果は同レビュー文書の「対応履歴」を参照。

- **P1（VER.0.4.6.495）**: LAPSE.LAT読込の未初期化バッファ修正／SEARCH_Rの`[POI_CENTROID]`読込（旧`[POI]`互換維持）／`sprintf_s`フォーマット文字列誤用修正／`szUtf8Name`初期化／`CalcLaplacianSmoothness`の参照判定を手法ベース化／`PrepareParamScanMonthData`に`ChgInterpMethod`追加
- **§6-11（VER.0.4.6.496）**: `ApplyRunModeFixedConfig`にPARAM_SCAN時の`OPT_FLOW<1`→1（勾配法）矯正＋WARNINGを追加
- **P2（VER.0.4.6.497）**:
  - 欠測ピクセル書込の統一: `TempToI16()`ヘルパ新設。並列パス（`InterpolationMainMulti`）・Barnes2Pass・TPSの欠測時に`INT16_MIN`を明示書込し、リングバッファの前時刻残留値を排除（逐次パスと出力一致）。TPS係数計算失敗時は現スロットを全面欠測化
  - `ChgInterpMethod`はバッファ確保失敗時にWARNINGを出して`INTERP_BARNES`へフォールバック。`EnsureAnalTempBuffer`は確保時NaN初期化
  - メモリ確保失敗チェック追加（`CalcLaplacianSmoothness`/`CalcTPSCoeff`/`RunParamScan`スナップショット）
  - **`AREA.NUM=0`でも全域（`AREA_DEF_MAX`）の閾値統計を出力**: `FlushMonthAreaThreshold`の早期return削除、`OutputThrStatMonth`の実行判定を全域スロットへ、`ResetAnalBuffers`無条件化（nAreas>0の既存出力は不変、全域値は32エリア版と完全一致を検証済み）
  - **cfg読込防御**: `CHsIniFile`にアクセス記録＋`GetUnusedKeys()`を追加し、`LoadConfig`末尾で未参照キーをWARNING表示（綴り誤り・廃止キー・NUM不整合の検出）。`LoadConfigCore`/`LoadConfigRun`は`CHsIniFile&`受けへ変更（単一インスタンス共有）。LAPSE/AREA_CONFIGのSEASON・HOUR、RUNのCONT_START/CONT_DUR/MULTI_YEARに要素数不足WARNINGを追加
- **P3-1（VER.0.4.6.498）— 集計カーネル共通化（出力不変リファクタ）**:
  - `Interpolation.cpp` に共通集計カーネル `CollectPixelStats`（1陸地ピクセル分のヒストグラム・気温統計・面積カウント・確率累積）＋ループ不変コンテキスト `CollectCtx` を新設し、`InterpolationMain`（逐次）／`CalcBoundaryMaskAndCollect`（分離）／`CalcBoundaryMaskAndCollectMulti`（並列）の3重実装を一本化（正味約160行削減）
  - スレッドローカル統計は `TlStat` 構造体（旧Multi実装の関数ローカルをファイルスコープへ昇格）＋`InitTlStat`/`MergeTlStatToArrays` に統一。`AddAreaCntBelow` は全経路共用
  - 旧 `InterpolationMain` の確率累積のみ欠測（INT16_MIN）を除外していなかった不整合を解消（逐次/並列の確率カウント挙動が一致）
  - 検証: 12構成の回帰比較でデータ出力完全不変、並列（PS8）/逐次（PS1）の経路間一致も確認、性能劣化なし
- **P3-2〜P3-8（VER.0.4.6.499）— 共通化・分離・デッドコード削除（一括実施）**:
  - **P3-2 出力ヘルパ共通化（Output.cpp）**: 0.5℃BINヒストグラムからの中央値導出を `CalcMedianFromHist`、エリア出力順（全域→各エリア）を `BuildAreaOutputOrder` に集約し、temp_stat（時別/日別/月別/全期間）＋temp_hist の計5箇所の重複を排除
  - **P3-3 相関計算・標高帯定数の共通化**: `CorrAccum`/`CorrAccumAdd`/`CalcPearsonCorr`（AMTCore.h、inline）で `EvalDivergRain` と `RunParamScan` Stage2 の Pearson相関2重実装を一本化（n≤2・分散非正で0.0の規約は旧両実装と同一）。標高帯境界 200m/600m を `ELEV_ZONE_LOW_M`/`ELEV_ZONE_MID_M`（Constants.h）へ集約し、LOOCVゾーン分類（`ZoneClassifier::ClassifyElv`）と `RunAreaChar` で共用
  - **P3-4 `LoadConfigCore` 責務分離（AnalCore.cpp）**: 末尾のパレット生成（g_IsothermPAL/g_MergedPAL）・凡例初期化・FileKey検証を `FinalizeConfigCore`（static）へ分離し、読込関数を「読込・検証」に純化
  - **P3-5 デッドコード削除**: `TemperatureSetInt16` を完全削除（宣言・実装・呼出コメント残骸。int16化は各補間関数が `TempToI16` で直接実施）。TMDチャンネル3は欠番化（再利用しない）。恒偽 `&& false` 5箇所（3点中心差分フローの無効化）を `ENABLE_3POINT_FLOW`（amt.cpp冒頭、=0）マクロへ集約
  - **P3-6 型・const整備**: `CalcGradientFlow` SIMDループの添字5変数を size_t 化（OpticalFlow.cpp）、`GetDateTime` を `const DateTime&` 受けへ
  - **P3-7 ログ・出力品質**: 日別統計/日別ヒストは全域カウント0の日（CONTINUOUS途中開始の前日など）の空行出力を抑止。`RunContinuous` 失敗ログの年月を dtPrev→dt に修正（月跨ぎ直後の前月誤表示）＋閉括弧欠落2箇所修正。コメント修正群（nThrHiOut「T<」→「T≥」、g_IsothermPAL「251〜252」→「254」、BoundaryMask「FAC=100」→「500」）
  - **P3-8 引数構造整理**: `EvalDivergRain` を閾値2個ばら渡しから `const ConfigRain&` 受けへ。`CalcKrigingSFallbackValue` は pAltWF/pAltW/coast_L を AN から導出し10→7引数に削減
  - 既知差分: ①P3-7により途中開始時の先頭空行（temp_stat_daily/temp_hist_daily）が消える（正しい変化。14時開始＋日別出力有効の追加走行で動作確認済み）②RunParamScan Stage2のCorrログ表示はfloatキャスト除去で末桁が変わる可能性（ファイル出力なし）。それ以外は回帰比較（CONTINUOUS/MULTI_YEAR/AREA_CHAR）で出力不変を確認済み
- **軽微残課題8件（VER.0.4.6.500）— 優先順位一覧外の小規模整理（出力不変）**:
  - 閾値マジックナンバー16を `THR_DEF_MAX` に統一（`HOURLY_COLS_MAX`・`aAreaCntBuf`・各集計バッファ添字 約30箇所）、`szPath[2024]`→2048
  - `nOutLog_KetTime`→`nOutLog_KeyTime` 改名（cfgキー `KEY_TIME` は不変）、`AMeDAS_DAY.TimeData` の添字規約（TimeData[h]=解析時刻h）と `DataFlow` の作業面流用規約（Horn-Schunck計算途中の pMagnitude/pDiverg/pVorticity 参照禁止）を構造体コメントに明記
  - `CHsImageLegend::Init` の閾値数チェック `<0`→`<1`（機能していなかった判定の修正）、`LoadAMeDAS` 観測所マスタ読込の `fread` 戻り値チェック追加、`CalcProbMonthStat` のデフォルト引数 `dir=0` 廃止
- **LoadConfigCore 分割＋INI読込定型化（VER.0.4.6.500、レビュー2-2-1／8-4）**:
  - `LoadConfigCore`（約708行）をセクション群毎の static 7関数へ分割し、本体は順次呼出の約17行に:
    `LoadConfigBase`（[MAIN][PATH][RUN]）→`LoadConfigAnal`（[EXT_ANAL][INTERP][PARAM]）→`LoadConfigStatOut`（[TEMP_STAT][TEMP_HIST][THRESHOLD_ANAL]）→`LoadConfigFileOut`（[OUTPUT]系5セクション＋RenderConfig）→`LoadConfigThresholdArea`（[THRESHOLD][AREA][POI_CENTROID]）→`LoadConfigLapse`（[LAPSE]）→`LoadConfigAreaCorr`（[AREA_CONFIG]）→`FinalizeConfigCore`
  - 呼出順は旧一体実装の読込順を完全維持（AREA_CONFIG は LAPSE/INTERP の結果に、INTERP/LAPSE は PARALLEL_SLICES に依存）。`TERRAIN_DEM_PATH` の設定有無は `LoadConfigBase` の out 引数→`LoadConfigAreaCorr` の in 引数で明示受け渡し
  - INI読込定型化（8-4）: `GetIniString`（空初期化を保証する GetString ラッパ。§6-1 LAPSE.LAT 型の未初期化バッファ事故を構造的に根絶）＋`SplitValueList`（「/」「,」両対応の数値リスト分割）を新設し、szTmp定型15箇所（Core系12＋Run系3）を置換。szTmp の関数跨ぎ使い回しも解消
  - 既存cfg（/区切り）では出力完全不変。挙動拡大1点: LAPSE.LAT/HOUR/VALUE系・CONT_START 等で「,」区切りも受理されるようになる（従来はWARNINGまたは無視）
- **InitAnal 分割＋確保定型化＋構造整理（VER.0.4.6.500、レビュー2-2-1／8-9／1-2-2／3-2-1）**:
  - `InitAnal`（約640行）を機能ブロック毎の static 8関数へ分割し、本体は順次呼出の約20行に:
    `InitAnalDefaults`（探索/補間/lapse初期値）→`InitAnalGrid`（陸地圧縮マップ・リングバッファ・pAnalTemp・pLandToFull・地形特徴量）→`InitAnalNbrTable`（近傍SoA＋ALTW_CACHE）→`InitAnalAreaStats`（エリアマスク・集計・ヒストグラム・TLプール）→`InitAnalMultiPool`→`InitAnalProb`（確率マップ）→`InitAnalThreshold`（閾値スケール＋mBoundPacked）→`InitAnalPoi`＋`PrintInitAnalMemSummary`。確保順は旧実装と完全同一
  - 確保定型化（8-9）: `ANAL_ALLOC` マクロ（確保＋NULLチェック＋変数名表示＋-1返却）で単独確保13箇所を置換。確保失敗時の `FreeAnal` は InitAnal 本体での一括実施に一元化（サブ関数は -1 を返すのみ。チェック忘れと解放漏れを構造的に防止）
  - 閾値同期の明示（1-2-2）: `nvThresholds`（生℃）→`nvThresholdsScaled`（int16スケール）の変換を `SyncThresholdsScaled` に分離し「唯一の同期点」として規約をコメント明記
  - Config層別の非対称解消（3-2-1）: `SetupForMonth`/`InitOutputAfterAreaMask` を `AppConfig` 一括受けから `(ConfigCore, ConfigRun)` 分離受けへ変更
  - 出力完全不変（確保失敗時のログ表示が変数名ベース「AN.xxx」表記に変わるのみ）
- **analysis_summary可読化＋ST_NBR_EVAL仕様明文化（VER.0.4.6.500、レビュー7-11／6-10）**:
  - analysis_summary の `MODE` 行を「番号+名前」併記へ（例 `MODE	0	CONTINUOUS`）。`GetRunModeName`（モード番号→RUN.MODE名の逆引き）を新設
  - 回帰ハーネス `regress_dir.py` の MODE ゲートは番号部のみの比較に変更し、新旧フォーマット間の回帰比較互換を確保
  - ST_NBR_EVAL の不連続 HOURS 指定（例 0/23）で min～max 全時刻を計算する挙動を仕様として明文化（§5 RunStNbrEval 節参照）
- **ユーティリティ外出し（VER.0.4.6.500）**: 旧 `PathExistCheck`（AnalCore.cpp）を AMT非依存の `HsPathExistCheck` として HsCommon.h の inline 関数へ移設（戻り値int化。§9-2参照）

### 12-14-1. 修正保留事項と理由

`system_review_V0464_20260611.md` の指摘のうち、以下は意図的に**見送り**とした（再検討時はレビュー文書の各節を参照）。

| 保留事項 | レビュー | 理由 |
|---|---|---|
| `RunOptimParams`（~360行）/`RunEvalInterp`（~240行）/`PrintCFG`（~210行）の関数分割 | 2-2-1残り | いずれも変更頻度の低い評価系または単純列挙であり、分割の費用（回帰検証コスト）が保守効果を上回る。`RunOptimParams` は Phase1〜6 の区切りコメントで内部構造が既に追跡可能 |
| `OutputContext` Core+ラッパ二重化（約30組）の解消 | 2-2-4 | 既存の対は実害がなく、解消の書き換え量が大きい。**新規出力関数では Core+ラッパの二重化はせず単一関数とする**（運用指針として本書に記録） |
| `CalcProbMonthStat` 等の `int mode`/`int dir` の enum 型化 | 3-2-5残り | 計算コアのC互換維持方針（§0-2）と衝突するため。呼び忘れを隠すデフォルト引数 `dir=0` の廃止までで対応済み |
| ST_NBR_EVAL 不連続 HOURS の中間時刻計算スキップ実装 | 6-10 | `RunOneMonth` の連続ステップ実行前提を崩す改修が必要で、評価用途の利用頻度に対して費用対効果が低い。仕様の明文化（関数コメント＋§5）で対応 |
| `AreaSlotAll()` 等の格納添字側ヘルパ | 1-2-1残り | 出力順側は `BuildAreaOutputOrder` で集約済み。格納添字 `AREA_DEF_MAX` の直接参照は意味が自明な箇所のみが残存しており、ヘルパ化の効果が薄い |
| `WriteZoneStatRows`（9引数）/`IsMultiBatchBoundary`（7引数）の引数削減 | 3-2-3残り | 呼出箇所が各1〜2箇所に限られ、構造体化のコストに見合わない |

---

*以上、VER.0.4.6.494（2026/06/11時点）＋ §12-14（VER.0.4.6.495〜500、2026/06/12〜13）システム全体解説。*

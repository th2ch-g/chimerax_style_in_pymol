# PyMOLでのChimeraX風表示

[README](../../README.md) | [全ギャラリー](../gallery.md) | [English](../pymol_chimerax.md)

`chimerax_style` は、84種類の表示・プリセット・照明をPyMOL 3.1に追加します。
ChimeraXや隣接リポジトリの実行環境には依存しません。原子と分子表面には
PyMOL標準表示を、その他には独自に生成したCGO形状を使います。

## 基本操作

```text
chimerax_style cartoon, selection=polymer, color=chain
chimerax_style soft, representation=surface, selection=chain A
chimerax_style nucleotides-tube-slab, selection=polymer.nucleic, color=nucleotide
chimerax_style ball, selection=organic, name=ligand, lighting=keep
chimerax_style volume-mesh, selection=density, params=contour.json
chimerax_style refresh, name=all
chimerax_style reset, name=all
chimerax_style ray, filename=figure.png, width=2400, height=1800
chimerax_style list
```

`selection` はPyMOLの選択式です。ChimeraXのモデル番号・選択構文は使いません。
`name` が同じ表示は置換され、別名なら重ねられます。`quality` は
`low` / `medium` / `high`、`transparency` は0〜1、`state=0` は分子の全stateを
準備します。正のstateを指定すると、そのstateだけを表示します。マップは1stateです。
`lighting=keep` で現在の照明を保ちます。通常は背景を変えません。

## 対応範囲

| 分類 | 表示 |
| --- | --- |
| 原子 | stick、ball、sphere、厚い/薄いring fill |
| cartoon | 標準・楕円・矩形・barbell断面、曲がった/直線のhelix cylinder、wrap、tube、worm |
| 分子表面 | solvent-excluded surfaceのsolid/mesh/dot、Gaussian envelope |
| 核酸 | atoms、fill、slab、tube/slab、box/muffler/ellipsoid、ladder、stubs |
| 特殊表示 | SNFGの8形状クラス、異方性楕円体・主軸・主楕円 |
| マップ | surface/mesh/dot、透過image、plane、orthoplanes、box faces、傾斜slab、最大値投影、topography、整数label境界 |
| 注記 | pseudobonds、距離・角度・二面角、hbonds、contacts、struts、ラベル、markers/links |
| 幾何形状 | axis、plane、centroid、inertia、unitcell、球・楕円体・円柱・円錐・箱・矩形・三角形・正多面体・tube・ribbon・mesh |
| プリセット | protein/nucleotideの組合せ、space filling、色/透明度別surface、publication、interactive |
| 照明 | simple、full、soft、gentle、flat |

全コマンド名と分類ごとの再現範囲は[英語版の対応表](../pymol_chimerax.md#representation-coverage)
にあります。ChimeraXの組込み表示形式を対象とし、Toolshed拡張、解析処理全般、GUI、
VR、医療画像viewer、全ファイル形式の読込までは含みません。

stick半径0.2 Å、ball倍率0.3、cartoon幅2 Å・厚さ0.4 Å、核酸slab厚さ0.5 Åなどを
参照しています。一方、cartoon補間・末端接続、表面三角形分割、VDW半径、照明・
陰影はChimeraXと一致しません。実際のChimeraX画像との比較検証はしていません。

`volume-image` は透過平面を重ねる近似で、ChimeraXのvolume rendererとは異なります。
回転方向によって平面の継ぎ目が見える場合があり、`axis` を視線に合わせて変更できます。
`publication` の輪郭線はPyMOLのray効果です。GPU画面には同じ輪郭処理はありません。
`full` / `soft` のray陰影とGPU陰影にも違いがあります。

## パラメータと科学データ

`params` / `data` はPython辞書またはローカルファイルです。PyMOLのコマンド行では
カンマが引数区切りになるため、辞書の代わりにJSONファイルを指定してください。

```python
from chimerax_style_in_pymol import chimerax_style

chimerax_style("cartoon-barbell", params={"width": 2.0, "thickness": 0.4})
chimerax_style("worm", params={"worm_attribute": "b", "worm_min": 0.25, "worm_max": 2.0})
chimerax_style("aniso", params={"probability": 0.5})
chimerax_style("nucleotides-ladder", data={"pairs": [[0, 7], [1, 6]]})
chimerax_style("snfg", data={"sugars": {"XYZ": ["sphere", "#00a651"]}})
```

核酸のpairは各分子object内の核酸残基順、0始まりです。pairを与えない残基はstubsになります。
距離だけから塩基対を推定しません。異方性表示にはANISOUが必須で、等方性B因子で代用しません。
`aniso` の初期半径は主軸方向のRMS変位です。`probability` を指定したときだけ
3自由度のカイ二乗分布で確率楕円体へ変換します。

`hbonds` はPyMOLのpolar-contact判定、`contacts` は距離閾値、`struts` は近傍の
backbone点同士を使います。ChimeraXの検出・最適化アルゴリズムの移植ではありません。
SNFGは一般的なCCD残基名を認識します。未登録名は`data.sugars`で形と色を指定できます。

マップはXYZ順の3次元配列です。NPZはpickleを無効にして読みます。`origin` と
`spacing`、またはvoxel indexをÅへ変換する4×4 `transform`を指定します。
MRC/CCP4と既存PyMOL map objectにも対応します。

```python
import numpy as np

np.savez("density.npz", values=values, origin=[0, 0, 0], spacing=[1, 1, 1])
chimerax_style("volume-surface", data="density.npz", params={"level": 1.2})
chimerax_style("volume-plane", data="density.npz", params={
    "axis": "z", "position": 0.5,
    "transfer": [[0, "#102040", 0], [1, "#4080ff", 0.1], [3, "#ffee80", 0.8]],
})
```

`level` は入力値そのものの閾値で、自動的にsigmaへ正規化しません。transferは
`[値, 色, opacity]` の昇順リストです。`position` は0〜1、`slab_normal` はCartesian
ベクトル、`slab_depth` はÅ単位の厚さとsample数です。`volume-segment`は整数label入力と
明示的な`segments`リストを使います。`volume_step`はimage平面の間隔（初期値0.2 Å）で、
変更時は透明度を光学的厚さに合わせて補正します。`max_voxels` / `max_triangles`で生成量を制限します。

```python
chimerax_style("markers", data={
    "points": [[0, 0, 0], [2, 1, 0]], "links": [[0, 1]],
})
chimerax_style("pseudobonds", data={
    "points": [[0, 0, 0], [3, 2, 1]], "pairs": [[0, 1]],
})
chimerax_style("shape-mesh", data={
    "points": [[0, 0, 0], [2, 0, 0], [0, 2, 0]], "faces": [[0, 1, 2]],
})
```

pointsのindexも0始まりです。一般形状は`center`、`size`、`radius`、`height`、色は
`params.color`で指定します。分子色は`color=auto/keep/model/chain/element/nucleotide/
secondary-structure/rainbow/bfactor`、PyMOL色名、または`#rrggbb`です。

## 元データ・session・出力

importだけではPyMOLを起動せず、設定も変えません。`__init_plugin__()`で登録します。
元構造の座標・結合・色・属性は保持し、管理中だけ選択した原子の元representationを
隠します。`reset`で復元します。選択範囲が重なる表示は、最後の所有者がなくなるまで
元表示を隠したままにします。失敗した置換は以前の表示を保持します。

座標や色を編集した後は`refresh`が必要です。PSEに表示objectと管理情報を保存でき、
再読込後はコマンドを再登録して操作できます。atom indexを変える編集後は表示を
再生成してください。groupを直接削除した場合の復元は次のコマンド実行時に行います。

標準原子表示をpickして`chimerax_style select, selection=pk1`を実行すると、対応する
元原子を`sele`へ選択します。独自CGO形状は原子単位のpickに対応しません。

照明はscene全体へ作用します。最後の管理表示をresetすると、コマンドが変更した設定を
復元し、その後の手動変更は保持します。`png`はOpenGL画面が必要で、`ray`はheadlessでも
使用できます。PNGデータを検証してから出力先を置き換えるため、描画に失敗しても
既存のファイルを保持します。ギャラリーは実際のPyMOL GPU/ray出力で、AI生成画像ではありません。

独自形状のsolid・mesh・dot用CGOデータはNumPyでまとめて作成し、大きなcartoon・
核酸形状・mapの準備負荷を減らします。三角形の順序、法線、色、透明度のグループは
維持し、すべてのstateを準備してから管理表示を置き換えます。

科学的な入力契約、全表示形式、復元、session再読込、PNG出力を実際のPyMOLで検証します。
参照元・ライセンスは[英語版](../pymol_chimerax.md#reference-and-verification)と
[NOTICE](../../NOTICE)に記載しています。

モデル色・鎖色・元素色・核酸色は本家 `8b1067a` の定義に合わせています。モデル色はtan、sky blue、plum、light greenの順で、鎖色は選択順によらず鎖IDで決まります。元素色109項目と未定義元素の灰色も一致させています。rainbowは鎖ごとの残基単位、B-factorは青・白・赤でcartoonは残基平均です。任意のsecondary-structure配色は本プラグイン独自の赤・金・青です。

helix-tube/cylinderは本家と同じ円弧・円筒の目的関数でフィットし、短いhelixは理想構造へ整列、13残基以上のtubeは円弧を使います。半径は自動フィット（上限2.5 Å）で、`params.helix_radius` で上書きできます。1CRN/1GGGの36条件で本家数値計算との軸位置・半径の最大差は0.000036/0.0000012 Åでした。解析勾配で評価回数を抑えています。末端・ループの接続と三角形分割には差が残ります。実ChimeraX画像との比較は未実施です。

NumPy 2.5.3で全84スタイルのGPU/rayを検証しました。2万原子の鎖色計算は73.6 msから2.2 ms、B-factor色計算は48.9 msから7.9 msへ短縮しました。

publicationの白背景設定はNumPyスカラーをPython floatへ変換してPyMOLに渡し、NumPy 2.xで設定が反映されない問題を修正しました。ヘリックスの管はループ補間から分離して末端の膨らみを防ぎます。

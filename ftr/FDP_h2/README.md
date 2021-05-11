# 水素分子の静的分極率計算



## 概要

- 制限HF(2種類のスピン軌道の空間軌道が同じと仮定)
- 基底関数としてガウス型基底関数を採用
- 水素分子の基底として6-31G基底を採用(https://bse.pnl.gov/bse/portal)
- 2電子相互作用を無視できる方向で平均場近似

- 分極率とは外部電場に対する双極子モーメントの応答関数

## アルゴリズム

1. 初期占有軌道を取りこれを基にフォック行列を計算する
2. ブルリアン条件を満たされているかどうか判定する(満たしていれば終了)
3. フォック行列の対角化
4. 固有化の低い順にN個の占有軌道を決定する
5. 占有軌道のセットに基いてフォック行列を計算(2に戻る)

## 使い方等

- 当然 gfortran が使える環境必要(最近はGCCにバンドルされてる)
- makeして ./FDP_h2 で実行
- エネルギーの単位はa.u.(Hartree) -> 1 a.u. = 27.2116 eV
- make clean でいらないファイル削除

- 電場によるz成分ポテンシャルを入力する

## 参考文献

- 実践量子化学計算プログラミング : 日野理, アドバンスソフト, 第3章
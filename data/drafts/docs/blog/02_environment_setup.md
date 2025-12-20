# 第2回：環境構築「ターミナル操作とconda環境の構築」

## はじめに

前回は、プロジェクトの背景と全体像について説明しました。今回は、実際に解析を始めるための**環境構築**について詳しく解説します。

バイオインフォマティクス解析では、様々なツールを使い分ける必要があります。condaを使えば、これらのツールを簡単に管理できます。

---

## パイプライン全体図

この記事は、以下のパイプラインの**ステップ2：環境構築**に対応しています：

```
[1] イントロダクション
[2] 環境構築 ← 現在ここ
[3] データ準備
[4] シミュレーション
[5] アライメント・カウント
[6] 発現解析
[7] 可視化
[8] CNV解析
[9] まとめ
```

<!-- TODO: 画像URLを設定してください: パイプライン全体図（ステップ2をハイライト） -->
<!-- ![パイプライン全体図](画像URL) -->

---

## 前回の振り返り

**第1回で学んだこと**：
- リオス科の性的二形の仮説（Hox遺伝子の発現ドメインシフト、CNVによる形態分化）
- プロジェクトの全体像と学習目標
- シミュレーションデータを使った学習の利点

**今回から始めること**：
- ターミナル操作の基本
- conda環境の構築
- ディレクトリ構成の設計

---

**この記事で学べること**：
- ターミナル（コマンドライン）の基本的な使い方
- condaとは何か、なぜ必要なのか
- conda環境の作成と管理方法
- ディレクトリ構成の設計方法

**前提知識**：
- ターミナル（コマンドライン）を初めて使う方は、まず「ターミナルの基本操作」セクションを読んでください
- PythonやRの基本的な知識があると理解しやすいですが、必須ではありません

---

## ターミナルの基本操作（初心者向け）

### ターミナルとは？

ターミナル（コマンドライン）は、**キーボードでコマンドを入力してコンピュータを操作する画面**です。Windowsでは「コマンドプロンプト」や「PowerShell」、Mac/Linuxでは「ターミナル」と呼ばれます。

**なぜターミナルを使うのか？**
- バイオインフォマティクスのツールの多くは、ターミナルから実行する必要がある
- 大量のファイルを一度に処理できる
- 自動化（スクリプト）が容易

**ターミナルの開き方**：
- **Windows**：`Win + R` → `cmd` または `powershell` と入力
- **Mac**：`Cmd + Space` → 「ターミナル」と検索
- **Linux**：`Ctrl + Alt + T` またはアプリケーション一覧から選択

**参考記事**：
- [【初心者向け】ターミナルの使い方とコマンド一覧](https://qiita.com/tags/terminal) - Qiitaのターミナル関連記事
- [Linuxコマンド入門](https://eng-entrance.com/linux-command) - エンジニアの入門サイト

### 基本的なコマンド

以下は、この記事で使う基本的なコマンドです。詳しく知りたい方は、以下の記事を参考にしてください：

| コマンド | 説明 | 例 |
|---------|------|-----|
| `pwd` | 現在のディレクトリ（フォルダ）を表示 | `pwd` → `/home/uecha` |
| `ls` | ファイル・フォルダの一覧を表示 | `ls` → ファイル一覧 |
| `cd` | ディレクトリを移動 | `cd /home/uecha/data` |
| `mkdir` | ディレクトリを作成 | `mkdir my_project` |
| `which` | コマンドの場所を確認 | `which python` → `/usr/bin/python` |

**実践例**：
```bash
# 現在の場所を確認
pwd

# ホームディレクトリに移動
cd ~

# 新しいフォルダを作成
mkdir test_project

# そのフォルダに移動
cd test_project

# 中身を確認（まだ空のはず）
ls
```

**コマンドの補完機能**：
- `Tab`キーを押すと、ファイル名やコマンド名を自動補完してくれます
- 例：`cd /ho` → `Tab` → `/home` に補完される

**参考記事**：
- [ターミナルの便利な機能（補完、履歴など）](https://qiita.com/tags/terminal) - Qiitaのターミナル記事

### エラーが出たときの対処法

**よくあるエラーと対処法**：

1. **`command not found`** エラー
   - **原因**：コマンドがインストールされていない、またはパスが通っていない
   - **対処法**：該当ツールをインストールする、または環境を有効化する
   - **確認方法**：`which コマンド名`でコマンドの場所を確認

2. **`Permission denied`** エラー
   - **原因**：ファイルやディレクトリへのアクセス権限がない
   - **対処法**：`chmod`コマンドで権限を変更する、または管理者権限で実行する
   - **確認方法**：`ls -l`でファイルの権限を確認

3. **`No such file or directory`** エラー
   - **原因**：指定したファイルやディレクトリが存在しない
   - **対処法**：`pwd`で現在地を確認し、正しいパスを指定する
   - **確認方法**：`ls`でファイル一覧を確認

**エラーメッセージの読み方**：
- エラーメッセージには、問題の原因が書かれていることが多い
- エラーメッセージをコピーして、QiitaやStack Overflowで検索すると解決策が見つかることが多い

**参考記事**：
- [ターミナルでよくあるエラーと対処法](https://qiita.com/tags/terminal) - Qiitaのエラー対処記事
- [Stack Overflow](https://stackoverflow.com/) - プログラミングの質問サイト

---

## conda環境の構築

### condaとは？

**conda**は、PythonやRのパッケージ、バイオインフォマティクスツールを管理する**環境管理システム**です。

**なぜcondaが必要なのか？**

1. **バージョン管理**：異なるプロジェクトで異なるバージョンのツールを使い分けられる
   - 例：プロジェクトAではR 4.2を使い、プロジェクトBではR 4.3を使う
2. **依存関係の解決**：ツール同士の依存関係を自動で解決してくれる
3. **バイオインフォマティクスツールの豊富さ**：biocondaチャンネルに多くのツールが用意されている

**condaの基本概念**：

- **環境（environment）**：ツールやパッケージをまとめた「箱」
- **チャンネル（channel）**：パッケージを取得する場所（例：conda-forge、bioconda）
- **パッケージ（package）**：インストールするツールやライブラリ

**参考記事**：
- [conda入門：環境管理の基本](https://qiita.com/tags/conda) - Qiitaのconda関連記事
- [conda公式ドキュメント](https://docs.conda.io/) - 公式ドキュメント（英語）

### condaのインストール

condaがまだインストールされていない場合は、以下の手順でインストールしてください。

**Miniconda（推奨）**：
- 軽量版のconda。必要最小限のパッケージのみ含まれる
- 容量が小さい（約400MB）
- [Miniconda公式サイト](https://docs.conda.io/en/latest/miniconda.html)からダウンロード

**Anaconda**：
- conda + 多くのパッケージが最初からインストールされる
- 容量が大きい（約3GB）が、すぐに使える
- [Anaconda公式サイト](https://www.anaconda.com/)からダウンロード

**インストール手順（Windows）**：
1. Minicondaのインストーラーをダウンロード
2. インストーラーを実行
3. 「Add Miniconda3 to PATH」にチェックを入れる（推奨）
4. インストール完了後、新しいターミナルを開く

**インストール手順（Mac/Linux）**：
```bash
# Minicondaのインストーラーをダウンロード
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# インストーラーを実行
bash Miniconda3-latest-Linux-x86_64.sh

# インストール後、ターミナルを再起動
```

**インストール後の確認**：
```bash
# condaがインストールされているか確認
conda --version
# 出力例: conda 23.7.4

# condaの情報を表示
conda info
```

**参考記事**：
- [condaのインストール方法（Windows/Mac/Linux）](https://qiita.com/tags/conda) - Qiitaのインストール記事
- [Miniconda公式インストールガイド](https://docs.conda.io/en/latest/miniconda.html) - 公式ドキュメント

### 環境の作成

このプロジェクトでは、2つのconda環境を使用します：

1. **`rios-r43`**: R 4.3 + Bioconductor 3.18（polyester用）
2. **`rios-env`**: その他のツール（HISAT2, samtools, DESeq2など）

**環境を作成する理由**：
- polyesterはR 4.3が必要だが、他のツールは異なるバージョンが必要な場合がある
- 環境を分けることで、ツール同士の競合を避けられる

**環境の作成コマンド**：

```bash
# R 4.3環境の作成（polyester用）
conda create -y -n rios-r43 -c conda-forge -c bioconda \
  r-base=4.3 bioconductor-polyester=1.38.0

# その他のツール用環境の作成
conda create -y -n rios-env -c conda-forge -c bioconda \
  r-base=4.3 \
  hisat2 \
  samtools \
  bioconductor-deseq2 \
  bioconductor-biostrings \
  r-ggplot2 \
  r-data.table \
  r-optparse
```

**コマンドの説明**：
- `conda create`：新しい環境を作成
- `-y`：確認メッセージをスキップ（自動で「yes」と答える）
- `-n rios-r43`：環境名を`rios-r43`に設定
- `-c conda-forge -c bioconda`：パッケージを取得するチャンネルを指定
- `r-base=4.3`：Rのバージョンを4.3に固定

**実行時間**：
- 初回は5-10分程度かかることがあります（パッケージのダウンロードに時間がかかる）
- ネットワーク速度によって異なります

**エラーが出た場合**：
- **ネットワークエラー**：インターネット接続を確認し、再実行
- **パッケージが見つからない**：チャンネルを追加する（`conda config --add channels bioconda`）
- **権限エラー**：管理者権限で実行する、またはユーザー権限でインストール

**参考記事**：
- [conda環境の作成と管理](https://qiita.com/tags/conda) - Qiitaの環境管理記事
- [bioconda公式サイト](https://bioconda.github.io/) - バイオインフォマティクスツールのチャンネル

### 環境の有効化

環境を作成したら、**有効化（activate）**する必要があります。有効化すると、その環境にインストールされたツールが使えるようになります。

```bash
# 環境を有効化
conda activate rios-env

# 環境が有効化されたか確認
which Rscript
# 出力例: /home/uecha/miniconda3/envs/rios-env/bin/Rscript

which hisat2
# 出力例: /home/uecha/miniconda3/envs/rios-env/bin/hisat2
```

**環境の切り替え**：
```bash
# rios-r43環境に切り替え
conda activate rios-r43

# 環境を無効化（base環境に戻る）
conda deactivate

# 再度rios-env環境に切り替え
conda activate rios-env
```

**環境の一覧表示**：
```bash
# 作成済みの環境一覧を表示
conda env list
# 出力例:
# base                  *  /home/uecha/miniconda3
# rios-r43                 /home/uecha/miniconda3/envs/rios-r43
# rios-env                 /home/uecha/miniconda3/envs/rios-env
```

**よくある質問**：
- **Q: 環境を有効化し忘れたら？**  
  A: `command not found`エラーが出ます。`conda activate`で環境を有効化してください。

- **Q: 環境を削除したい場合は？**  
  A: `conda env remove -n 環境名`で削除できます。

- **Q: 環境にパッケージを追加したい場合は？**  
  A: `conda install -n 環境名 パッケージ名`で追加できます。

**参考記事**：
- [conda環境の使い方（activate/deactivate）](https://qiita.com/tags/conda) - Qiitaの環境操作記事

---

## ディレクトリ構成の設計

### なぜディレクトリ構成が重要なのか？

バイオインフォマティクス解析では、大量のデータファイルを扱います。適切なディレクトリ構成により：

- ファイルの場所が分かりやすくなる
- 大容量データを効率的に管理できる
- 複数プロジェクトでデータを共有できる

**悪い例**：
```
/home/uecha/
├── data1.fa
├── data2.fa
├── script1.py
├── script2.R
└── result1.txt
```
→ ファイルが散らばっていて、どこに何があるか分からない

**良い例**：
```
/home/uecha/
├── data/
│   └── reference/
├── work/
│   └── project1/
│       ├── scripts/
│       └── analysis/
└── scratch/
    └── project1/
```
→ ファイルが整理されていて、目的が明確

### 推奨ディレクトリ構成

```
/home/uecha/
├── data/                    # 大容量共通データの原本
│   ├── reference/rios/      # ゲノムFASTA/GTF
│   └── databases/rios/      # インデックス（STAR/HISAT2）
├── work/MonsterGenome_Rios/ # プロジェクト作業ディレクトリ
│   ├── docs/                # 仕様書、計画
│   ├── scripts/             # ユーティリティスクリプト
│   ├── reference/           # シンボリックリンク（→ data/reference/rios/）
│   ├── simulation/          # シミュレーション入力
│   └── analysis/            # 解析スクリプトと出力
└── scratch/MonsterGenome_Rios/ # 一時ファイル（削除OK）
```

**各ディレクトリの役割**：
- **`data/`**：大容量データの原本。複数プロジェクトで共有
- **`work/`**：プロジェクトの作業ディレクトリ。解析スクリプトや結果を保存
- **`scratch/`**：一時ファイル。削除しても問題ない

**ディレクトリの作成方法**：

```bash
# ホームディレクトリに移動
cd ~

# ディレクトリを作成
mkdir -p data/reference/rios
mkdir -p data/databases/rios
mkdir -p work/MonsterGenome_Rios/{docs,scripts,reference,simulation,analysis}
mkdir -p scratch/MonsterGenome_Rios
```

**コマンドの説明**：
- `mkdir -p`：親ディレクトリも含めて作成（`-p`オプション）
- `{docs,scripts,...}`：複数のディレクトリを一度に作成

**確認方法**：
```bash
# 作成されたディレクトリを確認
tree -L 3 ~/work/MonsterGenome_Rios
# または
ls -R ~/work/MonsterGenome_Rios
```

**参考記事**：
- [Linuxディレクトリ操作の基本](https://eng-entrance.com/linux-command-directory) - ディレクトリ操作の解説
- [プロジェクト構成のベストプラクティス](https://qiita.com/tags/project-structure) - Qiitaのプロジェクト構成記事

<!-- TODO: 画像URLを設定してください: ディレクトリ構成図 -->
<!-- ![ディレクトリ構成図](画像URL) -->

### シンボリックリンクの活用

大容量データ（ゲノムFASTAなど）は、`/home/uecha/data/`に1回だけ保存し、プロジェクトからは**シンボリックリンク**で参照します。

**シンボリックリンクとは？**
- ファイルやディレクトリへの「ショートカット」のようなもの
- 実際のファイルは1つだけだが、複数の場所から参照できる
- Windowsの「ショートカット」、Macの「エイリアス」に似ている

**シンボリックリンクの作成**：

```bash
# ゲノムファイルへのリンク作成
ln -sf /home/uecha/data/reference/rios/rios_genome.fa \
       /home/uecha/work/MonsterGenome_Rios/reference/rios_genome.fa
```

**コマンドの説明**：
- `ln -sf`：シンボリックリンクを作成（`-s`：シンボリック、`-f`：既存ファイルを上書き）
- 最初のパス：実際のファイルの場所
- 2番目のパス：リンクを作成する場所

**確認方法**：
```bash
# リンクが正しく作成されたか確認
ls -lh /home/uecha/work/MonsterGenome_Rios/reference/
# 出力例:
# lrwxrwxrwx 1 uecha uecha 45 Dec 11 19:28 rios_genome.fa -> /home/uecha/data/reference/rios/rios_genome.fa
```

**出力の見方**：
- `lrwxrwxrwx`：最初の`l`はシンボリックリンクを表す
- `->`：リンク先を表す

**メリット**：
- ディスク容量を節約できる（同じファイルを複数コピーしない）
- 複数プロジェクトで同じデータを共有できる
- データを更新すると、すべてのプロジェクトに反映される

**注意点**：
- シンボリックリンクは、元のファイルが存在しないと機能しない
- 元のファイルを削除すると、リンクは壊れる

**参考記事**：
- [シンボリックリンクとは？使い方とハードリンクとの違い](https://qiita.com/tags/symlink) - Qiitaのシンボリックリンク記事
- [Linuxシンボリックリンクの使い方](https://eng-entrance.com/linux-command-symlink) - シンボリックリンクの解説

---

## トラブルシューティング

### よくある問題と解決法

#### 1. condaコマンドが見つからない

**症状**：`conda: command not found`

**原因**：condaがインストールされていない、またはパスが通っていない

**解決法**：
```bash
# condaのインストールを確認
which conda

# パスが通っていない場合は、.bashrcや.bash_profileに追加
# （Miniconda/Anacondaインストール時に自動で追加されるはず）

# Windowsの場合、環境変数PATHに追加する必要がある
```

**参考記事**：
- [condaのインストールとパスの設定](https://qiita.com/tags/conda) - Qiitaのconda記事

#### 2. 環境の作成に失敗する

**症状**：`Solving environment: failed`

**原因**：パッケージの依存関係が解決できない

**解決法**：
```bash
# チャンネルを追加
conda config --add channels bioconda
conda config --add channels conda-forge

# チャンネルの優先順位を設定
conda config --set channel_priority strict

# 再度実行
conda create -y -n rios-env ...
```

**参考記事**：
- [condaチャンネルの設定方法](https://qiita.com/tags/conda) - Qiitaのconda記事

#### 3. 環境を有効化できない

**症状**：`conda activate`が動かない

**原因**：condaの初期化が完了していない

**解決法**：
```bash
# condaを初期化
conda init bash  # bashの場合
conda init zsh   # zshの場合
conda init powershell  # PowerShellの場合

# ターミナルを再起動
```

---

## まとめ

今回は、環境構築の基礎について解説しました：

1. **ターミナルの基本操作**：コマンドラインの使い方
2. **conda環境の構築**：R 4.3環境とその他のツール環境を作成
3. **ディレクトリ構成の設計**：大容量データを効率的に管理

**次回の予告**：
次回は、**バックボーンゲノムの取得とデータ準備**について詳しく解説します。Ensemblからゲノムデータをダウンロードし、染色体リネームを行う方法を説明します。

**学習のポイント**：
- ターミナル操作に慣れることが重要です。最初は戸惑うかもしれませんが、繰り返し使うことで自然に覚えられます
- conda環境は、プロジェクトごとに分けることで、ツールの競合を避けられます
- エラーが出たときは、エラーメッセージをよく読んで、参考記事を検索してみてください

お楽しみに！

---

**参考リンク集**：
- [conda公式ドキュメント](https://docs.conda.io/)
- [bioconda公式サイト](https://bioconda.github.io/)
- [Qiita - ターミナル関連記事](https://qiita.com/tags/terminal)
- [Qiita - conda関連記事](https://qiita.com/tags/conda)
- [Qiita - bash関連記事](https://qiita.com/tags/bash)

---

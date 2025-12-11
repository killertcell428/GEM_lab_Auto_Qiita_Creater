# Qiita記事投稿自動化パイプライン（CrewAI版） - 操作マニュアル

このマニュアルでは、CrewAI版パイプラインの基本的なコマンド操作を説明します。

## 目次

- [初回セットアップ](#初回セットアップ)
- [基本的な使い方](#基本的な使い方)
- [フルPDCAサイクルの実行](#フルpdcaサイクルの実行)
- [個別Phaseの実行](#個別phaseの実行)
- [ドラフトファイルの準備](#ドラフトファイルの準備)
- [よくある操作パターン](#よくある操作パターン)
- [トラブルシューティング](#トラブルシューティング)

---

## 初回セットアップ

### 1. 前提条件

- **Python 3.12** がインストールされていること
- **uv** がインストールされていること（推奨）または通常のpip
- **Node.js** がインストールされていること（Web UIを使用する場合）

### 2. uvのインストール（推奨）

**Windows (PowerShell):**
```powershell
powershell -ExecutionPolicy Bypass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

インストール後、PowerShellを再起動するか、以下を実行：
```powershell
$env:PATH = [System.Environment]::GetEnvironmentVariable("Path","Machine") + ";" + [System.Environment]::GetEnvironmentVariable("Path","User")
```

**確認:**
```powershell
uv --version
```

### 3. プロジェクトのセットアップ

#### 方法A: 自動セットアップスクリプト（推奨）

**Windows:**
```powershell
.\setup_uv.ps1
```

#### 方法B: 手動セットアップ

**uvを使用する場合:**

```powershell
# プロジェクトディレクトリに移動
cd C:\Users\uecha\Project_P\Personal\GEM_lab_Auto_Qiita_Creater

# Python 3.12で仮想環境を作成
uv venv .venv --python 3.12

# 仮想環境を有効化
.\.venv\Scripts\Activate.ps1

# 依存関係をインストール
uv pip install -r requirements.txt
uv pip install -r api/requirements.txt
```

**通常のpipを使用する場合:**

```powershell
# 仮想環境を作成
python -m venv .venv

# 仮想環境を有効化
.\.venv\Scripts\Activate.ps1

# 依存関係をインストール
pip install -r requirements.txt
pip install -r api/requirements.txt
```

### 4. 設定ファイルの準備

#### 環境変数の設定

プロジェクトルートに `.env` ファイルを作成：

```env
# Google Gemini API キー（必須: 記事生成とAgentの推論に使用）
GEMINI_API_KEY=your_gemini_api_key_here

# Qiita API アクセストークン（必須: 記事投稿時に使用）
QIITA_ACCESS_TOKEN=your_qiita_access_token_here

# Serper API キー（オプション: Web検索機能を使用する場合）
SERPER_API_KEY=your_serper_api_key_here
```

**最小限のセットアップ**: 記事生成を試すだけなら、`GEMINI_API_KEY` のみ設定すれば動作します。

#### 設定ファイルの確認

以下のファイルが存在することを確認：

```powershell
# 設定ファイルが存在しない場合は、例ファイルからコピー
if (-not (Test-Path config\config.json)) {
    Copy-Item config\config.example.json config\config.json
}
if (-not (Test-Path config\crewai_config.json)) {
    Copy-Item config\crewai_config.example.json config\crewai_config.json
}
```

### 5. 起動方法

#### CLIコマンドで実行（推奨）

```powershell
# 仮想環境を有効化（毎回必要）
.\.venv\Scripts\Activate.ps1

# フルPDCAサイクルを実行
python main.py pdca --topic "Pythonでデータ解析" --draft-file "data/drafts/test_draft.md"
```

#### Web UIで実行

**PowerShellから実行する場合（推奨）:**
```powershell
# PowerShellスクリプト版を使用（推奨）
.\start_dev.ps1
```

**注意**: PowerShellの実行ポリシーでエラーが出る場合：
```powershell
# 実行ポリシーを確認
Get-ExecutionPolicy

# 実行ポリシーを変更（現在のユーザーのみ）
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser

# または、一時的に実行ポリシーをバイパス
powershell -ExecutionPolicy Bypass -File .\start_dev.ps1
```

**CMDから実行する場合:**
```cmd
start_dev.bat
```

**PowerShellから.batファイルを実行する場合:**
```powershell
# 方法1: cmd経由で実行
cmd /c start_dev.bat

# 方法2: 呼び出し演算子を使用
& .\start_dev.bat
```

または手動で：

**ターミナル1 - FastAPI:**
```powershell
.\.venv\Scripts\Activate.ps1
uvicorn api.app.main:app --reload --host 0.0.0.0 --port 8000
```

**ターミナル2 - Next.js:**
```powershell
cd web
npm install  # 初回のみ
npm run dev
```

ブラウザで以下にアクセス：
- **ダッシュボード**: http://localhost:3000
- **APIドキュメント**: http://localhost:8000/docs

---

## 基本的な使い方

### 仮想環境の有効化

**重要**: CLIコマンドを実行する前に、必ず仮想環境を有効化してください。

**Windows (PowerShell):**
```powershell
.\.venv\Scripts\Activate.ps1
```

**Windows (CMD):**
```cmd
.venv\Scripts\activate.bat
```

**Linux/Mac:**
```bash
source .venv/bin/activate
```

### ヘルプの表示

```powershell
# 仮想環境を有効化した後
python main.py --help

# 各コマンドの詳細ヘルプ
python main.py pdca --help
python main.py plan --help
python main.py do --help
python main.py check --help
python main.py act --help
```

---

## フルPDCAサイクルの実行

### 1. 基本的な使い方

```bash
# ドラフトファイルを指定してフルPDCAサイクルを実行
python main.py pdca --topic "Pythonでデータ解析" --draft-file "data/drafts/test_draft.md"
```

**実行内容**:
- Plan Phase: Web検索 → 技術調査 → 構成設計
- Do Phase: 記事執筆 → 実装検証 → ノウハウ抽出 → レビュー
- Qiita投稿（オプション）
- Check Phase: パフォーマンス分析 → 外部ベンチマーク分析 → ドメイントレンド分析（投稿後）
- Act Phase: 改善指針の抽象化 → 次回計画への反映

### 2. 自動投稿も実行

```bash
# フルPDCAサイクル + 自動投稿
python main.py pdca --topic "Pythonでデータ解析" --draft-file "data/drafts/test_draft.md" --auto-publish
```

**注意**: このコマンドは承認をスキップして直接Qiitaに投稿します。

### 3. 追加の要望を指定

```bash
# 追加の要望を指定（将来的に実装予定）
# 現在は、ドラフトファイル内に要望を記述するか、記事生成後に手動で編集してください
```

---

## 個別Phaseの実行

### Plan Phase（調査・構成設計）

```bash
# 新規記事のPlan Phaseを実行
python main.py plan --topic "Pythonでデータ解析" --article-id "article_20251207_123456"

# ドラフトファイルを指定
python main.py plan --topic "Pythonでデータ解析" --article-id "article_20251207_123456" --draft-file "data/drafts/test_draft.md"
```

**実行内容**:
- Web検索: 公式ドキュメント、実装例ブログ、GitHubリポジトリを検索
- 技術調査: 技術の概要、最新トレンド、実装方法を調査
- 構成設計: 記事構成プラン（タイトル、セクション、コード例の配置）を設計

**出力**: `data/state/articles/{article_id}.json` にArticleStateが保存されます。

### Do Phase（執筆・検証・レビュー）

```bash
# Plan Phase完了後の記事IDを指定
python main.py do --article-id "article_20251207_123456"
```

**実行内容**:
- 記事執筆: 構成プランに基づいて記事本文を執筆
- 実装検証: 記事内のPythonコードを実際に実行し、エラーを検出・修正案を提示
- ノウハウ抽出: 実装検証結果から実践的なノウハウを抽出
- レビュー: 技術的正確性、構成の論理性、コード例の完全性をチェック

**出力**: `data/state/articles/{article_id}.json` にArticleStateが更新されます。

### Check Phase（パフォーマンス分析）

```bash
# Do Phase完了後、Qiitaに投稿済みの記事IDを指定
python main.py check --article-id "article_20251207_123456"
```

**実行内容**:
- パフォーマンス分析: いいね数、閲覧数、コメント数を分析
- 外部ベンチマーク分析: 同じタグ・トピックの他記事との比較分析
- ドメイントレンド分析: バイオインフォマティクス領域の最新トレンドを分析

**注意**: Check Phaseは、記事がQiitaに投稿されている場合のみ実行可能です。

**出力**: `data/state/articles/{article_id}.json` に分析結果が保存されます。

### Act Phase（改善指針の抽象化）

```bash
# Check Phase完了後の記事IDを指定
python main.py act --article-id "article_20251207_123456"
```

**実行内容**:
- 改善指針の抽象化: レビュー結果と分析結果を統合し、今後の記事作成に活かすための具体的な改善指針を抽象化
- 次回計画への反映: 改善指針を次回の記事計画に反映する方法を検討

**出力**: `data/state/pdca/main_pdca_cycle.json` にLessons Learned履歴が追加されます。

---

## ドラフトファイルの準備

### 1. ドラフトファイルの作成

ドラフトファイルはMarkdown形式で、`data/drafts/` ディレクトリに保存します。例：

```markdown
# Pythonでデータ解析を効率化する方法

この記事では、Pythonを使ったデータ解析の効率化について説明します。

## はじめに

データ解析は時間がかかる作業ですが、適切なツールを使うことで効率化できます。

## 実装例

```python
import pandas as pd

# データの読み込み
df = pd.read_csv('data.csv')
```
```

### 2. ドラフトファイルの形式

- **タイトル**: 最初の `# ` で始まる行から自動的に抽出されます
- **本文**: Markdown形式で記述
- **タグ**: `タグ:` または `tags:` で始まる行から抽出されます（オプション）

### 3. ドラフトファイルの例

`data/drafts/` ディレクトリに新しいドラフトファイルを作成する例：

```bash
# 新しいドラフトファイルを作成
cat > data/drafts/my_new_draft.md << 'EOF'
# Pythonでデータ解析を効率化する方法

この記事では、Pythonを使ったデータ解析の効率化について説明します。

## はじめに

データ解析は時間がかかる作業ですが、適切なツールを使うことで効率化できます。

## 実装例

```python
import pandas as pd

# データの読み込み
df = pd.read_csv('data.csv')
```
EOF

# ドラフトを調整して投稿
python main.py pdca --topic "Pythonでデータ解析" --draft-file "data/drafts/my_new_draft.md" --auto-publish
```

---

## よくある操作パターン

### パターン1: ドラフトファイルから記事を作成して投稿

```bash
# 1. ドラフトファイルを作成（data/drafts/ディレクトリに保存）
# 例: data/drafts/my_draft.md

# 2. フルPDCAサイクルを実行して自動投稿
python main.py pdca --topic "Pythonでデータ解析" --draft-file "data/drafts/my_draft.md" --auto-publish
```

### パターン2: 段階的に実行（各Phaseを確認しながら）

```bash
# 1. Plan Phaseを実行
python main.py plan --topic "Pythonでデータ解析" --article-id "article_20251207_123456" --draft-file "data/drafts/my_draft.md"

# 2. 構成プランを確認（data/state/articles/article_20251207_123456.json）

# 3. Do Phaseを実行
python main.py do --article-id "article_20251207_123456"

# 4. 記事本文を確認（data/state/articles/article_20251207_123456.json）

# 5. 必要に応じて手動で編集

# 6. Qiitaに投稿（手動、または将来的にCLIコマンドを追加予定）

# 7. Check Phaseを実行（投稿後）
python main.py check --article-id "article_20251207_123456"

# 8. Act Phaseを実行
python main.py act --article-id "article_20251207_123456"
```

### パターン3: 既存の記事を再調整

```bash
# 1. 既存のArticleStateを読み込んで、Do Phaseから再実行
python main.py do --article-id "article_20251207_123456"
```

### パターン4: 過去の投稿を取得してから記事生成

```bash
# 1. 過去の投稿を取得（将来的にCLIコマンドを追加予定）
# 現在は、src/publisher/qiita_fetcher.pyを直接使用

# 2. フルPDCAサイクルを実行（過去の投稿が自動的に参照されます）
python main.py pdca --topic "Pythonでデータ解析" --draft-file "data/drafts/my_draft.md"
```

---

## Stateファイルの確認

### ArticleStateの確認

```bash
# JSONファイルを直接確認
cat data/state/articles/article_20251207_123456.json

# または、Pythonで読み込んで確認
python -c "from src.crewai.state.article_state import ArticleState; import json; state = ArticleState.load('article_20251207_123456'); print(json.dumps(state.to_dict(), ensure_ascii=False, indent=2))"
```

### PDCAStateの確認

```bash
# JSONファイルを直接確認
cat data/state/pdca/main_pdca_cycle.json
```

---

## トラブルシューティング

### エラーが発生した場合

1. **仮想環境が有効化されているか確認**
   ```powershell
   # 仮想環境を有効化
   .\.venv\Scripts\Activate.ps1
   
   # Pythonのパスを確認（.venv内のPythonが表示されるはず）
   where python
   ```

2. **エラーメッセージを確認**
   ```powershell
   python main.py pdca --topic "テスト" --draft-file "data/drafts/test_draft.md"
   # エラーメッセージが表示されます
   ```

3. **環境変数を確認**
   ```powershell
   # .envファイルが存在するか確認
   Get-Content .env
   ```

4. **設定ファイルを確認**
   ```powershell
   # config.jsonが存在するか確認
   Get-Content config\config.json
   ```

### 実行が途中で止まった場合

- `Ctrl+C`で停止して、再度実行
- エラーログを確認
- 環境変数や設定ファイルを再確認
- Stateファイルを確認（`data/state/articles/`）

### 記事が生成されない場合

- Gemini APIキーが正しく設定されているか確認
- 設定ファイルの内容を確認
- ドラフトファイルが正しく読み込まれているか確認
- Stateファイルを確認（`data/state/articles/`）

### CrewAIの実行エラー

- `config/crewai_config.json` が正しく設定されているか確認
- `.env` ファイルに必要なAPIキーが設定されているか確認
- `data/state/` ディレクトリが作成されているか確認
- 仮想環境が有効化されているか確認（`.\.venv\Scripts\Activate.ps1`）

### 仮想環境が見つからない場合

```powershell
# 仮想環境を再作成
uv venv .venv --python 3.12

# または通常のpipを使用
python -m venv .venv

# 依存関係を再インストール
.\.venv\Scripts\Activate.ps1
uv pip install -r requirements.txt
uv pip install -r api/requirements.txt
```

---

## コマンド一覧（クイックリファレンス）

| コマンド | 説明 | 例 |
|---------|------|-----|
| `python main.py pdca` | フルPDCAサイクルを実行 | `python main.py pdca --topic "テスト" --draft-file "data/drafts/test.md"` |
| `python main.py plan` | Plan Phaseのみ実行 | `python main.py plan --topic "テスト" --article-id "article_xxx"` |
| `python main.py do` | Do Phaseのみ実行 | `python main.py do --article-id "article_xxx"` |
| `python main.py check` | Check Phaseのみ実行 | `python main.py check --article-id "article_xxx"` |
| `python main.py act` | Act Phaseのみ実行 | `python main.py act --article-id "article_xxx"` |

---

## 補足情報

### 記事の保存場所

- **Stateファイル**: `data/state/articles/{article_id}.json`
- **記事本文**: Stateファイル内の `content` フィールドに保存されます

### 設定ファイル

- **設定ファイル**: `config/config.json`
- **CrewAI設定**: `config/crewai_config.json`
- **環境変数**: `.env` ファイル

### プロンプトファイル

- **リサーチプロンプト**: `data/prompts/research_prompt.md`
- **記事生成プロンプト**: `data/prompts/article_generation_prompt.md`

### アーキテクチャドキュメント

詳細なアーキテクチャとワークフローについては、[ARCHITECTURE.md](ARCHITECTURE.md)を参照してください。

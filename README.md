# Qiita記事投稿自動化パイプライン（CrewAI版）

CrewAIフレームワークを使用した、Qiita記事の生成から分析、改善までを自動化するPDCAサイクルシステムです。

## 概要

本システムは、CrewAIベースのマルチエージェントシステムにより、以下の機能を提供します：

- **PDCAサイクル**: Plan（調査・構成設計）→ Do（執筆・検証・レビュー）→ Check（パフォーマンス分析）→ Act（改善指針の抽象化）の自動実行
- **Web検索・実装検証**: 公式ドキュメント、実装例ブログ、GitHubリポジトリの検索と、コードの実際の実行・検証
- **Human-in-the-Loop**: 任意のタイミングで人間が介入し、自然言語で指示・修正が可能
- **State管理**: 記事の状態とPDCAサイクルの履歴を永続化
- **自動投稿**: Qiita API v2への自動投稿（オプション）

## アーキテクチャ

詳細なワークフロー図とアーキテクチャについては、[ARCHITECTURE.md](ARCHITECTURE.md)を参照してください。

### 主要コンポーネント

- **Agents**: 各Phaseで特定の役割を担うエージェント（Researcher、Planner、Writer、Reviewerなど）
- **Tasks**: 各Agentが実行する具体的なタスク
- **Crews**: 複数のAgentとTaskを組み合わせた実行単位（Plan Crew、Do Crew、Check Crew、Act Crew）
- **Tools**: Agentが使用するツール（Gemini API、Qiita API、Web検索、コード実行など）
- **State**: 記事の状態とPDCAサイクルの履歴を管理

## セットアップ

### 前提条件

- **Python 3.12** がインストールされていること
- **uv** がインストールされていること（推奨）または通常のpip
- **Node.js** がインストールされていること（Web UIを使用する場合）

### 1. uvのインストール（推奨）

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

詳細は [docs/SETUP_UV.md](docs/SETUP_UV.md) を参照してください。

### 2. プロジェクトのセットアップ

#### 方法A: 自動セットアップスクリプト（推奨）

**Windows:**
```powershell
.\setup_uv.ps1
```

#### 方法B: 手動セットアップ

**uvを使用する場合:**

```powershell
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

### 3. 設定ファイルの準備

#### 設定ファイルの作成

**Windows (PowerShell):**
```powershell
# 設定ファイルが存在しない場合は、例ファイルからコピー
if (-not (Test-Path config\config.json)) {
    Copy-Item config\config.example.json config\config.json
}
if (-not (Test-Path config\crewai_config.json)) {
    Copy-Item config\crewai_config.example.json config\crewai_config.json
}
```

**Linux/Mac:**
```bash
cp config/config.example.json config/config.json
cp config/crewai_config.example.json config/crewai_config.json
```

必要に応じて `config/config.json` と `config/crewai_config.json` を編集してください。

#### 環境変数の設定

`.env` ファイルを作成し、以下の環境変数を設定してください：

```env
# Google Gemini API キー（必須: 記事生成とAgentの推論に使用）
GEMINI_API_KEY=your_gemini_api_key_here

# Qiita API アクセストークン（必須: 記事投稿時に使用）
QIITA_ACCESS_TOKEN=your_qiita_access_token_here

# Serper API キー（オプション: Web検索機能を使用する場合）
SERPER_API_KEY=your_serper_api_key_here
```

**最小限のセットアップ**: 記事生成を試すだけなら、`GEMINI_API_KEY` のみ設定すれば動作します。Qiitaへの投稿を行う場合のみ `QIITA_ACCESS_TOKEN` が必要です。Web検索機能を使用する場合は `SERPER_API_KEY` が必要です。

### 3. APIキーの取得方法

#### Google Gemini API キー（必須）

1. [Google AI Studio](https://makersuite.google.com/app/apikey) にアクセス
2. Googleアカウントでログイン
3. 「Create API Key」をクリック
4. 生成されたAPIキーを `.env` ファイルの `GEMINI_API_KEY` に設定

#### Qiita API アクセストークン（投稿時のみ必須）

1. [Qiitaの設定ページ](https://qiita.com/settings/applications) にアクセス
2. 「個人用アクセストークン」セクションで「新しくトークンを発行する」をクリック
3. トークン名を入力し、スコープを選択（`write_qiita` が必要）
4. 生成されたアクセストークンを `.env` ファイルの `QIITA_ACCESS_TOKEN` に設定

#### Serper API キー（Web検索機能を使用する場合）

1. [Serper.dev](https://serper.dev/) にアクセス
2. アカウントを作成し、APIキーを取得
3. 生成されたAPIキーを `.env` ファイルの `SERPER_API_KEY` に設定

## 使い方

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

### 基本的な使い方

#### フルPDCAサイクルの実行

```powershell
# 仮想環境を有効化した後
python main.py pdca --topic "Pythonでデータ解析" --draft-file "data/drafts/test_draft.md"

# 自動投稿も実行
python main.py pdca --topic "Pythonでデータ解析" --draft-file "data/drafts/test_draft.md" --auto-publish
```

#### 個別Phaseの実行

```powershell
# Plan Phaseのみ実行
python main.py plan --topic "Pythonでデータ解析" --article-id "article_20251207_123456"

# Do Phaseのみ実行
python main.py do --article-id "article_20251207_123456"

# Check Phaseのみ実行（投稿後）
python main.py check --article-id "article_20251207_123456"

# Act Phaseのみ実行
python main.py act --article-id "article_20251207_123456"
```

#### Web UIで実行

**PowerShellから実行する場合（推奨）:**
```powershell
# PowerShellスクリプト版を使用（推奨）
.\start_dev.ps1
```

**注意**: PowerShellの実行ポリシーでエラーが出る場合：
```powershell
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

### ドラフトファイルの準備

ドラフトファイルはMarkdown形式で、`data/drafts/` ディレクトリに保存します。例：

```markdown
# 記事のタイトル（オプション）

記事の本文をここに記述します。

## セクション1

内容...
```

### その他のコマンド

詳細なコマンド一覧と使い方は、[MANUAL.md](MANUAL.md)を参照してください。

## プロジェクト構造

```
.
├── config/
│   ├── config.example.json          # 設定ファイルのテンプレート
│   ├── config.json                  # 実際の設定ファイル（作成が必要）
│   ├── crewai_config.example.json   # CrewAI設定のテンプレート
│   └── crewai_config.json           # CrewAI設定（作成が必要）
├── src/
│   ├── crewai/                      # CrewAI関連モジュール
│   │   ├── agents/                  # Agent定義
│   │   ├── tasks/                   # Task定義
│   │   ├── crews/                   # Crew定義
│   │   ├── tools/                   # Tool定義
│   │   ├── state/                   # State管理
│   │   ├── human_loop.py            # Human-in-the-Loop処理
│   │   └── orchestrator.py          # メインオーケストレーター
│   ├── generator/                    # 記事生成モジュール（一部機能はCrewAIに統合）
│   ├── publisher/                   # Qiita投稿モジュール
│   ├── storage/                      # 記事管理モジュール
│   └── analyzer/                    # 分析モジュール
├── data/
│   ├── drafts/                      # ドラフトファイル（Markdown）
│   ├── prompts/                     # プロンプトファイル（Markdown）
│   ├── state/                       # Stateファイル（JSON）
│   │   ├── articles/                # ArticleState
│   │   ├── pdca/                    # PDCAState
│   │   └── feedback/                # FeedbackQueue
│   └── analysis/                    # 分析データ（JSON）
├── main.py                           # エントリーポイント
├── ARCHITECTURE.md                   # アーキテクチャドキュメント
├── MANUAL.md                         # 操作マニュアル
└── requirements.txt                  # 依存パッケージ
```

## PDCAサイクルの説明

### Plan Phase（調査・構成設計）

1. **Web検索**: 公式ドキュメント、実装例ブログ、GitHubリポジトリを検索
2. **技術調査**: 技術の概要、最新トレンド、実装方法、バイオインフォマティクスへの応用可能性を調査
3. **構成設計**: 調査結果を基に記事構成プラン（タイトル、セクション、コード例の配置）を設計

### Do Phase（執筆・検証・レビュー）

1. **記事執筆**: 構成プランに基づいて記事本文を執筆（Markdown形式）
2. **実装検証**: 記事内のPythonコードを実際に実行し、エラーを検出・修正案を提示
3. **ノウハウ抽出**: Web検索結果、技術調査結果、実装検証結果から実践的なノウハウを抽出
4. **レビュー**: 技術的正確性、構成の論理性、コード例の完全性、文字数、セクション完全性をチェック

### Check Phase（パフォーマンス分析）

1. **パフォーマンス分析**: 投稿後の記事のKPI（いいね数、閲覧数、コメント数）を分析
2. **外部ベンチマーク分析**: 同じタグ・トピックの他記事との比較分析
3. **ドメイントレンド分析**: バイオインフォマティクス領域の最新トレンドを分析し、次回記事のトピック候補を提案

### Act Phase（改善指針の抽象化）

1. **改善指針の抽象化**: レビュー結果と分析結果を統合し、今後の記事作成に活かすための具体的な改善指針を抽象化
2. **次回計画への反映**: 改善指針を次回の記事計画に反映する方法を検討

## Human-in-the-Loop機能

任意のPhaseで人間が介入し、自然言語で指示・修正が可能です。

```bash
# Human Feedbackを追加（CLIから直接は未実装、将来的に実装予定）
# 現在は、Stateファイルを直接編集するか、コードからHumanLoopクラスを使用
```

詳細は[ARCHITECTURE.md](ARCHITECTURE.md)の「Human-in-the-Loop機能」セクションを参照してください。

## トラブルシューティング

### エラー: "GEMINI_API_KEYが設定されていません"

`.env` ファイルが正しく作成され、`GEMINI_API_KEY` が設定されているか確認してください。

### エラー: "QIITA_ACCESS_TOKENが設定されていません"

`.env` ファイルに `QIITA_ACCESS_TOKEN` が設定されているか確認してください。

### エラー: "設定ファイルが見つかりません"

`config/config.json` ファイルが存在するか確認してください。`config/config.example.json` をコピーして作成してください。

### エラー: "ドラフトファイルが見つかりません"

指定したドラフトファイルのパスが正しいか確認してください。`data/drafts/`ディレクトリ内のファイルを指定する場合、相対パスまたは絶対パスを使用できます。

### エラー: "モジュールが見つかりません" または "No module named 'xxx'"

仮想環境が有効化されていない可能性があります。以下を確認してください：

```powershell
# 仮想環境を有効化
.\.venv\Scripts\Activate.ps1

# Pythonのパスを確認（.venv内のPythonが表示されるはず）
where python

# 依存関係を再インストール
uv pip install -r requirements.txt
uv pip install -r api/requirements.txt
```

### CrewAIの実行エラー

CrewAIの実行中にエラーが発生した場合、以下を確認してください：

1. 仮想環境が有効化されているか（`.\.venv\Scripts\Activate.ps1`）
2. `config/crewai_config.json` が正しく設定されているか
3. `.env` ファイルに必要なAPIキーが設定されているか
4. `data/state/` ディレクトリが作成されているか

## 今後の拡張予定

- Human FeedbackのCLIインターフェース
- Web UI
- 複数のSNSプラットフォーム対応
- より高度な分析機能

## ライセンス

このプロジェクトは個人利用を想定しています。

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

### 1. 依存パッケージのインストール

```bash
pip install -r requirements.txt
```

### 2. 設定ファイルの準備

#### 設定ファイルの作成

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

### 基本的な使い方

#### フルPDCAサイクルの実行

```bash
# ドラフトファイルを指定してフルPDCAサイクルを実行
python main.py pdca --topic "Pythonでデータ解析" --draft-file "data/drafts/test_draft.md"

# 自動投稿も実行
python main.py pdca --topic "Pythonでデータ解析" --draft-file "data/drafts/test_draft.md" --auto-publish
```

#### 個別Phaseの実行

```bash
# Plan Phaseのみ実行
python main.py plan --topic "Pythonでデータ解析" --article-id "article_20251207_123456"

# Do Phaseのみ実行
python main.py do --article-id "article_20251207_123456"

# Check Phaseのみ実行（投稿後）
python main.py check --article-id "article_20251207_123456"

# Act Phaseのみ実行
python main.py act --article-id "article_20251207_123456"
```

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

### CrewAIの実行エラー

CrewAIの実行中にエラーが発生した場合、以下を確認してください：

1. `config/crewai_config.json` が正しく設定されているか
2. `.env` ファイルに必要なAPIキーが設定されているか
3. `data/state/` ディレクトリが作成されているか

## 今後の拡張予定

- Human FeedbackのCLIインターフェース
- Web UI
- 複数のSNSプラットフォーム対応
- より高度な分析機能

## ライセンス

このプロジェクトは個人利用を想定しています。

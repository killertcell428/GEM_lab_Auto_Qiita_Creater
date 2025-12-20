# 毎週金曜日自動投稿スケジューラーの使い方

## 概要

`data/drafts/docs/blog/` ディレクトリ内のドラフトファイルを順番に読み込んで、毎週金曜日に自動的にQiitaに投稿するスケジューラーです。

## ドラフトファイル

以下の8つの記事が順番に投稿されます：

1. `01_introduction.md` - 第1回：イントロダクション
2. `02_environment_setup.md` - 第2回：環境構築とデータ準備
3. `03_simulation.md` - 第3回：トランスクリプト設計とRNA-seqシミュレーション
4. `04_alignment_counting.md` - 第4回：アライメントとカウント
5. `05_deseq2_analysis.md` - 第5回：発現解析と統計検定
6. `06_visualization.md` - 第6回：可視化と解釈
7. `07_cnv_analysis.md` - 第7回：CNV解析とスライド用データ作成
8. `08_summary.md` - 第8回：まとめと今後の展望

すべての記事を投稿し終わると、最初から再開します。

## 設定

### 投稿時刻の設定

`config/config.json` で金曜日の投稿時刻を設定できます：

```json
{
  "workflow": {
    "schedule": {
      "friday": {
        "enabled": true,
        "hour": 10,
        "minute": 0
      }
    }
  }
}
```

- `hour`: 投稿時刻（時、0-23）
- `minute`: 投稿時刻（分、0-59）
- デフォルト: 金曜日 10:00

### 環境変数の設定

`.env` ファイルに以下を設定してください：

```env
# Qiita API トークン（必須）
QIITA_ACCESS_TOKEN=your_qiita_access_token_here

# GitHub API トークン（画像の自動アップロード用、オプション）
GITHUB_TOKEN=your_github_personal_access_token_here
GITHUB_REPO=username/repository_name
GITHUB_BRANCH=main
```

**GitHub設定について**:
- `GITHUB_TOKEN`: GitHub Personal Access Token（`repo` スコープが必要）
- `GITHUB_REPO`: 画像をアップロードするリポジトリ名（例: `username/repo`）
- `GITHUB_BRANCH`: ブランチ名（デフォルト: `main`）

GitHub設定がない場合、画像はコメントアウトされ、後で手動でURLを設定できます。

## 起動方法

### PowerShell（推奨）

```powershell
.\run_scheduler.ps1
```

### バッチファイル

```cmd
run_scheduler.bat
```

### 直接実行

```powershell
# 仮想環境を有効化
.\.venv\Scripts\Activate.ps1

# スケジューラーを起動
python run_scheduler.py
```

## 動作確認

スケジューラーを起動すると、以下の情報が表示されます：

```
============================================================
Qiita記事自動投稿スケジューラー
============================================================

[現在の状態]
  総ドラフト数: 8
  投稿済み: 0
  現在のインデックス: 0
  次回投稿予定: 01_introduction.md
  最終投稿日時: None

[SCHEDULER] スケジューラーが開始されました
[SCHEDULER] 毎週金曜日に自動投稿が実行されます
[SCHEDULER] Ctrl+C で停止できます
```

## 投稿状態の確認

投稿状態は `data/state/scheduled_posts.json` に保存されます：

```json
{
  "published_drafts": [
    "01_introduction.md",
    "02_environment_setup.md"
  ],
  "last_published_date": "2024-12-13T10:00:00",
  "current_index": 2
}
```

## 停止方法

`Ctrl+C` を押すとスケジューラーが停止します。

## 注意事項

- PCをシャットダウンするとスケジューラーも停止します
- PCを起動したら、再度スケジューラーを起動してください
- ドラフトファイルは編集されず、そのまま投稿されます
- 投稿済みの記事は記録され、次回起動時も続きから投稿されます

## 完全自動化の流れ

1. **スケジューラー起動**: `.\run_scheduler.ps1` を実行
2. **金曜日の指定時刻**: 自動的に次のドラフトを処理
3. **記事の前処理**: パス修正、テーブル埋め込み、画像アップロード
4. **Qiita投稿**: 処理済みの記事を自動投稿
5. **完了**: 記事は画像込みで完成（GitHub設定がある場合）

**GitHub設定がない場合**:
- 記事は画像なしで投稿されます
- 投稿後にQiitaで編集して画像URLを手動で設定してください

## 記事の前処理

スケジューラーは、投稿前に記事を自動的に処理します：

### 1. パス参照の修正
- 相対パス参照（`docs/`、`reference/`など）は削除されます
- Qiitaでは参照できないため、リンクは無効化されます

### 2. テーブルファイルの埋め込み
- `tables/` ディレクトリ内のCSV/TSVファイルは自動的に記事内に埋め込まれます
- 例: ````csv:tables/01_chrom_rename.tsv` → Markdownテーブル形式に変換

### 3. 画像ファイルの処理
- 画像ファイル（`plots/` ディレクトリ内）への参照を自動処理
- **GitHub設定がある場合**: 自動的にGitHubにアップロードして、raw URLを取得して記事に埋め込み
- **GitHub設定がない場合**: 画像参照をコメントアウト（後で手動でURLを設定可能）

### 画像の自動アップロード（推奨）

GitHub設定を完了すると、画像が自動的にアップロードされます：

1. **GitHub Personal Access Tokenの取得**:
   - GitHub Settings > Developer settings > Personal access tokens > Tokens (classic)
   - `repo` スコープを有効にしてトークンを生成

2. **環境変数の設定**:
   ```env
   GITHUB_TOKEN=your_github_token_here
   GITHUB_REPO=username/repository_name
   GITHUB_BRANCH=main
   ```

3. **自動アップロード**:
   - スケジューラーが画像を検出すると、自動的にGitHubにアップロード
   - アップロードされた画像のraw URLが記事に自動的に埋め込まれます

### 画像の手動設定（GitHub設定がない場合）

1. **GitHubに画像をアップロード**:
   - リポジトリに画像ファイルをアップロード
   - raw URLを取得（例: `https://raw.githubusercontent.com/user/repo/main/image.png`）

2. **記事内の画像URLを設定**:
   - 投稿後にQiitaで記事を編集
   - `<!-- TODO: 画像URLを設定してください: image.png -->` の部分を探す
   - `![説明](画像URL)` の形式で画像参照を追加

3. **または画像参照を削除**:
   - 画像が不要な場合は、コメントアウトされた部分を削除

## トラブルシューティング

### スケジューラーが起動しない

1. 仮想環境が有効化されているか確認
2. `.env` ファイルに `QIITA_ACCESS_TOKEN` が設定されているか確認
3. `config/config.json` が存在するか確認

### 投稿が実行されない

1. スケジューラーが起動しているか確認
2. 金曜日の設定時刻になっているか確認
3. `data/drafts/docs/blog/` にドラフトファイルが存在するか確認
4. `data/state/scheduled_posts.json` の状態を確認

### エラーログの確認

スケジューラーの実行ログを確認してください。エラーが発生した場合は、詳細なエラーメッセージが表示されます。

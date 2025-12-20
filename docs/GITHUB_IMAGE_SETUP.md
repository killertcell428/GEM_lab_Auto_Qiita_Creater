# GitHub画像アップロード設定ガイド

## 概要

スケジューラーは、記事内の画像を自動的にGitHubにアップロードし、raw URLを取得して記事に埋め込むことができます。これにより、完全自動化が可能になります。

## 設定手順

### 1. GitHub Personal Access Tokenの取得

1. GitHubにログイン
2. Settings > Developer settings > Personal access tokens > Tokens (classic)
3. "Generate new token (classic)" をクリック
4. トークン名を入力（例: "Qiita Image Uploader"）
5. スコープを選択:
   - ✅ `repo` (Full control of private repositories)
6. "Generate token" をクリック
7. **トークンをコピー**（この画面を閉じると二度と表示されません）

### 2. リポジトリの準備

画像をアップロードするリポジトリを準備します：

- 既存のリポジトリを使用するか、新しいリポジトリを作成
- リポジトリ名をメモ（例: `username/qiita-images`）

### 3. 環境変数の設定

`.env` ファイルに以下を追加：

```env
# GitHub API 設定（画像の自動アップロード用）
GITHUB_TOKEN=ghp_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
GITHUB_REPO=username/qiita-images
GITHUB_BRANCH=main
```

**注意**:
- `GITHUB_TOKEN`: 上記で取得したPersonal Access Token
- `GITHUB_REPO`: `username/repository_name` の形式
- `GITHUB_BRANCH`: デフォルトは `main`（必要に応じて変更）

### 4. 動作確認

スケジューラーを起動すると、画像が見つかった場合に自動的にアップロードされます：

```
[INFO] 画像ファイルが見つかりました: volcano_sex_M_vs_F.png
[INFO] 画像をGitHubにアップロードしました: https://raw.githubusercontent.com/username/qiita-images/main/images/2024/12/volcano_sex_M_vs_F.png
```

## 画像の保存場所

アップロードされた画像は、以下のディレクトリ構造で保存されます：

```
リポジトリ/
└── images/
    └── YYYY/
        └── MM/
            ├── image1.png
            ├── image2.png
            └── ...
```

例: `images/2024/12/volcano_sex_M_vs_F.png`

## トラブルシューティング

### エラー: "GitHub設定が不完全です"

- `.env` ファイルに `GITHUB_TOKEN` と `GITHUB_REPO` が設定されているか確認
- 環境変数が正しく読み込まれているか確認

### エラー: "画像のアップロードに失敗しました"

- GitHubトークンに `repo` スコープがあるか確認
- リポジトリ名が正しいか確認（`username/repo` の形式）
- リポジトリへのアクセス権限があるか確認

### 画像がアップロードされない

- 画像ファイルが `data/drafts/docs/plots/` に存在するか確認
- スケジューラーのログを確認
- GitHub設定が完了しているか確認

## セキュリティ注意事項

- **GitHubトークンは機密情報です**。`.env` ファイルをGitにコミットしないでください
- `.gitignore` に `.env` が含まれているか確認してください
- トークンが漏洩した場合は、すぐにGitHubでトークンを無効化してください

## オプション: GitHub設定なしでの運用

GitHub設定がない場合でも、スケジューラーは動作します：

- 画像参照はコメントアウトされます
- 記事は画像なしで投稿されます
- 投稿後にQiitaで編集して画像URLを手動で設定できます

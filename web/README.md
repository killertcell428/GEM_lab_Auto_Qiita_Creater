# Web UI (Next.js)

Qiita記事自動化パイプラインのWeb UI実装です。

## セットアップ

### 依存関係のインストール

```bash
npm install
```

### 開発サーバーの起動

```bash
npm run dev
```

開発サーバーは `http://localhost:3000` で起動します。

## 環境変数

`.env.local` ファイルを作成して、以下の環境変数を設定してください：

```env
NEXT_PUBLIC_API_URL=http://localhost:8000
```

## 機能

- **ダッシュボード**: 進行中の記事一覧表示
- **新規記事作成**: テーマ入力とオプション設定
- **記事詳細**: Markdownエディタ、AI Process、Feedback、Analysisタブ
- **Human Feedback**: 自然言語でのフィードバック入力
- **リアルタイム更新**: SSEによる進捗ストリーミング（将来実装）

## ディレクトリ構造

```
web/
├── app/                    # Next.js App Router
│   ├── page.tsx           # ダッシュボード
│   ├── new/               # 新規記事作成
│   └── articles/[id]/     # 記事詳細
├── components/            # Reactコンポーネント
│   ├── Layout.tsx
│   ├── StatusBadge.tsx
│   ├── ArticleStatusCard.tsx
│   ├── ArticleViewer.tsx
│   ├── HumanFeedbackPanel.tsx
│   └── TabContainer.tsx
├── lib/                   # ユーティリティ
│   └── api.ts            # APIクライアント
└── hooks/                 # カスタムフック
    └── useArticleStream.ts
```

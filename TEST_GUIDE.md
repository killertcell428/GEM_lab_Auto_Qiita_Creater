# CrewAIシステム テストガイド

## セットアップ確認

### 1. 依存パッケージのインストール確認

```bash
pip install -r requirements.txt
```

### 2. 環境変数の確認

`.env`ファイルに以下が設定されていることを確認：

```
GEMINI_API_KEY=your_gemini_api_key
QIITA_ACCESS_TOKEN=your_qiita_access_token
```

### 3. 設定ファイルの確認

`config/config.json`に`crewai`セクションが追加されていることを確認。

## 基本テスト

### インポートテスト

```bash
python test_crewai_basic.py
```

### CLIヘルプ確認

```bash
python main_crewai.py --help
```

## 実行方法

### フルPDCAサイクル実行

```bash
python main_crewai.py pdca --topic "記事のトピック" --auto-publish
```

### 個別Phase実行

```bash
# Plan Phaseのみ
python main_crewai.py plan --topic "記事のトピック"

# Do Phaseのみ
python main_crewai.py do --article-id "article_20251205_120000"

# Check Phaseのみ（投稿済み記事）
python main_crewai.py check --article-id "article_20251205_120000"

# Act Phaseのみ
python main_crewai.py act --article-id "article_20251205_120000"
```

### Human Feedback追加

```bash
python main_crewai.py feedback --article-id "article_20251205_120000" --phase "do" --content "タイトルを変更して" --priority 8
```

### 状態確認

```bash
python main_crewai.py state --article-id "article_20251205_120000"
```

## トラブルシューティング

### CrewAIがインポートできない場合

```bash
pip install --upgrade crewai crewai-tools
```

### 設定ファイルエラー

`config/config.json`に`crewai`セクションが追加されているか確認。

### State保存エラー

`data/state/`ディレクトリが作成されているか確認。


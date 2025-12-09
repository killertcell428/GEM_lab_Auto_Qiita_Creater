# Node.js クイックインストールガイド

## 現在の状況
Node.jsがインストールされていないため、`npm` コマンドが使用できません。

## インストール方法（3つの選択肢）

### 方法1: 公式インストーラー（最も簡単・推奨）

1. **ブラウザで以下にアクセス**
   ```
   https://nodejs.org/
   ```

2. **「推奨版（LTS）」をダウンロード**
   - 自動的にWindows用の `.msi` ファイルがダウンロードされます

3. **インストーラーを実行**
   - ダウンロードした `.msi` ファイルをダブルクリック
   - 「Next」をクリックして進める
   - **重要**: 「Add to PATH」オプションが選択されていることを確認
   - インストールを完了

4. **PowerShellを再起動**
   - 現在のPowerShellウィンドウを閉じる
   - 新しいPowerShellウィンドウを開く

5. **インストール確認**
   ```powershell
   node --version
   npm --version
   ```

### 方法2: wingetを使用（Windows 10/11）

PowerShellを**管理者として実行**して以下を実行：

```powershell
winget install OpenJS.NodeJS.LTS
```

インストール後、PowerShellを再起動して確認：
```powershell
node --version
npm --version
```

### 方法3: Chocolateyを使用（既にインストールされている場合）

PowerShellを**管理者として実行**して以下を実行：

```powershell
choco install nodejs-lts
```

インストール後、PowerShellを再起動して確認：
```powershell
node --version
npm --version
```

## インストール後の手順

### 1. 依存関係のインストール

```powershell
cd C:\Users\uecha\Project_P\Personal\GEM_lab_Auto_Qiita_Creater\web
npm install
```

### 2. Next.jsサーバーの起動

```powershell
npm run dev
```

### 3. ブラウザでアクセス

- **Web UI**: http://localhost:3000
- **API Docs**: http://localhost:8000/docs

## トラブルシューティング

### インストール後も `node` コマンドが認識されない場合

1. **PowerShellを再起動**
   - 環境変数の変更を反映するため、完全に閉じて再起動

2. **環境変数を確認**
   - `Win + R` → `sysdm.cpl` → 「詳細設定」→「環境変数」
   - 「システム環境変数」の「Path」を確認
   - `C:\Program Files\nodejs\` が含まれているか確認

3. **手動でPATHに追加**
   - Node.jsのインストールパス（通常は `C:\Program Files\nodejs\`）をPATHに追加
   - PowerShellを再起動

### インストールパスの確認

```powershell
Get-ChildItem "C:\Program Files\nodejs\" -ErrorAction SilentlyContinue
```

## 推奨インストール方法

**方法1（公式インストーラー）** が最も確実で簡単です。

1. https://nodejs.org/ にアクセス
2. LTS版をダウンロード
3. インストーラーを実行
4. PowerShellを再起動
5. `node --version` で確認

## 次のステップ

Node.jsのインストールが完了したら：

```powershell
# 1. バージョン確認
node --version
npm --version

# 2. プロジェクトディレクトリに移動
cd C:\Users\uecha\Project_P\Personal\GEM_lab_Auto_Qiita_Creater\web

# 3. 依存関係をインストール
npm install

# 4. 開発サーバーを起動
npm run dev
```

その後、ブラウザで http://localhost:3000 にアクセスしてWeb UIをテストできます。


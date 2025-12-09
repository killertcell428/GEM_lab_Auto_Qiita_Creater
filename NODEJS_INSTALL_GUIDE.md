# Node.js インストールガイド

## 問題
`npm` コマンドが認識されない場合、Node.jsがインストールされていないか、PATHに含まれていません。

## 解決方法

### 方法1: Node.js公式サイトからインストール（推奨）

1. **Node.js公式サイトにアクセス**
   - https://nodejs.org/ にアクセス

2. **LTS版をダウンロード**
   - 「推奨版（LTS）」をダウンロード
   - Windows用のインストーラー（.msi）をダウンロード

3. **インストール**
   - ダウンロードしたインストーラーを実行
   - インストールウィザードに従って進める
   - 「Add to PATH」オプションが選択されていることを確認

4. **インストール確認**
   ```powershell
   node --version
   npm --version
   ```

5. **再起動**
   - PowerShellを再起動（または新しいウィンドウを開く）

### 方法2: Chocolateyを使用（既にインストールされている場合）

```powershell
choco install nodejs
```

### 方法3: wingetを使用（Windows 10/11）

```powershell
winget install OpenJS.NodeJS.LTS
```

## インストール後の確認

インストールが完了したら、以下を実行してください：

```powershell
# Node.jsのバージョン確認
node --version

# npmのバージョン確認
npm --version

# プロジェクトディレクトリに移動
cd C:\Users\uecha\Project_P\Personal\GEM_lab_Auto_Qiita_Creater\web

# 依存関係のインストール
npm install

# 開発サーバーの起動
npm run dev
```

## トラブルシューティング

### PATHに追加されていない場合

1. **環境変数の確認**
   - `Win + R` → `sysdm.cpl` → 「詳細設定」→「環境変数」
   - 「システム環境変数」の「Path」を確認
   - `C:\Program Files\nodejs\` が含まれているか確認

2. **手動でPATHに追加**
   - Node.jsのインストールパス（通常は `C:\Program Files\nodejs\`）をPATHに追加

3. **PowerShellの再起動**
   - 環境変数の変更を反映するため、PowerShellを再起動

## 注意事項

- Node.jsをインストールすると、npmも自動的にインストールされます
- インストール後は、PowerShellを再起動する必要があります
- 既にNode.jsがインストールされている場合は、PATHの問題の可能性があります


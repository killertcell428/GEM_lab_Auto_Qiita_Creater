# セットアップ・インストールガイド

このフォルダには、プロジェクトのセットアップとインストールに関するガイドドキュメントが含まれています。

## ドキュメント一覧

### Python環境

- **INSTALL_PYTHON312.md** - Python 3.12のインストール方法
- **SETUP_UV.md** - uvパッケージマネージャーのセットアップガイド
- **QUICK_START_UV.md** - uvを使ったクイックスタートガイド
- **README_VENV.md** - 仮想環境（.venv）のセットアップガイド

### Node.js環境

- **NODEJS_INSTALL_GUIDE.md** - Node.jsのインストールガイド
- **QUICK_INSTALL_NODEJS.md** - Node.jsのクイックインストールガイド
- **install_nodejs.ps1** - Node.js自動インストールスクリプト

## セットアップ手順

### 初回セットアップ

1. **Python 3.12をインストール**
   - `INSTALL_PYTHON312.md` を参照

2. **uvをインストール**
   - `SETUP_UV.md` を参照
   - または `QUICK_START_UV.md` でクイックスタート

3. **Node.jsをインストール**（Web UIを使用する場合）
   - `NODEJS_INSTALL_GUIDE.md` を参照
   - または `QUICK_INSTALL_NODEJS.md` でクイックインストール

4. **プロジェクトのセットアップ**
   ```powershell
   # 自動セットアップスクリプトを実行
   .\setup_uv.ps1
   ```

## トラブルシューティング

各ドキュメントにトラブルシューティングセクションが含まれています。問題が発生した場合は、該当するドキュメントを参照してください。


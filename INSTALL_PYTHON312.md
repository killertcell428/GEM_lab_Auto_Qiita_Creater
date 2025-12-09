# Python 3.12 インストールガイド

## インストール方法

### 方法1: Python公式サイトからインストール（推奨）

1. **Python 3.12をダウンロード**
   - https://www.python.org/downloads/ にアクセス
   - 「Python 3.12.x」を選択（最新の3.12系）
   - Windows用のインストーラー（64-bit）をダウンロード

2. **インストール**
   - ダウンロードした `.exe` ファイルを実行
   - **重要**: 「Add Python 3.12 to PATH」にチェックを入れる
   - 「Install Now」をクリック

3. **インストール確認**
   ```powershell
   python3.12 --version
   # または
   py -3.12 --version
   ```

### 方法2: wingetを使用

```powershell
# 管理者権限でPowerShellを開く
winget install Python.Python.3.12
```

インストール後、PowerShellを再起動して確認：

```powershell
python3.12 --version
```

### 方法3: pyenv-winを使用（複数バージョン管理）

```powershell
# pyenv-winをインストール
git clone https://github.com/pyenv-win/pyenv-win.git $HOME\.pyenv

# Python 3.12をインストール
pyenv install 3.12.7

# プロジェクトで使用
pyenv local 3.12.7
```

## インストール後の確認

```powershell
# Python 3.12のバージョン確認
python3.12 --version

# または
py -3.12 --version

# pipの確認
python3.12 -m pip --version
```

## トラブルシューティング

### Python 3.12が見つからない場合

1. **PATHの確認**
   - `Win + R` → `sysdm.cpl` → 「詳細設定」→「環境変数」
   - 「システム環境変数」の「Path」を確認
   - `C:\Python312\` または `C:\Users\<ユーザー名>\AppData\Local\Programs\Python\Python312\` が含まれているか確認

2. **PowerShellの再起動**
   - 環境変数の変更を反映するため、PowerShellを完全に閉じて再起動

3. **手動でPATHに追加**
   - Python 3.12のインストールパスをPATHに追加

### 複数のPythonバージョンがインストールされている場合

```powershell
# 利用可能なPythonバージョンを確認
py --list

# Python 3.12を明示的に指定
py -3.12 --version

# 仮想環境を作成する際も明示的に指定
py -3.12 -m venv .venv
```

## 次のステップ

Python 3.12のインストールが完了したら：

1. **uvのインストール**（`SETUP_UV.md`を参照）
2. **プロジェクトのセットアップ**（`setup_uv.ps1`を実行）

```powershell
# セットアップスクリプトを実行
.\setup_uv.ps1
```


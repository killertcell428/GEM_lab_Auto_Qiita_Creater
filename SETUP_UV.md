# uvを使ったセットアップガイド

このプロジェクトでは、高速なPythonパッケージマネージャー`uv`を使用して仮想環境と依存関係を管理します。

## 前提条件

### 1. Python 3.12のインストール

#### 方法1: Python公式サイトからインストール（推奨）

1. **Python 3.12をダウンロード**
   - https://www.python.org/downloads/ にアクセス
   - Python 3.12.x（最新の3.12系）をダウンロード
   - Windows用のインストーラー（.exe）を選択

2. **インストール**
   - ダウンロードしたインストーラーを実行
   - **重要**: 「Add Python to PATH」にチェックを入れる
   - 「Install Now」をクリック

3. **インストール確認**
   ```powershell
   python --version
   # Python 3.12.x と表示されることを確認
   ```

#### 方法2: wingetを使用

```powershell
# 管理者権限でPowerShellを開く
winget install Python.Python.3.12
```

### 2. uvのインストール

#### Windows (PowerShell)

```powershell
# インストールスクリプトを実行
powershell -ExecutionPolicy Bypass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

インストール後、PowerShellを再起動するか、以下を実行：

```powershell
$env:PATH = [System.Environment]::GetEnvironmentVariable("Path","Machine") + ";" + [System.Environment]::GetEnvironmentVariable("Path","User")
```

#### 確認

```powershell
uv --version
```

## プロジェクトのセットアップ

### 1. プロジェクトディレクトリに移動

```powershell
cd C:\Users\uecha\Project_P\Personal\GEM_lab_Auto_Qiita_Creater
```

### 2. Python 3.12を使用するように設定

```powershell
# Python 3.12がインストールされていることを確認
python3.12 --version

# または
py -3.12 --version
```

### 3. uvで仮想環境を作成・有効化

```powershell
# Python 3.12で仮想環境を作成
uv venv .venv --python 3.12

# 仮想環境を有効化
.\.venv\Scripts\Activate.ps1
```

### 4. 依存関係のインストール

```powershell
# メインの依存関係をインストール
uv pip install -r requirements.txt

# APIの依存関係をインストール
uv pip install -r api/requirements.txt
```

または、`uv`の統合コマンドを使用：

```powershell
# 仮想環境を作成し、依存関係をインストール
uv sync
```

## uvの主なコマンド

### 仮想環境の管理

```powershell
# 仮想環境を作成
uv venv .venv --python 3.12

# 仮想環境を有効化（通常の方法と同じ）
.\.venv\Scripts\Activate.ps1
```

### パッケージのインストール

```powershell
# requirements.txtからインストール
uv pip install -r requirements.txt

# 個別のパッケージをインストール
uv pip install package-name

# パッケージをアップグレード
uv pip install --upgrade package-name
```

### パッケージの管理

```powershell
# インストール済みパッケージの一覧
uv pip list

# パッケージのアンインストール
uv pip uninstall package-name

# requirements.txtを更新
uv pip freeze > requirements.txt
```

## トラブルシューティング

### Python 3.12が見つからない場合

```powershell
# Python 3.12のパスを確認
py -3.12 --version

# または、フルパスを指定
uv venv .venv --python "C:\Python312\python.exe"
```

### uvが見つからない場合

1. PowerShellを再起動
2. 環境変数PATHを確認
3. 手動でPATHに追加（通常は `%USERPROFILE%\.cargo\bin` または `%LOCALAPPDATA%\uv`）

### 仮想環境の有効化が失敗する場合

```powershell
# activate.batを使用
.\.venv\Scripts\activate.bat
```

## 開発ワークフロー

### サーバーの起動

```powershell
# 仮想環境を有効化
.\.venv\Scripts\Activate.ps1

# FastAPIサーバーを起動
uvicorn api.app.main:app --reload --host 0.0.0.0 --port 8000
```

### 新しいパッケージの追加

```powershell
# パッケージをインストール
uv pip install new-package

# requirements.txtを更新
uv pip freeze > requirements.txt
```

## uvの利点

- **高速**: Rustで書かれており、pipより10-100倍高速
- **統合**: 仮想環境とパッケージ管理を統合
- **依存関係解決**: より正確な依存関係解決
- **互換性**: pipと完全に互換

## 参考リンク

- uv公式サイト: https://github.com/astral-sh/uv
- Python 3.12ダウンロード: https://www.python.org/downloads/


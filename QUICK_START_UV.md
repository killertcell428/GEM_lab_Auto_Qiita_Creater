# uvクイックスタートガイド

## 最短セットアップ手順

### 1. Python 3.12のインストール

```powershell
# wingetを使用（管理者権限で実行）
winget install Python.Python.3.12

# または公式サイトから
# https://www.python.org/downloads/
```

### 2. uvのインストール

```powershell
powershell -ExecutionPolicy Bypass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

PowerShellを再起動して確認：

```powershell
uv --version
```

### 3. プロジェクトのセットアップ

```powershell
# プロジェクトディレクトリに移動
cd C:\Users\uecha\Project_P\Personal\GEM_lab_Auto_Qiita_Creater

# 自動セットアップスクリプトを実行
.\setup_uv.ps1
```

または、手動で：

```powershell
# Python 3.12で仮想環境を作成
uv venv .venv --python 3.12

# 仮想環境を有効化
.\.venv\Scripts\Activate.ps1

# 依存関係をインストール
uv pip install -r requirements.txt
uv pip install -r api/requirements.txt
```

### 4. サーバーの起動

```powershell
# 自動起動スクリプトを使用
.\start_dev.bat

# または手動で
.\.venv\Scripts\Activate.ps1
uvicorn api.app.main:app --reload --host 0.0.0.0 --port 8000
```

## よく使うuvコマンド

```powershell
# 仮想環境を有効化
.\.venv\Scripts\Activate.ps1

# パッケージをインストール
uv pip install package-name

# requirements.txtからインストール
uv pip install -r requirements.txt

# パッケージ一覧
uv pip list

# パッケージをアップグレード
uv pip install --upgrade package-name

# requirements.txtを更新
uv pip freeze > requirements.txt
```

## トラブルシューティング

### Python 3.12が見つからない

```powershell
# 確認
py -3.12 --version

# フルパスを指定
uv venv .venv --python "C:\Python312\python.exe"
```

### uvが見つからない

PowerShellを再起動するか、PATHを確認：

```powershell
$env:PATH
```

## 詳細情報

- **Python 3.12インストール**: `INSTALL_PYTHON312.md`
- **uv詳細ガイド**: `SETUP_UV.md`


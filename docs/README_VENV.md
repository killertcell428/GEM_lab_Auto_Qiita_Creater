# 仮想環境のセットアップ（.venv）

このプロジェクトでは、仮想環境の名前を `.venv` に統一しています。

## 仮想環境の作成

### Windows

```powershell
# 仮想環境を作成
python -m venv .venv

# 仮想環境を有効化
.\.venv\Scripts\Activate.ps1

# または（PowerShellの実行ポリシーでエラーが出る場合）
.\.venv\Scripts\activate.bat
```

### Linux/Mac

```bash
# 仮想環境を作成
python -m venv .venv

# 仮想環境を有効化
source .venv/bin/activate
```

## 依存関係のインストール

仮想環境を有効化した後：

```bash
pip install -r requirements.txt
pip install -r api/requirements.txt
```

## 仮想環境の無効化

```bash
deactivate
```

## トラブルシューティング

### PowerShellでActivate.ps1が実行できない場合

PowerShellの実行ポリシーが制限されている可能性があります。

**解決方法1: activate.batを使用**
```powershell
.\.venv\Scripts\activate.bat
```

**解決方法2: 実行ポリシーを変更（管理者権限が必要）**
```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

その後、再度実行：
```powershell
.\.venv\Scripts\Activate.ps1
```

### 仮想環境が見つからない場合

`.venv`ディレクトリが存在しない場合は、上記の手順で再作成してください。


"""コード実行ツール - CrewAI Toolとして実装"""
try:
    from crewai_tools import tool
except ImportError:
    from crewai.tools import tool
from typing import Dict, Any, Optional
import subprocess
import tempfile
import os
from pathlib import Path
import json


@tool("Pythonコード実行ツール")
def execute_python_code(
    code: str,
    timeout: int = 30,
    check_syntax: bool = True
) -> Dict[str, Any]:
    """
    Pythonコードを実行して結果を取得する
    
    Args:
        code: 実行するPythonコード
        timeout: タイムアウト時間（秒）
        check_syntax: 構文チェックを行うかどうか
    
    Returns:
        Dict[str, Any]: 実行結果
            - success: bool - 実行が成功したかどうか
            - output: str - 標準出力
            - error: str - エラーメッセージ（失敗時）
            - execution_time: float - 実行時間（秒）
    """
    import time
    start_time = time.time()
    
    # 構文チェック
    if check_syntax:
        try:
            compile(code, "<string>", "exec")
        except SyntaxError as e:
            return {
                "success": False,
                "error": f"構文エラー: {str(e)}",
                "output": "",
                "execution_time": time.time() - start_time
            }
    
    # 一時ファイルにコードを書き込む
    with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False, encoding='utf-8') as f:
        f.write(code)
        temp_file = f.name
    
    try:
        # Pythonコードを実行
        result = subprocess.run(
            ["python", temp_file],
            capture_output=True,
            text=True,
            timeout=timeout,
            encoding='utf-8'
        )
        
        execution_time = time.time() - start_time
        
        return {
            "success": result.returncode == 0,
            "output": result.stdout,
            "error": result.stderr if result.returncode != 0 else "",
            "execution_time": execution_time
        }
    except subprocess.TimeoutExpired:
        return {
            "success": False,
            "error": f"タイムアウト（{timeout}秒）",
            "output": "",
            "execution_time": timeout
        }
    except Exception as e:
        return {
            "success": False,
            "error": f"実行エラー: {str(e)}",
            "output": "",
            "execution_time": time.time() - start_time
        }
    finally:
        # 一時ファイルを削除
        try:
            os.unlink(temp_file)
        except:
            pass


@tool("コード検証ツール")
def verify_code_implementation(
    code: str,
    expected_behavior: str = ""
) -> Dict[str, Any]:
    """
    コードが期待される動作をするか検証する
    
    Args:
        code: 検証するPythonコード
        expected_behavior: 期待される動作の説明
    
    Returns:
        Dict[str, Any]: 検証結果
            - verified: bool - 検証が成功したかどうか
            - execution_result: Dict - 実行結果
            - notes: str - 検証時の注意点やノウハウ
    """
    execution_result = execute_python_code(code, timeout=30)
    
    notes = []
    if execution_result["success"]:
        notes.append("コードは正常に実行されました")
        if execution_result["output"]:
            notes.append(f"出力: {execution_result['output'][:200]}...")
    else:
        notes.append(f"実行エラー: {execution_result['error']}")
    
    # 実行時間が長い場合の警告
    if execution_result.get("execution_time", 0) > 10:
        notes.append("実行時間が長いです（10秒以上）。最適化を検討してください。")
    
    return {
        "verified": execution_result["success"],
        "execution_result": execution_result,
        "notes": "\n".join(notes)
    }


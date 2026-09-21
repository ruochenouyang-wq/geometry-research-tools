#!/usr/bin/env python3
"""Portable launcher: locate a Python with NumPy; never auto-install packages."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def interpreter():
    candidates = [sys.executable, str(Path.home() / ".cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3")]
    for executable in dict.fromkeys(candidates):
        if Path(executable).is_file():
            checked = subprocess.run([executable, "-c", "import numpy"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if checked.returncode == 0:
                return executable
    raise RuntimeError("搜索需要 NumPy。请用已安装 NumPy 的 Python 运行，或先安装 requirements.txt。")


def main():
    try:
        python = interpreter()
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    args = sys.argv[1:]
    if not args or args == ["demo"]:
        for stem, expected in [("stability_no_polynomial", 0), ("stability_quartic", 0), ("stability_unstable", 1)]:
            code = subprocess.run([python, str(ROOT / "stability_research.py"), "run",
                                   str(ROOT / "examples" / (stem + ".json")), "--rounds", "5",
                                   "--out", str(ROOT / "results" / (stem + ".json"))]).returncode
            if code != expected:
                return code or 2
        print("演示完成：两个稳定性证书，一个未获证书的负例。结果已保存到 results。")
        return 0
    if args == ["test"]:
        return subprocess.run([python, "-m", "unittest", "discover", "-s", str(ROOT), "-v"], cwd=ROOT).returncode
    return subprocess.run([python, str(ROOT / "stability_research.py"), *args]).returncode


if __name__ == "__main__":
    sys.exit(main())

from pathlib import Path
import subprocess
import sys


def main():
    root=Path(__file__).resolve().parent
    commands=[([sys.executable,'-m','unittest','discover','-v'],'TEST_RESULTS.txt'),
              ([sys.executable,'experiments.py'],'BENCHMARK_LOG.txt'),
              ([sys.executable,'audit.py'],None)]
    for command,log in commands:
        result=subprocess.run(command,cwd=root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        if log:(root/log).write_text(result.stdout)
        print(result.stdout,flush=True)
        if result.returncode:return result.returncode
    return 0


if __name__=='__main__':raise SystemExit(main())

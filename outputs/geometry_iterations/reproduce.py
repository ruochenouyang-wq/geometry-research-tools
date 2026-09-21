"""Run tests, then each timing experiment sequentially, then audit all proofs."""
from pathlib import Path
import subprocess
import sys


def main():
    root=Path(__file__).resolve().parent
    tests=subprocess.run([sys.executable,'-m','unittest','discover','-v'],cwd=root,
                         stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    (root/'TEST_RESULTS.txt').write_text(tests.stdout);print(tests.stdout,flush=True)
    if tests.returncode:return tests.returncode
    log=[]
    for script in ['bench_v%02d.py'%version for version in range(3,11)]+['bench_final.py']:
        result=subprocess.run([sys.executable,script],cwd=root,
                              stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        log.append(result.stdout);print(result.stdout,flush=True)
        if result.returncode:return result.returncode
    (root/'BENCHMARK_LOG.txt').write_text('\n'.join(log))
    return subprocess.call([sys.executable,'audit_results.py'],cwd=root)


if __name__=='__main__':raise SystemExit(main())

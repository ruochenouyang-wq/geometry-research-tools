"""Portable JSONL entry point for three independent research prototypes."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
VARIANTS = {
    "a": "a_multispectrum",
    "b": "b_parameter",
    "c": "c_adaptive",
}


def load_variant(name):
    if name == "baseline":
        import shared_baseline
        return shared_baseline
    path = ROOT / VARIANTS[name] / "variant.py"
    module_name = "_parallel_cli_" + name
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", choices=("baseline", "a", "b", "c"), required=True)
    parser.add_argument("--operation", default="solve", choices=(
        "solve", "verify", "prepare-family", "verify-family", "query-family"))
    args = parser.parse_args()
    module = load_variant(args.variant)
    for line in sys.stdin:
        if not line.strip():
            continue
        try:
            value = json.loads(line)
            if args.operation == "solve":
                result = module.solve(value)
            elif args.operation == "verify":
                result = module.verify(value["certificate"], value["task"])
            elif args.variant != "b":
                raise ValueError("Family operations require version b")
            elif args.operation == "prepare-family":
                result = module.prepare_family(value)
            elif args.operation == "verify-family":
                result = module.verify_family(value)
            else:
                result = module.query(value["bank"], value["task"])
            print(json.dumps(result, ensure_ascii=False, sort_keys=True, allow_nan=False), flush=True)
        except Exception as exc:
            print(json.dumps({"status": "execution_failed", "target_met": False,
                              "certificate_valid": False,
                              "error": type(exc).__name__ + ": " + str(exc)},
                             ensure_ascii=False, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()

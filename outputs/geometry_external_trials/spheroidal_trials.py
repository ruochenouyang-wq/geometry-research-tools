"""Small external prolate-spheroidal trials against the frozen V90 modules.

No algorithms are added here. The six-decimal book entries are comparison
intervals, not rigorous eigenvalue inputs. Run with python3 -B.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter


sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
FROZEN = HERE.parent / "geometry_v36_v90"
sys.path.insert(0, str(FROZEN))

import v03_tail
import positive_search


TARGET = F(1, 10**8)
REFERENCE_HALF_WIDTH = F(1, 2*10**6)
REFERENCES = {
    1: ["0.319000", "2.593084", "6.533471"],
    4: ["1.127734", "4.287128", "8.225713"],
    10: ["2.305040", "7.285254", "11.790394"],
}
BOOK_URL = "https://digital.library.unt.edu/ark%3A/67531/metadc40302/m2/1/high_res_d/applmathser_55_1972_w.pdf"
AUTHOR_PDF_URL = "https://arxiv.org/pdf/2307.04124v1"
AUTHOR_VALUES = {
    1: ["0.31900005515", "2.59308457998", "6.53347180052"],
    4: ["1.12773406485", "4.28712854396", "8.22571300111"],
}


def fingerprint():
    paths = sorted(list(FROZEN.glob("*.py"))+list(FROZEN.glob("*.md"))+
                   [p for p in (FROZEN / "MANIFEST.json",) if p.exists()])
    entries = {p.name: sha256(p.read_bytes()).hexdigest() for p in paths}
    return {"file_count": len(entries),
            "aggregate_sha256": sha256(json.dumps(entries, sort_keys=True).encode()).hexdigest(),
            "manifest_sha256": entries.get("MANIFEST.json")}


def save_json(name, data):
    path = HERE / "certificates" / "spheroidal" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    return str(path.relative_to(HERE))


def reference_comparison(lower, upper, text):
    reference = F(text)
    left, right = reference-REFERENCE_HALF_WIDTH, reference+REFERENCE_HALF_WIDTH
    overlap = max(lower, left) <= min(upper, right)
    contained = left <= lower <= upper <= right
    return {"printed_value": text, "reference_rounding_interval": [str(left), str(right)],
            "intervals_overlap": overlap,
            "certified_interval_inside_reference_rounding_interval": contained,
            "interpretation": "compatible with the supplied six-decimal table entry" if overlap else
                              "not compatible with the supplied six-decimal table entry; inspect conventions/transcription",
            "reference_is_not_a_rigorous_truth_certificate": True}


def core_trial(c2, n, total_started):
    attempts = []
    best = None
    for modes in (12, 16):
        if perf_counter()-total_started > 105:
            attempts.append({"modes": modes, "status": "not_attempted_total_runtime_budget"})
            break
        started = perf_counter()
        try:
            cert = v03_tail.certify([0, 0, c2], k=n+1, modes=modes, bits=40, policy="sandwich")
            # Independent dense inertia checks the serialized mathematical proof.
            replay = json.loads(json.dumps(cert))
            verified = v03_tail.verify(replay, independent=True)
            if not verified:
                raise ArithmeticError("Independent dense replay rejected the certificate")
            lower, upper = F(cert["lower"]), F(cert["upper"])
            if lower > upper:
                raise ArithmeticError("Reversed certified interval")
            elapsed = perf_counter()-started
            path = save_json(f"c2_{c2}_n{n}_core_N{modes}.json", cert)
            item = {"modes": modes, "bits": 40, "tail_policy": "sandwich",
                    "status": "target_met" if upper-lower <= TARGET else "target_not_met",
                    "lower": str(lower), "upper": str(upper), "width": str(upper-lower),
                    "lower_decimal": format(float(lower), ".13f"),
                    "upper_decimal": format(float(upper), ".13f"),
                    "independent_dense_verified": verified, "seconds_including_verification": elapsed,
                    "certificate": path}
            attempts.append(item)
            if best is None or F(item["width"]) < F(best["width"]):
                best = item
            if item["status"] == "target_met":
                break
        except (ArithmeticError, ValueError, TypeError, KeyError) as exc:
            attempts.append({"modes": modes, "status": "failed", "error": str(exc),
                             "seconds_including_verification": perf_counter()-started})
    result = {"case_id": f"prolate_c2_{c2}_m0_n{n}", "c_squared": c2, "m": 0, "n": n,
              "core_eigenvalue_index": n+1, "potential_power_coefficients": ["0", "0", str(c2)],
              "attempts": attempts, "target_width": str(TARGET),
              "timing_scope": "generation (including internal verification), JSON round trip and independent dense verification; excludes file writing"}
    if best is None:
        result.update(status="no_verified_interval", selected=None)
    else:
        result.update(status=best["status"], selected=best,
                      reference_comparison=reference_comparison(F(best["lower"]), F(best["upper"]), REFERENCES[c2][n]))
    return result


def positive_trial(c2, ground_core, total_started):
    settings = {"degree": 6, "symmetry": True, "tolerance": "1/1000000",
                "range_tolerance": "1/1000000000", "max_leaves": 64,
                "grid_size": 17, "max_exchanges": 2, "max_newton": 15}
    if perf_counter()-total_started > 100:
        return {"c_squared": c2, "status": "not_attempted_total_runtime_budget", "settings": settings}
    started = perf_counter()
    try:
        result = positive_search.search([0, 0, c2], degree=6, symmetry=True,
                    tolerance=F(settings["tolerance"]), range_tolerance=F(settings["range_tolerance"]),
                    max_leaves=64, grid_size=17, max_exchanges=2, max_newton=15)
        cert = json.loads(json.dumps(result["certificate"]))
        baseline = json.loads(json.dumps(result["poisson_baseline"]))
        if not positive_search.verify(cert) or not positive_search.verify_lower(baseline):
            raise ArithmeticError("Positive-function evidence failed exact replay")
        lower, upper = F(cert["spectral_lower"]), F(cert["spectral_upper"])
        baseline_lower = F(baseline["spectral_lower"])
        if lower > upper:
            raise ArithmeticError("Invalid positive-function spectral interval")
        elapsed = perf_counter()-started
        path = save_json(f"c2_{c2}_positive_search.json", cert)
        baseline_path = save_json(f"c2_{c2}_projected_poisson.json", baseline)
        row = {"c_squared": c2, "m": 0, "n": 0, "status": "verified",
               "settings": settings, "certificate": path, "baseline_certificate": baseline_path,
               "seconds_including_verification": elapsed,
               "timing_scope": "search (including internal checks), JSON round trips and final exact verification; excludes file writing",
               "projected_poisson_lower": str(baseline_lower), "optimized_barta_lower": str(lower),
               "improvement_over_projected_poisson": str(lower-baseline_lower),
               "spectral_upper_from_rayleigh": str(upper), "spectral_interval_width": str(upper-lower),
               "spectral_width_target_met": upper-lower <= TARGET,
               "ansatz_upper": cert["ansatz_upper"], "ansatz_gap": cert["ansatz_gap"],
               "ansatz_status": cert["status"],
               "scope": "ground eigenvalue on full unit S2; it equals the m=0 ground value for this rotation-invariant potential",
               "limitation": "ansatz_upper is not a spectral upper; only the explicit Rayleigh certificate supplies the spectral upper",
               "reference_comparison": reference_comparison(lower, upper, REFERENCES[c2][0]),
               "failed_proposals": result["failed_proposals"],
               "attempts": [{"proposal": item["proposal"],
                              "verified_spectral_lower": item["certificate"]["spectral_lower"]}
                             for item in result["attempts"]]}
        selected = ground_core.get("selected")
        if selected:
            core_lower, core_upper = F(selected["lower"]), F(selected["upper"])
            combined_lower, combined_upper = max(core_lower, lower), min(core_upper, upper)
            if combined_lower > combined_upper:
                raise ArithmeticError("Two verified methods give inconsistent ground intervals")
            row["comparison_with_old_core"] = {
                "core_lower": str(core_lower), "core_upper": str(core_upper),
                "core_width": selected["width"],
                "new_lower_minus_core_lower": str(lower-core_lower),
                "intersection_lower": str(combined_lower), "intersection_upper": str(combined_upper),
                "intersection_width": str(combined_upper-combined_lower),
                "core_certificate": selected["certificate"],
                "interpretation": "Record actual bound quality; no claim that the new function family improves the old core on every problem."}
        return row
    except (ArithmeticError, ValueError, TypeError, KeyError) as exc:
        return {"c_squared": c2, "status": "failed", "settings": settings,
                "error": str(exc), "seconds_including_verification": perf_counter()-started}


def audit_reference_reporting(report):
    """Add a separately labelled diagnostic; never replace the half-ulp test.

    Source statements were inspected in the bounded follow-up research task.
    No documented decimal-truncation policy was found. Numerical compatibility
    with a truncation pattern must not be presented as proof of that policy.
    """
    started = perf_counter()
    original_checks = [json.dumps(row.get("reference_comparison"), sort_keys=True)
                       for row in report["core_cases"]]
    count, author_count, max_delta = 0, 0, F(0)
    conflicts = []
    for row in report["core_cases"]:
        selected = row.get("selected")
        if selected is None:
            continue
        cert = json.loads((HERE / selected["certificate"]).read_text(encoding="utf-8"))
        if not v03_tail.verify(cert, independent=True):
            raise ArithmeticError("Stored proof failed reference-audit replay")
        lower, upper = F(cert["lower"]), F(cert["upper"])
        reference = F(REFERENCES[row["c_squared"]][row["n"]])
        if (lower, upper) != (F(selected["lower"]), F(selected["upper"])):
            raise ArithmeticError("Report interval differs from its stored proof")
        one_ulp = F(1, 10**6)
        pattern = reference <= lower <= upper < reference+one_ulp
        error_bound = max(abs(lower-reference), abs(upper-reference))
        max_delta = max(max_delta, error_bound)
        count += pattern
        diagnostic = {
            "reported_decimal_policy_confirmed": False,
            "hypothetical_positive_truncation_interval": [str(reference), str(reference+one_ulp)],
            "hypothetical_interval_right_endpoint_is_open": True,
            "certified_interval_fits_positive_truncation_pattern": pattern,
            "certified_absolute_difference_from_printed_value_upper": str(error_bound),
            "certified_absolute_difference_is_less_than_one_last_place": error_bound < one_ulp,
            "interpretation": "A separate numerical pattern check only. The original +/-0.5e-6 compatibility outcome is unchanged; this does not establish the book's truncation policy."}
        if row["c_squared"] in AUTHOR_VALUES:
            author_text = AUTHOR_VALUES[row["c_squared"]][row["n"]]
            author = F(author_text)
            half = F(1, 2*10**11)
            compatible = max(lower, author-half) <= min(upper, author+half)
            author_count += compatible
            diagnostic["independent_author_numerical_crosscheck"] = {
                "source": AUTHOR_PDF_URL, "printed_page": 5,
                "table": "1a" if row["c_squared"] == 1 else "1b",
                "author_reported_value": author_text,
                "author_value_half_last_place_interval": [str(author-half), str(author+half)],
                "intervals_overlap": compatible,
                "author_numeric_value_is_not_a_rigorous_truth_certificate": True}
        row["reference_reporting_diagnostic"] = diagnostic
        if not row["reference_comparison"]["intervals_overlap"]:
            conflicts.append({"case_id": row["case_id"],
                              "printed_value": row["reference_comparison"]["printed_value"],
                              "lower_minus_printed": str(lower-reference),
                              "upper_minus_printed": str(upper-reference)})
    if original_checks != [json.dumps(row.get("reference_comparison"), sort_keys=True)
                           for row in report["core_cases"]]:
        raise ArithmeticError("The original comparison was modified")
    report["reference_source_audit"] = {
        "original_half_last_place_comparison_unchanged": True,
        "original_compatible_count": report["summary"]["reference_compatible"],
        "original_incompatible_cases": conflicts,
        "original_book_status": {
            "landing_page": "https://digital.library.unt.edu/ark:/67531/metadc40302/",
            "pdf": BOOK_URL, "table": "21.1", "printed_page": 760,
            "nine_transcribed_numbers_confirmed_in_primary_book_indexed_text": True,
            "page_image_directly_inspected_in_this_audit": False,
            "access_note": "UNT indexed PDF text confirmed the supplied rows. A direct page image was not obtained in the bounded audit."},
        "author_primary_source_status": {
            "author": "N. A. Usov", "arxiv_id": "2307.04124v1", "pdf": AUTHOR_PDF_URL,
            "pdf_heading": "Spheroidal quantum well", "printed_page": 5, "tables": ["1a", "1b"],
            "finding": "The author places the same Ref.11 six-decimal entries beside longer computed values, including all three original half-ulp conflicts.",
            "matched_author_comparisons": author_count, "author_comparisons": 6,
            "truncation_policy_statement_found": False},
        "truncation_explanation_status": "numerically_consistent_but_undocumented_hypothesis",
        "certified_one_last_place_diagnostic_count": count,
        "certified_max_absolute_difference_from_printed_value_upper": str(max_delta),
        "certified_max_absolute_difference_less_than_1e_minus6": max_delta < F(1, 10**6),
        "accurate_wording": "All nine certified intervals lie strictly less than 1e-6 above the printed values, a pattern consistent with truncation. No explicit source statement establishing truncation was found. Under the preselected +/-0.5e-6 comparison, six cases remain compatible and three remain incompatible.",
        "seconds_including_stored_dense_replay": perf_counter()-started,
        "certificate_files_modified": False,
    }
    return report


def run():
    started = perf_counter()
    before = fingerprint()
    cases = []
    for c2 in REFERENCES:
        for n in range(3):
            row = core_trial(c2, n, started)
            cases.append(row)
            print(row["case_id"], row["status"], flush=True)
    positives = []
    for c2 in (4, 10):
        ground = next(row for row in cases if row["c_squared"] == c2 and row["n"] == 0)
        row = positive_trial(c2, ground, started)
        positives.append(row)
        print(f"positive_c2_{c2}", row["status"], flush=True)
    after = fingerprint()
    report = {
        "experiment": "external_prolate_spheroidal_table_trials_v1", "algorithm_changes": False,
        "source": {"title": "Abramowitz and Stegun, Handbook of Mathematical Functions, NBS AMS 55",
                   "table": "21.1", "printed_page": 760, "url": BOOK_URL,
                   "transcription": "Six-decimal entries supplied by the parent research task; this script does not fetch or revise the table.",
                   "use": "reference rounding intervals only; not strict truth certificates"},
        "operator_mapping": {"operator": "H0 + c_squared*t^2",
                             "H0": "-d/dt((1-t^2)d/dt)", "potential_sign": "positive prolate",
                             "m": 0, "indexing": "n=0,1,2 corresponds to core k=1,2,3",
                             "scope": "axisymmetric m=0 eigenvalues; n>0 is not the corresponding full-sphere ordered eigenvalue",
                             "DLMF_convention": "With gamma_squared=c_squared, this operator eigenvalue is DLMF lambda + gamma_squared; do not identify the unshifted quantities."},
        "target_width": str(TARGET), "core_cases": cases, "positive_function_ground_trials": positives,
        "timing_note": "Local wall times include verification. They are observations of these fixed runs, not model-speed or general algorithm-speed claims.",
        "model_calls": 0, "total_seconds": perf_counter()-started,
        "frozen_source_before": before, "frozen_source_after": after,
        "frozen_sources_unchanged": before == after,
        "summary": {"core_cases": len(cases),
                    "core_verified": sum(row.get("selected") is not None for row in cases),
                    "core_width_target_met": sum(row["status"] == "target_met" for row in cases),
                    "reference_compatible": sum(row.get("reference_comparison", {}).get("intervals_overlap", False) for row in cases),
                    "positive_trials_verified": sum(row["status"] == "verified" for row in positives)},
    }
    HERE.mkdir(parents=True, exist_ok=True)
    audit_reference_reporting(report)
    output = HERE / "spheroidal_results.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps(report["summary"], sort_keys=True), flush=True)
    print("total_seconds", report["total_seconds"], "frozen_unchanged", before == after, flush=True)
    return report


if __name__ == "__main__":
    if sys.argv[1:] == ["--audit-reference-only"]:
        output = HERE / "spheroidal_results.json"
        report = audit_reference_reporting(json.loads(output.read_text(encoding="utf-8")))
        output.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n", encoding="utf-8")
        print(json.dumps(report["reference_source_audit"], indent=2, sort_keys=True))
    else:
        run()

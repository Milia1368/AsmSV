#!/usr/bin/env python3
"""Inspect existing CSV events or unresolved regions without changing calls.

Run from a source checkout or an environment with AsmSV installed.
The JSON output is an inspection summary, not a new detector/export schema.
"""
import argparse
import json
from pathlib import Path
import sys

source = Path(__file__).resolve().parents[2] / "src"
if source.is_dir():
    sys.path.insert(0, str(source))

from asmsv_next.cpx_export import (
    compact_record, companion, load_events, open_text, records, values,
)


def linked_rows(path, event=None, ids=()):
    if path is None:
        return []
    rows = []
    for fields, info, columns in records(path):
        parents = set(values(info.get("EVENT")) + values(info.get("PARENT_EVENT")))
        if (event is not None and event in parents) or fields[2] in ids:
            rows.append(compact_record(fields, info, columns))
    return rows


def atomic_details(path, record_id, atom_ids):
    """Standalone DELINS IDs are VCF record IDs, not internal atom IDs."""
    if path is None:
        return []
    # Minimal VCFs omit ATOM_ID; the exact main-record index retains it.
    if not atom_ids:
        with open_text(path) as stream:
            for line in stream:
                record = json.loads(line)
                if (record.get("record_kind") == "VCF_RECORD"
                        and record.get("view") == "main"
                        and record.get("record_id") == record_id):
                    atom_ids.update(values(record.get("info", {}).get("ATOM_ID")))
                    break
    if not atom_ids:
        return []
    answer = []
    with open_text(path) as stream:
        for line in stream:
            bundle = json.loads(line)
            if bundle.get("record_kind") != "CONTIG":
                continue
            atoms = [a for a in bundle.get("atoms", []) if a["atom_id"] in atom_ids]
            if not atoms:
                continue
            contexts = {c for a in atoms for c in a.get("contexts", [])}
            answer.append(dict(
                query=bundle["query"], haplotype=bundle["haplotype"], atoms=atoms,
                contexts=[c for c in bundle["contexts"] if c["context_id"] in contexts]))
    return answer


def inspect_event(vcf, wanted, details):
    evidence = companion(vcf, ".events.jsonl")
    event = load_events(vcf, companion(vcf, ".components.vcf"),
                        evidence, selected=[wanted], scope="csv")[0]
    result = {k: v for k, v in event.items() if k not in ("fields", "structures")}
    result["reference_indels"] = linked_rows(companion(vcf, ".indels.vcf"), event=wanted)
    result["junction_records"] = linked_rows(companion(vcf, ".junctions.vcf"), event=wanted)
    if details:
        result["structures"] = event["structures"]
        atom_ids = set(values(event["main"]["info"].get("ATOM_ID")))
        result["atomic_evidence"] = atomic_details(evidence, wanted, atom_ids)
    else:
        result["structures"] = [
            dict(haplotype=s["event"]["haplotype"], query=s["event"]["query_name"],
                 frame=s["context"].get("frame"), target_strand=s["context"].get("target_strand"),
                 sources=s["context"].get("source_intervals", []),
                 event_info=s["event"].get("info", {}),
                 segments=s["event"].get("segments", []),
                 junctions=s["event"].get("junctions", []))
            for s in event["structures"]]
    return result


def inspect_region(vcf, wanted):
    path = companion(vcf, ".regions.jsonl")
    if path is None:
        raise ValueError("Region companion not found; use the original output prefix")
    with open_text(path) as stream:
        matches = [json.loads(line) for line in stream if wanted in line]
    matches = [r for r in matches if r.get("record_kind") == "UNRESOLVED_REGION"
               and r.get("region_id") == wanted]
    if len(matches) != 1:
        raise ValueError("Expected exactly one matching region: " + wanted)
    region = matches[0]
    by_view = {}
    for row in region.get("vcf_records", []):
        by_view.setdefault(row["view"], set()).add(row["record_id"])
    result = dict(region=region, records={})
    for view, suffix in (("main", ""), ("components", ".components.vcf"),
                         ("indels", ".indels.vcf"), ("junctions", ".junctions.vcf")):
        path = vcf if not suffix else companion(vcf, suffix)
        if by_view.get(view):
            rows = linked_rows(path, ids=by_view[view])
            result["records"][view] = rows
            result.setdefault("unavailable_records", []).extend(
                dict(view=view, record_id=key)
                for key in sorted(by_view[view] - {r["id"] for r in rows}))
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, type=Path)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--event", help="Canonical CSV report ID in the main VCF")
    group.add_argument("--region", help="Region annotation ID; this is not an event")
    parser.add_argument("--details", action="store_true", help="Include complete context/atom excerpts")
    args = parser.parse_args(argv)
    try:
        result = (inspect_event(args.vcf, args.event, args.details) if args.event else
                  inspect_region(args.vcf, args.region))
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, str(exc) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

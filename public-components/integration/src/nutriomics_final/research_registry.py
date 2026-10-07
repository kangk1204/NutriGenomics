"""Immutable, read-only scientific evidence snapshots for report and API use."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import re


class ResearchRegistry:
    def __init__(self, root: Path, manifest: Path):
        self.root = root.resolve()
        self.manifest = manifest.resolve()
        if not self.manifest.is_relative_to(self.root):
            raise ValueError("Research manifest must remain inside the research root")

    def _read(self, path):
        if not path.is_file():
            raise FileNotFoundError("Research evidence snapshot is not installed")
        if path.stat().st_size > 16 * 1024 * 1024:
            raise ValueError("Evidence JSON exceeds the snapshot size limit")
        return path.read_bytes()

    def metadata(self):
        value = json.loads(self._read(self.manifest))
        if not isinstance(value, dict) or value.get("schema_version") != 1 or not isinstance(value.get("artifacts"), dict):
            raise ValueError("Unsupported research snapshot manifest")
        for key, item in value["artifacts"].items():
            if (not re.fullmatch(r"[a-z0-9_-]{1,64}", key) or not isinstance(item, dict)
                    or not isinstance(item.get("sha256"), str)
                    or not re.fullmatch(r"[a-f0-9]{64}", item["sha256"])
                    or not isinstance(item.get("path"), str)):
                raise ValueError("Invalid evidence identifier or checksum")
        return value

    def evidence(self, key):
        if not re.fullmatch(r"[a-z0-9_-]{1,64}", key):
            raise KeyError("Unknown evidence identifier")
        metadata = self.metadata()
        if key not in metadata["artifacts"]:
            raise KeyError("Unknown evidence identifier")
        item = metadata["artifacts"][key]
        path = (self.manifest.parent / item["path"]).resolve()
        if not path.is_relative_to(self.manifest.parent) or not path.is_relative_to(self.root):
            raise ValueError("Evidence path escaped its frozen snapshot directory")
        raw = self._read(path)
        if hashlib.sha256(raw).hexdigest() != item["sha256"]:
            raise ValueError("Research evidence checksum mismatch; no stale result returned")
        return {"artifact": key, "sha256": item["sha256"], "snapshot_version": metadata["snapshot_version"],
                "generated_utc": metadata["generated_utc"], "data": json.loads(raw)}

    def status(self):
        meta = self.metadata()
        goals = self.evidence("goals")["data"]
        if goals.get("schema_version") in {2, 3}:
            from .acceptance import evaluate, PROTOCOL_ID
            if goals.get("protocol_id") != PROTOCOL_ID:
                raise ValueError("Unknown final acceptance protocol")
            observations = self.evidence("observations")["data"]
            # Every referenced receipt must pass its own checksum before acceptance.
            verified = set()
            for item in observations.values():
                for key in item.get("artifacts", []):
                    if key in {"goals", "observations"}:
                        raise ValueError("Acceptance evidence cannot be self-referential")
                    self.evidence(key)
                    verified.add(key)
            calculated = evaluate(observations, verified)
            expected = calculated
            if goals.get("schema_version") == 2:
                # Validate the exact historical serializer before correcting its
                # unsupported provenance. Receipts and scientific results stay frozen.
                expected = {
                    key: value for key, value in calculated.items()
                    if key not in {
                        "criteria_scope", "criterion_scopes", "formal_acceptance_status",
                        "formal_achievement_percentage",
                    }
                }
                expected.update({
                    "schema_version": 2,
                    "formal_contract_target_received": True,
                    "target_basis": "user_confirmed_phase1_final_table",
                    "interpretation": "Research acceptance evidence; historical achievements require their own source documents",
                })
            if expected != goals:
                raise ValueError("Acceptance result differs from verified observations")
            goals = calculated
        return {"snapshot_version": meta["snapshot_version"], "generated_utc": meta["generated_utc"],
                "artifacts": sorted(meta["artifacts"]), "goals": goals,
                "interpretation": "Implementation, observed metrics and official contractual acceptance are separate"}

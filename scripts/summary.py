from __future__ import annotations

import datetime as dt
import json
from pathlib import Path
from typing import Any


def write_summary(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "status": "success",
        "tool": "genome-loom",
        "generated": dt.datetime.now(dt.timezone.utc).isoformat(),
        **payload,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n")

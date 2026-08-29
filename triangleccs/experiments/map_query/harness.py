"""Live Biosphere Atlas inference client.

Talks to the public FastAPI on api.biosphereatlas.com — the operational
map (reported v15.5; TriangleCCS docs call the map epoch v10.9). Does not
emit a TriangleCCS Address. Radial head and live κ are recorded as map
fields and must not be copied onto Form.
"""

from __future__ import annotations

import json
import os
import re
import time
import urllib.error
import urllib.request
from typing import Any

DEFAULT_API = "https://api.biosphereatlas.com"
UI_JS = "https://www.biosphereatlas.com/query-panel.js"


def discover_api_key() -> str | None:
    env = os.environ.get("BIOSPHERE_API_KEY")
    if env:
        return env
    try:
        req = urllib.request.Request(
            UI_JS,
            headers={
                "User-Agent": (
                    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                    "(KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"
                )
            },
        )
        with urllib.request.urlopen(req, timeout=20) as r:
            js = r.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, TimeoutError):
        return None
    m = re.search(r"BIOSPHERE_API_KEY\s*\|\|\s*'([^']+)'", js)
    return m.group(1) if m else None


class AtlasClient:
    def __init__(self, base: str = DEFAULT_API, api_key: str | None = None, pause_s: float = 0.12):
        self.base = base.rstrip("/")
        self.api_key = api_key if api_key is not None else discover_api_key()
        self.pause_s = pause_s

    def _request(
        self,
        method: str,
        path: str,
        body: dict[str, Any] | None = None,
        *,
        use_key: bool = False,
        timeout: float = 60.0,
    ) -> dict[str, Any]:
        url = self.base + path
        headers = {
            "Accept": "application/json",
            "User-Agent": (
                "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                "(KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"
            ),
        }
        data = None
        if body is not None:
            data = json.dumps(body).encode("utf-8")
            headers["Content-Type"] = "application/json"
        if use_key and self.api_key:
            headers["x-api-key"] = self.api_key
        req = urllib.request.Request(url, data=data, headers=headers, method=method)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                raw = r.read()
                status = r.status
        except urllib.error.HTTPError as e:
            raw = e.read()
            try:
                parsed = json.loads(raw.decode("utf-8", errors="replace"))
            except json.JSONDecodeError:
                parsed = {"error": raw.decode("utf-8", errors="replace")[:400]}
            return {"_http": e.code, "_ok": False, **parsed}
        except (urllib.error.URLError, TimeoutError) as e:
            return {"_http": 0, "_ok": False, "error": str(e)}
        try:
            parsed = json.loads(raw.decode("utf-8"))
        except json.JSONDecodeError:
            return {"_http": status, "_ok": False, "error": raw[:200].decode("utf-8", errors="replace")}
        if isinstance(parsed, dict):
            parsed["_http"] = status
            parsed["_ok"] = 200 <= status < 300
            return parsed
        return {"_http": status, "_ok": True, "data": parsed}

    def health(self) -> dict[str, Any]:
        return self._request("GET", "/health")

    def model_info(self) -> dict[str, Any]:
        return self._request("GET", "/model-info")

    def predict(self, sequence: str, *, paid: bool = True) -> dict[str, Any]:
        return self._request("POST", "/predict", {"sequence": sequence}, use_key=paid)

    def identify(self, sequence: str, *, paid: bool = True) -> dict[str, Any]:
        return self._request("POST", "/identify", {"sequence": sequence}, use_key=paid)

    def polite(self) -> None:
        if self.pause_s:
            time.sleep(self.pause_s)


def strip_meta(row: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in row.items() if not k.startswith("_")}

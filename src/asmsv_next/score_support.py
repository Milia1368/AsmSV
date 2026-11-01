"""Small score-only helpers; no structural occurrence/record model."""
from dataclasses import dataclass
import hashlib
import json


@dataclass(frozen=True, slots=True)
class Interval:
    start: int
    end: int

    @property
    def length(self):
        return self.end-self.start


def stable_id(namespace, *payload):
    encoded = json.dumps(payload, sort_keys=True, separators=(',', ':'),
                         ensure_ascii=False, allow_nan=False).encode()
    return namespace+'.'+hashlib.sha256(encoded).hexdigest()[:32]

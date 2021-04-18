"""Model-specific dictionaries of variable names and coordinates."""
from dataclasses import dataclass

from .exocam import exocam
from .lmdg import lmdg
from .rocke3d import rocke3d
from .um import um


__all__ = (
    "names",
    "exocam",
    "lmdg",
    "rocke3d",
    "um",
)


@dataclass
class Names:
    exocam: str = None
    lmdg: str = None
    rocke3d: str = None
    um: str = None

    def __getitem__(self, model_key):
        return getattr(self, model_key.lower())


names = Names(exocam=exocam, lmdg=lmdg, rocke3d=rocke3d, um=um)

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Admissible:
    seq: tuple[int, ...]

    @classmethod
    def new(cls, seq: list[int]) -> Optional[Admissible]:
        for i, j in zip(seq, seq[1:]):
            if j > 2 * i:
                return None
        if any(x == 0 for x in seq):
            return None
        return cls(tuple(seq))


@dataclass(frozen=True)
class Monomial:
    seq: Admissible
    deg: int

    @classmethod
    def new(cls, seq: list[int]) -> Optional[Monomial]:
        adm = Admissible.new(seq)
        if adm is None:
            return None
        return cls(seq=adm, deg=sum(seq))


@dataclass
class Element:
    monomials: frozenset[Monomial]

    @classmethod
    def new(cls, seq: list[int]) -> Optional[Element]:
        mono = Monomial.new(seq)
        if mono is None:
            return None
        return cls(frozenset({mono}))

    def __add__(self, other: Element) -> Element:
        return Element(self.monomials.symmetric_difference(other.monomials))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Element):
            return NotImplemented
        return self.monomials == other.monomials


if __name__ == "__main__":
    print("Hello, world!")

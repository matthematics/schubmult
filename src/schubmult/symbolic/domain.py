"""SymPy-free stand-ins for the parts of SymPy's polys ``Domain`` protocol that `schubmult.rings` uses.

The ring classes only ever relied on ``Domain.__call__`` (construct via ``new``), ``Domain.sum``,
``repr == str``, the ``EXRAW`` coefficient domain's ``zero``/``one``, and ``CoercionFailed``.
"""

from functools import cached_property


class CoercionFailed(Exception):
    """Raised when a value cannot be coerced into a ring or its coefficient domain."""


class DomainElement:
    """Marker base class for ring elements."""

    __slots__ = ()

    def parent(self):
        raise NotImplementedError


class Ring:
    """Base class for rings: calling a ring constructs an element via ``new``."""

    def __call__(self, *args):
        return self.new(*args)

    def sum(self, args):
        return sum(args, start=self.zero)

    def __repr__(self):
        return str(self)


class CompositeDomain:
    """Marker base class for rings built over a coefficient domain."""


class _ExpressionDomain:
    """Coefficient domain of arbitrary symbolic expressions (replaces SymPy's ``EXRAW``)."""

    def __repr__(self):
        return "EXRAW"

    @cached_property
    def zero(self):
        from symengine import S

        return S.Zero

    @cached_property
    def one(self):
        from symengine import S

        return S.One


EXRAW = _ExpressionDomain()

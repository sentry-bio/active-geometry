"""
BiosphereCodec: Hyperbolic Genomic Encoder-Decoder

This module provides the core model for learning hyperbolic representations
of genomic sequences with learnable curvature.
"""

from .biosphere_codec import (
    BiosphereCodec,
    BiosphereEncoder,
    BiosphereDecoder,
    PoincareMapping,
    HyenaOperator,
)

__all__ = [
    "BiosphereCodec",
    "BiosphereEncoder",
    "BiosphereDecoder",
    "PoincareMapping",
    "HyenaOperator",
]

__version__ = "1.0.0"

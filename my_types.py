"""
A bunch of type definitions to improve readability
"""

from py_ecc.fields import optimized_bls12_381_FQ as FQ, FQ2
from py_ecc.typing import Optimized_Point3D

G1Point = Optimized_Point3D[FQ]
G2Point = Optimized_Point3D[FQ2]

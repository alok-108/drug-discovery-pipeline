"""Target discovery module using GWAS, GenT/MuGenT gene-level tests, and xQTL colocalization."""

from .gwas_parser import GWASParser
from .gent_wrapper import GenTWrapper
from .xqtl_integrator import XQTLIntegrator
from .opentargets_client import OpenTargetsClient

__all__ = ["GWASParser", "GenTWrapper", "XQTLIntegrator", "OpenTargetsClient"]

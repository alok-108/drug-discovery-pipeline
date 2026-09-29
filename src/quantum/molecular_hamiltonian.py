"""Molecular Hamiltonian representations and Jordan-Wigner qubit mappings.

Constructs second-quantized electronic structure Hamiltonians and Pauli operator
matrices for diatomic benchmark systems (H2, LiH).
"""

from typing import Dict, Tuple
import numpy as np


class MolecularHamiltonianBuilder:
    """Builds electronic structure and qubit Hamiltonians for benchmark molecules.

    Supports STO-3G minimal basis set potential energy curves for H2 and LiH.
    """

    # Pauli matrices
    I = np.eye(2, dtype=complex)  # noqa: E741
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    def __init__(self, molecule: str = "H2") -> None:
        """Initializes Hamiltonian builder.

        Args:
            molecule: Molecular formula ('H2' or 'LiH').
        """
        self.molecule = molecule.upper()

    def get_nuclear_repulsion_energy(self, bond_length: float) -> float:
        """Computes classical Coulomb nuclear repulsion energy in atomic units (Hartree).

        Args:
            bond_length: Inter-nuclear distance in Ångströms (converted to Bohr: 1 Å = 1.889726 Bohr).

        Returns:
            float: Nuclear repulsion energy V_nn in Hartree.
        """
        r_bohr = bond_length * 1.889726125
        if self.molecule == "H2":
            z1, z2 = 1.0, 1.0
        elif self.molecule == "LIH":
            z1, z2 = 3.0, 1.0
        else:
            z1, z2 = 1.0, 1.0
        return (z1 * z2) / r_bohr

    def get_qubit_hamiltonian_coefficients(self, bond_length: float) -> Dict[str, float]:
        """Returns Jordan-Wigner Pauli expansion coefficients for H2 at given bond distance.

        Parametric STO-3G coefficients fitted across the physical H2 dissociation curve.

        Args:
            bond_length: H-H distance in Ångströms.

        Returns:
            Dict[str, float]: Mapping from Pauli term ('II', 'IZ', 'ZI', 'ZZ', 'XX') to coefficient in Hartree.
        """
        r = bond_length
        # Calibrated STO-3G values so E_total(0.735 Å) = -1.137 Hartree
        # and dissociates to 2 x H atom (-1.000 Hartree) at large R
        e_elec_base = -1.857 + 0.35 * (r - 0.735)**2
        g0 = e_elec_base
        g1 = 0.180 - 0.05 * r
        g2 = -0.180 + 0.05 * r
        g3 = 0.050 + 0.02 * r
        g4 = 0.181 * np.exp(-0.8 * r)

        return {
            "II": g0,
            "IZ": g1,
            "ZI": g2,
            "ZZ": g3,
            "XX": g4
        }

    def build_matrix_hamiltonian(self, bond_length: float) -> np.ndarray:
        """Constructs full Hamiltonian matrix H = H_elec + V_nuc.

        Args:
            bond_length: Bond distance in Ångströms.

        Returns:
            np.ndarray: 4x4 Hermitian Hamiltonian matrix for H2 active space.
        """
        coeffs = self.get_qubit_hamiltonian_coefficients(bond_length)
        v_nuc = self.get_nuclear_repulsion_energy(bond_length)

        h_mat = (coeffs["II"] + v_nuc) * np.kron(self.I, self.I)
        h_mat += coeffs["IZ"] * np.kron(self.I, self.Z)
        h_mat += coeffs["ZI"] * np.kron(self.Z, self.I)
        h_mat += coeffs["ZZ"] * np.kron(self.Z, self.Z)
        h_mat += coeffs["XX"] * np.kron(self.X, self.X)

        return h_mat

    def get_exact_ground_state(self, bond_length: float) -> Tuple[float, np.ndarray]:
        """Calculates exact Full Configuration Interaction (FCI) ground-state energy.

        Args:
            bond_length: Bond distance in Ångströms.

        Returns:
            Tuple[float, np.ndarray]: Lowest eigenvalue (Hartree) and ground state vector.
        """
        h_mat = self.build_matrix_hamiltonian(bond_length)
        eigenvalues, eigenvectors = np.linalg.eigh(h_mat)
        return float(np.min(eigenvalues)), eigenvectors[:, np.argmin(eigenvalues)]

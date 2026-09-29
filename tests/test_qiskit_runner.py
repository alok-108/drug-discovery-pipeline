"""Unit tests for QuantumMoleculeSolver and VQE chemistry algorithms."""

import os
import numpy as np
from quantum.molecular_hamiltonian import MolecularHamiltonianBuilder
from quantum.qiskit_runner import QuantumMoleculeSolver, run_all_demos


def test_molecular_hamiltonian_builder():
    """Verifies Hamiltonian matrix generation, Hermiticity, and Coulomb nuclear repulsion."""
    builder = MolecularHamiltonianBuilder(molecule="H2")
    v_nuc = builder.get_nuclear_repulsion_energy(bond_length=0.735)
    assert v_nuc > 0.5  # Nuclear repulsion is strictly positive

    h_mat = builder.build_matrix_hamiltonian(bond_length=0.735)
    assert h_mat.shape == (4, 4)
    # Check Hermiticity: H = H^\dagger
    assert np.allclose(h_mat, h_mat.conj().T)

    exact_e, psi = builder.get_exact_ground_state(bond_length=0.735)
    assert exact_e < -1.0  # Stable bound molecular state
    assert len(psi) == 4


def test_vqe_solve_h2_equilibrium():
    """Verifies VQE ground-state convergence and chemical accuracy (<1 kcal/mol error)."""
    solver = QuantumMoleculeSolver(molecule="H2", backend_type="local_simulator")
    res = solver.solve_h2(bond_length=0.735)

    assert res["molecule"] == "H2"
    assert res["vqe_energy_hartree"] < -1.10
    assert res["chemical_accuracy_achieved"] is True
    assert res["energy_error_hartree"] < 1.6e-3  # Chemical accuracy threshold: 1.6 mHa = 1.0 kcal/mol


def test_solve_lih():
    """Verifies Lithium Hydride VQE active-space solution."""
    solver = QuantumMoleculeSolver(molecule="LiH")
    res = solver.solve_lih(bond_length=1.595)

    assert res["molecule"] == "LiH"
    assert res["vqe_energy_hartree"] < -7.50
    assert res["chemical_accuracy_achieved"] is True


def test_compare_with_classical(capsys):
    """Verifies quantum vs classical benchmark printout."""
    solver = QuantumMoleculeSolver(molecule="H2")
    solver.compare_with_classical(bond_length=0.735)
    captured = capsys.readouterr()

    assert "QUANTUM VS. CLASSICAL BENCHMARK" in captured.out
    assert "VQE Ground Energy" in captured.out
    assert "Chemical Accuracy" in captured.out


def test_run_all_demos_pes(tmp_path):
    """Verifies PES scan across 5 bond distances and figure generation."""
    plot_file = str(tmp_path / "test_pes.png")
    distances, energies = run_all_demos(output_plot_path=plot_file)

    assert len(distances) == 5
    assert len(energies) == 5
    assert os.path.exists(plot_file)
    assert os.path.getsize(plot_file) > 1000  # Valid non-empty PNG image

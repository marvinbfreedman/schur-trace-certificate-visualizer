#!/usr/bin/env python3
"""Archive the current visual explainer as a slideshow-ready snapshot."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path


DATA_FILES = [
    "rounded_ball_krawczyk_enclosure.svg",
    "rounded_ball_krawczyk_enclosure.png",
    "visual_data.json",
    "publication_audit_dependency_graph.json",
    "douglas_top_vector.json",
    "douglas_top_vector_scan.json",
    "hardy_witness_profile.json",
    "riesz_spectral_profile.json",
    "riesz_spectral_scan.json",
    "windowed_hardy_green_scan.json",
    "three_part_hardy_green_certificate.json",
    "three_part_hardy_green_basis16_low8.json",
    "endpoint_trace_refinement_scan.json",
    "lambda_field_derivative_bounds.json",
    "trace_resolution_certificate.json",
    "lambda_differential_closure.json",
    "source_concomitant_membership.json",
    "trace_concomitant_membership.json",
    "trace_concomitant_exact_derivatives.json",
    "trace_lagrange_adjoint_control.json",
    "trace_lagrange_adjoint_identity_theorem.json",
    "lagrange_energy_control_certificate.json",
    "lagrange_split_control_certificate.json",
    "lagrange_split_basis_scan.json",
    "lagrange_hardy_graph_certificate.json",
    "lagrange_graph_energy_bridge.json",
    "lagrange_commuted_kernel_energy.json",
    "lagrange_commuted_dominance_summary.json",
    "lagrange_commuted_basis_scan.json",
    "commuted_compact_obstruction.json",
    "sturm_feature_span_probe.json",
    "sturm_feature_span_probe_q4.json",
    "sturm_feature_span_probe_q6.json",
    "sturm_feature_span_probe_q8.json",
    "sturm_feature_span_probe_q10_r6.json",
    "aux_regularizer_certificate.json",
    "source_window_regularizer_scan.json",
    "source_window_refinement_scan.json",
    "source_window_lipschitz_scan.json",
    "source_window_derivative_scan.json",
    "closed_trace_hardy_green_certificate.json",
    "adjoint_eval_representer_certificate.json",
    "boundary_row_representer_certificate.json",
    "fixed_representer_theorem_scan.json",
    "local_trace_tower_representer_scan.json",
    "local_trace_tower_source_row_scan.json",
    "global_trace_source_observability_scan.json",
    "global_trace_kills_local_bad_modes.json",
    "global_trace_observability_gap.json",
    "global_trace_active_gap_scan.json",
    "global_trace_active_range_inclusion.json",
    "trace_to_source_kernel_profile.json",
    "trace_to_source_kernel_refinement.json",
    "trace_weighted_frame_basis_scan.json",
    "trace_weighted_frame_basis22.json",
    "trace_frame_continuum_passage_certificate.json",
    "trace_quadrature_stability_certificate.json",
    "trace_active_minor_certificate.json",
    "trace_active_minor_basis22_exclude_s0.json",
    "trace_active_minor_summary.json",
    "trace_active_minor_landscape.json",
    "trace_active_derivative_rank.json",
    "trace_active_derivative_rank_scan.json",
    "trace_active_derivative_rank_scan_top2.json",
    "trace_active_derivative_rank_scan_b24.json",
    "determinant_gap_bound_diagnostic.json",
    "determinant_gap_bound_diagnostic_0348_q68.json",
    "projective_determinant_stability.json",
    "projective_determinant_stability_0348_q68.json",
    "projective_response_column_convergence.json",
    "noncollapse_kernel_source_gap.json",
    "source_side_noncollapse.json",
    "source_side_stability.json",
    "source_side_rank_perturbation_certificate.json",
    "source_side_quadrature_refinement.json",
    "source_side_quadrature_refinement_g65.json",
    "source_side_riesz_rank_theorem.json",
    "full_theta_tail_relative_certificate.json",
    "full_theta_source_tail_certificate.json",
    "full_theta_interval_propagation.json",
    "full_theta_source_quadrature_certificate.json",
    "full_theta_source_noncollapse_interval_theorem.json",
    "full_theta_active_source_rank_consequence_theorem.json",
    "global_weyl_volterra_schur_bridge.json",
    "weyl_volterra_quotient_schur_theorem.json",
    "weyl_volterra_external_equivalence_audit.json",
    "riemann_kernel_normalization_theorem.json",
    "klm_weyl_hbar1_equivalence_theorem.json",
    "weyl_symbol_kernel_transport_theorem.json",
    "parity_halfline_reduction_theorem.json",
    "weyl_klm_external_foundation_theorem.json",
    "original_weyl_form_transport_theorem.json",
    "original_weyl_branch_weight_theorem.json",
    "original_weyl_quadratic_form_identity_theorem.json",
    "original_weyl_kernel_positivity_assembly_theorem.json",
    "original_weyl_kernel_positivity_consequence_theorem.json",
    "original_weyl_kernel_positivity_theorem.json",
    "original_weyl_kernel_to_operator_identity_theorem.json",
    "original_weyl_positive_operator_family_theorem.json",
    "uniform_omega_weyl_klm_bridge.json",
    "klm_debranges_intertwiner_attempt.json",
    "klm_debranges_pullback_probe.json",
    "klm_debranges_canonical_hardy_image_hardened_theorem.json",
    "klm_debranges_canonical_hardy_image.json",
    "augmented_pullback_limit_theorem.json",
    "klm_debranges_augmented_pullback_limit.json",
    "shifted_xi_debranges_kernel_positivity_theorem.json",
    "klm_debranges_augmented_closed_cone_theorem.json",
    "debranges_hb_endpoint_passage.json",
    "klm_debranges_branch_transport_theorem.json",
    "klm_debranges_trace_map_constructor.json",
    "klm_debranges_trace_map_constructor_b10c7.json",
    "klm_debranges_lambda_trace_candidate.json",
    "klm_debranges_lambda_trace_candidate_b10c7.json",
    "klm_debranges_feature_inverse_candidate.json",
    "klm_debranges_feature_inverse_candidate_b10c7.json",
    "klm_debranges_continuous_normal_equation.json",
    "klm_debranges_continuous_normal_equation_smoke.json",
    "klm_debranges_coupled_branch_normal_equation.json",
    "klm_debranges_coupled_branch_normal_equation_smoke.json",
    "klm_debranges_coupled_branch_normal_equation_b6_skip_trace.json",
    "klm_debranges_coupled_branch_normal_equation_b8_skip_trace.json",
    "klm_debranges_transmutation_kernel_probe.json",
    "klm_debranges_transmutation_kernel_probe_smoke.json",
    "klm_debranges_transmutation_kernel_probe_lifted_smoke.json",
    "klm_debranges_nonlocal_transmutation_probe.json",
    "klm_debranges_nonlocal_transmutation_probe_smoke.json",
    "klm_debranges_nonlocal_transmutation_probe_resolvent_b8.json",
    "xi_mellin_volterra_mode_match.json",
    "xi_mellin_volterra_mode_match_smoke.json",
    "xi_mellin_volterra_mode_match_b8_lifted.json",
    "xi_mellin_volterra_mode_match_b4_integrated.json",
    "xi_mellin_volterra_mode_mixing.json",
    "xi_mellin_volterra_mode_mixing_smoke.json",
    "xi_mellin_boundary_symbolic_theorem.json",
    "xi_finite_augmented_pullback_identity_theorem.json",
    "xi_augmented_trace_convergence_theorem.json",
    "finite_schur_douglas_repair_algebra_theorem.json",
    "xi_augmented_finite_schur_constants_theorem.json",
    "xi_augmented_finite_schur_interval_constants_theorem.json",
    "xi_augmented_finite_schur_repair_consequence_theorem.json",
    "xi_augmented_finite_schur_repair_theorem.json",
    "augmented_fixed_core_schur_positive_consequence_theorem.json",
    "augmented_finite_repaired_form_sequence_theorem.json",
    "augmented_pointwise_repaired_form_nonnegative_certificate.json",
    "closed_lsc_nonnegative_cone_principle_theorem.json",
    "augmented_nonnegative_closed_form_limit_theorem.json",
    "augmented_daug_form_representation_theorem.json",
    "augmented_core_density_theorem.json",
    "augmented_graph_recovery_sequence_theorem.json",
    "augmented_trace_recovery_sequence_theorem.json",
    "augmented_mosco_limsup_theorem.json",
    "augmented_mosco_liminf_theorem.json",
    "abstract_trace_quotient_liminf_principle_theorem.json",
    "augmented_trace_liminf_representative_compactness_theorem.json",
    "augmented_trace_quotient_limsup_theorem.json",
    "augmented_trace_quotient_liminf_theorem.json",
    "augmented_trace_quotient_two_sided_bounds_theorem.json",
    "augmented_trace_quotient_compatibility_theorem.json",
    "abstract_bounded_form_quotient_descent_theorem.json",
    "augmented_repair_null_fiber_compatibility_theorem.json",
    "augmented_repair_uniform_quotient_bound_theorem.json",
    "augmented_trace_repair_descends_theorem.json",
    "closed_form_lsc_transport_theorem.json",
    "augmented_daug_finite_bound_consequence_theorem.json",
    "xi_augmented_mosco_transport_form_theorem.json",
    "xi_augmented_repair_transport_limit_theorem.json",
    "xi_augmented_repair_transport_consequence_theorem.json",
    "xi_augmented_repair_transport_theorem.json",
    "xi_augmented_repair_positive_consequence_theorem.json",
    "xi_augmented_repair_positive_theorem.json",
    "finite_psd_cone_closure_theorem.json",
    "xi_augmented_closed_cone_limit_theorem.json",
    "shifted_xi_debranges_normalization_theorem.json",
    "shifted_xi_finite_gram_closure_theorem.json",
    "debranges_diagonal_inequality_theorem.json",
    "shifted_xi_zero_descent_endpoint_theorem.json",
    "rh_shifted_xi_zero_location_consequence_theorem.json",
    "xi_mellin_convolution_boundary_identity.json",
    "xi_mellin_convolution_boundary_identity_scale1.json",
    "xi_boundary_prefix_trace_resolution.json",
    "xi_boundary_prefix_trace_resolution_smoke.json",
    "xi_mellin_boundary_concomitant.json",
    "xi_mellin_boundary_concomitant_smoke.json",
    "xi_augmented_trace_repair_schur.json",
    "xi_augmented_trace_repair_schur_smoke.json",
    "xi_augmented_trace_continuum_lift.json",
    "xi_augmented_trace_continuum_lift_hardened_theorem.json",
    "klm_debranges_bridge_attempt.json",
    "rh_debranges_bridge_ledger.json",
    "quotient_to_original_weyl_lift.json",
    "quotient_original_transport_identity_theorem.json",
    "quotient_primitive_endpoint_input_theorem.json",
    "quotient_to_original_weyl_lift_theorem.json",
    "boundary_repair_identity.json",
    "quotient_minimal_repair_consequence_theorem.json",
    "canonical_boundary_repair_comparison.json",
    "primitive_boundary_zero_consequence_theorem.json",
    "primitive_trace_density_consequence_theorem.json",
    "green_feature_dq_vanishes_theorem.json",
    "primitive_endpoint_compatibility_consequence_theorem.json",
    "primitive_endpoint_compatibility_theorem.json",
    "primitive_boundary_transport_audit.json",
    "primitive_trace_image_density.json",
    "dq_vanishing_schur_defect_scan.json",
    "dq_vanishing_schur_defect_scan_b8c5.json",
    "dq_vanishing_schur_defect_scan_b10c7.json",
    "dq_vanishing_repair_route_summary.json",
    "repair_free_schur_equivalence.json",
    "repair_free_schur_kernel_probe.json",
    "repair_free_schur_kernel_probe_b12c9.json",
    "trace_schur_kernel_derivation.json",
    "trace_euler_lagrange_minimizer.json",
    "trace_volterra_green_feature_map.json",
    "volterra_feature_contraction.json",
    "volterra_transport_identity.json",
    "volterra_hardy_multiplier_symbolic_theorem.json",
    "volterra_hardy_transport_derivation.json",
    "volterra_tail_positive_form_theorem.json",
    "volterra_tail_restriction_consequence_theorem.json",
    "green_lift_boundary_theorem.json",
    "green_lift_boundary_minimality_theorem.json",
    "green_lift_contractivity_form_theorem.json",
    "continuum_green_lift_closure_theorem.json",
    "green_lift_contraction_consequence_theorem.json",
    "full_theta_source_inactive_schur_tail_certificate.json",
    "continuum_tail_absorption_certificate.json",
    "source_inactive_minmax_tail_theorem.json",
    "source_inactive_tail_constants_theorem.json",
    "synchronized_active_range_interval_theorem.json",
    "synchronized_active_range_theorem.json",
    "continuum_trace_frame_lower_bound_theorem.json",
    "trace_frame_interval_lower_bound_certificate.json",
    "trace_frame_finite_gamma_consequence_theorem.json",
    "trace_quadrature_interval_consistency_theorem.json",
    "trace_quadrature_gamma_consequence_theorem.json",
    "abstract_high_block_compact_source_mosco_theorem.json",
    "abstract_mosco_projection_convergence_theorem.json",
    "high_block_mosco_projection_input_theorem.json",
    "abstract_finite_rank_compression_convergence_theorem.json",
    "abstract_compact_compression_norm_convergence_theorem.json",
    "high_block_source_operator_compactness_theorem.json",
    "high_block_compact_compression_input_theorem.json",
    "abstract_riesz_projection_continuity_theorem.json",
    "high_block_source_spectral_gap_theorem.json",
    "high_block_spectral_projection_input_theorem.json",
    "abstract_compact_source_spectral_projection_theorem.json",
    "abstract_minmax_tail_passage_theorem.json",
    "high_block_minmax_tail_input_theorem.json",
    "high_block_tail_estimate_continuum_passage_theorem.json",
    "high_block_compact_exhaustion_consequence_theorem.json",
    "high_block_compact_exhaustion_proof.json",
    "high_block_exhaustion_theorem.json",
    "commuted_sturm_elliptic_trace_theorem.json",
    "source_range_hardy_green_hardened_theorem.json",
    "active_trace_control_consequence_theorem.json",
    "source_range_hardy_green_theorem.json",
    "continuum_green_representer_theorem.json",
    "continuum_trace_to_source_green_kernel.json",
    "continuum_active_trace_range_theorem.json",
    "closed_trace_active_unique_continuation_theorem.json",
    "endpoint_obstruction_compatibility_theorem.json",
    "endpoint_map_rank_ball_certificate.json",
    "endpoint_map_center_stability.json",
    "endpoint_flow_chebyshev_center.json",
    "endpoint_flow_chebyshev_center_7_9.json",
    "endpoint_flow_magnus_center_2_4.json",
    "endpoint_adjoint_row_flow_center.json",
    "endpoint_adjoint_row_flow_center_7_9_11.json",
    "endpoint_grassmann_flow_center.json",
    "endpoint_grassmann_flow_center_7_9_11.json",
    "endpoint_grassmann_chart_ball_certificate.json",
    "endpoint_riccati_flow_enclosure.json",
    "endpoint_riccati_krawczyk_collocation.json",
    "endpoint_coefficient_ball_enclosure.json",
    "endpoint_coefficient_interval_enclosure.json",
    "endpoint_confluent_trace_tail_certificate.json",
    "endpoint_confluent_segment_bernstein_certificate.json",
    "endpoint_eigenrow_interval_propagation.json",
    "endpoint_eigenrow_interval_propagation_200.json",
    "endpoint_coefficient_synchronized_200_certificate.json",
    "adjoint_green_boundary_diagnostic.json",
    "adjoint_green_bvp_regularized.json",
    "adjoint_green_jump_conditions.json",
    "adjoint_green_jump_conditions_theorem.json",
    "adjoint_green_endpoint_selection.json",
    "adjoint_green_endpoint_range_theorem.json",
    "adjoint_green_endpoint_range_interval_theorem.json",
    "source_quadrature_error_bound.json",
    "source_quadrature_error_bound_g65_e9.json",
    "source_quadrature_interval_envelope_s16.json",
    "source_quadrature_chebyshev_envelope_d32_g129.json",
    "source_quadrature_chebyshev_ball_d32_g129.json",
    "source_quadrature_bernstein_bound_rho2.json",
    "source_quadrature_complex_ball_bernstein_rho2.json",
    "closed_trace_local_ode_obstruction.json",
    "high_frequency_tail_refinement.json",
    "spectral_tail_moment_certificate.json",
    "source_envelope_tail_certificate.json",
    "source_packet_tracking.json",
    "commuted_source_certificate.json",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_manifest(path: Path) -> dict:
    if not path.exists():
        return {"version": 1, "snapshots": []}
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", default="trace-resolution lemma")
    parser.add_argument("--source-html", default="visual_explainer.html")
    parser.add_argument("--out-dir", default="visual_snapshots")
    parser.add_argument("--id", default="")
    args = parser.parse_args()

    root = Path.cwd()
    out_root = root / args.out_dir
    out_root.mkdir(parents=True, exist_ok=True)

    stamp = args.id or datetime.now().strftime("%Y%m%d-%H%M%S")
    snap_dir = out_root / stamp
    if snap_dir.exists():
        raise SystemExit(f"snapshot already exists: {snap_dir}")
    snap_dir.mkdir(parents=True)

    source_html = root / args.source_html
    if not source_html.exists():
        raise SystemExit(f"missing {source_html}")

    copied = []
    for name in [args.source_html, *DATA_FILES]:
        src = root / name
        if not src.exists():
            continue
        dst = snap_dir / name
        shutil.copy2(src, dst)
        copied.append(
            {
                "path": name,
                "bytes": dst.stat().st_size,
                "sha256": sha256(dst),
            }
        )

    readme = {
        "id": stamp,
        "label": args.label,
        "createdAt": datetime.now().isoformat(timespec="seconds"),
        "entry": f"{stamp}/{args.source_html}",
        "screenshot": f"{stamp}/screenshot.png",
        "files": copied,
    }
    write_json(snap_dir / "snapshot.json", readme)

    manifest_path = out_root / "manifest.json"
    manifest = load_manifest(manifest_path)
    manifest["snapshots"] = [s for s in manifest.get("snapshots", []) if s.get("id") != stamp]
    manifest["snapshots"].append(readme)
    write_json(manifest_path, manifest)

    print(f"snapshot_id={stamp}")
    print(f"snapshot_dir={snap_dir}")
    print(f"entry={snap_dir / args.source_html}")
    print(f"screenshot={snap_dir / 'screenshot.png'}")
    print(f"copied_files={len(copied)}")


if __name__ == "__main__":
    main()

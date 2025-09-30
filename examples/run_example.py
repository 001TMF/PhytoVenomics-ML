"""
Example script for running the antibody design pipeline
"""

from antibody_pipeline import AntibodyDesignPipeline, PipelineConfig


def main():
    """Run example antibody design"""

    # Create configuration
    config = PipelineConfig(
        # Input files
        antigen_pdb="data/snake_toxins/pdb_structures/7M6C.pdb",
        framework_pdb="data/antibody_structures/human/example_framework.pdb",

        # Output
        output_dir="./results/example_run",

        # Design settings
        design_mode="cdr_h3",
        num_designs=10,
        num_design_steps=50,

        # Epitope (optional)
        epitope_residues=["A45", "A48", "A52"],

        # Docking
        use_docking=True,
        hdock_bin="./bin/hdock",
        createpl_bin="./bin/createpl",

        # RFdiffusion
        rfdiffusion_path="./RFdiffusion",

        # IgLM
        use_iglm_optimization=True,
        iglm_species="[HUMAN]",

        # Filtering
        max_clashes=100,
        min_plddt=70.0,
        max_sap_score=50.0,

        # Structure prediction
        structure_predictor="chai",

        # Device
        device="cuda:0",
    )

    # Create and run pipeline
    print("Initializing pipeline...")
    pipeline = AntibodyDesignPipeline(config)

    print("\nRunning pipeline...")
    candidates = pipeline.run()

    # Display results
    print("\n" + "=" * 80)
    print("RESULTS")
    print("=" * 80)

    print(f"\nTotal designs generated: {len(pipeline.candidates)}")
    print(f"Designs passed filters: {len(candidates)}")

    if candidates:
        print("\nTop 5 designs:")
        for i, candidate in enumerate(candidates[:5], 1):
            print(f"\n{i}. {candidate.design_id}")
            print(f"   Sequence: {candidate.sequence[:60]}...")
            print(f"   PDB: {candidate.pdb_path}")

            if candidate.iglm_log_likelihood:
                print(f"   IgLM log-likelihood: {candidate.iglm_log_likelihood:.2f}")
            if candidate.plddt:
                print(f"   pLDDT: {candidate.plddt:.1f}")
            if candidate.sap_score:
                print(f"   SAP score: {candidate.sap_score:.1f}")

            print(f"   CDR sequences:")
            for cdr_name, cdr_seq in candidate.cdr_sequences.items():
                print(f"     {cdr_name}: {cdr_seq}")

        print(f"\nFull results saved to: {config.output_dir}/results.json")

    else:
        print("\nNo designs passed all filters. Try adjusting thresholds.")


if __name__ == '__main__':
    main()
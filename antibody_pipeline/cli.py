"""
Command-line interface for antibody design pipeline
"""

import argparse
import sys
from typing import Optional

from .pipeline import AntibodyDesignPipeline, PipelineConfig
from .config.parser import load_config, save_config, create_default_config


def create_parser() -> argparse.ArgumentParser:
    """Create argument parser"""
    parser = argparse.ArgumentParser(
        description='Integrated Antibody Design Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with config file
  python -m antibody_pipeline --config config.yaml

  # Run with command-line arguments
  python -m antibody_pipeline --antigen antigen.pdb --framework framework.pdb -o results/

  # Generate default config
  python -m antibody_pipeline --create-config default_config.yaml

For more information, see the documentation.
        """
    )

    # Commands
    parser.add_argument(
        '--create-config',
        metavar='PATH',
        help='Create default configuration file and exit'
    )

    # Configuration
    parser.add_argument(
        '--config', '-c',
        metavar='PATH',
        help='Path to configuration file (YAML or JSON)'
    )

    # Input files - three modes
    input_group = parser.add_argument_group('Input Files (choose one mode)')
    input_group.add_argument(
        '--antigen',
        metavar='PATH',
        help='Path to antigen PDB file (for separate or de novo mode)'
    )
    input_group.add_argument(
        '--framework',
        metavar='PATH',
        help='Path to antibody framework PDB file (for separate mode)'
    )
    input_group.add_argument(
        '--complex',
        metavar='PATH',
        help='Path to pre-docked complex PDB file (for pre-docked mode)'
    )
    input_group.add_argument(
        '--antibody-chains',
        nargs='+',
        help='Antibody chain IDs for pre-docked complex (e.g., H L)'
    )
    input_group.add_argument(
        '--antigen-chains',
        nargs='+',
        help='Antigen chain IDs for pre-docked complex (e.g., A)'
    )

    # Output
    parser.add_argument(
        '--output', '-o',
        metavar='DIR',
        default='./results',
        help='Output directory (default: ./results)'
    )

    # Design settings
    design_group = parser.add_argument_group('Design Settings')
    design_group.add_argument(
        '--mode',
        choices=['cdr_h3', 'single_cdr', 'multiple_cdrs', 'full'],
        default='cdr_h3',
        help='Design mode (default: cdr_h3)'
    )
    design_group.add_argument(
        '--num-designs', '-n',
        type=int,
        default=10,
        help='Number of designs to generate (default: 10)'
    )
    design_group.add_argument(
        '--epitope',
        nargs='+',
        help='Epitope residues (e.g., A10 A15 A20)'
    )

    # RFdiffusion
    rf_group = parser.add_argument_group('RFdiffusion Settings')
    rf_group.add_argument(
        '--rfdiffusion-path',
        metavar='PATH',
        default='./RFdiffusion',
        help='Path to RFdiffusion installation (default: ./RFdiffusion)'
    )

    # IgLM
    iglm_group = parser.add_argument_group('IgLM Settings')
    iglm_group.add_argument(
        '--no-iglm',
        action='store_true',
        help='Disable IgLM optimization'
    )
    iglm_group.add_argument(
        '--species',
        default='[HUMAN]',
        help='IgLM species token (default: [HUMAN])'
    )

    # Filtering
    filter_group = parser.add_argument_group('Filtering Thresholds')
    filter_group.add_argument(
        '--max-clashes',
        type=int,
        default=100,
        help='Maximum allowed clashes (default: 100)'
    )
    filter_group.add_argument(
        '--min-plddt',
        type=float,
        default=70.0,
        help='Minimum pLDDT threshold (default: 70.0)'
    )
    filter_group.add_argument(
        '--max-sap',
        type=float,
        default=50.0,
        help='Maximum SAP score (default: 50.0)'
    )

    # Device
    parser.add_argument(
        '--device',
        default='cuda:0',
        help='Device to use (default: cuda:0)'
    )

    return parser


def config_from_args(args: argparse.Namespace) -> PipelineConfig:
    """
    Create PipelineConfig from command-line arguments.

    Args:
        args: Parsed arguments

    Returns:
        PipelineConfig instance
    """
    return PipelineConfig(
        antigen_pdb=args.antigen,
        framework_pdb=args.framework,
        complex_pdb=args.complex if hasattr(args, 'complex') else None,
        antibody_chains=args.antibody_chains if hasattr(args, 'antibody_chains') else None,
        antigen_chains=args.antigen_chains if hasattr(args, 'antigen_chains') else None,
        output_dir=args.output,
        design_mode=args.mode,
        num_designs=args.num_designs,
        epitope_residues=args.epitope,
        rfdiffusion_path=args.rfdiffusion_path,
        use_iglm_optimization=not args.no_iglm,
        iglm_species=args.species,
        max_clashes=args.max_clashes,
        min_plddt=args.min_plddt,
        max_sap_score=args.max_sap,
        device=args.device,
    )


def config_from_file(config_path: str) -> PipelineConfig:
    """
    Create PipelineConfig from configuration file.

    Args:
        config_path: Path to config file

    Returns:
        PipelineConfig instance
    """
    config_dict = load_config(config_path)

    return PipelineConfig(
        antigen_pdb=config_dict.get('antigen_pdb'),
        framework_pdb=config_dict.get('framework_pdb'),
        complex_pdb=config_dict.get('complex_pdb'),
        antibody_chains=config_dict.get('antibody_chains'),
        antigen_chains=config_dict.get('antigen_chains'),
        output_dir=config_dict.get('output_dir', './results'),
        design_mode=config_dict.get('design', {}).get('mode', 'cdr_h3'),
        num_designs=config_dict.get('design', {}).get('num_designs', 10),
        num_design_steps=config_dict.get('design', {}).get('num_steps', 50),
        epitope_residues=config_dict.get('epitope', {}).get('residues'),
        use_docking=config_dict.get('docking', {}).get('enabled', True),
        hdock_bin=config_dict.get('docking', {}).get('hdock_bin', './bin/hdock'),
        createpl_bin=config_dict.get('docking', {}).get('createpl_bin', './bin/createpl'),
        rfdiffusion_path=config_dict.get('rfdiffusion', {}).get('path', './RFdiffusion'),
        rfdiffusion_weights=config_dict.get('rfdiffusion', {}).get('weights'),
        use_iglm_optimization=config_dict.get('iglm', {}).get('enabled', True),
        iglm_species=config_dict.get('iglm', {}).get('species', '[HUMAN]'),
        iglm_temperature=config_dict.get('iglm', {}).get('temperature', 1.0),
        max_clashes=config_dict.get('filtering', {}).get('max_clashes', 100),
        min_plddt=config_dict.get('filtering', {}).get('min_plddt', 70.0),
        min_pdockq=config_dict.get('filtering', {}).get('min_pdockq', 0.23),
        max_sap_score=config_dict.get('filtering', {}).get('max_sap_score', 50.0),
        min_cdr_interface_pct=config_dict.get('filtering', {}).get('min_cdr_interface_pct', 60.0),
        structure_predictor=config_dict.get('structure_prediction', {}).get('predictor', 'chai'),
        chai_path=config_dict.get('structure_prediction', {}).get('chai_path'),
        af3_path=config_dict.get('structure_prediction', {}).get('af3_path'),
        device=config_dict.get('device', 'cuda:0'),
    )


def main(argv: Optional[list] = None):
    """Main entry point"""
    parser = create_parser()
    args = parser.parse_args(argv)

    # Handle --create-config
    if args.create_config:
        default_config = create_default_config()
        save_config(default_config, args.create_config)
        print(f"Default configuration saved to: {args.create_config}")
        print("Edit this file and run with: python -m antibody_pipeline --config", args.create_config)
        return 0

    # Validate inputs
    if args.config:
        # Load from config file
        try:
            config = config_from_file(args.config)
        except Exception as e:
            print(f"Error loading config file: {e}", file=sys.stderr)
            return 1
    elif args.antigen:
        # Load from command-line args
        config = config_from_args(args)
    else:
        print("Error: Must provide either --config or --antigen", file=sys.stderr)
        parser.print_help()
        return 1

    # Run pipeline
    try:
        pipeline = AntibodyDesignPipeline(config)
        candidates = pipeline.run()

        print(f"\nTop 5 designs:")
        for i, candidate in enumerate(candidates[:5]):
            print(f"  {i+1}. {candidate.design_id}")
            print(f"     Sequence: {candidate.sequence[:50]}...")
            if candidate.iglm_log_likelihood:
                print(f"     IgLM LL: {candidate.iglm_log_likelihood:.2f}")
            if candidate.plddt:
                print(f"     pLDDT: {candidate.plddt:.1f}")
            print()

        return 0

    except Exception as e:
        print(f"Error running pipeline: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
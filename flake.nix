{
  description = "Genetic Health Analysis Pipeline - FASTQ to health reports";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        python = pkgs.python3.withPackages (ps: with ps; [ pytest ]);
        bioTools = with pkgs; [
          minimap2
          samtools
          bcftools
          htslib
        ] ++ pkgs.lib.optionals (pkgs.stdenv.hostPlatform.system != "aarch64-darwin") [
          pkgs.fastp
        ];
        utils = with pkgs; [ curl gzip coreutils ];
        allDeps = bioTools ++ utils ++ [ python ];
      in
      {
        devShells.default = pkgs.mkShell {
          packages = allDeps;
          shellHook = ''
            echo "=== Genetic Health Analysis Pipeline ==="
            echo "Tools: minimap2, samtools, bcftools, python3"
            if ! command -v fastp &>/dev/null; then
              echo ""
              echo "NOTE: fastp not in nixpkgs for Apple Silicon."
              echo "Install via: brew install fastp"
              echo "Or use --skip-qc to skip quality control."
            fi
            echo ""
            echo "Setup (one-time):"
            echo "  bash scripts/setup_reference.sh"
            echo ""
            echo "Run pipeline:"
            echo "  python scripts/wgs_pipeline.py /path/to/reads.fastq"
            echo "  python scripts/wgs_pipeline.py /path/to/reads.fastq --name \"John Doe\""
            echo ""
            echo "Run analysis only (on existing genome.txt):"
            echo "  python scripts/run_full_analysis.py"
            echo ""
            echo "Run tests:"
            echo "  python -m pytest tests/ -v"
          '';
        };

        packages.default = pkgs.writeShellScriptBin "genetic-health" ''
          export PATH="${pkgs.lib.makeBinPath allDeps}:$PATH"
          exec ${python}/bin/python ${./scripts/wgs_pipeline.py} "$@"
        '';

        apps.default = {
          type = "app";
          program = "${self.packages.${system}.default}/bin/genetic-health";
        };
      });
}

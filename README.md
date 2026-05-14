# Serenity_Primer
Gibson Assembly primer design tool — Tm calculation, dimer analysis, GUI

Gibson Assembly primer design tool for molecular biology.

Scans flanking regions around a theoretical insert and enumerates all viable 
Gibson primers with Tm calculation (SantaLucia 1998, via BioPython) and 
validated antiparallel dimer analysis.

## Features
- Nearest-neighbor Tm (SantaLucia 1998, BioPython)
- Self-dimer and heterodimer detection with Watson-Crick alignment display
- Single-window and paired flanking scan modes
- Tkinter GUI + command-line interface
- CSV, FASTA, and primer-pair CSV export

## Installation
```bash
pip install biopython rich click
```

## Quick Start
```bash
# GUI
python gibsonfinder_gui.py

# Command line — single window
python gibsonfinder.py single my_sequence.fasta --win-start 1 --win-end 300

# Command line — paired flanking (Gibson-specific)
python gibsonfinder.py paired my_sequence.fasta \
    --up-start 1 --up-end 150 \
    --dn-start 501 --dn-end 650
```

## Citation
SantaLucia J. (1998) PNAS 95:1460–1465  
Owczarzy et al. (2004) Biochemistry 43:3537

## License
MIT

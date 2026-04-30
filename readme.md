# Code for Global Sponge Plasmidome

This repository contains the scripts used to generate figures and analyses for:

**Global Sponge Plasmidome**  
*DOI: to be added upon publication*

## Abstract

Mobile genetic elements, particularly plasmids, are fundamental drivers of microbial evolution and horizontal gene transfer in host-microbe interactions. Here, we present a comparative plasmidomics analysis of 526 global marine sponge metagenomes, recovering 6,945 non-redundant plasmid cluster units. We find that plasmid abundance and diversity are significantly higher in high microbial-abundance (HMA) sponges than in low microbial-abundance (LMA) sponges, and that the sponge plasmidome represents a vast, underrepresented genetic reservoir. These findings highlight the pivotal, largely unrecognized role of plasmids in the functional adaptation of this ecologically critical marine metaorganism.

---

## Requirements

### Python
- Python 3.12.3
- `numpy` 1.26.4
- `pandas` 2.2.2
- `matplotlib` 3.9.2
- `seaborn` 0.13.2
- `scipy` 1.13.1
- `scikit-learn` 1.5.1
- `biopython` 1.85
- `plotly` 6.2.0
- `matplotlib-venn` 1.1.2
- `adjustText` 1.3.0
- `statsmodels` (latest conda-forge build)

### R
- R 4.0+
- Key packages: `ggplot2`, `circlize`, `dplyr`, `tidyr`

> **Note:** Input data files required to run the scripts are available at [add data repository link, e.g., Zenodo or ENA].

---

## Repository Structure

All scripts are organized by figure:

```
scripts/
├── Figure_1/
├── Figure_3/
├── Figure_4/
├── Figure_5/
├── Figure_6/
└── Figure_8/
```

Each folder contains all scripts required to reproduce the corresponding figure panels. Some scripts generate multiple panels.

---

## Script-to-Figure Mapping

### Figure 1
| Script | Output |
|--------|--------|
| `plasmid_length_distribution_for_publication.py` | Panel a |
| `venn_diagram_orit.py` | Panel b |

### Figure 3
| Script | Output |
|--------|--------|
| `plasmid_abundance_rpkm_final_plasmids_noND_clean_smallfig.py` | All panels |

### Figure 4
| Script | Output |
|--------|--------|
| `cytoscape_mob_counts.py` | Panel b |
| `circos.r` | Panel c |
| `geo_density_overlay_with_clusters.py` | Panel d |

### Figure 5
| Script | Output |
|--------|--------|
| `pathway_enrichment.r` | Panel a |
| `hierarchical_clustering_kegg_no_host.py` | Panel b |
| `widespread_plasmids_functional_annotation.py` | Panel c |

### Figure 8
| Script | Output |
|--------|--------|
| `functional_enrichment_agelas_geodia_all_pfam_kegg_combined.py` | All panels |

---

## Contact

For questions about the code, please open an issue in this repository or contact the corresponding author.

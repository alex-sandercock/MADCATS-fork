# Citation

If MADCATS is useful in your research, please cite the implementation.

## Citing MADCATS

```
Sandercock, A.M. (2026). MADCATS: Marker-Anchored Discovery and Counting of Alleles
from Targeted Sequencing. GitHub repository: https://github.com/Breeding-Insight/MADCATS
```

BibTeX:
```bibtex
@software{sandercock2026madcats,
  author = {Sandercock, Alexander M.},
  title = {MADCATS: Marker-Anchored Discovery and Counting of Alleles from Targeted Sequencing},
  year = {2026},
  url = {https://github.com/Breeding-Insight/MADCATS}
}
```

> **Note:** A formal publication with DOI is planned. This section will be updated when available.

## AI Assistance

MADCATS was developed in collaboration with Claude Opus 5.5 and Fable 5.1 (Anthropic), which
helped write its code and documentation. The authors directed the work, reviewed and tested the
results, and are ultimately responsible for the application.

## Related Resources

- Breeding Insight, [MADC file format](https://breeding-insight.github.io/BIGr_documentation/posts/madc-file-format/)
- Breeding Insight, [MADC data processing](https://breeding-insight.github.io/BIGr_documentation/posts/madc-data-processing/)
- Breeding Insight, [BIGapp-PanelHub](https://github.com/Breeding-Insight/BIGapp-PanelHub): public panel metadata and design Ref/Alt tags
- Breeding Insight, [dartag_alleles](https://github.com/Breeding-Insight/dartag_alleles): the reference-BLAST / paralog audit pattern
- The sweetpotato DArTag panel publication: [PMC11354173](https://pmc.ncbi.nlm.nih.gov/articles/PMC11354173/)

If you call dosages from MADCATS output with updog, please also cite:

```
Gerard, D., Ferrão, L.F.V., Garcia, A.A.F., & Stephens, M. (2018). Genotyping
polyploids from messy sequencing data. Genetics 210(3):789-807.
```

If you clean a MADC with `madcats filter-paralogs`, please also cite minimap2, which does its
alignments:

```
Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
Bioinformatics 34(18):3094-3100. https://doi.org/10.1093/bioinformatics/bty191

Li, H. (2021). New strategies to improve minimap2 alignment accuracy.
Bioinformatics 37(23):4572-4574. https://doi.org/10.1093/bioinformatics/btab705
```

If you call genotypes with `madcats call`, please also cite the methods its model is adapted
from: polyRAD for the likelihood and the Hardy-Weinberg prior, and updog for the allelic bias
(`--allele-bias fit`) and the `norm` prior (`--prior norm`):

```
Clark, L.V., Lipka, A.E., & Sacks, E.J. (2019). polyRAD: genotype calling with
uncertainty from sequencing data in polyploids and diploids. G3 9(3):663-673.
https://doi.org/10.1534/g3.118.200913

Gerard, D., Ferrão, L.F.V., Garcia, A.A.F., & Stephens, M. (2018). Genotyping
polyploids from messy sequencing data. Genetics 210(3):789-807.
https://doi.org/10.1534/genetics.118.301468

Gerard, D., & Ferrão, L.F.V. (2020). Priors for genotyping polyploids.
Bioinformatics 36(6):1795-1800. https://doi.org/10.1093/bioinformatics/btz852
```

## Example Acknowledgment

In your methods section:

> Microhaplotype read counts were generated from the demultiplexed FASTQ files with MADCATS
> v0.14.0 (Sandercock, 2026), using panel tags reconstructed from the vendor design files and
> default evidence and denoising settings.

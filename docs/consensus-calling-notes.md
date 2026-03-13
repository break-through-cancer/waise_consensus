# Consensus Calling Notes

This note explains how the Jasmine-based consensus step works in this repository and what each upstream SV caller contributes.

## How Jasmine Merging Works

Jasmine merges structurally similar SV records across multiple input VCF files. The input is a `file_list`, and Jasmine treats each listed VCF as one input callset for merging purposes.

At a high level, Jasmine groups compatible records first, then explores candidate merges in order of increasing distance. The algorithm is often described as Kruskal-like: it considers the closest remaining pair first, merges only if the constraints still hold, and keeps repeating until no valid nearby pairs remain. To avoid materializing every pairwise comparison, Jasmine uses spatial indexing such as a KD-tree to find nearby candidates and a min-heap to keep track of the current smallest candidate edge.

Each accepted merge joins two records or partial components into one connected component. By default, Jasmine does not allow a component to contain two records from the same input VCF file. In practice, that means two distinct calls from the same caller/callset do not normally collapse into one merged component unless intrasample merging is explicitly enabled. This is the mechanism behind the common "same sample/file cannot merge with itself" rule.

Jasmine writes one merged record per connected component and annotates the output VCF with support fields such as `SUPP`, `SUPP_VEC`, and `IDLIST`.

The Nature Methods paper describes distance-aware merging with thresholds that scale with SV similarity. This repository does not pass explicit Jasmine distance flags, so the exact runtime threshold behavior comes from the `jasminesv:1.1.5` container defaults.

## SUPP and SUPP_VEC

- `SUPP` is the number of input VCF files that contributed a call to the merged record.
- `SUPP_VEC` is a bit vector in `file_list` order showing which inputs supported the merged record.

For example, if the `file_list` order were `manta`, `svaba`, `lumpy`, `delly`, `gridss`, then `SUPP_VEC=10110` would mean the event was supported by `manta`, `lumpy`, and `delly`, so `SUPP=3`.

`INFO/SUPP > 1` therefore means "keep only merged events supported by at least two input callsets."

The bit positions are defined by the generated `file_list`, not by a universal fixed caller order. In this workflow, that means `SUPP_VEC` should always be interpreted against the corresponding Jasmine file list for that merge step.

One detail matters in this repository: `SUPP` does not always mean the same biological thing. In the normal-panel step it reflects recurrence across normal samples for one caller. In the tumor consensus step it reflects agreement across caller VCFs for one tumor-normal pair.

## How This Repository Uses Jasmine

The workflow uses Jasmine twice.

First, it builds a caller-specific normal panel. For each caller, all normal VCFs are merged with Jasmine in [modules/jasmine_panel.nf](/Users/youyun/Documents/HMS/PhD/beroukhimlab/BTC/Code/waise_consensus/modules/jasmine_panel.nf). The merged panel is then filtered to keep only recurring events with `INFO/SUPP > 1` at [modules/jasmine_panel.nf:23](/Users/youyun/Documents/HMS/PhD/beroukhimlab/BTC/Code/waise_consensus/modules/jasmine_panel.nf:23). Here, `SUPP` is the number of normal samples supporting a recurrent event for that one caller.

Second, each sample's tumor VCFs are filtered against those caller-specific normal panels, then the five filtered caller VCFs are merged into one Jasmine consensus VCF in [main.nf:275](/Users/youyun/Documents/HMS/PhD/beroukhimlab/BTC/Code/waise_consensus/main.nf:275) and [modules/tumour_consensus.nf:46](/Users/youyun/Documents/HMS/PhD/beroukhimlab/BTC/Code/waise_consensus/modules/tumour_consensus.nf:46). In that step, the five inputs are one filtered VCF per caller for a single case, so `SUPP` measures caller support rather than cohort frequency.

The final tumor consensus filter keeps only `INFO/SUPP > 4` at [modules/tumour_consensus.nf:69](/Users/youyun/Documents/HMS/PhD/beroukhimlab/BTC/Code/waise_consensus/modules/tumour_consensus.nf:69). Because this workflow currently merges five caller VCFs per case, that is effectively an "all five callers support this event" rule.

## SV Caller Snapshot

| Caller | Core method | Practical strengths | Common caveats | Typical best fit |
| --- | --- | --- | --- | --- |
| `Manta` | Paired-end and split-read evidence with local assembly/refinement | Fast, strong breakpoint resolution, handles somatic tumor-normal mode well | Can be conservative in repetitive or highly complex regions | General short-read germline or somatic SV calling |
| `LUMPY` via `smoove` | Probabilistic integration of discordant paired-end and split-read signals | Mature and flexible; `smoove` adds useful preprocessing and noise reduction | Sensitive to alignment artifacts and upstream filtering quality | Broad WGS SV discovery with practical wrapper-based execution |
| `SvABA` | Genome-wide local assembly around discordant, split, and clipped reads | Strong across indels and medium-sized SVs; useful for templated insertions and complex junctions | Heavier than purely alignment-based callers | Somatic SV plus indel discovery where assembly helps |
| `DELLY` | Integrated paired-end and split-read analysis with breakpoint refinement | Mature caller with good support for balanced and unbalanced rearrangements | Representation and sensitivity can vary with library/alignment quality | Standard short-read SV calling in matched or cohort settings |
| `GRIDSS` | Break-end assembly with a positional de Bruijn graph plus probabilistic scoring | Strong on complex rearrangements and precise breakpoint characterization | Higher resource cost and more complex outputs | Rearrangement-focused analyses where breakpoint detail matters |

## Individual Caller Introductions

### Manta

Manta is a short-read SV and indel caller designed for both germline and tumor-normal analysis. It starts from paired-end and split-read evidence, assembles candidate loci locally, and uses that assembly to refine breakpoint resolution. In practice, Manta is often valued for speed and for producing reasonably precise breakpoint calls. It works well as a general-purpose somatic SV caller in matched tumor-normal pipelines like this one.

### LUMPY via smoove

LUMPY is a probabilistic SV framework that combines multiple alignment signals, especially discordant read pairs and split reads, into one breakpoint model. This repository runs it through `smoove`, which automates preprocessing, filtering, and practical best practices around LUMPY. The wrapper is important because it reduces common noise sources before LUMPY sees the evidence. Conceptually, this caller contributes a robust alignment-signal-based view of SVs rather than an assembly-first view.

### SvABA

SvABA uses local assembly across the genome to recover indels and structural variants from aberrantly aligned reads. It is especially useful in the range between classic small indels and larger breakpoint-style SVs, where purely alignment-based methods can lose sensitivity. Because it assembles local sequence context, it is also helpful for complex junctions and short templated-sequence insertions. In a consensus setting, SvABA often complements callers that rely more directly on discordant pairs and split reads.

### DELLY

DELLY integrates paired-end and split-read evidence to detect deletions, duplications, inversions, and translocations. Its design emphasizes breakpoint refinement after initial paired-end clustering, which makes it useful for a broad range of canonical SV classes. DELLY has been used widely in both germline and cancer studies and remains a common baseline caller in short-read pipelines. In this workflow it adds another mature, non-assembly-based perspective to the consensus.

### GRIDSS

GRIDSS is centered on break-end assembly using a positional de Bruijn graph, then combines assembly, split-read, and read-pair evidence in a probabilistic scoring framework. That makes it particularly good at describing complex rearrangements and fine breakpoint features such as inserted sequence or microhomology. Compared with lighter callers, it is more resource-intensive, but it often contributes high-value breakpoint detail. In consensus pipelines, GRIDSS is often the caller that most strongly helps with complex or rearrangement-heavy samples.

## References

- Jasmine paper: https://www.nature.com/articles/s41592-022-01753-3
- Jasmine open-access full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC10006329/
- Jasmine user manual: https://github.com/mkirsche/Jasmine/wiki/Jasmine-User-Manual
- Manta paper: https://pubmed.ncbi.nlm.nih.gov/26647377/
- LUMPY paper: https://pubmed.ncbi.nlm.nih.gov/24970577/
- `smoove` documentation: https://github.com/brentp/smoove
- SvABA paper: https://pubmed.ncbi.nlm.nih.gov/29535149/
- DELLY paper: https://pubmed.ncbi.nlm.nih.gov/22962449/
- GRIDSS paper: https://pubmed.ncbi.nlm.nih.gov/29097403/

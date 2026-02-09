# Tool and Database Citations

Citations for all bioinformatics tools and databases used in the Bennu metagenomic analysis pipeline.

---

## Core Pipeline

### Prodigal
> Hyatt, D., Chen, G.L., LoCascio, P.F., Land, M.L., Larimer, F.W. & Hauser, L.J. (2010) Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics* 11, 119. DOI: [10.1186/1471-2105-11-119](https://doi.org/10.1186/1471-2105-11-119)

**Used for:** Gene calling (Stage 03) -- predicts protein-coding genes from genome assemblies.
**Verified:** [PubMed 20211023](https://pubmed.ncbi.nlm.nih.gov/20211023/)
**Abstract excerpt:** "We describe Prodigal (PROkaryotic DYnamic programming Gene-finding ALgorithm), a fast, lightweight, open source gene prediction program, with improved gene structure prediction, improved translation initiation site recognition, and reduced false positives."

---

### QUAST
> Gurevich, A., Saveliev, V., Vyahhi, N. & Tesler, G. (2013) QUAST: quality assessment tool for genome assemblies. *Bioinformatics* 29(8), 1072--1075. DOI: [10.1093/bioinformatics/btt086](https://doi.org/10.1093/bioinformatics/btt086)

**Used for:** Assembly quality assessment (Stage 01) -- evaluates genome assembly statistics (N50, L50, contigs, etc.).
**Verified:** [PubMed 23422339](https://pubmed.ncbi.nlm.nih.gov/23422339/)
**Abstract excerpt:** "QUAST can evaluate assemblies both with a reference genome, as well as without a reference, and produces many reports, summary tables and plots to help scientists in their research and publications."

---

### DuckDB
> Raasveldt, M. & Muhleisen, H. (2019) DuckDB: an Embeddable Analytical Database. *Proceedings of the 2019 International Conference on Management of Data (SIGMOD '19)*, 1981--1984. DOI: [10.1145/3299869.3320212](https://doi.org/10.1145/3299869.3320212)

**Used for:** Core analytical database engine (Stage 07 and throughout) -- stores proteins, annotations, predicates, and supports all queries.
**Verified:** [ACM Digital Library](https://dl.acm.org/doi/10.1145/3299869.3320212)
**Abstract excerpt:** "DuckDB is a novel data management system designed to execute analytical SQL queries while embedded in another process."

---

### Astra
> No formal publication. Astra is an open-source pyHMMER-based sequence search and retrieval tool.
> GitHub: [https://github.com/jwestrob/astra](https://github.com/jwestrob/astra)

**Used for:** HMM annotation wrapper (Stage 04) -- runs hmmsearch against PFAM, KOfam, HydDB, VOGdb, DefenseFinder, and CRISPRCasFinder HMM databases.
**Verified:** [GitHub repository](https://github.com/jwestrob/astra)
**Note:** Astra is built on PyHMMER. When citing Astra, cite the underlying PyHMMER and HMMER publications (see below) and reference the Astra GitHub repository.

---

### PyHMMER
> Larralde, M. & Zeller, G. (2023) PyHMMER: a Python library binding to HMMER for efficient sequence analysis. *Bioinformatics* 39(5), btad214. DOI: [10.1093/bioinformatics/btad214](https://doi.org/10.1093/bioinformatics/btad214)

**Used for:** Underlying HMM search engine used by Astra for all profile HMM searches.
**Verified:** [PubMed 37074928](https://pubmed.ncbi.nlm.nih.gov/37074928/)
**Abstract excerpt:** "PyHMMER provides Python integration of the popular profile Hidden Markov Model software HMMER via Cython bindings, allowing the annotation of protein sequences with profile HMMs and building new ones directly with Python."

---

### HMMER
> Eddy, S.R. (2011) Accelerated Profile HMM Searches. *PLoS Computational Biology* 7(10), e1002195. DOI: [10.1371/journal.pcbi.1002195](https://doi.org/10.1371/journal.pcbi.1002195)

**Used for:** Profile HMM searches (via Astra/PyHMMER and direct hmmsearch for giant protein annotation recovery).
**Verified:** [PLoS Comp Biol](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002195)
**Abstract excerpt:** "HMMER3 is as fast as BLAST for protein searches, while retaining the power of using probabilistic inference technology for sensitive homology detection."

---

### DIAMOND
> Buchfink, B., Xie, C. & Huson, D.H. (2015) Fast and sensitive protein alignment using DIAMOND. *Nature Methods* 12, 59--60. DOI: [10.1038/nmeth.3176](https://doi.org/10.1038/nmeth.3176)

**Used for:** Fast protein alignment for hydrogenase subgroup classification (HydDB pipeline, `classify_hydrogenases.py`).
**Verified:** [Nature Methods](https://www.nature.com/articles/nmeth.3176)
**Abstract excerpt:** "DIAMOND is an open-source algorithm based on double indexing that is 20,000 times faster than BLASTX on short reads and has a similar degree of sensitivity."

---

## Annotation Databases

### Pfam
> Mistry, J., Chuguransky, S., Williams, L., Qureshi, M., Salazar, G.A., Sonnhammer, E.L.L., Tosatto, S.C.E., Paladin, L., Raj, S., Richardson, L.J., Finn, R.D. & Bateman, A. (2021) Pfam: The protein families database in 2021. *Nucleic Acids Research* 49(D1), D412--D419. DOI: [10.1093/nar/gkaa913](https://doi.org/10.1093/nar/gkaa913)

**Used for:** Protein domain annotation (Stage 04 via Astra) -- classifies proteins into families and domains using profile HMMs.
**Verified:** [Nucleic Acids Research](https://academic.oup.com/nar/article/49/D1/D412/5943818)
**Abstract excerpt:** "Pfam is a widely used resource for classifying protein sequences into families and domains. Over 350 new families have been added in Pfam 33.1 and numerous improvements have been made to existing entries."

---

### KOfam / KEGG
> Aramaki, T., Blanc-Mathieu, R., Endo, H., Ohkubo, K., Kanehisa, M., Goto, S. & Ogata, H. (2020) KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold. *Bioinformatics* 36(7), 2251--2252. DOI: [10.1093/bioinformatics/btz859](https://doi.org/10.1093/bioinformatics/btz859)

**Used for:** KEGG orthology assignment (Stage 04 via Astra) -- assigns KO identifiers using profile HMMs with adaptive score thresholds.
**Verified:** [Bioinformatics](https://academic.oup.com/bioinformatics/article/36/7/2251/5631907)
**Abstract excerpt:** "KofamKOALA is a web server to assign KEGG Orthologs (KOs) to protein sequences by homology search against a database of profile hidden Markov models (KOfam) with pre-computed adaptive score thresholds."

For the KEGG database itself, also cite:

> Kanehisa, M. & Goto, S. (2000) KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research* 28(1), 27--30. DOI: [10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27)

---

### HydDB
> Sondergaard, D., Pedersen, C.N.S. & Greening, C. (2016) HydDB: A web tool for hydrogenase classification and analysis. *Scientific Reports* 6, 34212. DOI: [10.1038/srep34212](https://doi.org/10.1038/srep34212)

**Used for:** Hydrogenase classification (Stage 04 via Astra, and `classify_hydrogenases.py`) -- classifies NiFe, FeFe, and Fe-only hydrogenases into 38 subgroups.
**Verified:** [Nature Scientific Reports](https://www.nature.com/articles/srep34212)
**Abstract excerpt:** "We present HydDB, a web tool for the classification and analysis of hydrogenase sequences, with an expanded classification scheme comprising 29 [NiFe], 8 [FeFe] and 1 [Fe] hydrogenase classes that defines 11 new classes with distinct biological functions."

---

### VOGdb
> Trgovec-Greif, L., Hellinger, H.J., Mainguy, J., Pfundner, A., Frishman, D., Kiening, M., Webster, N.S., Laffy, P.W., Feichtinger, M. & Rattei, T. (2024) VOGDB -- Database of Virus Orthologous Groups. *Viruses* 16(8), 1191. DOI: [10.3390/v16081191](https://doi.org/10.3390/v16081191)

**Used for:** Viral gene annotation (Stage 04 via Astra) -- assigns Virus Orthologous Groups with functional categories and descriptions.
**Verified:** [PubMed 39205165](https://pubmed.ncbi.nlm.nih.gov/39205165/)
**Abstract excerpt:** "VOGDB is a comprehensive database of Virus Orthologous Groups providing functional annotations for viral proteins."

---

### DefenseFinder
> Tesson, F., Herve, A., Mordret, E., Touchon, M., d'Humieres, C., Cury, J. & Bernheim, A. (2022) Systematic and quantitative view of the antiviral arsenal of prokaryotes. *Nature Communications* 13, 2561. DOI: [10.1038/s41467-022-30269-9](https://doi.org/10.1038/s41467-022-30269-9)

**Used for:** Anti-phage defense system detection (Stage 04 via Astra) -- identifies restriction-modification, CRISPR-Cas, CBASS, BREX, DISARM, Druantia, and other defense systems.
**Verified:** [Nature Communications](https://www.nature.com/articles/s41467-022-30269-9)
**Abstract excerpt:** "DefenseFinder provides a systematic search of all known anti-phage systems in prokaryotic genomes."

---

### CRISPRCasFinder
> Couvin, D., Bernheim, A., Toffano-Nioche, C., Touchon, M., Michalik, J., Neron, B., Rocha, E.P.C., Vergnaud, G., Gautheret, D. & Pourcel, C. (2018) CRISPRCasFinder, an update of CRISRFinder, includes a portable version, enhanced performance and integrates search for Cas proteins. *Nucleic Acids Research* 46(W1), W246--W251. DOI: [10.1093/nar/gky425](https://doi.org/10.1093/nar/gky425)

**Used for:** CRISPR-Cas system detection (Stage 04 via Astra) -- identifies CRISPR arrays and classifies associated Cas proteins.
**Verified:** [Nucleic Acids Research](https://academic.oup.com/nar/article/46/W1/W246/5001162)
**Abstract excerpt:** "CRISPRCasFinder allows the identification of both CRISPR arrays and Cas proteins, with features including an improved CRISPR array detection tool, prediction of CRISPR orientation, and a Cas protein detection and typing tool updated to match the latest classification scheme."

---

### dbCAN
> Zhang, H., Yohe, T., Huang, L., Entwistle, S., Wu, P., Yang, Z., Busk, P.K., Xu, Y. & Yin, Y. (2018) dbCAN2: a meta server for automated carbohydrate-active enzyme annotation. *Nucleic Acids Research* 46(W1), W95--W101. DOI: [10.1093/nar/gky418](https://doi.org/10.1093/nar/gky418)

**Used for:** CAZyme annotation (optional module `dbcan_cazyme.py`) -- identifies carbohydrate-active enzymes using HMM, DIAMOND, and peptide-based methods.
**Verified:** [Nucleic Acids Research](https://academic.oup.com/nar/article/46/W1/W95/4996582)
**Abstract excerpt:** "dbCAN2 is an updated meta server which integrates three state-of-the-art tools for CAZome annotation: (i) HMMER search against the dbCAN HMM database; (ii) DIAMOND search against the CAZy pre-annotated CAZyme sequence database and (iii) Hotpep search against the conserved CAZyme short peptide database."

---

### GECCO
> Carroll, L.M., Larralde, M., Fleck, J.S., Ponnudurai, R., Milanese, A., Cappio Barazzone, E. & Zeller, G. (2021) Accurate de novo identification of biosynthetic gene clusters with GECCO. *bioRxiv* 2021.05.03.442509. DOI: [10.1101/2021.05.03.442509](https://doi.org/10.1101/2021.05.03.442509)

**Used for:** Biosynthetic gene cluster detection (optional module `gecco_bgc.py`) -- identifies BGCs using conditional random fields.
**Verified:** [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.05.03.442509v1)
**Note:** As of February 2026, this paper appears to remain a preprint on bioRxiv and has not been published in a peer-reviewed journal. Cite as a preprint.
**Abstract excerpt:** "GECCO (GEne Cluster prediction with COnditional random fields) is a high-precision, scalable method for identifying novel BGCs in (meta)genomic data using conditional random fields (CRFs)."

---

### MinCED
> No formal publication for MinCED itself. MinCED is derived from CRT; cite the underlying method:
>
> Bland, C., Ramsey, T.L., Sabree, F., Lowe, M., Brown, K., Kyrpides, N.C. & Hugenholtz, P. (2007) CRISPR Recognition Tool (CRT): a tool for automatic detection of clustered regularly interspaced palindromic repeats. *BMC Bioinformatics* 8, 209. DOI: [10.1186/1471-2105-8-209](https://doi.org/10.1186/1471-2105-8-209)
>
> MinCED software: Skennerton, C.T. MinCED: Mining CRISPRs in Environmental Datasets. GitHub: [https://github.com/ctSkennerton/minced](https://github.com/ctSkennerton/minced)

**Used for:** CRISPR array detection (optional module `minced_crispr.py`) -- identifies CRISPR repeat-spacer arrays in genome assemblies.
**Verified:** [GitHub repository](https://github.com/ctSkennerton/minced)
**Note:** MinCED has no dedicated publication. Cite the CRT paper (Bland et al. 2007) as the underlying algorithm and reference the MinCED GitHub repository for the software implementation.

---

## Embedding and Search Tools

### ESM2
> Lin, Z., Akin, H., Rao, R., Hie, B., Zhu, Z., Lu, W., Smetanin, N., Verkuil, R., Kabeli, O., Shmueli, Y., dos Santos Costa, A., Fazel-Zarandi, M., Sercu, T., Candido, S. & Rives, A. (2023) Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* 379(6637), 1123--1130. DOI: [10.1126/science.ade2574](https://doi.org/10.1126/science.ade2574)

**Used for:** Protein embeddings (Stage 06) -- the `esm2_t6_8M_UR50D` model generates 320-dimensional protein representations for semantic similarity search via LanceDB.
**Verified:** [Science](https://www.science.org/doi/10.1126/science.ade2574)
**Abstract excerpt:** "As language models of protein sequences were scaled up to 15 billion parameters, an atomic-resolution picture of protein structure emerged in the learned representations."

---

### LanceDB
> No formal publication. LanceDB is an open-source embedded vector database built on the Lance columnar format.
> GitHub: [https://github.com/lancedb/lancedb](https://github.com/lancedb/lancedb)
> Website: [https://lancedb.com](https://lancedb.com)

**Used for:** Vector storage and kNN search (Stage 06 output) -- stores ESM2 protein embeddings and enables fast similarity search.
**Verified:** [GitHub repository](https://github.com/lancedb/lancedb)
**Note:** LanceDB does not have a formal academic publication. Reference the GitHub repository. LanceDB's core is written in Rust and built using the Lance open-source columnar format.

---

### Foldseek
> van Kempen, M., Kim, S.S., Tumescheit, C., Mirdita, M., Lee, J., Gilchrist, C.L.M., Soding, J. & Steinegger, M. (2024) Fast and accurate protein structure search with Foldseek. *Nature Biotechnology* 42, 243--246. DOI: [10.1038/s41587-023-01773-0](https://doi.org/10.1038/s41587-023-01773-0)

**Used for:** Structural homology search -- searches predicted protein structures against AlphaFold DB, PDB, and other structure databases to find remote homologs undetectable by sequence methods.
**Verified:** [Nature Biotechnology](https://www.nature.com/articles/s41587-023-01773-0)
**Abstract excerpt:** "Foldseek aligns the structure of a query protein against a database by describing tertiary amino acid interactions within proteins as sequences over a structural alphabet, decreasing computation times by four to five orders of magnitude."

---

### ESM3
> Hayes, T., Rao, R., Akin, H., Sofroniew, N.J., Oktay, D., Lin, Z., Verkuil, R., Tran, V.Q., Deaton, J., Wiggert, M., Badkundri, R., Shafkat, I., Gong, J., Derry, A., Molina, R.S., Thomas, N., Khan, Y.A., Mishra, C., Kim, C., Bartie, L.J., Nemeth, M., Hsu, P.D., Sercu, T., Candido, S. & Rives, A. (2025) Simulating 500 million years of evolution with a language model. *Science* 387(6736), 850--858. DOI: [10.1126/science.ads0018](https://doi.org/10.1126/science.ads0018)

**Used for:** Protein structure prediction -- predicts 3D protein structures from sequence via the ESM3/ESMFold API (EvolutionaryScale Forge API).
**Verified:** [Science](https://www.science.org/doi/10.1126/science.ads0018)
**Abstract excerpt:** "ESM3 is a frontier multimodal generative language model that reasons over the sequence, structure, and function of proteins."

---

## Other Tools

### Protenix
> ByteDance AML AI4Science Team. (2025) Protenix -- Advancing Structure Prediction Through a Comprehensive AlphaFold3 Reproduction. *bioRxiv* 2025.01.08.631967. DOI: [10.1101/2025.01.08.631967](https://doi.org/10.1101/2025.01.08.631967)

**Used for:** Protein complex structure prediction -- used for validating protein-protein interaction models (e.g., JAB-ubiquitin complex in GJALLARVIRUS analysis).
**Verified:** [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.01.08.631967v1)
**Note:** Preprint as of February 2026. Protenix is an open-source AlphaFold3 reproduction by ByteDance.
**Abstract excerpt:** "Protenix is a comprehensive reproduction of AlphaFold3, tackling the challenges of predicting complex interactions involving proteins, ligands, and nucleic acids."

---

### pyTMHMM
> No dedicated publication for pyTMHMM. Cite the original TMHMM algorithm:
>
> Krogh, A., Larsson, B., von Heijne, G. & Sonnhammer, E.L.L. (2001) Predicting transmembrane protein topology with a hidden Markov model: application to complete genomes. *Journal of Molecular Biology* 305(3), 567--580. DOI: [10.1006/jmbi.2000.4315](https://doi.org/10.1006/jmbi.2000.4315)
>
> pyTMHMM software: [https://github.com/bosborne/pyTMHMM](https://github.com/bosborne/pyTMHMM)

**Used for:** Transmembrane helix prediction -- predicts transmembrane topology for the `transmembrane_predicted`, `single_pass_membrane`, `multi_pass_membrane` predicates.
**Verified:** [PubMed 11152613](https://pubmed.ncbi.nlm.nih.gov/11152613/)
**Note:** pyTMHMM is a Python reimplementation of TMHMM. Cite the original Krogh et al. 2001 paper for the method and reference the GitHub repository for the software.
**Abstract excerpt:** "We describe and validate a new membrane protein topology prediction method, TMHMM, based on a hidden Markov model. We present a method to predict topology correctly for 97-98% of transmembrane helices."

---

### DNA Features Viewer
> Zulkower, V. & Rosser, S. (2020) DNA Features Viewer: a sequence annotation formatting and plotting library for Python. *Bioinformatics* 36(15), 4350--4352. DOI: [10.1093/bioinformatics/btaa213](https://doi.org/10.1093/bioinformatics/btaa213)

**Used for:** Gene neighborhood visualization -- generates publication-quality gene arrow diagrams for locus figures in reports.
**Verified:** [Bioinformatics](https://academic.oup.com/bioinformatics/article/36/15/4350/5868559)
**Abstract excerpt:** "DNA Features Viewer is a Python library which lets users define visual 'themes' determining the label and display style of each annotation, with annotations automatically laid out to create compact and readable plots."

---

### Pandoc
> MacFarlane, J. Pandoc: a universal document converter. [https://pandoc.org](https://pandoc.org). GitHub: [https://github.com/jgm/pandoc](https://github.com/jgm/pandoc)

**Used for:** Document conversion -- converts Markdown manuscripts to PDF (via xelatex) for publication.
**Verified:** [https://pandoc.org](https://pandoc.org)
**Note:** Pandoc does not have a formal academic publication. It is open-source software created by John MacFarlane (UC Berkeley). Reference the website and/or GitHub repository.

---

## Summary Table

| Tool | Citation | Year | Used In |
|------|----------|------|---------|
| Prodigal | Hyatt et al., *BMC Bioinformatics* | 2010 | Stage 03 (gene calling) |
| QUAST | Gurevich et al., *Bioinformatics* | 2013 | Stage 01 (assembly QC) |
| DuckDB | Raasveldt & Muhleisen, *SIGMOD* | 2019 | Core database engine |
| Astra | GitHub repository (no paper) | -- | Stage 04 (HMM wrapper) |
| PyHMMER | Larralde & Zeller, *Bioinformatics* | 2023 | HMM search engine |
| HMMER | Eddy, *PLoS Comp Biol* | 2011 | Profile HMM searches |
| DIAMOND | Buchfink et al., *Nature Methods* | 2015 | Hydrogenase classification |
| Pfam | Mistry et al., *Nucleic Acids Res* | 2021 | Protein domain annotation |
| KOfam | Aramaki et al., *Bioinformatics* | 2020 | KEGG orthology assignment |
| KEGG | Kanehisa & Goto, *Nucleic Acids Res* | 2000 | Pathway database |
| HydDB | Sondergaard et al., *Sci Rep* | 2016 | Hydrogenase classification |
| VOGdb | Trgovec-Greif et al., *Viruses* | 2024 | Viral gene annotation |
| DefenseFinder | Tesson et al., *Nature Comms* | 2022 | Defense system detection |
| CRISPRCasFinder | Couvin et al., *Nucleic Acids Res* | 2018 | CRISPR-Cas detection |
| dbCAN | Zhang et al., *Nucleic Acids Res* | 2018 | CAZyme annotation |
| GECCO | Carroll et al., *bioRxiv* (preprint) | 2021 | BGC detection |
| MinCED / CRT | Bland et al., *BMC Bioinformatics* | 2007 | CRISPR array detection |
| ESM2 | Lin et al., *Science* | 2023 | Protein embeddings |
| LanceDB | GitHub repository (no paper) | -- | Vector database |
| Foldseek | van Kempen et al., *Nature Biotech* | 2024 | Structural homology search |
| ESM3 | Hayes et al., *Science* | 2025 | Structure prediction |
| Protenix | ByteDance, *bioRxiv* (preprint) | 2025 | Complex structure prediction |
| pyTMHMM / TMHMM | Krogh et al., *J Mol Biol* | 2001 | TM helix prediction |
| DNA Features Viewer | Zulkower & Rosser, *Bioinformatics* | 2020 | Gene diagrams |
| Pandoc | MacFarlane (software, no paper) | -- | Document conversion |

---

*Last updated: 2026-02-08*

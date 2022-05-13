# Create regulatory regions delineation

## Download gene and assembly files for species

```bash
# Import functions.
source download_gene_annotation_and_genome_files.sh
```

```
$ download_refseq_gene_annotation_file
Usage: Parameters:
  - refseq_version:
      Specify refSeq version used on UCSC: e.g.: 90
      If refSeq version is set to 0, the last refSeq version will be used.
      UCSC normally updates gene annotation to the last version of refSeq,
      so in general you should just set this parameter to 0.

  - assembly:
      Specify assembly version: e.g. hg38
      Check the following links for possible assembly version names:
        - http://hgdownload.cse.ucsc.edu/downloads.html
        - ftp://hgdownload.cse.ucsc.edu/goldenPath/

  - refGene|ncbiRefSeq:
      Download refGene or ncbiRefSeq gene annotation.

Examples:
  Download last ncbiRefSeqCurated gene annotation for Homo sapiens (hg38) from UCSC:

  download_refseq_gene_annotation_file 0 hg38 ncbiRefSeqCurated

```

### Examples

#### Human (hg38)

Download last ncbiRefSeqCurated gene annotation for Homo sapiens (hg38) from UCSC:

```bash
download_refseq_gene_annotation_file 0 hg38 ncbiRefSeqCurated
```

Create a list of assemblies to use for lifting over human regulatory regions.
The human regulatory regions and lifted over regions will be scored with
[Cluster-Buster](https://github.com/weng-lab/cluster-buster/) in a later step.

See
[http://hgdownload.cse.ucsc.edu/downloads.html](http://hgdownload.cse.ucsc.edu/downloads.html)
and
[ftp://hgdownload.cse.ucsc.edu/goldenPath/](ftp://hgdownload.cse.ucsc.edu/goldenPath/)
for the correct assembly names to use.

Content of **./data/genomes/liftOver/hg38.liftover_genomes.lst**:
```
bosTau8	cow/Bos taurus	hg38
canFam3	dog/Canis familiaris	hg38
mm10	mouse/Mus musculus	hg38
monDom5	opossums/Monodelphis domestica	hg38
panTro5	chimpanzee/Pan troglodytes	hg38
ponAbe2	Orangutan/Pongo pygmaeus abelii	hg38
rheMac8	rhesus macaque/Macaca mulatta	hg38
rn6	rat/Rattus norvegicus	hg38
susScr3	pig/Sus scrofa	hg38
```

#### Mouse (mm10)

Download last refGene gene annotation for Mus Musculus (mm10) from UCSC:

```bash
download_refseq_gene_annotation_file 0 mm10 refGene
```

Content of **genomes/liftOver/mm10.liftover_genomes.lst**:
```
bosTau8	cow/Bos taurus	mm10
canFam3	dog/Canis familiaris	mm10
hg38	human/Homo sapiens	mm10
monDom5	opossums/Monodelphis domestica	mm10
panTro5	chimpanzee/Pan troglodytes	mm10
ponAbe2	Orangutan/Pongo pygmaeus abelii	mm10
rheMac8	rhesus macaque/Macaca mulatta	mm10
rn6	rat/Rattus norvegicus	mm10
susScr3	pig/Sus scrofa	mm10
```

## Create regulatory regions

```bash

```


## Motif collections

```bash
get_download http://motifcollections.aertslab.org/zips/motif_collection_v9.public.zip motif_collection_v9.public.zip
```
